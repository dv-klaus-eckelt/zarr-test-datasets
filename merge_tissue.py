#!/usr/bin/env python3
"""
Merge multiple S-BIAD2146 SpatialData slides into one tissue-level
SpatialData Zarr v3 store.

Renames image/shape/points keys per slide, adds tissue/disease/patient_status
metadata, updates region columns, concatenates tables, and writes a single
merged SpatialData store.

Usage:
    python merge_tissue.py --input /tmp/slides --metadata sbiad2146_metadata.json \\
        --tissue skin --output test-data/skin.spatialdata.zarr

    python merge_tissue.py --input /tmp/slides --metadata sbiad2146_metadata.json \\
        --slides skin_s1,skin_s2 --output test-data/custom.zarr

    python merge_tissue.py --input /tmp/slides --metadata sbiad2146_metadata.json \\
        --tissue skin --output test-data/skin.spatialdata.zarr --clean
"""

from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
import warnings
from pathlib import Path
from typing import Any, Dict, List

warnings.filterwarnings("ignore")

import anndata as ad
import numpy as np
import pandas as pd
import spatialdata as sd
from spatialdata import SpatialData
from spatialdata.models import TableModel


# ---------------------------------------------------------------------------
# Core merge logic
# ---------------------------------------------------------------------------


def merge_slides(
    slide_dir: Path,
    slide_ids: List[str],
    metadata: Dict[str, Dict[str, str]],
    output_path: Path,
    clean: bool = False,
) -> Path:
    """
    Merge multiple SpatialData slides into one tissue-level SpatialData Zarr v3 store.

    Strategy:
    1. Load each slide, rename coordinate systems to avoid conflicts.
    2. Rename all element keys: ``{original_key}_{slide_id}``.
    3. Add metadata columns (tissue, disease, slide_id, patient_status) to all table obs.
    4. Update region column in linked tables to match renamed shape keys.
    5. Concatenate tables of the same type across slides.
    6. Create a new SpatialData with all elements and write as Zarr v3.

    Args:
        slide_dir: Directory containing unzipped slide directories
                   (e.g., ``sdata_skin_s1.zarr/``).
        slide_ids: List of slide IDs to merge (e.g., ``["skin_s1", "skin_s2"]``).
        metadata: Slide metadata dictionary keyed by slide ID.
        output_path: Output path for the merged Zarr v3 store.
        clean: If True, delete input slide directories after a successful merge.

    Returns:
        The output path.
    """
    total = len(slide_ids)

    # ------------------------------------------------------------------
    # Validate inputs
    # ------------------------------------------------------------------
    for sid in slide_ids:
        slide_path = slide_dir / f"sdata_{sid}.zarr"
        if not slide_path.exists():
            print(
                f"[ERROR] Slide directory not found: {slide_path}",
                file=sys.stderr,
            )
            sys.exit(1)
        if sid not in metadata:
            print(
                f"[ERROR] Slide '{sid}' not found in metadata",
                file=sys.stderr,
            )
            sys.exit(1)

    print(f"Merging {total} slide(s) into: {output_path}")
    print(f"Source directory: {slide_dir}")
    print()

    # ------------------------------------------------------------------
    # Collect elements from all slides
    # ------------------------------------------------------------------
    all_images: Dict[str, Any] = {}
    all_labels: Dict[str, Any] = {}
    all_shapes: Dict[str, Any] = {}
    all_points: Dict[str, Any] = {}
    tables_by_type: Dict[str, List[ad.AnnData]] = {}

    for i, slide_id in enumerate(slide_ids):
        slide_path = slide_dir / f"sdata_{slide_id}.zarr"
        print(f"Processing slide {i + 1}/{total}: {slide_id} ...")

        sdata = sd.read_zarr(str(slide_path))
        slide_info = metadata[slide_id]

        # --------------------------------------------------------------
        # Rename coordinate systems to avoid collisions across slides
        # --------------------------------------------------------------
        cs_rename = {cs: f"{cs}_{slide_id}" for cs in sdata.coordinate_systems}
        if cs_rename:
            sdata.rename_coordinate_systems(cs_rename)

        # --------------------------------------------------------------
        # Rename images, shapes, points, labels — prefix with slide_id
        # --------------------------------------------------------------
        shape_rename_map: Dict[str, str] = {}

        for img_name, img in sdata.images.items():
            new_name = f"{img_name}_{slide_id}"
            all_images[new_name] = img

        for shape_name, shape in sdata.shapes.items():
            new_name = f"{shape_name}_{slide_id}"
            all_shapes[new_name] = shape
            shape_rename_map[shape_name] = new_name

        for pt_name, pt in sdata.points.items():
            new_name = f"{pt_name}_{slide_id}"
            all_points[new_name] = pt

        for label_name, label in sdata.labels.items():
            new_name = f"{label_name}_{slide_id}"
            all_labels[new_name] = label

        # --------------------------------------------------------------
        # Process tables: add metadata, update region, collect by type
        # --------------------------------------------------------------
        for table_name, table in sdata.tables.items():
            # Add tissue-level metadata columns to each cell
            table.obs["tissue"] = slide_info["tissue"]
            table.obs["disease"] = slide_info["disease"]
            table.obs["slide_id"] = slide_id
            table.obs["patient_status"] = slide_info["patient_status"]

            # Update region column to match renamed shape keys
            if "spatialdata_attrs" in table.uns:
                attrs = table.uns["spatialdata_attrs"]
                region_key = attrs.get("region_key", "region")
                if region_key in table.obs.columns:
                    old_regions = table.obs[region_key].astype(str)
                    new_regions = old_regions.map(
                        lambda r: shape_rename_map.get(r, r)
                    )
                    table.obs[region_key] = new_regions.astype("category")

            tables_by_type.setdefault(table_name, []).append(table)

        print(f"  → images: {len(sdata.images)}, shapes: {len(sdata.shapes)}, "
              f"points: {len(sdata.points)}, tables: {len(sdata.tables)}")

    # ------------------------------------------------------------------
    # Concatenate tables of the same type
    # ------------------------------------------------------------------
    print()
    merged_tables: Dict[str, ad.AnnData] = {}
    all_shape_names = list(all_shapes.keys())

    for tname, tables in tables_by_type.items():
        print(f"Merging table '{tname}' from {len(tables)} slide(s) ...")

        if len(tables) == 1:
            merged = tables[0]
        else:
            # Remove spatialdata_attrs before concatenation to avoid uns
            # merge conflicts (region values differ per slide).
            has_spatial_attrs = any("spatialdata_attrs" in t.uns for t in tables)
            saved_attrs: Dict[str, Any] = {}
            if has_spatial_attrs:
                # Save attrs from the first table as reference
                first_table_with_attrs = next(
                    (t for t in tables if "spatialdata_attrs" in t.uns), None
                )
                if first_table_with_attrs is not None:
                    saved_attrs = dict(first_table_with_attrs.uns["spatialdata_attrs"])
                for t in tables:
                    t.uns.pop("spatialdata_attrs", None)

            merged = ad.concat(tables, axis=0, join="outer", merge="first")

            # Reset obs index after concatenation (duplicate indices may exist)
            merged.obs = merged.obs.reset_index(drop=True)

            # Re-parse with TableModel if this is a linked table
            if saved_attrs:
                region_key = saved_attrs.get("region_key", "region")
                instance_key = saved_attrs.get("instance_key", "cell_id")

                # Determine which shapes this table references
                if region_key in merged.obs.columns:
                    region_values = sorted(
                        set(merged.obs[region_key].astype(str).unique())
                    )
                else:
                    region_values = []

                merged_tables[tname] = TableModel.parse(
                    merged,
                    region=region_values,
                    region_key=region_key,
                    instance_key=instance_key,
                    overwrite_metadata=True,
                )
            else:
                merged_tables[tname] = merged

        # Handle single-table case (also needs re-parsing)
        if tname not in merged_tables:
            if "spatialdata_attrs" in merged.uns:
                attrs = merged.uns["spatialdata_attrs"]
                region_key = attrs.get("region_key", "region")
                instance_key = attrs.get("instance_key", "cell_id")
                if region_key in merged.obs.columns:
                    region_values = sorted(
                        set(merged.obs[region_key].astype(str).unique())
                    )
                else:
                    region_values = []
                merged_tables[tname] = TableModel.parse(
                    merged,
                    region=region_values,
                    region_key=region_key,
                    instance_key=instance_key,
                    overwrite_metadata=True,
                )
            else:
                merged_tables[tname] = merged

        n_cells = merged_tables[tname].n_obs
        print(f"  → {tname}: {n_cells} cells")

    # ------------------------------------------------------------------
    # Build merged SpatialData
    # ------------------------------------------------------------------
    print()
    print("Building merged SpatialData ...")

    # Sanitize obsm column names before building SpatialData (Zarr v3 forbids "/")
    merged_tables = _sanitize_obsm_columns(merged_tables)

    all_elements: Dict[str, Any] = {}
    all_elements.update(all_images)
    all_elements.update(all_labels)
    all_elements.update(all_shapes)
    all_elements.update(all_points)
    all_elements.update(merged_tables)

    merged_sdata = SpatialData.init_from_elements(all_elements)

    # ------------------------------------------------------------------
    # Write as Zarr v3
    # ------------------------------------------------------------------
    print(f"Writing Zarr v3 store to: {output_path}")
    merged_sdata.write(str(output_path), overwrite=True, sdata_formats=[])

    print(f"[DONE] Merged {total} slide(s) into: {output_path}")
    print(f"  Images:   {len(all_images)}")
    print(f"  Labels:   {len(all_labels)}")
    print(f"  Shapes:   {len(all_shapes)}")
    print(f"  Points:   {len(all_points)}")
    print(f"  Tables:   {len(merged_tables)}")

    # ------------------------------------------------------------------
    # Clean up input slide directories if requested
    # ------------------------------------------------------------------
    if clean:
        print()
        print("Cleaning up input slide directories ...")
        for sid in slide_ids:
            slide_path = slide_dir / f"sdata_{sid}.zarr"
            if slide_path.exists():
                shutil.rmtree(slide_path)
                print(f"  Removed: {slide_path}")

            zip_path = slide_dir / f"sdata_{sid}.zarr.zip"
            if zip_path.exists():
                zip_path.unlink()
                print(f"  Removed: {zip_path}")

    return output_path


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _sanitize_obsm_columns(tables_dict):
    """Replace forward slashes in obsm DataFrame column names (Zarr v3 limitation)."""
    replaced = 0
    for tname, table in tables_dict.items():
        for key in list(table.obsm.keys()):
            val = table.obsm[key]
            if isinstance(val, pd.DataFrame):
                cols = val.columns.tolist()
                new_cols = [str(c).replace("/", "_org_") for c in cols]
                if new_cols != cols:
                    table.obsm[key].columns = new_cols
                    replaced += sum(1 for oc, nc in zip(cols, new_cols) if oc != nc)
    if replaced:
        print(f"  Sanitized {replaced} obsm column name(s) with '/' → '_org_'")
    return tables_dict


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=(
            "Merge multiple S-BIAD2146 SpatialData slides into one tissue-level "
            "SpatialData Zarr v3 store."
        ),
    )
    parser.add_argument(
        "--input",
        type=str,
        required=True,
        metavar="DIR",
        help="Directory containing the unzipped slide directories "
        "(e.g., /tmp/sbiad2146_skin/)",
    )
    parser.add_argument(
        "--metadata",
        type=str,
        required=True,
        metavar="JSON",
        help="Path to sbiad2146_metadata.json",
    )
    parser.add_argument(
        "--tissue",
        type=str,
        metavar="TISSUE",
        help="Tissue name to merge all slides for that tissue (e.g., 'skin')",
    )
    parser.add_argument(
        "--slides",
        type=str,
        metavar="SLIDES",
        help="Comma-separated list of specific slide IDs "
        "(e.g., 'skin_s1,skin_s2') — alternative to --tissue",
    )
    parser.add_argument(
        "--output",
        type=str,
        required=True,
        metavar="PATH",
        help="Output path for the merged Zarr v3 store "
        "(e.g., test-data/skin.spatialdata.zarr)",
    )
    parser.add_argument(
        "--clean",
        action="store_true",
        help="Delete input slide directories after successful merge",
    )
    return parser


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    # ------------------------------------------------------------------
    # Resolve input arguments
    # ------------------------------------------------------------------
    if not args.tissue and not args.slides:
        parser.error("Either --tissue or --slides must be specified")
    if args.tissue and args.slides:
        parser.error("Use --tissue or --slides, not both")

    input_dir = Path(args.input).resolve()
    if not input_dir.is_dir():
        print(f"[ERROR] Input directory not found: {input_dir}", file=sys.stderr)
        sys.exit(1)

    metadata_path = Path(args.metadata).resolve()
    if not metadata_path.is_file():
        print(f"[ERROR] Metadata file not found: {metadata_path}", file=sys.stderr)
        sys.exit(1)

    with open(metadata_path) as f:
        all_metadata: Dict[str, Dict[str, str]] = json.load(f).get("slides", {})

    if args.tissue:
        slide_ids = sorted(
            sid for sid, info in all_metadata.items()
            if info["tissue"] == args.tissue
        )
        if not slide_ids:
            available = sorted({info["tissue"] for info in all_metadata.values()})
            print(
                f"[ERROR] No slides found for tissue '{args.tissue}'. "
                f"Available tissues: {', '.join(available)}",
                file=sys.stderr,
            )
            sys.exit(1)
    else:
        slide_ids = [s.strip() for s in args.slides.split(",") if s.strip()]
        # Validate all slide IDs exist in metadata
        for sid in slide_ids:
            if sid not in all_metadata:
                print(
                    f"[ERROR] Slide '{sid}' not found in metadata. "
                    f"Available slides: {', '.join(sorted(all_metadata.keys()))}",
                    file=sys.stderr,
                )
                sys.exit(1)

    output_path = Path(args.output).resolve()

    # ------------------------------------------------------------------
    # Run merge
    # ------------------------------------------------------------------
    merge_slides(
        slide_dir=input_dir,
        slide_ids=slide_ids,
        metadata=all_metadata,
        output_path=output_path,
        clean=args.clean,
    )


if __name__ == "__main__":
    main()
