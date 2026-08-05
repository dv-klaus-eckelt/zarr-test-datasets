#!/usr/bin/env python3
"""
Transform the raw Visium brain SpatialData store into the trimmed test dataset
in test-data/visium_brain.spatialdata.zarr.

Source dataset ("Visium HD Mouse Brain" / visium_associated_xenium_io, two
serial sections ST8059048 + ST8059050), downloadable from the SpatialData
example-datasets index:
    https://spatialdata.scverse.org/en/stable/tutorials/notebooks/datasets/README.html

The raw store holds, per section, a single-scale `*_hires_image` and a
single-scale `*_lowres_image`. This script:
  1. drops the `*_lowres_image` elements - they carry an empty
     coordinateTransformations list, so they belong to no coordinate system and
     nothing (no shapes, no table region) references them;
  2. re-parses each `*_hires_image` into a 4-level multiscale pyramid
     (scale factor 2 per level) with tiled chunks, preserving the original
     pixels at scale0, the identity transformation into the section's
     coordinate system, and the omero r/g/b channel labels;
  3. renames `<section>_hires_image` -> `<section>_image`. Note the bare
     section name is unavailable: SpatialData requires element names to be
     unique across element types, and `<section>` is already the shapes
     element that the table's `region` column points at. The `_image` suffix
     also keeps the image.startswith(region) convention used by the other
     datasets in this repo;
  4. writes a fresh store with the current spatialdata format, i.e. Zarr v3
     with OME-NGFF 0.5-dev-spatialdata metadata.

Shapes (the Visium spot circles) and the table are passed through untouched.

Usage:
    python transform_visium_brain.py --input /tmp/visium_brain_raw.zarr \\
        --output test-data/visium_brain.spatialdata.zarr

    python transform_visium_brain.py --input /tmp/visium_brain_raw.zarr \\
        --output test-data/visium_brain.spatialdata.zarr --overwrite
"""

from __future__ import annotations

import argparse
import sys
import warnings
from pathlib import Path

warnings.filterwarnings("ignore")

import spatialdata as sd
from spatialdata.models import Image2DModel
from spatialdata.transformations import get_transformation

# One extra pyramid level per entry, each halving y and x.
SCALE_FACTORS = [2, 2, 2]

# Tiled chunks so a viewer can fetch single tiles per level; the raw store uses
# one chunk per whole image, which defeats the point of a pyramid.
CHUNKS = (3, 512, 512)

LOWRES_SUFFIX = "_lowres_image"
HIRES_SUFFIX = "_hires_image"
NEW_SUFFIX = "_image"


def transform(input_path: Path, output_path: Path, overwrite: bool) -> None:
    src = sd.read_zarr(input_path)

    images = {}
    for name, image in src.images.items():
        if name.endswith(LOWRES_SUFFIX):
            print(f"dropping {name}")
            continue
        if not name.endswith(HIRES_SUFFIX):
            print(f"keeping {name} as-is (unexpected image name)")
            images[name] = image
            continue

        new_name = name[: -len(HIRES_SUFFIX)] + NEW_SUFFIX
        transformations = dict(get_transformation(image, get_all=True))
        c_coords = [str(c) for c in image.coords["c"].values]
        images[new_name] = Image2DModel.parse(
            image.data,
            dims=("c", "y", "x"),
            c_coords=c_coords,
            transformations=transformations,
            scale_factors=SCALE_FACTORS,
            chunks=CHUNKS,
        )
        print(
            f"{name} -> {new_name}: "
            f"{len(SCALE_FACTORS) + 1} levels, channels={c_coords}, "
            f"coordinate systems={list(transformations)}"
        )

    dst = sd.SpatialData(
        images=images,
        shapes=dict(src.shapes),
        tables=dict(src.tables),
    )
    dst.write(output_path, overwrite=overwrite)
    print(f"\nwrote {output_path}")
    print(dst)


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    parser.add_argument(
        "--input",
        required=True,
        type=Path,
        help="raw Visium brain SpatialData store (single-scale hires + lowres images)",
    )
    parser.add_argument(
        "--output",
        default=Path("test-data/visium_brain.spatialdata.zarr"),
        type=Path,
        help="destination store (default: %(default)s)",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="overwrite the destination store if it already exists",
    )
    args = parser.parse_args()

    if not args.input.exists():
        print(f"error: input store not found: {args.input}", file=sys.stderr)
        return 1
    if args.output.exists() and not args.overwrite:
        print(
            f"error: output store already exists: {args.output} (pass --overwrite)",
            file=sys.stderr,
        )
        return 1
    if args.input.resolve() == args.output.resolve():
        print("error: --input and --output must differ", file=sys.stderr)
        return 1

    transform(args.input, args.output, args.overwrite)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())



# analysis example: https://spatialdata.scverse.org/en/latest/tutorials/notebooks/notebooks/examples/technology_visium.html
# mapping of cell types: https://cell2location.readthedocs.io/en/latest/notebooks/cell2location_short_demo.html