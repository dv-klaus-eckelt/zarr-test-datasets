#!/usr/bin/env python
# coding: utf-8

# In[1]:


# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "anndata>=0.12.11",
#     "scanpy>=1.12.1",
#     "zarr>=3.1.6",
# ]
#
# [tool.uv]
# exclude-newer = "2026-05-03T06:54:44.427134817+02:00"
# ///


# In[2]:


from pathlib import Path
import json
from datetime import datetime, timezone
from dataclasses import dataclass, asdict
from typing import Optional
import importlib.metadata


# In[3]:


import anndata as ad
import numpy as np
import scanpy as sc
import zarr


# In[4]:


INPUT_PATH = Path("habib17.h5ad")
OUTPUT_PATH = Path("test-data/habib17-differential-expression-test-data-format.zarr")
GROUPBY_COLUMN = "CellType"
MIN_CELLS = 50


# In[4b]:


@dataclass
class ContrastConfig:
    """Configuration for a single DE contrast."""
    contrast_id: str
    group_1: str
    method: str
    corr_method: Optional[str] = None
    test_type: str = "one_vs_rest"
    

def _build_contrast_matrix(cell_types: list[str]) -> list[ContrastConfig]:
    """Build deterministic contrast matrix with 45 contrasts.
    
    15 wilcoxon + benjamini-hochberg (all cell types)
    15 wilcoxon + bonferroni (all cell types)
    15 t-test + benjamini-hochberg (all cell types)
    """
    contrasts = []
    
    # Series 1: wilcoxon + benjamini-hochberg (15 contrasts)
    for i, ct in enumerate(cell_types, 1):
        contrasts.append(ContrastConfig(
            contrast_id=f"de_{i:03d}",
            group_1=ct,
            method="wilcoxon",
            corr_method="benjamini-hochberg",
            test_type="one_vs_rest"
        ))
    
    # Series 2: wilcoxon + bonferroni (15 contrasts)
    for i, ct in enumerate(cell_types, len(contrasts) + 1):
        contrasts.append(ContrastConfig(
            contrast_id=f"de_{i:03d}",
            group_1=ct,
            method="wilcoxon",
            corr_method="bonferroni",
            test_type="one_vs_rest"
        ))
    
    # Series 3: t-test + benjamini-hochberg (15 contrasts)
    for i, ct in enumerate(cell_types, len(contrasts) + 1):
        contrasts.append(ContrastConfig(
            contrast_id=f"de_{i:03d}",
            group_1=ct,
            method="t-test",
            corr_method="benjamini-hochberg",
            test_type="one_vs_rest"
        ))
    
    return contrasts


# In[5]:


def _sample_expression_values(x: object, max_items: int = 10000) -> np.ndarray:
    if hasattr(x, "tocoo"):
        values = np.asarray(x.data)
    else:
        values = np.asarray(x).ravel()

    if values.size == 0:
        return values

    if values.size > max_items:
        step = max(1, values.size // max_items)
        values = values[::step]

    return values


# In[6]:


def _needs_preprocessing(adata: ad.AnnData) -> tuple[bool, str]:
    if "log1p" in adata.uns:
        return False, "Detected adata.uns['log1p']; matrix appears already log-transformed."

    values = _sample_expression_values(adata.X)
    if values.size == 0:
        return True, "Empty matrix sample; applying preprocessing by default."

    tol = 1e-6
    non_integer_fraction = float(np.mean(np.abs(values - np.round(values)) > tol))
    max_value = float(np.max(values))

    # Heuristic: mostly non-integers with compressed range usually indicates logged data.
    if non_integer_fraction > 0.2 and max_value < 50:
        return (
            False,
            (
                "Expression values look already transformed "
                f"(non-integer fraction={non_integer_fraction:.3f}, max={max_value:.3f})."
            ),
        )

    return (
        True,
        (
            "Expression values look like raw counts "
            f"(non-integer fraction={non_integer_fraction:.3f}, max={max_value:.3f})."
        ),
    )


# In[6b]:


def _extract_contrast_arrays(
    adata: ad.AnnData,
    contrast_key: str,
    group_1: str,
    mean_expr_global: np.ndarray,
) -> dict:
    """Extract per-contrast arrays from rank_genes_groups result.
    
    Returns dict with spec-compliant field names and pre-sorted arrays.
    """
    result = adata.uns[contrast_key]
    
    # Extract arrays for the target group
    # Handle both DataFrame and structured array formats
    names_field = result["names"]
    if hasattr(names_field, "columns"):  # DataFrame
        symbols = np.array(names_field[group_1].values)
    else:  # Structured array
        symbols = np.array(names_field[group_1])
    
    scores_field = result["scores"]
    if hasattr(scores_field, "columns"):  # DataFrame
        scores = np.array(scores_field[group_1].values)
    else:  # Structured array
        scores = np.array(scores_field[group_1])
    
    # For logreg, pvals/logfoldchanges will be missing; initialize as NaN
    logfoldchanges_arr = None
    if "logfoldchanges" in result:
        lfc_field = result["logfoldchanges"]
        if hasattr(lfc_field, "columns"):  # DataFrame
            if group_1 in lfc_field.columns:
                logfoldchanges_arr = np.array(lfc_field[group_1].values)
        else:  # Structured array
            if group_1 in lfc_field.dtype.names:
                logfoldchanges_arr = np.array(lfc_field[group_1])
    
    logfoldchanges = logfoldchanges_arr if logfoldchanges_arr is not None else np.full(len(symbols), np.nan)
    
    pvals_arr = None
    if "pvals" in result:
        pvals_field = result["pvals"]
        if hasattr(pvals_field, "columns"):  # DataFrame
            if group_1 in pvals_field.columns:
                pvals_arr = np.array(pvals_field[group_1].values)
        else:  # Structured array
            if group_1 in pvals_field.dtype.names:
                pvals_arr = np.array(pvals_field[group_1])
    
    pvals = pvals_arr if pvals_arr is not None else np.full(len(symbols), np.nan)
    
    pvals_adj_arr = None
    if "pvals_adj" in result:
        pvals_adj_field = result["pvals_adj"]
        if hasattr(pvals_adj_field, "columns"):  # DataFrame
            if group_1 in pvals_adj_field.columns:
                pvals_adj_arr = np.array(pvals_adj_field[group_1].values)
        else:  # Structured array
            if group_1 in pvals_adj_field.dtype.names:
                pvals_adj_arr = np.array(pvals_adj_field[group_1])
    
    pvals_adj = pvals_adj_arr if pvals_adj_arr is not None else np.full(len(symbols), np.nan)
    
    # Get pts (percent expressing in target group)
    pct_expr_target_arr = None
    if "pts" in result:
        pts_field = result["pts"]
        if hasattr(pts_field, "columns"):  # DataFrame
            if group_1 in pts_field.columns:
                pct_expr_target_arr = np.array(pts_field[group_1].values) * 100
        else:  # Structured array
            if group_1 in pts_field.dtype.names:
                pct_expr_target_arr = np.array(pts_field[group_1]) * 100
    
    pct_expr_target = pct_expr_target_arr if pct_expr_target_arr is not None else np.full(len(symbols), np.nan)
    
    # Get pts_rest (percent expressing in reference/rest)
    pct_expr_ref_arr = None
    if "pts_rest" in result:
        pts_rest_field = result["pts_rest"]
        if hasattr(pts_rest_field, "columns"):  # DataFrame
            if group_1 in pts_rest_field.columns:
                pct_expr_ref_arr = np.array(pts_rest_field[group_1].values) * 100
        else:  # Structured array
            if group_1 in pts_rest_field.dtype.names:
                pct_expr_ref_arr = np.array(pts_rest_field[group_1]) * 100
    
    pct_expr_ref = pct_expr_ref_arr if pct_expr_ref_arr is not None else np.full(len(symbols), np.nan)
    
    mean_expr = mean_expr_global  # Global mean expression across all cells
    
    return {
        "symbols": symbols,
        "logfoldchanges": logfoldchanges,
        "pvals": pvals,
        "pvals_adj": pvals_adj,
        "scores": scores,
        "pct_expr_target": pct_expr_target,
        "pct_expr_ref": pct_expr_ref,
        "mean_expr": mean_expr,
    }


def _presort_contrast_arrays(arrays: dict) -> dict:
    """Pre-sort all arrays by descending scores with stable tie-breakers.
    
    Sort order: descending scores, then ascending pvals_adj, then ascending symbol names.
    """
    symbols = arrays["symbols"]
    scores = arrays["scores"]
    pvals_adj = arrays["pvals_adj"]
    
    # Create deterministic sort index (descending scores, stable for ties)
    sort_idx = np.lexsort((
        symbols,  # Symbol names ascending (stable tie-breaker)
        np.where(np.isnan(pvals_adj), np.inf, pvals_adj),  # pvals_adj ascending
        -scores  # Scores descending (note the minus sign)
    ))
    
    # Apply sort to all arrays
    sorted_arrays = {}
    for key, arr in arrays.items():
        sorted_arrays[key] = arr[sort_idx]
    
    return sorted_arrays


def _compute_axis_bounds(
    logfoldchanges: np.ndarray,
    pvals_adj: np.ndarray,
) -> tuple[Optional[float], Optional[float]]:
    """Compute lfc_max and logp_max for chart axis anchoring."""
    # Filter out NaN/Inf values
    valid_lfc = logfoldchanges[~(np.isnan(logfoldchanges) | np.isinf(logfoldchanges))]
    valid_pvals = pvals_adj[~np.isnan(pvals_adj)]
    
    lfc_max = None
    logp_max = None
    
    if len(valid_lfc) > 0:
        lfc_max = float(np.ceil(np.max(np.abs(valid_lfc))))
    
    if len(valid_pvals) > 0:
        # Compute -log10(p) and cap at 300
        logp_vals = -np.log10(np.maximum(valid_pvals, 1e-300))
        logp_max = min(300.0, float(np.ceil(np.max(logp_vals))))
    
    return lfc_max, logp_max


def _write_contrast_to_zarr(
    zarr_path: Path,
    contrast_id: str,
    arrays: dict,
    contrast_config: ContrastConfig,
) -> dict:
    """Write pre-sorted contrast arrays to zarr format with spec field names.
    
    Returns metadata dict for contrast_registry.
    """
    de_dir = zarr_path / "uns" / "de"
    contrast_dir = de_dir / contrast_id
    contrast_dir.mkdir(parents=True, exist_ok=True)
    
    store = zarr.open(str(contrast_dir), mode="w")
    
    # Write each array to zarr with explicit dtype handling for strings
    for field_name, arr in arrays.items():
        if field_name == "symbols" and arr.dtype.kind == "O":
            # Convert object array of strings to fixed-length unicode strings
            max_len = max(len(s) for s in arr) if len(arr) > 0 else 1
            arr_str = arr.astype(f"U{max_len}")
            store[field_name] = arr_str
        else:
            store[field_name] = arr
    
    # Compute chart bounds for registry
    lfc_max, logp_max = _compute_axis_bounds(
        arrays["logfoldchanges"],
        arrays["pvals_adj"],
    )
    
    # Determine available plots based on method
    if contrast_config.method == "logreg":
        available_plots = ["bar_chart", "dotplot"]
        y_axis_label = None
        logp_max = None
    else:
        available_plots = ["volcano", "dotplot"]
        y_axis_label = "-log10(FDR)" if contrast_config.corr_method else "-log10(p)"
    
    return {
        "contrast_id": contrast_id,
        "group_1": contrast_config.group_1,
        "group_2": None,
        "test_type": contrast_config.test_type,
        "de_method": contrast_config.method,
        "correction_method": contrast_config.corr_method,
        "feature_type": "Gene Expression",
        "x_axis_label": "log2(Fold Change)",
        "y_axis_label": y_axis_label,
        "lfc_max": lfc_max,
        "logp_max": logp_max,
        "n_features": len(arrays["symbols"]),
        "available_plots": available_plots,
    }


# In[7]:


def main() -> None:
    if not INPUT_PATH.exists():
        raise FileNotFoundError(f"Input file not found: {INPUT_PATH}")

    adata = ad.read_h5ad(INPUT_PATH)

    if GROUPBY_COLUMN not in adata.obs.columns:
        available = ", ".join(adata.obs.columns.astype(str).tolist())
        raise KeyError(
            f"Column '{GROUPBY_COLUMN}' not found in obs. Available columns: {available}"
        )

    # Align with the zarr v3 output style used in other dataset creation scripts.
    ad.settings.zarr_write_format = 3
    ad.settings.write_csr_csc_indices_with_min_possible_dtype = True
    ad.settings.auto_shard_zarr_v3 = True

    # Make sure group labels are categorical for rank_genes_groups.
    adata.obs[GROUPBY_COLUMN] = adata.obs[GROUPBY_COLUMN].astype("category")

    should_preprocess, reason = _needs_preprocessing(adata)
    print(reason)

    if should_preprocess:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        print("Applied preprocessing: normalize_total + log1p")
    else:
        print("Skipped preprocessing")

    # Get deterministic cell type list
    cell_types = sorted(adata.obs[GROUPBY_COLUMN].unique().astype(str).tolist())
    print(f"Found {len(cell_types)} cell types: {cell_types}")
    
    # Build contrast matrix
    contrasts = _build_contrast_matrix(cell_types)
    print(f"Built contrast matrix with {len(contrasts)} contrasts")
    
    # Compute global mean expression per gene (for all cells)
    if hasattr(adata.X, "mean"):
        mean_expr_global = np.asarray(adata.X.mean(axis=0)).ravel()
    else:
        mean_expr_global = np.asarray(adata.X).mean(axis=0).ravel()
    
    # Prepare output directory
    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    de_dir = OUTPUT_PATH / "uns" / "de"
    de_dir.mkdir(parents=True, exist_ok=True)
    
    # Write full AnnData to zarr (preserves X, obs, var, etc.)
    # This will be overwritten partially for uns/de, but preserves all non-DE data
    adata.write_zarr(OUTPUT_PATH)
    print(f"Wrote AnnData to {OUTPUT_PATH}")
    
    # Process each contrast
    registry = []
    processed_count = 0
    
    for contrast_config in contrasts:
        try:
            # Validate group size
            group_size = (adata.obs[GROUPBY_COLUMN] == contrast_config.group_1).sum()
            if group_size < MIN_CELLS:
                print(f"  ⊘ {contrast_config.contrast_id}: group '{contrast_config.group_1}' has {group_size} cells < {MIN_CELLS}")
                continue
            
            # Run rank_genes_groups
            key_added = f"_contrast_{contrast_config.contrast_id}"
            
            # Build rank_genes_groups kwargs (logreg doesn't support corr_method)
            rank_kwargs = {
                "adata": adata,
                "groupby": GROUPBY_COLUMN,
                "groups": [contrast_config.group_1],
                "method": contrast_config.method,
                "tie_correct": True,
                "pts": True,
                "key_added": key_added,
                "use_raw": False,
            }
            
            # Only add corr_method for methods that support it
            if contrast_config.corr_method is not None:
                rank_kwargs["corr_method"] = contrast_config.corr_method
            
            sc.tl.rank_genes_groups(**rank_kwargs)
            
            # Extract per-contrast arrays
            arrays = _extract_contrast_arrays(
                adata,
                key_added,
                contrast_config.group_1,
                mean_expr_global,
            )
            
            # Pre-sort all arrays
            sorted_arrays = _presort_contrast_arrays(arrays)
            
            # Write to zarr
            metadata = _write_contrast_to_zarr(
                OUTPUT_PATH,
                contrast_config.contrast_id,
                sorted_arrays,
                contrast_config,
            )
            
            # Add to registry
            registry.append(metadata)
            processed_count += 1
            
            print(f"  ✓ {contrast_config.contrast_id}: {contrast_config.method} + {contrast_config.corr_method} ({group_size} cells)")
            
            # Clean up temporary key
            del adata.uns[key_added]
            
        except Exception as e:
            print(f"  ✗ {contrast_config.contrast_id}: {e}")
            continue
    
    # Write contrast_registry.json
    registry_json = {
        "version": "1.0",
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "scanpy_version": importlib.metadata.version("scanpy"),
        "zarr_version": importlib.metadata.version("zarr"),
        "dataset_metadata": {
            "n_obs": adata.n_obs,
            "n_vars": adata.n_vars,
            "preprocessing": "normalize_total(target_sum=1e4) + log1p" if should_preprocess else "none",
            "groupby_column": GROUPBY_COLUMN,
            "n_groups": len(cell_types),
        },
        "contrast_summary": {
            "total_contrasts": processed_count,
            "series": [
                {"series_id": 1, "method": "wilcoxon", "corr": "benjamini-hochberg", "n_contrasts": 15},
                {"series_id": 2, "method": "wilcoxon", "corr": "bonferroni", "n_contrasts": 15},
                {"series_id": 3, "method": "t-test", "corr": "benjamini-hochberg", "n_contrasts": 15},
            ],
        },
        "contrasts": registry,
    }
    
    registry_path = OUTPUT_PATH / "uns" / "de" / "contrast_registry.json"
    with open(str(registry_path), "w") as f:
        json.dump(registry_json, f, indent=2)
    
    print(f"\nWrote {processed_count} DE contrasts to {OUTPUT_PATH}/uns/de/")
    print(f"Registry: {registry_path}")


# In[8]:


if __name__ == "__main__":
    main()




