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
from dataclasses import dataclass
from typing import Optional
import importlib.metadata


# In[3]:


import anndata as ad
import numpy as np
import scanpy as sc
import zarr
from zarr.core.dtype import VariableLengthUTF8


# In[4]:


INPUT_PATH = Path("habib17.h5ad")
OUTPUT_PATH = Path("test-data/habib17-differential-expression-test-data-format.zarr")
GROUPBY_COLUMN = "CellType"
MIN_CELLS = 5  # Minimum cells required in group_1 to run DE contrast (relaxed for testing)


# In[4b]:


@dataclass
class ContrastConfig:
    """Configuration for a single DE contrast."""
    contrast_id: str
    group_1: str
    method: str
    corr_method: Optional[str] = None
    test_type: str = "one_vs_rest"
    subset_column: Optional[str] = None
    subset_value: Optional[str] = None
    group_2: Optional[str] = None
    

def _build_contrast_matrix(
    cell_types: list[str],
    hippocampus_cell_types: list[str],
) -> list[ContrastConfig]:
    """Build deterministic contrast matrix with base and subset contrasts.
    
    15 wilcoxon + benjamini-hochberg (all cell types)
    15 t-test + bonferroni (all cell types)
    15 logreg multiclass entries (all cell types)
    N wilcoxon + benjamini-hochberg on subset region=Hippocampus
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
    
    # Series 2: t-test + bonferroni (15 contrasts)
    for i, ct in enumerate(cell_types, len(contrasts) + 1):
        contrasts.append(ContrastConfig(
            contrast_id=f"de_{i:03d}",
            group_1=ct,
            method="t-test",
            corr_method="bonferroni",
            test_type="one_vs_rest"
        ))
    
    # Series 3: multiclass logreg (15 entries, one folder per group)
    for i, ct in enumerate(cell_types, len(contrasts) + 1):
        contrasts.append(ContrastConfig(
            contrast_id=f"de_{i:03d}",
            group_1=ct,
            method="logreg",
            corr_method=None,
            test_type="multiclass"
        ))

    # Series 4: subset region=Hippocampus, one-vs-rest Wilcoxon.
    for i, ct in enumerate(hippocampus_cell_types, len(contrasts) + 1):
        contrasts.append(ContrastConfig(
            contrast_id=f"de_{i:03d}",
            group_1=ct,
            method="wilcoxon",
            corr_method="benjamini-hochberg",
            test_type="one_vs_rest",
            subset_column="region",
            subset_value="Hippocampus",
        ))
    
    return contrasts


def _map_region(cell_id: str) -> str:
    cell_id = cell_id.lower()
    if "pfc" in cell_id:
        return "PFC"
    if "hp" in cell_id:
        return "Hippocampus"
    return "Other"


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
    group_2: Optional[str] = None,
) -> dict:
    """Extract per-contrast arrays from rank_genes_groups result.
    
    Returns dict with spec-compliant field names and pre-sorted arrays.
    """
    result = adata.uns[contrast_key]
    
    # Extract arrays for the target group.
    # Scanpy rank_genes_groups names are based on adata.var_names.
    # Handle both DataFrame and structured array formats
    names_field = result["names"]
    if hasattr(names_field, "columns"):  # DataFrame
        gene_id = np.array(names_field[group_1].values)
    else:  # Structured array
        gene_id = np.array(names_field[group_1])
    
    scores_field = result["scores"]
    if hasattr(scores_field, "columns"):  # DataFrame
        scores = np.array(scores_field[group_1].values)
    else:  # Structured array
        scores = np.array(scores_field[group_1])
    
    # For logreg, pvals/effect-size source values may be missing; initialize as NaN
    effect_size_arr = None
    if "logfoldchanges" in result:
        lfc_field = result["logfoldchanges"]
        if hasattr(lfc_field, "columns"):  # DataFrame
            if group_1 in lfc_field.columns:
                effect_size_arr = np.array(lfc_field[group_1].values)
        else:  # Structured array
            if group_1 in lfc_field.dtype.names:
                effect_size_arr = np.array(lfc_field[group_1])
    else:
        print(f"Warning: 'logfoldchanges' field not found in result for contrast {contrast_key}. logfoldchanges will be NaN.")

    effect_size = effect_size_arr if effect_size_arr is not None else np.full(len(gene_id), np.nan)
    
    pvals_arr = None
    if "pvals" in result:
        pvals_field = result["pvals"]
        if hasattr(pvals_field, "columns"):  # DataFrame
            if group_1 in pvals_field.columns:
                pvals_arr = np.array(pvals_field[group_1].values)
        else:  # Structured array
            if group_1 in pvals_field.dtype.names:
                pvals_arr = np.array(pvals_field[group_1])
    else:
        print(f"Warning: 'pvals' field not found in result for contrast {contrast_key}. pvals will be NaN.")
    
    pvals = pvals_arr if pvals_arr is not None else np.full(len(gene_id), np.nan)
    
    pvals_adj_arr = None
    if "pvals_adj" in result:
        pvals_adj_field = result["pvals_adj"]
        if hasattr(pvals_adj_field, "columns"):  # DataFrame
            if group_1 in pvals_adj_field.columns:
                pvals_adj_arr = np.array(pvals_adj_field[group_1].values)
        else:  # Structured array
            if group_1 in pvals_adj_field.dtype.names:
                pvals_adj_arr = np.array(pvals_adj_field[group_1])
    else:
        print(f"Warning: 'pvals_adj' field not found in result for contrast {contrast_key}. pvals_adj will be NaN.")
    
    pvals_adj = pvals_adj_arr if pvals_adj_arr is not None else np.full(len(gene_id), np.nan)
    
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
    else:
        print(f"Warning: 'pts' field not found in result for contrast {contrast_key}. pct_expr_target will be NaN.")

    
    pct_expr_target = pct_expr_target_arr if pct_expr_target_arr is not None else np.full(len(gene_id), np.nan)
    
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
    else:
        print(f"Warning: 'pts_rest' field not found in result for contrast {contrast_key}. pct_expr_ref will be NaN.")
    
    pct_expr_ref = pct_expr_ref_arr if pct_expr_ref_arr is not None else np.full(len(gene_id), np.nan)
    
    mean_expr = mean_expr_global  # Global mean expression across all cells

    # Compute mean expression for target group (cells where GROUPBY_COLUMN == group_1)
    target_mask = np.asarray(adata.obs[GROUPBY_COLUMN] == group_1)
    n_target = int(target_mask.sum())
    if n_target > 0:
        X_target = adata.X[target_mask]
        if hasattr(X_target, "mean"):
            mean_expr_target = np.asarray(X_target.mean(axis=0)).ravel()
        else:
            mean_expr_target = np.asarray(X_target).mean(axis=0).ravel()
    else:
        mean_expr_target = np.full(adata.n_vars, np.nan)

    # Compute mean expression for reference/rest group.
    # When group_2 is specified (pairwise contrast), restrict to those cells only.
    # When group_2 is None (one_vs_rest), use all cells NOT in the target group.
    if group_2 is not None:
        ref_mask = np.asarray(adata.obs[GROUPBY_COLUMN] == group_2)
    else:
        ref_mask = ~target_mask
    n_ref = int(ref_mask.sum())
    if n_ref > 0:
        X_ref = adata.X[ref_mask]
        if hasattr(X_ref, "mean"):
            mean_expr_ref = np.asarray(X_ref.mean(axis=0)).ravel()
        else:
            mean_expr_ref = np.asarray(X_ref).mean(axis=0).ravel()
    else:
        mean_expr_ref = np.full(adata.n_vars, np.nan)

    return {
        "gene_id": gene_id,
        "effect_size": effect_size,
        "pvals": pvals,
        "pvals_adj": pvals_adj,
        "scores": scores,
        "pct_expr_target": pct_expr_target,
        "pct_expr_ref": pct_expr_ref,
        "mean_expr": mean_expr,
        "mean_expr_target": mean_expr_target,
        "mean_expr_ref": mean_expr_ref,
    }


def _presort_contrast_arrays(arrays: dict) -> dict:
    """Pre-sort all arrays by descending scores with stable tie-breakers.
    
    Sort order: descending scores, then ascending pvals_adj, then ascending gene IDs.
    """
    gene_id = arrays["gene_id"]
    scores = arrays["scores"]
    pvals_adj = arrays["pvals_adj"]
    
    # Create deterministic sort index (descending scores, stable for ties)
    sort_idx = np.lexsort((
        gene_id,  # Gene IDs ascending (stable tie-breaker)
        np.where(np.isnan(pvals_adj), np.inf, pvals_adj),  # pvals_adj ascending
        -scores  # Scores descending (note the minus sign)
    ))
    
    # Apply sort to all arrays
    sorted_arrays = {}
    for key, arr in arrays.items():
        sorted_arrays[key] = arr[sort_idx]
    
    return sorted_arrays


def _write_array(group: zarr.Group, field_name: str, arr: np.ndarray) -> None:
    if arr.dtype.kind in {"O", "U"}:
        string_array = group.create_array(
            field_name,
            shape=arr.shape,
            dtype=VariableLengthUTF8(),
            fill_value="",
        )
        string_array[:] = np.asarray(arr, dtype=object)
        return

    group[field_name] = arr


def _compute_axis_bounds(
    effect_size: np.ndarray,
    pvals_adj: np.ndarray,
) -> tuple[Optional[float], Optional[float]]:
    """Compute effect_size_max and significance_max for metadata anchoring."""
    # Filter out NaN/Inf values
    valid_effect_size = effect_size[~(np.isnan(effect_size) | np.isinf(effect_size))]
    valid_pvals = pvals_adj[~np.isnan(pvals_adj)]
    
    effect_size_max = None
    significance_max = None
    
    if len(valid_effect_size) > 0:
        effect_size_max = float(np.ceil(np.max(np.abs(valid_effect_size))))
    
    if len(valid_pvals) > 0:
        # Compute -log10(p) and cap at 300
        logp_vals = -np.log10(np.maximum(valid_pvals, 1e-300))
        significance_max = min(300.0, float(np.ceil(np.max(logp_vals))))
    
    return effect_size_max, significance_max


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
    
    # Write string arrays with explicit V3 UTF-8 metadata for cross-library compatibility.
    for field_name, arr in arrays.items():
        _write_array(store, field_name, arr)
    
    # Compute metadata bounds for registry
    effect_size_max, significance_max = _compute_axis_bounds(
        arrays["effect_size"],
        arrays["pvals_adj"],
    )
    
    # Persist metric labels only when the corresponding metric is available.
    if contrast_config.method == "logreg":
        significance_label = None
        significance_max = None
    else:
        significance_label = (
            "-log10(FDR)" if contrast_config.corr_method else "-log10(p)"
        ) if significance_max is not None else None

    effect_size_label = "log2(Fold Change)" if effect_size_max is not None else None

    has_effect_size = effect_size_max is not None
    has_significance = significance_max is not None
    has_pct_expr_target = bool(np.any(~np.isnan(arrays["pct_expr_target"])))
    has_pct_expr_ref = bool(np.any(~np.isnan(arrays["pct_expr_ref"])))
    has_mean_expr_target = bool(np.any(~np.isnan(arrays["mean_expr_target"])))
    has_mean_expr_ref = bool(np.any(~np.isnan(arrays["mean_expr_ref"])))
    correction_method = (
        contrast_config.corr_method if significance_max is not None else None
    )
    if significance_max is not None and correction_method is None:
        correction_method = "uncorrected"
    
    # Extract top 10 feature IDs from sorted gene_id values.
    gene_id_list = arrays["gene_id"].tolist()
    top_10_gene_ids = gene_id_list[:10]
    
    return {
        "contrast_id": contrast_id,
        "group_1": contrast_config.group_1,
        "group_2": contrast_config.group_2,
        "test_type": contrast_config.test_type,
        "subset_column": contrast_config.subset_column,
        "subset_value": contrast_config.subset_value,
        "contrast_column": GROUPBY_COLUMN,
        "de_method": contrast_config.method,
        "correction_method": correction_method,
        "has_effect_size": has_effect_size,
        "has_significance": has_significance,
        "has_pct_expr_target": has_pct_expr_target,
        "has_pct_expr_ref": has_pct_expr_ref,
        "has_mean_expr_target": has_mean_expr_target,
        "has_mean_expr_ref": has_mean_expr_ref,
        "feature_type": "Gene Expression",
        "effect_size_label": effect_size_label,
        "significance_label": significance_label,
        "effect_size_max": effect_size_max,
        "significance_max": significance_max,
        "top_10_gene_ids": top_10_gene_ids,
    }


# In[7]:


def main() -> None:
    if not INPUT_PATH.exists():
        raise FileNotFoundError(f"Input file not found: {INPUT_PATH}")

    adata = ad.read_h5ad(INPUT_PATH)

    if not adata.var_names.is_unique:
        duplicate_count = int(adata.var_names.duplicated().sum())
        print(
            "Detected duplicate var_names "
            f"({duplicate_count}). Applying var_names_make_unique() before DE."
        )
        adata.var_names_make_unique()

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
    adata.obs["region"] = [_map_region(idx) for idx in adata.obs_names]

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

    hippocampus_subset = adata[adata.obs["region"] == "Hippocampus"].copy()
    hippocampus_subset.obs[GROUPBY_COLUMN] = hippocampus_subset.obs[GROUPBY_COLUMN].astype("category")
    hippocampus_cell_types = sorted(
        hippocampus_subset.obs[GROUPBY_COLUMN].unique().astype(str).tolist()
    )
    print(
        f"Found {len(hippocampus_cell_types)} hippocampus subset cell types: "
        f"{hippocampus_cell_types}"
    )
    
    # Build contrast matrix
    contrasts = _build_contrast_matrix(cell_types, hippocampus_cell_types)
    print(f"Built contrast matrix with {len(contrasts)} contrasts")
    
    # Compute global mean expression per gene (for all cells)
    if hasattr(adata.X, "mean"):
        mean_expr_global = np.asarray(adata.X.mean(axis=0)).ravel()
    else:
        mean_expr_global = np.asarray(adata.X).mean(axis=0).ravel()

    if hasattr(hippocampus_subset.X, "mean"):
        mean_expr_hippocampus = np.asarray(hippocampus_subset.X.mean(axis=0)).ravel()
    else:
        mean_expr_hippocampus = np.asarray(hippocampus_subset.X).mean(axis=0).ravel()
    
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
    logreg_key = "_contrast_logreg_multiclass"
    logreg_is_ready = False
    
    for contrast_config in contrasts:
        try:
            is_hippocampus_subset = (
                contrast_config.subset_column == "region"
                and contrast_config.subset_value == "Hippocampus"
            )
            active_adata = hippocampus_subset if is_hippocampus_subset else adata
            active_mean_expr = mean_expr_hippocampus if is_hippocampus_subset else mean_expr_global

            # Validate group size
            group_size = (active_adata.obs[GROUPBY_COLUMN] == contrast_config.group_1).sum()
            if group_size < MIN_CELLS:
                print(f"  ⊘ {contrast_config.contrast_id}: group '{contrast_config.group_1}' has {group_size} cells < {MIN_CELLS}")
                continue
            
            # Run rank_genes_groups
            if contrast_config.method == "logreg":
                # Fit the multiclass model exactly once, then reuse per-group coefficients.
                if not logreg_is_ready:
                    sc.tl.rank_genes_groups(
                        adata,
                        groupby=GROUPBY_COLUMN,
                        method="logreg",
                        pts=True,
                        key_added=logreg_key,
                        use_raw=False,
                    )
                    logreg_is_ready = True
                key_added = logreg_key
            else:
                key_added = f"_contrast_{contrast_config.contrast_id}"
                sc.tl.rank_genes_groups(
                    active_adata,
                    groupby=GROUPBY_COLUMN,
                    groups=[contrast_config.group_1],
                    method=contrast_config.method,
                    corr_method=contrast_config.corr_method,
                    tie_correct=True,
                    pts=True,
                    key_added=key_added,
                    use_raw=False,
                )
            
            # Extract per-contrast arrays
            arrays = _extract_contrast_arrays(
                active_adata,
                key_added,
                contrast_config.group_1,
                active_mean_expr,
                contrast_config.group_2,
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
            
            # Clean up per-contrast temporary keys. Keep multiclass logreg key until all logreg entries are processed.
            if contrast_config.method != "logreg" and key_added in active_adata.uns:
                del active_adata.uns[key_added]
            
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
        },
        "contrast_summary": {
            "total_contrasts": processed_count,
            "series": [
                {
                    "series_id": 1,
                    "method": "wilcoxon",
                    "corr": "benjamini-hochberg",
                    "n_contrasts": len(cell_types),
                },
                {
                    "series_id": 2,
                    "method": "t-test",
                    "corr": "bonferroni",
                    "n_contrasts": len(cell_types),
                },
                {
                    "series_id": 3,
                    "method": "logreg",
                    "corr": None,
                    "n_contrasts": len(cell_types),
                },
                {
                    "series_id": 4,
                    "method": "wilcoxon",
                    "corr": "benjamini-hochberg",
                    "subset_column": "region",
                    "subset_value": "Hippocampus",
                    "n_contrasts": len(hippocampus_cell_types),
                },
            ],
        },
        "contrasts": registry,
    }

    if logreg_key in adata.uns:
        del adata.uns[logreg_key]
    
    registry_path = OUTPUT_PATH / "uns" / "de" / "contrast_registry.json"
    with open(str(registry_path), "w") as f:
        json.dump(registry_json, f, indent=2)
    
    print(f"\nWrote {processed_count} DE contrasts to {OUTPUT_PATH}/uns/de/")
    print(f"Registry: {registry_path}")


# In[8]:


if __name__ == "__main__":
    main()




