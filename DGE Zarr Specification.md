# **Zarr Structure for Differential Expression**

## **1\. The Zarr Directory Tree**

```
dataset.zarr/
└── uns/
 └── de/
 ├── contrast_registry.json \# The manifest for UI discovery
 ├── de_001/ \# Atomic folder for a specific contrast
 │ ├── scores/ \# Test scores and primary sort key (descending)
 │ ├── gene_id/ \# Feature IDs from var index (e.g., gene symbols, Ensembl IDs, numeric IDs)
 │ ├── effect_size/ \# Effect sizes (NaN for logreg)
 │ ├── pvals_adj/ \# FDR values (NaN for logreg)
 │ ├── pvals/ \# Raw p-values (NaN for logreg)
 │ ├── significance/ \# Precomputed -log10(adjusted p-value), data-driven cap (NaN for logreg)
 │ ├── pct_expr_target/ \# % in group_1 (NaN for logreg)
 │ ├── pct_expr_ref/ \# % in reference/rest (NaN for logreg)
 │ ├── mean_expr/ \# Global average expression
 │ ├── mean_expr_target/ \# Mean expression in group_1
 │ └── mean_expr_ref/ \# Mean expression in reference/rest
 └── de_002/ \# Parallel folder for the next test
```

Note: all contrast folders currently contain the same array keys across methods. `gene_id` values are sourced from AnnData var index labels (`var_names`). For multiclass `logreg`, arrays that Scanpy does not provide are emitted as `NaN` placeholders to keep a stable on-disk shape.

### The significance transform

`significance` is the **plot-ready** volcano y-value, precomputed at data-generation time so the frontend is a dumb consumer that plots the column directly (no client-side `-log10`).

- Value: `-log10(pvals_adj)` per gene.
- Zero / underflowed p-values (which give infinite `-log10`) are pinned to a **data-driven cap** `= ceil(max_finite * 1.5)`, i.e. the largest real `-log10` in the contrast scaled up by a factor (default `1.5`, configurable via `SIGNIFICANCE_GAP_FACTOR`). This leaves a clear, readable gap above the real values instead of the old arbitrary flat cap. If no zero-p exists the cap is just `ceil(max_finite)`.
- `NaN` where `pvals_adj` is `NaN` (e.g. `logreg`).
- The cap equals `nanmax(significance)` and is exposed in the registry as `significance_max`, so array and bound are always consistent and nothing exceeds the bound.

## **2\. The contrast_registry.json Schema**

This registry allows the frontend to descirbe the available contrasts, determine availablity of volcano plot and set chart labels.

> **Note — not preserved by an AnnData round-trip.** `contrast_registry.json` is a plain JSON file placed inside the store directory, **not** a Zarr node. This is intentional: it is a manifest developed for our apps, fetched directly by path over HTTP, and is not meant to be consumed by a generic Zarr/AnnData reader. As a consequence it is skipped during metadata consolidation, and a full AnnData round-trip (`ad.read_zarr(store)` → `adata.write_zarr(...)`) **drops the registry file**, because AnnData only serializes nodes it recognizes. Regenerate it by re-running the notebook if a store has been rewritten through AnnData.

| Column                  | Example Value        | Description                                                                                                       |
| :---------------------- | :------------------- | :---------------------------------------------------------------------------------------------------------------- |
| **contrast_id**         | "de_001"             | Matches the Zarr folder name.                                                                                     |
| **test_type**           | "multiclass"         | Frontend grouping label for contrast semantics (for example one_vs_rest, pairwise, multiclass, or custom labels). |
| **de_method**           | "wilcoxon"           | Free-form DE algorithm identifier (for example `wilcoxon`, `logreg`, `deseq2`, `mast`, `edger`).                  |
| **has_effect_size**     | true                 | Whether effect-size metrics are available for this contrast.                                                      |
| **has_significance**    | true                 | Whether significance metrics are available for this contrast.                                                     |
| **has_pct_expr_target** | true                 | Whether `% expressing` values for the target group are available.                                                 |
| **has_pct_expr_ref**    | true                 | Whether `% expressing` values for the reference/rest group are available.                                         |
| **has_mean_expr_target** | true                | Whether mean expression values for the target group are available.                                                |
| **has_mean_expr_ref**   | true                 | Whether mean expression values for the reference/rest group are available.                                        |
| **correction_method**   | "benjamini-hochberg" | P-value adjustment method when significance metrics are available; otherwise `null`.                              |
| **feature_type**        | "Gene Expression"    | Type of features for display/filtering. Currently not used.                                                       |
| **effect_size_label**   | "log2(Fold Change)"  | Eeffect-size label in the frontend.                                                                               |
| **significance_label**  | "-log10(FDR)"        | Significance label in the frontend when available; otherwise `null`.                                              |
| **effect_size_max**     | 8.5                  | Max absolute effect-size bound when available; otherwise `null`.                                                  |
| **significance_max**    | 30.0                 | Data-driven significance cap = `nanmax(significance)` (the volcano y-max) when available; otherwise `null`.        |
| **subset_column**       | "region"             | Metadata category defining test boundary (or `null`).                                                             |
| **subset_value**        | "Hippocampus"        | Subset value within `subset_column` (or `null`).                                                                  |
| **contrast_column**     | "CellType"           | Metadata key being tested.                                                                                        |
| **group_1**             | "Sst"                | Target group.                                                                                                     |
| **group_2**             | null                 | Reference group (`null` for one-vs-rest and current multiclass output).                                           |
| **top_10_gene_ids**     | ["A", "B", ...]      | First 10 feature IDs after score-based pre-sorting.                                                               |

## **3\. The Core Architectural Rules**

1. **The "Pre-sort" Rule:** Every array inside de_XXX/ **must** be sorted descending by the scores array before writing. This allows Zarrita.js to use slice(0, 50\) to get the top markers instantly.
2. **The "No Magic Strings" Rule:** For One-vs-Rest, group_2 must be null. The frontend handles the "Rest" label translation.
3. **The "No Overlap" Rule:** Comparisons must happen within a single column. For complex multi-column logic, use a composite column in .obs first.
4. **The "Stable Arrays" Rule:** Keep the same per-contrast array keys for all methods. If a method does not provide a metric (for current Scanpy `logreg`: `effect_size`, `pvals`, `pvals_adj`, `significance`, `pts`, `pts_rest`), emit `NaN` values in the corresponding arrays.
