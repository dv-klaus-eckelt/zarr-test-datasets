# **Production Zarr Structure for Differential Expression (v2.1)**

This structure is optimized for React/Zarrita.js frontends, ensuring sub-second rendering by pre-sorting data and providing all necessary metadata for dynamic visualization.

## **1\. The Zarr Directory Tree**

dataset.zarr/  
└── uns/  
 └── de/  
 ├── contrast_registry.json \# The manifest for UI discovery  
 ├── de_001/ \# Atomic folder for a specific contrast  
 │ ├── symbols/ \# \[Pre-sorted\] Gene symbols (e.g., "CD3D")  
 │ ├── logfoldchanges/ \# \[Pre-sorted\] Effect sizes (NaN for logreg)  
 │ ├── pvals_adj/ \# \[Pre-sorted\] FDR values (NaN for logreg)  
 │ ├── pvals/ \# \[Pre-sorted\] Raw p-values (NaN for logreg)  
 │ ├── scores/ \# \[Pre-sorted\] Primary sort key (descending)  
 │ ├── pct_expr_target/ \# \[Pre-sorted\] % in group_1 (NaN for logreg)  
 │ ├── pct_expr_ref/ \# \[Pre-sorted\] % in reference/rest (NaN for logreg)  
 │ └── mean_expr/ \# \[Pre-sorted\] Global average expression  
 └── de_002/ \# Parallel folder for the next test

Note: all contrast folders currently contain the same array keys across methods. For multiclass `logreg`, arrays that Scanpy does not provide are emitted as `NaN` placeholders to keep a stable on-disk shape.

## **2\. The contrast_registry.json Schema**

This registry allows the frontend to populate dropdowns, set chart labels, and anchor axis scales without downloading any heavy gene arrays.

| Column                | Example Value          | Description                                                                |
| :-------------------- | :--------------------- | :------------------------------------------------------------------------- |
| **contrast_id**       | "de_001"               | Matches the Zarr folder name.                                              |
| **test_type**         | "multiclass"           | UI category (one_vs_rest, pairwise, multiclass).                           |
| **de_method**         | "wilcoxon"             | Algorithm (`wilcoxon`, `t-test`, `logreg`).                                |
| **correction_method** | "benjamini-hochberg"   | P-value adjustment for non-logreg methods, else `null`.                    |
| **feature_type**      | "Gene Expression"      | Type of features for display/filtering.                                    |
| **x_axis_label**      | "log2(Fold Change)"    | Current generator output for all methods.                                  |
| **y_axis_label**      | "-log10(FDR)"          | Non-logreg significance axis label; `null` for `logreg`.                   |
| **available_plots**   | ["volcano", "dotplot"] | `wilcoxon`/`t-test`: volcano + dotplot; `logreg`: bar_chart + dotplot.     |
| **lfc_max**           | 8.5                    | Max absolute LFC bound for non-logreg volcano x-axis; `null` for `logreg`. |
| **logp_max**          | 300.0                  | Max -log10(p) bound (capped at 300) for non-logreg; `null` for `logreg`.   |
| **subset_column**     | "region"               | Metadata category defining test boundary (or `null`).                      |
| **subset_value**      | "Hippocampus"          | Subset value within `subset_column` (or `null`).                           |
| **contrast_column**   | "CellType"             | Metadata key being tested.                                                 |
| **group_1**           | "Sst"                  | Target group.                                                              |
| **group_2**           | null                   | Reference group (`null` for one-vs-rest and current multiclass output).    |
| **n_features**        | 15422                  | Number of ranked features in the corresponding `de_xxx` arrays.            |
| **top_10_genes**      | ["A", "B", ...]        | First 10 symbols after score-based pre-sorting.                            |

## **3\. The Core Architectural Rules**

1. **The "Pre-sort" Rule:** Every array inside de_XXX/ **must** be sorted descending by the scores array before writing. This allows Zarrita.js to use slice(0, 50\) to get the top markers instantly.
2. **The "No Magic Strings" Rule:** For One-vs-Rest, group_2 must be null. The frontend handles the "Rest" label translation.
3. **The "No Overlap" Rule:** Comparisons must happen within a single column. For complex multi-column logic, use a composite column in .obs first.
4. **The "Stable Arrays" Rule:** Keep the same per-contrast array keys for all methods. If a method does not provide a metric (for current Scanpy `logreg`: `logfoldchanges`, `pvals`, `pvals_adj`, `pts`, `pts_rest`), emit `NaN` values in the corresponding arrays.
5. **The "Layout Anchor" Rule:** For non-logreg contrasts, use `lfc_max` to set a symmetric X-axis (from -lfc_max to +lfc_max) and `logp_max` to set the Y-axis height. For `logreg`, both are `null`.
6. **The "Multiclass LogReg" Rule:** A single multiclass logistic regression model may be fit once and emitted as one contrast folder per group; each folder uses group-specific coefficients in `scores`, with `test_type` set to "multiclass".
