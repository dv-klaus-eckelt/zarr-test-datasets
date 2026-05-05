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
 │ ├── pct_expr_target/ \# % in group_1 (NaN for logreg)
 │ ├── pct_expr_ref/ \# % in reference/rest (NaN for logreg)
 │ └── mean_expr/ \# Global average expression
 └── de_002/ \# Parallel folder for the next test
```

Note: all contrast folders currently contain the same array keys across methods. `gene_id` values are sourced from AnnData var index labels (`var_names`). For multiclass `logreg`, arrays that Scanpy does not provide are emitted as `NaN` placeholders to keep a stable on-disk shape.

## **2\. The contrast_registry.json Schema**

This registry allows the frontend to descirbe the available contrasts, determine availablity of volcano plot and set chart labels.

| Column                  | Example Value        | Description                                                                                                       |
| :---------------------- | :------------------- | :---------------------------------------------------------------------------------------------------------------- |
| **contrast_id**         | "de_001"             | Matches the Zarr folder name.                                                                                     |
| **test_type**           | "multiclass"         | Frontend grouping label for contrast semantics (for example one_vs_rest, pairwise, multiclass, or custom labels). |
| **de_method**           | "wilcoxon"           | Free-form DE algorithm identifier (for example `wilcoxon`, `logreg`, `deseq2`, `mast`, `edger`).                  |
| **has_effect_size**     | true                 | Whether effect-size metrics are available for this contrast.                                                      |
| **has_significance**    | true                 | Whether significance metrics are available for this contrast.                                                     |
| **has_pct_expr_target** | true                 | Whether `% expressing` values for the target group are available.                                                 |
| **has_pct_expr_ref**    | true                 | Whether `% expressing` values for the reference/rest group are available.                                         |
| **correction_method**   | "benjamini-hochberg" | P-value adjustment method when significance metrics are available; otherwise `null`.                              |
| **feature_type**        | "Gene Expression"    | Type of features for display/filtering. Currently not used.                                                       |
| **effect_size_label**   | "log2(Fold Change)"  | Eeffect-size label in the frontend.                                                                               |
| **significance_label**  | "-log10(FDR)"        | Significance label in the frontend when available; otherwise `null`.                                              |
| **effect_size_max**     | 8.5                  | Max absolute effect-size bound when available; otherwise `null`.                                                  |
| **significance_max**    | 300.0                | Max significance bound when available; otherwise `null`.                                                          |
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
4. **The "Stable Arrays" Rule:** Keep the same per-contrast array keys for all methods. If a method does not provide a metric (for current Scanpy `logreg`: `effect_size`, `pvals`, `pvals_adj`, `pts`, `pts_rest`), emit `NaN` values in the corresponding arrays.
