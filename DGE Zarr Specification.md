# **Production Zarr Structure for Differential Expression (v2.0)**

This structure is optimized for React/Zarrita.js frontends, ensuring sub-second rendering by pre-sorting data and providing all necessary metadata for dynamic visualization.

## **1\. The Zarr Directory Tree**

dataset.zarr/  
└── uns/  
 └── de/  
 ├── contrast_registry.json \# The manifest for UI discovery  
 ├── de_001/ \# Atomic folder for a specific contrast  
 │ ├── symbols/ \# \[Pre-sorted\] Gene symbols (e.g., "CD3D")  
 │ ├── logfoldchanges/ \# \[Pre-sorted\] Effect sizes (x-axis data)  
 │ ├── pvals_adj/ \# \[Pre-sorted\] FDR values (y-axis data)  
 │ ├── pvals/ \# \[Pre-sorted\] Raw p-values  
 │ ├── scores/ \# \[Pre-sorted\] Primary sort key (Descending)  
 │ ├── pct_expr_target/ \# \[Pre-sorted\] % in group_1  
 │ ├── pct_expr_ref/ \# \[Pre-sorted\] % in group_2  
 │ └── mean_expr/ \# \[Pre-sorted\] Global average expression  
 └── de_002/ \# Parallel folder for the next test

## **2\. The contrast_registry.json Schema**

This registry allows the frontend to populate dropdowns, set chart labels, and anchor axis scales without downloading any heavy gene arrays.

| Column                | Example Value           | Description                                               |
| :-------------------- | :---------------------- | :-------------------------------------------------------- |
| **comparison_id**     | "de_001"                | Matches the Zarr folder name.                             |
| **test_type**         | "pairwise"              | UI category (e.g., one_vs_rest, pairwise).                |
| **de_method**         | "wilcoxon"              | Algorithm (e.g., wilcoxon, deseq2, mast).                 |
| **correction_method** | "benjamini-hochberg"    | The p-value adjustment used.                              |
| **model_formula**     | "\~ treatment \+ batch" | The statistical formula used (or null).                   |
| **feature_type**      | "Gene Expression"       | Type of features (Genes, ADTs, Peaks).                    |
| **x_axis_label**      | "log2(Fold Change)"     | Literal string for chart X-axis.                          |
| **y_axis_label**      | "-log10(FDR)"           | Literal string for chart Y-axis.                          |
| **lfc_max**           | 8.5                     | Global max absolute LFC for symmetric X-axis.             |
| **logp_max**          | 300.0                   | Global max \-log10(p) for Y-axis ceiling. Capped at 300\. |
| **subset_column**     | "tissue"                | Metadata category defining the test boundary (or null).   |
| **subset_value**      | "Lung"                  | The specific group the test is restricted to (or null).   |
| **contrast_column**   | "treatment"             | The metadata key being tested.                            |
| **group_1**           | "Drug_A"                | The Target group.                                         |
| **group_2**           | "Vehicle"               | The Reference group (or null for Rest).                   |

## **3\. The Core Architectural Rules**

1. **The "Pre-sort" Rule:** Every array inside de_XXX/ **must** be sorted descending by the scores array before writing. This allows Zarrita.js to use slice(0, 50\) to get the top markers instantly.
2. **The "No Magic Strings" Rule:** For One-vs-Rest, group_2 must be null. The frontend handles the "Rest" label translation.
3. **The "No Overlap" Rule:** Comparisons must happen within a single column. For complex multi-column logic, use a composite column in .obs first.
4. **The "Honest Labeling" Rule:** Always populate x_axis_label and y_axis_label based on the specific de_method (e.g., use "Beta" for MAST, "log2FC" for Scanpy).
5. **The "Layout Anchor" Rule:** Use lfc_max to set a symmetric X-axis (from \-lfc_max to \+lfc_max) and logp_max to set the Y-axis height. This prevents "jumping" axes when flipping between clusters.
