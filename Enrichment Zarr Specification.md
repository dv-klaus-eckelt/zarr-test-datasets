# **Zarr Structure for Gene Set Enrichment**

This is the enrichment analogue of `DGE Zarr Specification.md`. The store is
self-contained: a verbatim copy of the differential-expression (DE) store (`X`,
`obs`, `var`, `uns/de`) plus a new `uns/enrichment` group and per-cell activity
matrices in `obsm/`. Each enrichment is computed **from** a DE contrast and links
back to it, so the same frontend loaders used for DE (`loadColumnValues` over a
per-folder layout, plus a registry parse) apply here with the feature axis being
gene sets instead of genes.

## **1\. The Zarr Directory Tree**

```
dataset.zarr/
├── uns/
│   ├── de/                              # reused verbatim from the DE store
│   │   ├── contrast_registry.json
│   │   └── de_XXX/ …
│   └── enrichment/
│       ├── enrichment_registry.json     # manifest for UI discovery
│       ├── es_001/                      # atomic folder for one (contrast, method, collection)
│       │   ├── scores/                  # primary method statistic and sort key (descending)
│       │   ├── gene_set_id/             # gene set ids (e.g. HALLMARK_HYPOXIA)
│       │   ├── effect_size/             # signed effect (NES / coefficient); NaN for unsigned methods
│       │   ├── pvals/                   # raw p-values (NaN; decoupler returns adjusted only)
│       │   ├── pvals_adj/               # BH-adjusted p-values (NaN if the method gives none)
│       │   ├── significance/            # precomputed -log10(adjusted p-value), data-driven cap
│       │   │                            #   (NaN when no p-values)
│       │   ├── set_size/                # genes in the set (from the collection)
│       │   └── overlap_size/            # set genes present/tested in the data
│       └── es_002/ …                    # parallel folder for the next enrichment
└── obsm/
    ├── X_activity_mlm_progeny/          # per-cell activity matrix (cells × gene sets)
    ├── X_activity_ulm_collectri/
    └── X_activity_aucell_hallmark/
```

Note: all `es_XXX` folders contain the same array keys across methods. Values a
method does not produce are emitted as `NaN` placeholders to keep a stable
on-disk shape (e.g. `effect_size` is NaN for over-representation, `pvals` is NaN
because decoupler returns only adjusted p-values). `gene_set_id` values are the
gene set source names from the collection. Per-cell activity matrices store the
score matrix only; the gene set id for each column is recorded (in column order)
in the registry's `cell_activities[].gene_set_ids`.

## **2\. The per-`es_XXX` Arrays**

Every array in an `es_XXX/` folder has the same length (`n_gene_sets`) and is
pre-sorted descending by `scores`.

| Array              | Type   | Description                                                                       |
| :----------------- | :----- | :-------------------------------------------------------------------------------- |
| **gene_set_id**    | string | Gene set source name (for example `HALLMARK_HYPOXIA`). Mirrors DE `gene_id`.      |
| **scores**         | float  | Primary method statistic and descending sort key (NES, t-value, ln odds ratio).  |
| **effect_size**    | float  | Signed effect (NES / coefficient / ln odds). `NaN` for unsigned methods (ORA).   |
| **pvals**          | float  | Raw p-values. `NaN` — decoupler returns adjusted p-values only.                   |
| **pvals_adj**      | float  | Benjamini-Hochberg adjusted p-values. `NaN` if the method produces none.          |
| **significance**   | float  | Plot-ready `-log10(pvals_adj)` with a data-driven cap (see below). `NaN` if no p.  |
| **set_size**       | float  | Number of genes in the set, from the collection.                                  |
| **overlap_size**   | float  | Number of set genes present/tested in the data.                                   |

### The significance transform

`significance` is the **plot-ready** significance value, precomputed at
data-generation time (identical policy to the DE spec) so the frontend is a dumb
consumer that plots the column directly instead of computing `-log10`:

- Value: `-log10(pvals_adj)` per gene set.
- Zero / underflowed p-values (infinite `-log10`) are pinned to a **data-driven
  cap** `= ceil(max_finite * 1.5)` (factor configurable via
  `SIGNIFICANCE_GAP_FACTOR`, default `1.5`), leaving a readable gap above the
  largest real value; if no zero-p exists the cap is `ceil(max_finite)`.
- `NaN` where `pvals_adj` is `NaN` (a method that produces no p-values).
- The cap equals `nanmax(significance)` and is exposed as the registry's
  `significance_max`, so the array and the bound stay consistent.

## **3\. The enrichment_registry.json Schema**

The registry lets the frontend describe available enrichments, determine which
charts/labels apply, and link each enrichment back to its DE contrast — without
loading any per-gene-set arrays. It has two arrays: `contrast_enrichments[]` (one
per `es_XXX`) and `cell_activities[]` (one per `obsm/` activity matrix), plus
top-level `version`, `generated_at`, `decoupler_version`, `zarr_version`,
`dataset_metadata`, and `enrichment_summary`.

> **Note — not preserved by an AnnData round-trip.** `enrichment_registry.json`
> (like the DE `contrast_registry.json`) is a plain JSON file placed inside the
> store directory, **not** a Zarr node. This is intentional: it is a manifest
> developed for our apps, fetched directly by path over HTTP, and is not meant to
> be consumed by a generic Zarr/AnnData reader. As a consequence it is skipped
> during metadata consolidation (a harmless `ZarrUserWarning`), and a full
> AnnData round-trip (`ad.read_zarr(store)` → `adata.write_zarr(...)`) **drops
> the registry files**, because AnnData only serializes nodes it recognizes. The
> `uns/enrichment/es_XXX` arrays and `obsm/` matrices survive such a round-trip;
> the `*.json` manifests do not. Regenerate them by re-running the notebook if a
> store has been rewritten through AnnData.

### `contrast_enrichments[]`

| Column                  | Example Value                 | Description                                                                     |
| :---------------------- | :---------------------------- | :------------------------------------------------------------------------------ |
| **enrichment_id**       | "es_001"                      | Matches the Zarr folder name.                                                    |
| **source_contrast_id**  | "de_001"                      | DE contrast this enrichment was computed from; links to `uns/de`.               |
| **enrichment_method**   | "gsea"                        | Free-form enrichment algorithm identifier (`gsea`, `ora`, `mlm`, `ulm`, …).     |
| **gene_set_collection** | "MSigDB Hallmark"             | Collection label (`MSigDB Hallmark`, `GO:BP`, `OmniPath`, `PROGENy`, `CollecTRI`). |
| **n_gene_sets**         | 40                            | Number of gene sets (array length) in the folder.                               |
| **input_statistic**     | "de_scores"                   | Per-gene input to the method (`de_scores`, `top_50_up_genes`).                  |
| **feature_type**        | "Gene Set"                    | Always `Gene Set` for enrichment.                                               |
| **has_effect_size**     | true                          | Whether a signed effect-size metric is available.                               |
| **has_significance**    | true                          | Whether significance metrics are available.                                     |
| **has_set_size**        | true                          | Whether gene set sizes are available.                                           |
| **has_overlap**         | true                          | Whether overlap sizes are available.                                            |
| **has_leading_edge**    | false                         | Whether leading-edge / driver gene ids are available.                           |
| **correction_method**   | "benjamini-hochberg"          | P-value adjustment method when significance is available; otherwise `null`.     |
| **effect_size_label**   | "Normalized Enrichment Score" | Effect-size axis label when available; otherwise `null`.                        |
| **significance_label**  | "-log10(FDR)"                 | Significance axis label when available; otherwise `null`.                       |
| **effect_size_max**     | 3.0                           | Max absolute effect-size bound when available; otherwise `null`.                |
| **significance_max**    | 12.0                          | Data-driven significance cap = `nanmax(significance)` (plot y-max) when available; otherwise `null`. |
| **group_1**             | "ASC1"                        | Target group inherited from the source contrast.                                |
| **group_2**             | null                          | Reference group inherited from the source contrast (`null` for one-vs-rest).    |
| **test_type**           | "one_vs_rest"                 | Comparison semantics inherited from the source contrast.                        |
| **subset_column**       | null                          | Metadata category restricting the source contrast (or `null`).                  |
| **subset_value**        | null                          | Subset value within `subset_column` (or `null`).                                |
| **contrast_column**     | "CellType"                    | Metadata key tested in the source contrast.                                     |
| **top_10_gene_set_ids** | ["HALLMARK_HYPOXIA", …]       | First 10 gene set ids after score-based pre-sorting.                            |

### `cell_activities[]`

| Column                  | Example Value              | Description                                                          |
| :---------------------- | :------------------------- | :------------------------------------------------------------------ |
| **obsm_key**            | "X_activity_mlm_progeny"   | obsm key of the cells × gene sets activity matrix.                  |
| **method**              | "mlm"                      | Free-form enrichment algorithm identifier.                          |
| **gene_set_collection** | "PROGENy"                  | Collection used to compute the activities.                          |
| **value_label**         | "MLM t-value"              | Human-readable label for the activity values.                       |
| **gene_set_ids**        | ["Androgen", "EGFR", …]    | Gene set ids in column order of the obsm matrix.                    |
| **n_gene_sets**         | 14                         | Number of columns; equals `gene_set_ids` length.                   |

## **4\. The Core Architectural Rules**

1. **The "Pre-sort" Rule:** Every array inside `es_XXX/` **must** be sorted
   descending by the `scores` array before writing (ties broken by ascending
   `pvals_adj` then `gene_set_id`; `NaN` scores sort last). This lets the frontend
   use `slice(0, N)` to get the top gene sets instantly.
2. **The "No Magic Strings" Rule:** For one-vs-rest source contrasts, `group_2`
   is `null`. The frontend handles the "Rest" label translation. Availability is
   read from `has_*` flags, never inferred from `enrichment_method`.
3. **The "Stable Arrays" Rule:** Keep the same per-folder array keys for all
   methods. If a method does not provide a metric, emit `NaN` values in the
   corresponding array (`effect_size` for unsigned methods such as ORA, `pvals`
   for every method since decoupler returns only adjusted p-values, and
   `significance` when a method produces no p-values).
4. **The "Feature Axis" Rule:** The feature axis is **gene sets**, not genes.
   `gene_set_id` is the per-folder feature id (mirroring DE's `gene_id`), and
   `feature_type` is fixed to `"Gene Set"`. Each enrichment links back to the DE
   contrast it was derived from via `source_contrast_id`.
