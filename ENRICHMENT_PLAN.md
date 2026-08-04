# Gene Set Enrichment Test Dataset + Universal Schema — Implementation Plan

> Status: planned, not yet implemented. Pick up work from this file.

## Context

The repo already produces a differential-expression (DE) test dataset
(`calculate-differential-expression-data.ipynb` →
`test-data/habib17-differential-expression-test-data-format.zarr`) whose
`uns/de/` layout — one folder of pre-sorted arrays per contrast plus a
`contrast_registry.json` manifest — is a deliberately method-agnostic schema
consumed by the single_cell app (`DifferentialExpression.tsx`, loaded via
`loadColumnValues` + `contrast-registry.schema.ts`).

We want the analogous artifact for **gene set enrichment**: a second Python
script and a second, self-contained test dataset that adds enrichment results,
computed with the **decoupler** library across a diverse set of methods and gene
set collections, to define a universal enrichment schema the way the DE schema
was defined. Decisions confirmed with the user:

- **New store is self-contained**: base AnnData + `uns/de` (reused verbatim from
  the existing DE output) + new `uns/enrichment`.
- **Both representations**: contrast-based enrichment (mirrors the DE
  folder/registry pattern, arrays over gene sets) **and** per-cell activity
  (obsm matrices, cells × gene sets).
- **Real gene set collections** downloaded via decoupler's `dc.op` (OmniPath /
  MSigDB) — requires network access at run time.

Note: the single_cell app has **no** enrichment feature yet, so this dataset
*defines* the schema a future frontend will consume. Mirroring the DE structure
means that frontend can reuse the existing per-folder `loadColumnValues` loader
and a registry parse identical in shape to the DE path.

## Deliverables

1. `calculate-gene-set-enrichment-data.ipynb` — the new script.
2. `test-data/habib17-enrichment-test-data-format.zarr` — the new dataset.
3. `Enrichment Zarr Specification.md` — spec doc, parallel to `DGE Zarr Specification.md`.
4. `enrichment-registry.schema.ts` — Zod schema, parallel to `contrast-registry.schema.ts`.
5. `pyproject.toml` — add `decoupler` dependency; README row for the new dataset.

## Approach

### Base store: reuse existing DE output
The notebook builds the new store without re-running scanpy DE:
- `shutil.copytree` the existing
  `test-data/habib17-differential-expression-test-data-format.zarr` →
  `test-data/habib17-enrichment-test-data-format.zarr` (guarantees identical
  `uns/de`, `X`, `obs`, `var`). Fall back to a clear error if the DE store is
  absent (instruct user to run the DE notebook first).
- Read the AnnData back with `ad.read_zarr(new_path)` — gives the
  log-normalized `X` (for per-cell activity) and `obs`/`var` (human gene symbols,
  confirmed, e.g. `ISG15`, `SAMD11` — compatible with MSigDB Hallmark).

### Gene set collections (decoupler `dc.op`, decoupler 2.x)
Download the collections the user specified:
- `dc.op.hallmark(organism="human")` — MSigDB Hallmark (unweighted sets).
- **GO:BP** — MSigDB C5 GO:Biological Process via `dc.op.resource(...)` (filter to
  the GO:BP subcollection).
- **OmniPath** annotations via `dc.op.resource(...)` as gene sets.
- `dc.op.progeny(organism="human")` — PROGENy pathways (signed weights).
- `dc.op.collectri(organism="human")` — TF→target regulons (signed weights).

Each `net` is long-format `source, target, weight`. Record `n_gene_sets` and
`set_size` per source (`net.groupby("source").size()`).

### Contrast-based enrichment (mirrors DE) → `uns/enrichment/es_XXX/`
Input signature matrix: rows = selected DE contrasts, cols = genes, values = the
per-gene DE statistic read from `uns/de/de_XXX/` arrays (`scores`, reindexed to
full `var_names`). Use the one-vs-rest wilcoxon series (`de_001`…`de_015`, full
dataset, signed `effect_size` + `scores` present).

Run the **user-specified method × collection matrix**. decoupler 2.x `dc.mt.*`
methods take a `samples × features` DataFrame + `net` and return
`(score_df, padj_df)` shaped `contrasts × sources`; split per contrast into
folders (each (contrast, method, collection) triple = one `es_XXX` folder):

| method | collection | contrasts | input stat | signed? | notes |
|---|---|---|---|---|---|
| `dc.mt.gsea` | **Hallmark** | **all 15** one-vs-rest | ranked DE scores | yes | NES + pval (+ leading edge if available) |
| `dc.mt.ora` | **Hallmark** | **all 15** one-vs-rest | top-N up genes/contrast | no (odds) | over-representation; `effect_size` = ln odds ratio, or NaN |
| `dc.mt.gsea` | GO:BP | 3–5 contrasts | ranked DE scores | yes | large DB |
| `dc.mt.ora` | GO:BP | 3–5 contrasts | top-N up genes | no (odds) | |
| `dc.mt.ora` | OmniPath | 3–5 contrasts | top-N up genes | no (odds) | |
| `dc.mt.mlm` | PROGENy | 3–5 contrasts | DE scores | yes | multivariate, weighted pathways |
| `dc.mt.ulm` | CollecTRI | 3–5 contrasts | DE scores | yes | TF activity via regulons |

Total `es_XXX` ≈ 15 + 15 + 5×5 ≈ 45 folders. Use the same 3–5 contrast subset
across the non-Hallmark rows for comparability. Log any contrast/method skips
explicitly.

**Per-folder arrays** (one value per gene set, pre-sorted descending by `scores`,
same architectural rules as DE — pre-sort, stable arrays, NaN placeholders for
metrics a method doesn't produce):
- `gene_set_id` (source name, e.g. `HALLMARK_HYPOXIA`) — mirrors `gene_id`
- `scores` — primary method statistic + sort key
- `effect_size` — signed effect (NES / coefficient / ln odds); NaN for unsigned
- `pvals`, `pvals_adj` — raw + BH-adjusted (NaN if method gives none)
- `set_size` — genes in the set (from collection)
- `overlap_size` — set genes present/tested in the data
- `leading_edge_gene_ids` (optional) — comma-joined driver genes (GSEA leading
  edge / ORA overlap); omit if the method/version doesn't expose it

Reuse the DE writer helpers' conventions: `VariableLengthUTF8` for string arrays
(`_write_array`), lexsort pre-sort (`_presort_contrast_arrays`), and
`_compute_axis_bounds` for `effect_size_max` / `significance_max`.

### Per-cell activity → `obsm/` + registry section
Run a few `dc.mt.*` methods on the full AnnData (`X`), which populate
`adata.obsm['score_<method>']` (+ `padj_<method>`) as DataFrames (cells × gene
sets, column labels = gene set ids). Bounded combos aligned with the chosen
collections:
- `dc.mt.mlm` × PROGENy
- `dc.mt.ulm` × CollecTRI (TF activity)
- `dc.mt.aucell` × Hallmark (per-cell footprint method, for variety)

Persist by writing the updated AnnData's obsm to the store (distinct obsm keys so
DataFrame column labels — the gene set ids — survive). Record each in the
registry's `cell_activities` list with `obsm_key`, `method`,
`gene_set_collection`, `value_label`, and `gene_set_ids`.

### Registry: `uns/enrichment/enrichment_registry.json`
Mirror `contrast_registry.json`, with two arrays:

`contrast_enrichments[]` (one per `es_XXX`):
- `enrichment_id` ("es_001"), `source_contrast_id` ("de_001") — links back to DE
- inherited: `group_1`, `group_2`, `test_type`, `subset_column`, `subset_value`, `contrast_column`
- `enrichment_method` (free-form: "gsea"/"ora"/"mlm"/"ulm")
- `gene_set_collection` ("MSigDB Hallmark"/"GO:BP"/"OmniPath"/"PROGENy"/"CollecTRI"), `n_gene_sets`
- `input_statistic` ("de_scores"/"top_50_up_genes")
- `feature_type` = "Gene Set"
- `has_effect_size`, `has_significance`, `has_set_size`, `has_overlap`, `has_leading_edge`
- `correction_method`, `effect_size_label` (e.g. "Normalized Enrichment Score"),
  `significance_label` ("-log10(FDR)"), `effect_size_max`, `significance_max`
- `top_10_gene_set_ids`

`cell_activities[]` (one per obsm activity matrix): `obsm_key`, `method`,
`gene_set_collection`, `value_label`, `gene_set_ids`, `n_gene_sets`.

Plus top-level `version`, `generated_at`, `decoupler_version`, `zarr_version`,
`dataset_metadata`, and an `enrichment_summary` (series list) mirroring the DE
registry.

### Zod schema `enrichment-registry.schema.ts`
Parallel to `contrast-registry.schema.ts`: shared enrichment fields + the same
`has_effect_size`/`has_significance` discriminated union for label/bound nullability,
a `feature_type` fixed to gene-set semantics, plus a `CellActivitySchema`. Export
`EnrichmentSchema`, `EnrichmentsArraySchema`, `CellActivitiesArraySchema` and their
inferred TS types. Keep it method-agnostic (no discrimination on `enrichment_method`).

### Spec doc `Enrichment Zarr Specification.md`
Parallel to `DGE Zarr Specification.md`: directory tree
(`uns/enrichment/{enrichment_registry.json, es_XXX/…}` + `obsm/…`), the per-folder
array table, the registry table (both arrays), and the same four architectural
rules (Pre-sort, No Magic Strings, Stable Arrays, and a Feature-axis rule noting
the feature axis is gene sets, not genes).

### Dependency + README
- `pyproject.toml`: add `decoupler` to `dependencies` (decoupler 2.x resolves
  under the existing `exclude-newer = 2026-04-30` + `prerelease = "allow"`; it
  pulls `omnipath`). Then `uv sync`.
- README: add a row for `habib17-enrichment-test-data-format.zarr` in the
  Single Cell table (localhost + public links, following the DE row format).

## Critical files
- Create: `calculate-gene-set-enrichment-data.ipynb`, `Enrichment Zarr Specification.md`, `enrichment-registry.schema.ts`
- Modify: `pyproject.toml`, `README.md`
- Reuse (read/copy, do not modify): existing DE store `uns/de/*`, and the DE
  notebook's writer helper conventions (`_write_array`, `_presort_contrast_arrays`,
  `_compute_axis_bounds`, `VariableLengthUTF8`, zarr v3 settings)
- Reference (frontend contract, in /home/klaus/ws/aevidence, do not modify):
  `src/powernapps/single_cell/components/differential-expression/contrast-registry.schema.ts`;
  loader `src/powernapps/single_cell/store/zarr-utils.ts:354` (`loadColumnValues`)
  and `src/powernapps/single_cell/store/store.ts:1584-1640` (registry fetch/parse)

## Verification
1. `uv sync` succeeds with `decoupler` added; `python -c "import decoupler"` works.
2. Run the notebook end-to-end (needs network for `dc.op`); it writes
   `test-data/habib17-enrichment-test-data-format.zarr`.
3. Data-level checks (Python): `uns/de` still present and intact; open a few
   `uns/enrichment/es_XXX/` folders — every array length equals `gene_set_id`
   length; `scores` is non-increasing (pre-sort holds); NaN placeholders present
   for unsigned methods; `top_10_gene_set_ids` matches the first 10 `gene_set_id`.
4. Registry validates against `enrichment-registry.schema.ts` (a small ts/zod
   check or a mirrored Python assertion), and `cell_activities` obsm keys exist
   with matching `gene_set_ids` column labels.
5. Serve `python -m local_file_server ./test-data` and confirm the new store's
   registry + one `es_XXX` folder are fetchable over HTTP (frontend enrichment UI
   does not exist yet, so end-to-end stops at data/schema validation).
