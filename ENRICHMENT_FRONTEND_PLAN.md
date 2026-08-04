# Enrichment Tab — single_cell app implementation plan (future session)

> Status: planned, not started. The **data** side is done (this repo). This plan is
> for the **frontend** work in the `aevidence` repo. Re-verify frontend file
> paths/line numbers before editing — they drift.

## Context

The `zarr-test-datasets` repo now ships a gene-set-enrichment test store,
`test-data/habib17-enrichment-test-data-format.zarr`, whose schema deliberately
mirrors the differential-expression (DE) store so the app can reuse the DE
loading/rendering patterns. The single_cell app has **no** enrichment feature
yet; this plan defines the tab that consumes the store.

Authoritative references (this repo):
- `Enrichment Zarr Specification.md` — directory tree, per-folder arrays, registry tables, rules.
- `enrichment-registry.schema.ts` — Zod schema (copy into the frontend).
- `DGE Zarr Specification.md` / `contrast-registry.schema.ts` — the DE analog the app already implements.

## The data contract (what the app reads)

Store layout (self-contained: DE store + enrichment):
```
habib17-enrichment-test-data-format.zarr/
├── uns/
│   ├── de/…                                  # unchanged DE contrasts (existing app feature)
│   └── enrichment/
│       ├── enrichment_registry.json          # manifest: contrast_enrichments[] + cell_activities[]
│       └── es_001/ … es_055/                 # one folder per (contrast × method × collection)
└── obsm/
    ├── X_activity_mlm_progeny/               # per-cell activity matrices (cells × gene sets)
    ├── X_activity_ulm_collectri/
    └── X_activity_aucell_hallmark/
```

There are **two complementary representations**, both keyed on gene sets:

### A. `contrast_enrichments[]` → `es_XXX/` folders (group-level, volcano-style)
One entry per `(contrast × method × collection)`. Feature axis = gene sets.
Directly analogous to a DE contrast; links back to DE via `source_contrast_id`.

- Per-folder arrays (all equal length, pre-sorted descending by `scores`):
  `gene_set_id`, `scores`, `effect_size`, `pvals`, `pvals_adj`, `significance`,
  `set_size`, `overlap_size`.
- `significance` is the **precomputed** volcano y-value (`-log10(pvals_adj)` with a
  data-driven cap) — plot it directly, do **not** recompute `-log10`.
- Registry fields per entry: `enrichment_id`, `source_contrast_id`, `group_1`,
  `group_2`, `test_type`, `subset_column`, `subset_value`, `contrast_column`,
  `enrichment_method`, `gene_set_collection`, `n_gene_sets`, `input_statistic`,
  `has_effect_size`, `has_significance`, `has_set_size`, `has_overlap`,
  `has_leading_edge`, `correction_method`, `effect_size_label`,
  `significance_label`, `top_10_gene_set_ids` (+ optional `effect_size_max`,
  `significance_max`, `feature_type`).
- Current data: 55 folders — gsea×Hallmark (15), ora×Hallmark (15), gsea/ora×GO:BP
  (5 each), ora×OmniPath (5), mlm×PROGENy (5), ulm×CollecTRI (5). Two branches
  exercised: 30 effect+significance, 25 significance-only (ORA has no signed effect).

### B. `cell_activities[]` → `obsm/X_activity_*` matrices (single-cell, UMAP overlay)
One entry per `(method × collection)` run over the full expression matrix.
Shape = cells × gene sets. Row axis is identical to `X` / `obsm/X_umap`, so a
column is a drop-in replacement for a gene when coloring the UMAP.

- Registry fields per entry: `obsm_key`, `method`, `gene_set_collection`,
  `value_label`, `gene_set_ids` (column order), `n_gene_sets`.
- The matrix is a bare 2-D float array; the column for a pathway is
  `gene_set_ids.indexOf(name)`. Fetch that column and map onto `X_umap`.
- Current data: `X_activity_mlm_progeny` (14), `X_activity_ulm_collectri` (481),
  `X_activity_aucell_hallmark` (50).

**Relationship:** same collections/methods, different granularity and input.
`contrast_enrichments` answers "which gene sets are enriched in this cell
type/contrast?" (bar/volcano over gene sets); `cell_activities` answers "what is
this pathway's activity in each cell?" (continuous UMAP overlay). They are not
mathematically derived from each other and have no per-entry foreign key; they
share the gene-set-id namespace per collection.

## Proposed frontend implementation (two sub-features)

### 1. Contrast-based enrichment view — mirror the DE tab
Almost a structural clone of the DE feature, with genes → gene sets.

- **Schema/loader:** copy `enrichment-registry.schema.ts` into the app; fetch
  `/uns/enrichment/enrichment_registry.json` and `safeParse` with
  `EnrichmentsArraySchema` (+ `CellActivitiesArraySchema`), exactly like
  `store.ts` does for the DE registry (~`store.ts:1592-1613`). Store in Redux next
  to `deContrasts`.
- **Per-folder arrays:** reuse `loadColumnValues` (`store/zarr-utils.ts:354`) over
  `/uns/enrichment/<enrichment_id>/`, reading `gene_set_id`, `scores`,
  `effect_size`, `significance`, `pvals_adj`, `set_size`, `overlap_size`
  (mirror the hard-coded `Promise.all` in `DifferentialExpression.tsx:121-134`).
- **Selector:** an enrichment picker analogous to `ContrastSelector.tsx`, grouped
  by `enrichment_method` + `gene_set_collection`; labels from `group_1`/`group_2`,
  `top_10_gene_set_ids` for previews.
- **Plot:** a volcano-style `ScatterPlotFromTable` — x = `effect_size` (label
  `effect_size_label`), y = `significance` (label `significance_label`), points =
  gene sets. For ORA (unsigned, `has_effect_size=false`) fall back to a ranked
  bar/dot plot of `significance` or `scores`. Gate on `has_significance` /
  `has_effect_size` exactly like `useDEColumnVisibility.ts` / `DifferentialExpression.tsx:182`.
- **Table:** analog of `useDETableColumns.tsx` with columns `gene_set_id`,
  `scores`, `effect_size`, `significance`, `pvals_adj`, `set_size`, `overlap_size`
  (show/hide via the `has_*` flags).
- Because arrays are pre-sorted by `scores`, "top N gene sets" is a `slice(0, N)`.

### 2. Per-cell activity overlay — UMAP coloring
A new, smaller capability (no DE analog):

- List `cell_activities[]` in a selector (by `method` + `gene_set_collection`).
- On pathway pick: `col = gene_set_ids.indexOf(name)`, fetch that column from
  `obsm/<obsm_key>`, color the existing UMAP scatter (`obsm/X_umap`) by it —
  same code path as coloring by a gene, since rows align to cells.
- Colormap by `value_label`: `AUCell enrichment score` is ≥0 → sequential;
  `MLM t-value` / `ULM t-value` are signed → diverging, centered at 0.
- Naming gotcha: Hallmark ids drop the `HALLMARK_` prefix (e.g. `APOPTOSIS`, not
  `HALLMARK_APOPTOSIS`). Always resolve names from `gene_set_ids`.

## Reuse map (existing DE code to model on, `aevidence`)
`src/powernapps/single_cell/…`
- `store/zarr-utils.ts:354` `loadColumnValues` (generic array loader — reuse as-is).
- `store/store.ts:1592-1613` registry fetch/parse; `:2644` current-contrast selector.
- `components/differential-expression/`: `DifferentialExpression.tsx`,
  `ContrastSelector.tsx`, `useDERowsAndStats.ts`, `useDEColumnVisibility.ts`,
  `useDETableColumns.tsx`, `QCHeader.tsx`, `VolcanoData.ts`,
  `DEAnchorScatterPlot.tsx`, `contrast-registry.schema.ts`.
- `ScatterPlotFromTable` (auto-derives axis domains from data — see note below).
- Find where DE is registered as a tab (StudyPage / powernapp routing) and add an
  "Enrichment" tab beside it.

## Prerequisite / coupled change: consume the precomputed `significance` (DE first)
The data now ships a `significance` array in **both** DE and enrichment folders,
so the app should stop computing `-log10` itself. Doing this in the DE tab first
de-risks the enrichment tab (same pattern):
- `useDERowsAndStats.ts:50-56` — replace `Math.min(-Math.log10(pvalsAdj[i]), 300)`
  with a direct read of the `significance` value (map to the existing
  `negLog10PvalsAdj` field); delete the cap/`TODO`.
- `DifferentialExpression.tsx:123-146` — add `loadColumnValues(…, 'significance')`
  to the `Promise.all`, destructure, `lengths` check, return object.
- Optionally consume `significance_max` for the y-domain (currently a **dead
  field** — `ScatterPlotFromTable` auto-derives the domain from data, which already
  respects the baked-in cap, so this is cosmetic).

## Schema optionality note
In this repo, `effect_size_max`, `significance_max`, and `feature_type` were made
`.optional()` (they're unused by the app; the cap is baked into the `significance`
array and domains auto-infer). The **frontend's own copy** of the DE schema still
*requires* the two `_max` fields, so if you ever want producers to omit them,
mirror the optionality into the app's `contrast-registry.schema.ts` in the same
change. The new enrichment schema is already optional-tolerant.

## Open decisions for the implementation session
- ORA (unsigned) rendering: ranked bar/dot vs. volcano with x hidden.
- Whether the enrichment tab is one tab with a method/collection switch, or nested
  views; and whether to surface the DE↔enrichment link (`source_contrast_id`).
- Whether per-cell activity overlay lives in the enrichment tab or in the existing
  UMAP/expression coloring control (it's just another color source).
- `set_size`/`overlap_size` presentation (dot size? table only?).

## Verification (that session)
1. Enrichment registry `safeParse`s; tab lists 55 enrichments + 3 activities.
2. Selecting an enrichment loads its `es_XXX` arrays (lengths equal `gene_set_id`),
   renders the plot with correct axis labels, and the table shows the right columns.
3. `significance` plotted directly (no client-side `-log10`); zero-p gene sets sit
   at the capped top with a visible gap.
4. Picking a pathway (e.g. `APOPTOSIS` from `X_activity_aucell_hallmark`) colors the
   UMAP; signed activities use a diverging colormap.
5. Serve `python -m local_file_server ./test-data` (in this repo) and point the app
   at `…/habib17-enrichment-test-data-format.zarr/` (see README link row).
