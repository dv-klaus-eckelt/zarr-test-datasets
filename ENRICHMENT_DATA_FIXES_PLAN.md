# Enrichment test-data fixes plan

> Companion to the frontend "Gene Set Enrichment Tab" work in the `aevidence` repo.
> Building the frontend surfaced two data gaps to address at the source. This plan covers
> `calculate-gene-set-enrichment-data.ipynb`, `enrichment-registry.schema.ts`, and
> `Enrichment Zarr Specification.md`.

## Context

The generated store `test-data/habib17-enrichment-test-data-format.zarr` currently:

1. **Emits `has_leading_edge = false` for all 55 enrichments** and writes no leading-edge /
   driver-gene array (`has_leading_edge=False` is hardcoded in `write_enrichment_folder`).
   This blocks the "expand a pathway → show its driver genes" view that the frontend would
   otherwise offer; shipping leading-edge data re-enables it later.
2. **Hard-ties every enrichment to a DE *folder*** — `main()` builds each entry from
   `de_by_id[cid]` and always sets `source_contrast_id`. An enrichment that isn't derived
   from a precomputed DE contrast cannot be represented today.
   - Decision: **only the DE-folder link (`source_contrast_id`) is optional.** Every
     enrichment must still carry `contrast_column` + `group_1` (a real `.obs` grouping) so
     the single-cell app can always offer its analysis options (anchor-plot comparison +
     per-cell activity overlay). Do **not** make `contrast_column`/`group_1` nullable.

## D1. Emit leading-edge / driver genes

- Add a `leading_edge_genes` per-folder array in each `es_XXX/`: a **string array of length
  `n_gene_sets`**, each element a delimiter-joined (e.g. `;`) list of driver gene ids, empty
  string where N/A. This preserves the "Stable Arrays" (equal length) rule; the frontend
  splits on the delimiter. `_write_array` already handles string arrays (see `gene_set_id`).
- **ORA:** drivers = intersection of the ORA input gene list (`top_{ORA_N_UP}_up_genes`)
  with each set's members — cheap; compute in `assemble_arrays` (pass in the input gene list
  and `collection["net"]` membership).
- **GSEA:** the leading-edge subset (genes contributing up to the running-ES peak). First
  check whether `dc.mt.gsea` returns leading-edge membership; if not, compute it fgsea-style
  from the ranked signature + set membership, or (fallback) leave GSEA empty with
  `has_leading_edge=False` for GSEA only.
- **MLM / ULM / AUCell:** no leading-edge concept → empty array + `has_leading_edge=False`.
- In `write_enrichment_folder`, replace the hardcoded `has_leading_edge=False` with the
  per-entry populated flag.
- Update the `Enrichment Zarr Specification.md` array table with the new `leading_edge_genes`
  row and its string encoding.

## D2. Support enrichments not derived from a DE folder

- **Schema** (`enrichment-registry.schema.ts`): relax **only** `source_contrast_id` in
  `SharedEnrichmentFields` to `.nullable()` (or `.nullable().optional()`). **Keep
  `contrast_column` and `group_1` required** (`group_2` is already nullable).
  - The frontend copy in `aevidence`
    (`src/powernapps/single_cell/components/enrichment/enrichment-registry.schema.ts`)
    already applies this exact relaxation — keep the two in sync.
- **Notebook:** extend `EnrichJob` with an optional "standalone" descriptor: an explicit
  input gene list + a **still-required** `contrast_column` + `group_1` (i.e. the `.obs`
  grouping the gene sets pertain to), and no `contrasts` / `de_by_id` link. Add ≥1 such job
  to `build_jobs()` — e.g. ORA of a curated marker set for a specific cell type,
  `input_statistic="marker_list"`.
- **`write_enrichment_folder` / `main()`:** for a standalone job, set
  `source_contrast_id=None` but still populate `contrast_column`, `group_1`
  (and `group_2` / `subset_*` / `test_type` from the job); otherwise inherit from
  `de_by_id[cid]` as today.
- This exercises the frontend path where `source_contrast_id` is null yet the anchor-plot
  comparison and activity overlay still work.

## D3. Minor cleanups (optional; coordinate with the frontend loader/table)

- `pvals` is always `NaN` (decoupler returns adjusted p-values only). Either drop the
  `pvals` array + its frontend load/column, or keep and document. **Recommendation: drop**
  it to avoid an always-empty column. (If dropped, remove `loadColumnValues(..., 'pvals')`
  and the `pvals` table column in the `aevidence` enrichment tab.)
- Leave `significance_max` / `effect_size_max` / `feature_type` as-is (already optional and
  harmless; not consumed by the frontend).

## D4. Docs

Update `Enrichment Zarr Specification.md`: add the `leading_edge_genes` array row + its
string encoding, document that `source_contrast_id` may be null (standalone enrichments,
with `contrast_column`/`group_1` still required), and — if D3 is applied — remove `pvals`.

## Cross-repo sequencing

Part 2 is independent of the frontend work but **enables** an optional frontend follow-up:
once `leading_edge_genes` ships, the enrichment table can add an expandable row listing a
pathway's driver genes (and optionally link them into the DE gene view). Recommended order:

1. Ship the frontend enrichment tab against the **current** data (done in `aevidence`).
2. Apply this plan (D1–D4) and **regenerate** the store.
3. Optionally add the expandable driver-gene row in the frontend.

## Verification

Re-run `calculate-gene-set-enrichment-data.ipynb`, then:

- Assert every `es_XXX/` has a `leading_edge_genes` array of length `n_gene_sets`.
- Assert ≥1 registry entry has `source_contrast_id == None` (and still has
  `contrast_column` + `group_1`).
- `safeParse` the registry with the relaxed schema (both the repo copy and the `aevidence`
  copy).
- Spot-check that ORA driver genes ⊆ (set members ∩ input gene list).
