# **Authoring Differential Expression Data for the Single Cell App**

A producer-side guide. It describes the **contract** the Single Cell app consumes,
so you can expose DE results from any tool — Scanpy, DESeq2, edgeR, limma, MAST,
or a bespoke pipeline — without app changes. Nothing here is specific to the
reference datasets in this repository; those are just one implementation of this
contract.

You produce exactly two things inside an existing AnnData Zarr store:

1. **`uns/de/contrast_registry.json`** — a manifest the app reads once at load
   time to discover comparisons and configure labels and chart availability.
2. **One folder of 1-D arrays per comparison** (`uns/de/<contrast_id>/`) — the
   per-gene numbers, loaded lazily when a user selects that comparison.

The app is a deliberately dumb consumer: it plots and tabulates what you write.
It does no statistics, no correction, no `-log10`, and no sorting of its own.
Every derived value it needs is either precomputed by you or declared in the
registry.

See `Enrichment Zarr Specification.md` for the gene-set analogue, which follows
the same layout with gene sets on the feature axis.

---

## **1\. What the app actually loads**

| Step | Path | Behaviour |
| :--- | :--- | :--- |
| Registry fetch | `/uns/de/contrast_registry.json` | **404 / 403 → silent.** The store layer maps both to "absent", so the DE tab is simply empty with only a console note. A store with no DE data is a normal case, not an error. |
| | | Any *other* fetch failure (5xx, CORS, network) or malformed JSON → red notification. |
| Registry validation | (Zod schema, see `contrast-registry.schema.ts`) | Invalid → red notification *"Differential expression contrasts are not supported"* + issue summary in console. **All contrasts are rejected, not just the bad one.** |
| Per-contrast arrays | `/uns/de/<contrast_id>/<array>` | A missing **required** array throws `Failed to load column "<name>": …`; unequal lengths throw `Column length mismatch: <name> (n ≠ m)`. |

The practical consequence of row one: **a typo'd path fails the same way as no DE
data at all** — silently. If the tab is empty, check the URL before suspecting
your arrays.

Both Zarr **v2 and v3** work; the app normalises dtypes across them.
Consolidated metadata is used when present but is a pure performance
optimisation — the app also works without it.

> **The registry is not a Zarr node.** It is a plain JSON file fetched by path,
> so it must be a real file inside the store directory. Consequences: metadata
> consolidation skips it (a harmless `ZarrUserWarning`), and a full AnnData
> round-trip (`ad.read_zarr(store)` → `adata.write_zarr(...)`) **drops it**,
> because AnnData only serialises nodes it recognises. Re-emit the registry after
> any such rewrite.

---

## **2\. Directory layout to produce**

```
dataset.zarr/
└── uns/
    └── de/
        ├── contrast_registry.json   # manifest — one entry per contrast folder
        ├── de_001/                  # one folder per comparison; id is yours to choose
        │   ├── gene_id/             # feature ids — must match the var index
        │   ├── scores/              # primary statistic + descending sort key
        │   ├── effect_size/         # signed effect (volcano x-axis)
        │   ├── pvals/               # raw p-values
        │   ├── pvals_adj/           # corrected p-values
        │   ├── significance/        # recommended: plot-ready -log10(p) (volcano y-axis)
        │   ├── pct_expr_target/     # % of target-group cells expressing
        │   ├── pct_expr_ref/        # % of reference/rest cells expressing
        │   ├── mean_expr_target/    # mean expression in the target group
        │   ├── mean_expr_ref/       # mean expression in the reference/rest group
        │   └── mean_expr/           # mean expression across all cells
        └── de_002/ …
```

Folder names are arbitrary identifiers — `de_001`, `wilcoxon_Sst_vs_rest`,
anything stable and URL-safe. The only requirement is that the folder name
equals the entry's `contrast_id`.

---

## **3\. The per-contrast arrays**

| Array | Type | Required | What the app does with it |
| :--- | :--- | :--- | :--- |
| **gene_id** | string | **yes** | Row identity; defines the expected length; links a table row to the expression matrix. |
| **scores** | float | **yes** | Table column; the sort key you must pre-sort by. |
| **effect_size** | float | **yes** | Volcano **x-axis** and, via a derived robust z-score, its **point colour**; drives the up/down split in the QC header. |
| **pvals** | float | **yes** | Table column. |
| **pvals_adj** | float | **yes** | Table column; fallback source for the volcano y-value. |
| **significance** | float | *recommended* | Volcano y-axis, read directly. Only array the app tolerates as absent. |
| **pct_expr_target** | float | **yes** | Table column + dot-plot size. |
| **pct_expr_ref** | float | **yes** | Table column + dot-plot size. |
| **mean_expr_target** | float | **yes** | Table column + dot-plot colour. |
| **mean_expr_ref** | float | **yes** | Table column + dot-plot colour. |
| **mean_expr** | float | **yes** | Volcano **point size** + table column. Has no `has_*` flag — always write it. |

### The volcano plot's visual channels

Four arrays feed the plot, and two of them are easy to miss because no `has_*`
flag advertises them:

| Channel | Source | Detail |
| :--- | :--- | :--- |
| **x** | `effect_size` | Direct. |
| **y** | `significance` | Direct (falls back to `-log10(pvals_adj)`, §4). |
| **size** | `mean_expr` | Radius 2–9 px, scaled across the contrast. |
| **colour** | `effect_size` → robust z-score | Diverging `RdBu`, domain pinned to −3 … 0 … +3. |

**Size** is why `mean_expr` is required and flagless: it is not an optional
statistic but a display channel the volcano is wired to by default. A NaN-filled
`mean_expr` is not an error — every point simply collapses to the same radius and
the "Mean Expression" column is blank.

**Colour** is computed in the frontend, not read from disk: it is the
median/MAD robust z-score of `effect_size`
(`(value − median) / (MAD × 1.4826)`), evaluated per contrast over that
contrast's genes. The scale's domain is fixed at ±3, so anything beyond three
robust deviations saturates rather than stretching the ramp. Two consequences for
you: the colouring is *relative to the contrast*, so the same gene can be coloured
differently in two contrasts; and it degrades to a single colour when
`effect_size` is constant or all-`NaN` (MAD of 0 yields z = 0 everywhere).

**Why colour is derived rather than stored.** It is a display encoding, not a
result: no label field describes it, no registry entry declares it, and it means
nothing outside "make the volcano readable" — and it is a cheap, deterministic
function of an array the app already holds in full to draw the x-axis. Contrast
`significance`, which *is* stored, precisely because it cannot be reconstructed:
`-log10(0) = ∞` forces a capping policy that is a judgment call, and the cap has
to agree with `significance_max`. That is the dividing line for this contract —
**precompute what the app cannot reconstruct, or cannot reconstruct without extra
I/O; let it derive what is cheap and unambiguous.** So do not ship a robust
z-score array; it would be ignored.

### Rule 1 — Stable arrays: never omit, always NaN-fill

All eleven arrays except `significance` are fetched **unguarded and in
parallel**. Omitting one aborts the whole contrast with
`Failed to load column "<name>"`, even if the registry says the metric is
unavailable. If your method does not produce a metric, write a full-length
float array of `NaN` and set the matching `has_*` flag to `false`.

The flags control **table column visibility and chart availability only** — they
never make an array optional:

| Flag `false` | Hides |
| :--- | :--- |
| `has_effect_size` | effect-size column, robust z-score column, and (with significance) the volcano plot |
| `has_significance` | p-value, adjusted-p and significance columns, QC header, and (with effect size) the volcano plot |
| `has_pct_expr_target` / `has_pct_expr_ref` | that % column and its dot-plot half |
| `has_mean_expr_target` / `has_mean_expr_ref` | that mean column and its dot-plot half |

The **volcano plot** appears only when `has_effect_size` **and**
`has_significance` are both `true` — it needs both axes. The anchor plot beside
it is independent of these flags and always shown.

### Rule 2 — One dimension, equal lengths

Every array is 1-D and every array in a folder has the same length. The app
takes the expected length from `gene_id` and hard-fails on any mismatch with
`Column length mismatch: <name> (n ≠ m)`. Do not rely on Zarr fill values to
pad — write the real length.

### Rule 3 — Pre-sort descending by `scores`

Sort every array in the folder by `scores` descending (`NaN` last) **before
writing**, using one shared permutation. The app slices `[0, N)` to show top
markers, so the on-disk order *is* the ranking. Break ties deterministically
(ascending p-value, then `gene_id`) so repeated runs are byte-comparable.

### Rule 4 — Accepted dtypes

The loader is permissive; you rarely need to convert anything:

- **Numeric** — any int or float width. Integers larger than
  `Number.MAX_SAFE_INTEGER` lose precision (with a console warning), so prefer
  `float64`/`float32` for statistics.
- **String** — Zarr v3 variable-length UTF-8 (`string`), Zarr v2 `U`/`S`, or any
  array tagged `encoding-type: "string-array"`. Variable-length UTF-8 is the most
  portable choice.
- **Boolean** — read as `"true"`/`"false"` strings.
- **AnnData categorical group** (`categories` + `codes` children) — decoded to
  per-row labels. Valid for `gene_id`, though a plain string array is simpler.

Missing numbers must be `NaN`, not `null`, not a sentinel like `-1` or `0`.

---

## **4\. The `significance` array**

`significance` is the **plot-ready volcano y-value**: whatever number you want on
that axis, already transformed. The app plots it verbatim and never transforms it,
which is precisely why `significance_label` exists — the label, not the app, says
what the values mean.

`-log10(adjusted p)` labelled `"-log10(FDR)"` is the conventional choice, but
nothing in the contract requires it. `-log10` of a raw p-value, a `-log10(q)`, a
posterior probability, a local FDR, a signed test statistic, or a bespoke
confidence score are all valid as long as the label describes them and larger
means *more significant* (the axis and the QC header's 0.05 threshold both assume
that direction).

Whatever transform you choose, three properties keep plots readable:

- **No infinities.** If your transform can diverge — `-log10` of an underflowed
  or exactly-zero p-value is the usual culprit — clamp it. A data-driven cap such
  as `ceil(max_finite × 1.5)` pins the diverged points just above the largest real
  value, keeping them visible and distinguishable without letting one point set
  the axis range. A fixed cap works too.
- **`NaN` for "not computed"**, never `0` — `0` reads as "p = 1", a real and very
  unsurprising result.
- **Direction:** larger = more significant.

**If you omit the array**, the app falls back to computing
`min(-log10(pvals_adj), 300)` itself. Two consequences: you inherit `-log10(FDR)`
semantics whether or not they fit your method, and the fallback reads `pvals_adj`
**only** — so a contrast whose p-values live in `pvals` (uncorrected; see §5)
plots an empty y-axis. Write the array.

---

## **5\. `contrast_registry.json`**

Accepted at the top level as either a bare JSON array of entries or an object
with a `contrasts` key (the app reads `raw.contrasts ?? raw`). The object form is
recommended — it leaves room for your own provenance keys, which are ignored:

```json
{
  "version": "1.0",
  "generated_at": "2026-08-04T10:00:00Z",
  "contrasts": [ { "contrast_id": "de_001", "…": "…" } ]
}
```

### Entry fields

Every field below is **required by the validator** unless marked optional. There
is no partial acceptance: one malformed entry rejects the whole registry.

| Field | Examples | Required | Notes |
| :--- | :--- | :--- | :--- |
| **contrast_id** | `"de_001"`, `"wilcoxon_Sst_vs_rest"` | yes | Must equal the folder name. |
| **contrast_column** | `"CellType"`, `"condition"`, `"leiden"` | yes | An **obs** column. Must resolve as categorical. |
| **group_1** | `"Sst"`, `"Tumor"`, `"7"` | yes | Target group; exact-matched against obs values. Numeric-looking cluster labels are still strings. |
| **group_2** | `null`, `"Control"`, `"Astro"` | yes (nullable) | Reference group. `null` = one-vs-rest. |
| **test_type** | `"one_vs_rest"`, `"pairwise"`, `"multiclass"`, `"treatment_vs_control"` | yes | Free-form; groups the contrast selector. Underscores render as spaces, so `"one_vs_rest"` shows as *"One vs rest"*. |
| **de_method** | `"wilcoxon"`, `"t-test"`, `"logreg"`, `"deseq2"`, `"edger"`, `"limma-voom"`, `"mast"` | yes | Free-form; shown in the selector label. |
| **subset_column** | `null`, `"region"`, `"timepoint"` | yes (nullable) | Obs column restricting the domain, or `null` for a whole-dataset test. |
| **subset_value** | `null`, `"Hippocampus"`, `"day_7"` | yes (nullable) | Value within `subset_column`, or `null`. |
| **has_effect_size** | `true` | yes | See the flag table in §3. |
| **has_significance** | `true` | yes | |
| **has_pct_expr_target** | `true` | yes | |
| **has_pct_expr_ref** | `true` | yes | |
| **has_mean_expr_target** | `true` | yes | |
| **has_mean_expr_ref** | `true` | yes | |
| **correction_method** | `"benjamini-hochberg"`, `"bonferroni"`, `"storey-q"`, `null` | yes (nullable) | Lower-case; rendered capitalised as *"Benjamini-hochberg corrected"*, or *"Not corrected"* when `null`. |
| **effect_size_label** | `"log2(Fold Change)"`, `"Cohen's d"`, `"Coefficient"` | yes (nullable) | Volcano x-axis title. Name the actual quantity in `effect_size`. |
| **significance_label** | `"-log10(FDR)"`, `"-log10(p-value)"`, `"-log10(q)"` | yes (nullable) | Volcano y-axis title. Must describe whatever transform you put in `significance` (§4). |
| **effect_size_max** | `8.5`, `null` | yes | `max(abs(effect_size))`. A number when `has_effect_size`, `null` otherwise. |
| **significance_max** | `30.0`, `null` | yes | `nanmax(significance)`. A number when `has_significance`, `null` otherwise. |
| **top_10_gene_ids** | `["CD3D", "IL7R", …]`, `[]` | yes | Selector preview text; ids come from `gene_id`. Pass `[]` if you have nothing to show. |

### How entries are presented

Worth knowing when choosing values, since these strings are user-facing:

- **Selector group** — `test_type` verbatim, with underscores replaced by spaces
  and the first letter capitalised: `"one_vs_rest"` → *"One vs rest"*.
- **Selector label** — `"<Subset_column>: <subset_value> 🠚 "` (or `"Global 🠚"`)
  then `"<Contrast_column>: <group_1> vs <group_2 ?? 'rest'>"`, then `(<De_method>)`.
  So a fully-specified entry reads *"Region: Hippocampus 🠚 CellType: Sst vs rest
  (Deseq2)"*.
- **Description** — `"Top genes: "` plus `top_10_gene_ids`, resolved through the
  var display-name map where possible.

### The four valid flag/label combinations

The schema is a union of four shapes. A label must be a non-empty string when
its metric is available and `null` when it is not — mixing them is the most
common validation failure:

| `has_effect_size` | `has_significance` | `effect_size_label` / `_max` | `significance_label` / `_max` | `correction_method` |
| :--- | :--- | :--- | :--- | :--- |
| `true` | `true` | string / number | string / number | string **or** `null` |
| `true` | `false` | string / number | `null` / `null` | `null` |
| `false` | `true` | `null` / `null` | string / number | string **or** `null` |
| `false` | `false` | `null` / `null` | `null` / `null` | `null` |

`correction_method: null` with `has_significance: true` is deliberately legal: it
means *"significance available, no multiple-testing correction applied"* and the
app renders *"Not corrected"*. Use it for tests you report uncorrected — and then
put those p-values in `pvals`, leave `pvals_adj` as `NaN`, label significance
`-log10(p-value)`, and **write the `significance` array** (see §4).

---

## **6\. Linking to the rest of the store**

The DE tab is not standalone — it drives the cell-level views, and those links
are plain string equality. Get them wrong and the tab still loads while features
silently do nothing.

- **`contrast_column` must be an obs column** and must load as *categorical* —
  an AnnData categorical group, or a string/bool array. A numeric obs column is
  read as numerical data and will never match, leaving the comparison stuck
  loading.
- **`group_1`, `group_2`, `subset_value` are exact string matches** against that
  column's values. No trimming, no case folding. `group_2: null` means
  one-vs-rest: every other cell (within the subset, if any) becomes the
  reference.
- **`gene_id` values should match the `var` index** (`var_names`). Clicking a
  table row looks the gene up by raw id to drive the expression overlay; an
  unmatched id makes the click a no-op. Use the same identifier space as `var`,
  and let the app handle symbol display separately.
- **Comparisons live inside one column.** For multi-factor designs, materialise a
  composite obs column first (e.g. `obs["condition_x_celltype"]`) and point
  `contrast_column` at that.

---

## **7\. Reference implementation**

Rather than duplicate writer code here, use the notebook in this repository as the
worked example — it is the reference producer for this contract and stays in sync
with it by being run:

**`calculate-differential-expression-data.ipynb`**

Helpers worth reading (or lifting) in particular:

| Function | What it does |
| :--- | :--- |
| `_write_array` | Writes one 1-D array, using `VariableLengthUTF8` for strings so other libraries can read them. |
| `_presort_contrast_arrays` | Builds the single `scores`-descending permutation (`np.lexsort`) and applies it to every array — Rule 3. |
| `_compute_significance` | Produces the plot-ready `significance` array with a data-driven cap, and returns the cap for `significance_max`. |
| `_compute_axis_bounds` | Derives `effect_size_max` / `significance_max` from the arrays. |
| `_extract_contrast_arrays` | Pulls one contrast out of the upstream results into the eleven arrays. |
| `_write_contrast_to_zarr` | Writes one `de_XXX/` folder and its registry entry together. |

One step this notebook does *not* do, worth knowing if you write folders directly
rather than through AnnData: tagging the new nodes with AnnData encoding attributes
(`encoding-type: "dict"` on groups, `"array"` / `"string-array"` on arrays) and
rebuilding the root consolidated metadata. Without it the app still reads the store
fine, but `anndata` cannot round-trip it cleanly. See `finalize_metadata` in
`calculate-gene-set-enrichment-data.ipynb`, which does exactly that for both
`uns/de` and `uns/enrichment`.

That notebook is the enrichment counterpart and shares the same conventions.
