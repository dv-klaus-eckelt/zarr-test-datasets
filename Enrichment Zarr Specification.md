# **Authoring Gene Set Enrichment Data for the Single Cell App**

A producer-side guide, the gene-set analogue of `DGE Zarr Specification.md`. It
describes the **contract** the Single Cell app consumes, so you can expose
enrichment results from any engine — decoupler, fgsea, GSVA, clusterProfiler,
GSEApy, or a bespoke pipeline — without app changes. Nothing here is specific to
the reference datasets in this repository; those are one implementation of this
contract.

The structure mirrors DE with **gene sets on the feature axis instead of genes**.
If you have already produced DE data, everything you know transfers; the
differences are called out as they arise.

You produce up to three things inside an existing AnnData Zarr store:

1. **`uns/enrichment/enrichment_registry.json`** — a manifest read once at load
   time to discover enrichments and configure labels and chart availability.
2. **One folder of 1-D arrays per enrichment** (`uns/enrichment/<enrichment_id>/`) —
   the per-gene-set numbers, loaded lazily on selection.
3. **Optional per-cell activity matrices** in `obsm/` — cells × gene sets, which
   power the "score every cell for this gene set" map overlay.

The app is a deliberately dumb consumer: no statistics, no correction, no
`-log10`, no sorting.

---

## **1\. What the app actually loads**

| Step | Path | Behaviour |
| :--- | :--- | :--- |
| Registry fetch | `/uns/enrichment/enrichment_registry.json` | **404 / 403 → silent.** The store layer maps both to "absent", so the enrichment tab is simply not offered. |
| | | Any *other* fetch failure (5xx, CORS, network) or malformed JSON → red notification. |
| Registry validation | (Zod schema, see `enrichment-registry.schema.ts`) | Invalid → red notification *"Gene set enrichment is not supported"* + issue summary in console. **All entries are rejected, not just the bad one.** |
| Per-enrichment arrays | `/uns/enrichment/<enrichment_id>/<array>` | Missing **required** array → `Failed to load column "<name>": …`; unequal lengths → `Column length mismatch: <name> (n ≠ m)`. |
| Activity matrix | `/obsm/<obsm_key>` | Missing → error on overlay activation only. A group instead of an array → `obsm/<key> is a group; expected a 2-D activity array`. |

This matches DE exactly. The practical consequence of row one: **a typo'd path
fails the same way as no enrichment data at all** — silently, because 404 and 403
both mean "this store has no enrichment". If the tab never appears, check the URL
before suspecting your data.

Both Zarr **v2 and v3** work. Consolidated metadata is used when present but is a
pure performance optimisation.

> **The registry is not a Zarr node.** Like the DE `contrast_registry.json`, it is
> a plain JSON file fetched by path, so it must be a real file inside the store
> directory. Metadata consolidation skips it (a harmless `ZarrUserWarning`), and a
> full AnnData round-trip (`ad.read_zarr(store)` → `adata.write_zarr(...)`)
> **drops it** while keeping the `uns/enrichment/*` arrays and `obsm/` matrices.
> Re-emit the registry after any such rewrite.

---

## **2\. Directory layout to produce**

```
dataset.zarr/
└── uns/
    ├── de/                                  # optional: enrichment can link back to it
    │   └── …
    └── enrichment/
        ├── enrichment_registry.json         # manifest — two arrays, see §5
        ├── es_001/                          # one folder per (comparison, method, collection)
        │   ├── gene_set_id/                 # gene set ids, e.g. HALLMARK_HYPOXIA
        │   ├── scores/                      # primary statistic + descending sort key
        │   ├── effect_size/                 # signed effect (volcano x-axis)
        │   ├── pvals/                       # raw p-values
        │   ├── pvals_adj/                   # corrected p-values
        │   ├── significance/                # recommended: plot-ready -log10(p) (y-axis)
        │   ├── set_size/                    # genes in the set, per the collection
        │   └── overlap_size/                # set genes present/tested in the data
        └── es_002/ …
└── obsm/
    ├── X_activity_mlm_progeny/              # optional: cells × gene sets
    └── X_activity_aucell_hallmark/
```

Folder names are arbitrary stable identifiers; the only requirement is that a
folder name equals its entry's `enrichment_id`.

One folder represents **one comparison scored against one collection by one
method**. Running three methods over two collections for five comparisons is 30
folders — that is the intended granularity, and the selector groups them for the
user (§5).

---

## **3\. The per-enrichment arrays**

| Array | Type | Required | What the app does with it |
| :--- | :--- | :--- | :--- |
| **gene_set_id** | string | **yes** | Row identity; defines the expected length; links a row to its activity column. |
| **scores** | float | **yes** | Volcano **point colour** + table column; the sort key you must pre-sort by. |
| **effect_size** | float | **yes** | Volcano x-axis; drives the up/down split in the QC header. |
| **pvals** | float | **yes** | Table column. |
| **pvals_adj** | float | **yes** | Table column; fallback source for the volcano y-value. |
| **significance** | float | *recommended* | Volcano y-axis, read directly. Only array the app tolerates as absent. |
| **set_size** | float | **yes** | Table column. |
| **overlap_size** | float | **yes** | Table column. |

Eight arrays against DE's eleven — there is no expression-summary block, and
`gene_set_id` replaces `gene_id`.

### The volcano plot's visual channels

| Channel | Source | Detail |
| :--- | :--- | :--- |
| **x** | `effect_size` | Direct. |
| **y** | `significance` | Direct (falls back to `-log10(pvals_adj)`, §4). |
| **size** | *none* | Flat 6 px for every gene set. |
| **colour** | `scores` | Diverging `RdBu-reverse` centred at 0; domain from the data. |

Two differences from the DE tab worth knowing, since they change what your arrays
need to carry:

- **Size is constant.** DE sizes its points by `mean_expr`; enrichment has no
  expression-summary array and encodes nothing in radius. This is why the
  enrichment contract needs no `mean_expr` equivalent.
- **Colour comes straight from `scores`**, not from a derived robust z-score as in
  DE, and the scale is centred at 0 — so `scores` is assumed **signed** for
  colouring purposes. For a signed statistic (NES, t-value) that reads correctly:
  depleted sets blue, enriched red. For a strictly positive statistic (an odds
  ratio, an AUCell score) every point lands on one half of the ramp, which is
  legible but wastes contrast; put the signed quantity in `scores` where your
  method has one.

Users can re-map size and colour to any numeric column from the plot's settings
menu — these are just the defaults they land on.

### Rule 1 — Stable arrays: never omit, always NaN-fill

All seven arrays except `significance` are fetched **unguarded and in parallel**.
Omitting one aborts the whole enrichment with `Failed to load column "<name>"`,
even if the registry declares the metric unavailable. Write a full-length
`NaN` array and set the matching `has_*` flag to `false`.

The flags control **table column visibility and chart availability only**:

| Flag `false` | Hides |
| :--- | :--- |
| `has_effect_size` | effect-size column, and (with significance) the volcano plot |
| `has_significance` | p-value, adjusted-p and significance columns, QC header, and (with effect size) the volcano plot |
| `has_set_size` | set-size column |
| `has_overlap` | overlap-size column |

The **volcano plot** appears only when `has_effect_size` **and**
`has_significance` are both `true` — it needs both axes. (It is the same
volcano-shaped component the DE tab uses, with gene sets as points.) The anchor
plot beside it is independent of these flags and always shown.

This flag pattern maps cleanly onto real methods:

- **Over-representation / ORA / Fisher** — unsigned. `effect_size` is `NaN` (or
  put log-odds there if you have it) and `has_effect_size: false`. Table and
  ranking still work; the volcano plot is hidden.
- **GSEA / fgsea** — signed NES plus p-values. Everything on.
- **GSVA / AUCell / ssGSEA** — a score per set with no test. All p-value arrays
  `NaN`, `has_significance: false`. These often pair with an `obsm/` activity
  matrix (§6).

### Rule 2 — One dimension, equal lengths

Every array is 1-D and equal-length within a folder. The expected length comes
from `gene_set_id`; a mismatch hard-fails with
`Column length mismatch: <name> (n ≠ m)`.

### Rule 3 — Pre-sort descending by `scores`

Sort every array by `scores` descending (`NaN` last) with one shared permutation
before writing. The app slices `[0, N)`, so on-disk order *is* the ranking. Break
ties deterministically (ascending p-value, then `gene_set_id`).

### Rule 4 — Accepted dtypes

Identical to DE: any int/float width (prefer `float64`/`float32`); strings as
Zarr v3 variable-length UTF-8, Zarr v2 `U`/`S`, or anything tagged
`encoding-type: "string-array"`; booleans read as `"true"`/`"false"`; AnnData
categorical groups decoded to labels. `set_size` and `overlap_size` are counts
but are read as floats — write them as `float64` so `NaN` is representable.

Missing numbers must be `NaN`, never `null`, `-1` or `0`.

---

## **4\. p-values: adjusted, raw, or both**

The two p-value arrays exist so you can be explicit about what you computed
rather than making the app guess. Populate what you have and declare it:

| You have | `pvals` | `pvals_adj` | `correction_method` | `significance_label` |
| :--- | :--- | :--- | :--- | :--- |
| Raw + corrected | raw | corrected | `"benjamini-hochberg"` | `"-log10(FDR)"` |
| Corrected only | `NaN` | corrected | `"benjamini-hochberg"` | `"-log10(FDR)"` |
| Raw only (uncorrected) | raw | `NaN` | `null` | `"-log10(p-value)"` |
| No test | `NaN` | `NaN` | `null` | `null` |

The third row is worth knowing about, because several libraries return a single
p-value matrix whose correction status varies by method. decoupler ≥ 2, for
instance, overwrites its p-value matrix with BH-adjusted values in place for every
method **except** `mlm` (`decoupler/mt/_run.py`: `if name != "mlm"`) — so an `mlm`
result routed into `pvals_adj` and labelled `-log10(FDR)` would be silently
mislabelled. Check your engine's behaviour per method rather than assuming.

`correction_method: null` together with `has_significance: true` is a deliberately
legal combination meaning *"significance available, no multiple-testing correction
applied"*; the app renders *"Not corrected"*. It is not missing metadata.

### The `significance` array

`significance` is the **plot-ready volcano y-value**: whatever number you want on
that axis, already transformed. The app plots it verbatim and never transforms it,
which is precisely why `significance_label` exists — the label, not the app, says
what the values mean.

`-log10` of whichever p-value array you populated is the conventional choice, but
nothing in the contract requires it. A `-log10(q)`, a permutation-based FDR, a
posterior probability, or a bespoke confidence score are all valid as long as the
label describes them and larger means *more significant* (the axis and the QC
header's 0.05 threshold both assume that direction).

Whatever transform you choose, three properties keep plots readable:

- **No infinities.** If your transform can diverge — `-log10` of an underflowed or
  exactly-zero p-value is the usual culprit — clamp it. A data-driven cap such as
  `ceil(max_finite × 1.5)` pins the diverged points just above the largest real
  value, keeping them visible without letting one point set the axis range. A
  fixed cap works too.
- **`NaN` for "not computed"**, never `0` — `0` reads as "p = 1".
- **Direction:** larger = more significant.

**If you omit the array**, the app falls back to computing
`min(-log10(pvals_adj), 300)` itself. Two consequences: you inherit `-log10(FDR)`
semantics whether or not they fit your method, and the fallback reads `pvals_adj`
**only** — so an uncorrected result (row three above) plots an empty y-axis. Write
the array.

---

## **5\. `enrichment_registry.json`**

The file must contain a JSON **object** with exactly these two keys, each holding
an array:

- **`contrast_enrichments`** — one entry per `es_XXX` folder. The enrichments
  themselves.
- **`cell_activities`** — one entry per `obsm/` activity matrix. Required to be
  present, but `[]` is fine if you have no activity matrices.

Both keys must exist or validation fails. This is stricter than the DE registry,
which also accepts a bare top-level array (`[{…}, {…}]`) because it only ever
carries one list. Enrichment has two lists to carry, so it always needs the
wrapping object — writing just the array of enrichments is rejected.

Any *additional* top-level keys are ignored, so use them freely for provenance:

```json
{
  "version": "1.0",
  "generated_at": "2026-08-04T10:00:00Z",
  "contrast_enrichments": [ { "enrichment_id": "es_001", "…": "…" } ],
  "cell_activities": []
}
```

### `contrast_enrichments[]`

Every field is **required by the validator** unless marked optional. One
malformed entry rejects the whole registry.

| Field | Examples | Required | Notes |
| :--- | :--- | :--- | :--- |
| **enrichment_id** | `"es_001"`, `"gsea_hallmark_Sst"` | yes | Must equal the folder name. |
| **enrichment_method** | `"gsea"`, `"ora"`, `"gsva"`, `"aucell"`, `"mlm"`, `"ulm"`, `"fgsea"`, `"ssgsea"` | yes | Free-form. Part of the selector group label. |
| **gene_set_collection** | `"MSigDB Hallmark"`, `"GO:BP"`, `"Reactome"`, `"PROGENy"`, `"CollecTRI"` | yes | Free-form label, written the way you want it read. Part of the group label, **and the join key to `cell_activities`** (§6). |
| **contrast_column** | `"CellType"`, `"condition"`, `"leiden"` | yes | An **obs** column. Must resolve as categorical. |
| **group_1** | `"ASC1"`, `"Tumor"`, `"7"` | yes | Target group; exact-matched against obs values. Numeric-looking cluster labels are still strings. |
| **group_2** | `null`, `"Control"` | yes (nullable) | Reference group. `null` = one-vs-rest. |
| **test_type** | `"one_vs_rest"`, `"pairwise"`, `"multiclass"` | yes | Free-form comparison semantics. |
| **subset_column** | `null`, `"region"` | yes (nullable) | Obs column restricting the domain, or `null` for a whole-dataset run. |
| **subset_value** | `null`, `"Hippocampus"` | yes (nullable) | Value within `subset_column`, or `null`. |
| **n_gene_sets** | `50` | yes | Array length. Validated, not consumed. |
| **input_statistic** | `"de_scores"`, `"top_50_up_genes"`, `"log_normalized_counts"`, `"ranked_logfc"` | yes | What you fed the method — useful provenance when the same collection is scored several ways. Validated, not consumed. |
| **has_effect_size** | `true` | yes | See the flag table in §3. |
| **has_significance** | `true` | yes | |
| **has_set_size** | `true` | yes | |
| **has_overlap** | `true` | yes | |
| **has_leading_edge** | `false` | yes | Reserved; see below. Write `false`. |
| **correction_method** | `"benjamini-hochberg"`, `"bonferroni"`, `"permutation-fdr"`, `null` | yes (nullable) | Lower-case; rendered capitalised as *"Benjamini-hochberg corrected"*, or *"Not corrected"* when `null`. |
| **effect_size_label** | `"Normalized Enrichment Score"`, `"MLM t-value"`, `"GSVA enrichment score"`, `"ln(Odds Ratio)"` | yes (nullable) | Volcano x-axis title. Name the actual quantity in `effect_size`. |
| **significance_label** | `"-log10(FDR)"`, `"-log10(p-value)"`, `"-log10(q)"` | yes (nullable) | Volcano y-axis title. Must describe whatever transform you put in `significance` (§4). |
| **top_10_gene_set_ids** | `["HALLMARK_HYPOXIA", …]`, `[]` | yes | Selector preview (first three shown); ids come from `gene_set_id`. `[]` is fine. |
| **source_contrast_id** | `"de_001"`, `null` | optional | Links back to a `uns/de` folder. Omit or `null` when the enrichment did not come from one. |
| **effect_size_max** | `3.0`, `null` | yes | `max(abs(effect_size))`. A number when `has_effect_size`, `null` otherwise. |
| **significance_max** | `12.0`, `null` | yes | `nanmax(significance)`. A number when `has_significance`, `null` otherwise. |

Note that `contrast_column` + `group_1` are required even when the enrichment was
not derived from a DE contrast: they are how the app maps an enrichment to actual
cells for the comparison plot and activity overlay. Every entry must name a real
obs grouping.

#### `has_leading_edge` — reserved for driver genes

The **leading edge** is the subset of genes actually responsible for a set's
enrichment: for GSEA, the genes contributing up to the running-enrichment peak;
for over-representation, the intersection of your input gene list with the set's
members. It answers "*which* genes made this set come up?".

This is a **planned feature**, and the flag is its reserved slot. The intended
shape is a `leading_edge_genes` per-folder string array — one delimiter-joined gene
list per set, empty string where the method has no leading-edge concept — which
keeps the equal-length rule and lets the table expand a row to show drivers.

Until that array is specified and the app reads it, **write `false`**. If you
already have driver genes, hold them outside this contract rather than inventing a
column: the delimiter and array name need to be agreed once so every producer
writes the same thing, otherwise the feature arrives to inconsistent data.


### The four valid flag/label combinations

The schema is a union of four shapes. A label must be a non-empty string when its
metric is available and `null` when it is not — mixing them is the most common
validation failure:

| `has_effect_size` | `has_significance` | `effect_size_label` / `_max` | `significance_label` / `_max` | `correction_method` |
| :--- | :--- | :--- | :--- | :--- |
| `true` | `true` | string / number | string / number | string **or** `null` |
| `true` | `false` | string / number | `null` / `null` | `null` |
| `false` | `true` | `null` / `null` | string / number | string **or** `null` |
| `false` | `false` | `null` / `null` | `null` / `null` | `null` |

### How entries are presented

Worth knowing when choosing field values, since these strings are user-facing:

- **Selector group** — `"<Enrichment_method> · <gene_set_collection>"`, e.g.
  *"Gsea · MSigDB Hallmark"*. Only the first letter is capitalised, so write
  `gene_set_collection` the way you want it read.
- **Selector label** — `"<Subset_column>: <subset_value> 🠚 "` (or `"Global 🠚"`)
  followed by `"<Contrast_column>: <group_1> vs <group_2 ?? 'rest'>"`.
- **Description** — `"Top: "` plus the first three `top_10_gene_set_ids`.

---

## **6\. Per-cell activity matrices (`cell_activities[]`)**

Optional. Each entry describes one `obsm/` matrix of per-cell scores, which lets a
user click a gene set in the table and colour the embedding by that set's activity
in every cell.

| Field | Examples | Notes |
| :--- | :--- | :--- |
| **obsm_key** | `"X_activity_mlm_progeny"`, `"X_activity_aucell_hallmark"` | Key under `obsm/`. Must be a **2-D array**, not a group. |
| **gene_set_collection** | `"PROGENy"`, `"MSigDB Hallmark"` | **Must string-match** an enrichment entry's `gene_set_collection`, or the overlay never activates. |
| **method** | `"mlm"`, `"ulm"`, `"aucell"`, `"gsva"` | Free-form. |
| **value_label** | `"MLM t-value"`, `"AUCell score"`, `"GSVA enrichment score"` | Legend label for the overlay. Display only — the colour scale is inferred from the values, not this string. |
| **gene_set_ids** | `["Androgen", "EGFR", …]` | Column order of the matrix. Looked up by `indexOf`. |
| **n_gene_sets** | `14` | Must equal `gene_set_ids` length and the matrix's column count. |

Three requirements that are easy to get wrong:

1. **Shape** — `(n_obs, n_gene_sets)`, rows in the same order as `obs`. The app
   reads `shape[0]`/`shape[1]` and slices a column; a transposed matrix produces
   silently wrong colouring.
2. **`gene_set_ids` order must match the matrix columns exactly.** The gene set
   clicked in the table is resolved by `gene_set_ids.indexOf(id)`; a
   reordered list maps the user to the wrong column. Ids should also match the
   `gene_set_id` values in the corresponding enrichment folders, otherwise the
   lookup fails and the click is a no-op.
3. **Signedness is inferred from the values, not from your labels.** If the loaded
   column contains any negative value the overlay uses a diverging scale centred
   on zero; a strictly non-negative column gets a sequential scale. You do not
   declare this anywhere — write honest scores and the colouring follows. (Earlier
   builds pattern-matched `value_label` against `/t-value/i` for this, so a method
   whose label did not contain "t-value" was mis-coloured. Fixed; `value_label` is
   now purely a legend string.)

A matching entry for a matrix written as `obsm/X_activity_gsva_hallmark`:

```json
{
  "obsm_key": "X_activity_gsva_hallmark",
  "method": "gsva",
  "gene_set_collection": "MSigDB Hallmark",
  "value_label": "GSVA enrichment score",
  "gene_set_ids": ["HALLMARK_HYPOXIA", "HALLMARK_GLYCOLYSIS"],
  "n_gene_sets": 2
}
```

`float32` is plenty for the matrix itself. See `write_obsm_array` and
`run_cell_activities` in the notebook (§8) for the writing side.

---

## **7\. Linking to the rest of the store**

Plain string equality throughout. Get these wrong and the tab still loads while
features silently do nothing.

- **`contrast_column` must be an obs column** loading as *categorical* (AnnData
  categorical group, or string/bool array). A numeric column never matches and
  leaves the comparison stuck loading.
- **`group_1`, `group_2`, `subset_value` are exact string matches** against that
  column's values. `group_2: null` means one-vs-rest.
- **`gene_set_collection` joins enrichments to activities.** Same spelling,
  spacing and case in both arrays, or the overlay silently stays unavailable.
- **`gene_set_id` values should be consistent** between an enrichment folder and
  the matching `cell_activities[].gene_set_ids`.
- **`source_contrast_id`, when set, should name an existing `uns/de` folder.**
  Optional — omit it for enrichments computed straight from expression.

---

