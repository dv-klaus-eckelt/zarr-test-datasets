import { z } from "zod";

/**
 * Fields shared by every entry in
 * uns/enrichment/enrichment_registry.json -> contrast_enrichments[].
 *
 * An enrichment entry mirrors a DE contrast entry but the feature axis is gene
 * sets, not genes. Each entry links back to the DE contrast it was computed from
 * via `source_contrast_id`.
 */
const SharedEnrichmentFields = {
  enrichment_id: z
    .string()
    .min(1)
    .describe(
      "Stable enrichment identifier and folder key (for example es_001).",
    ),
  source_contrast_id: z
    .string()
    .min(1)
    .nullable()
    .optional()
    .describe(
      "DE contrast this enrichment was computed from (for example de_001); links to uns/de. Null/absent when not derived from a DE folder.",
    ),
  group_1: z
    .string()
    .min(1)
    .describe(
      "Target group label inherited from the source contrast (shown in UI selectors/titles).",
    ),
  group_2: z
    .string()
    .nullable()
    .describe(
      "Reference group label inherited from the source contrast. Null for one-vs-rest and current multiclass outputs.",
    ),
  test_type: z
    .string()
    .min(1)
    .describe(
      "Frontend grouping label for comparison semantics inherited from the source contrast.",
    ),
  subset_column: z
    .string()
    .nullable()
    .describe(
      "Optional obs column that restricted the source contrast domain. Null when no subset restriction was applied.",
    ),
  subset_value: z
    .string()
    .nullable()
    .describe(
      "Optional subset value within subset_column. Null when the source contrast ran on the full dataset.",
    ),
  contrast_column: z
    .string()
    .min(1)
    .describe(
      "Obs column that defined the compared groups in the source contrast (for example CellType).",
    ),
  enrichment_method: z
    .string()
    .min(1)
    .describe(
      "Enrichment algorithm identifier as a free-form string (for example gsea, ora, mlm, ulm).",
    ),
  gene_set_collection: z
    .string()
    .min(1)
    .describe(
      "Gene set collection label (for example MSigDB Hallmark, GO:BP, OmniPath, PROGENy, CollecTRI).",
    ),
  n_gene_sets: z
    .number()
    .int()
    .nonnegative()
    .describe("Number of gene sets (rows/arrays) written for this enrichment."),
  input_statistic: z
    .string()
    .min(1)
    .describe(
      "Per-gene input fed to the method (for example de_scores, top_50_up_genes).",
    ),
  has_effect_size: z
    .boolean()
    .describe("Whether a signed effect-size metric is available for this enrichment."),
  has_significance: z
    .boolean()
    .describe("Whether significance metrics are available for this enrichment."),
  has_set_size: z
    .boolean()
    .describe("Whether gene set sizes (genes per set from the collection) are available."),
  has_overlap: z
    .boolean()
    .describe("Whether overlap sizes (set genes present/tested in the data) are available."),
  has_leading_edge: z
    .boolean()
    .describe("Whether leading-edge / overlap driver gene ids are available."),
  top_10_gene_set_ids: z
    .array(z.string().min(1))
    .describe("Top 10 gene set ids ranked by score in descending order."),
};

const EffectAndSignificanceEnrichmentSchema = z.object({
  ...SharedEnrichmentFields,
  has_effect_size: z.literal(true),
  has_significance: z.literal(true),
  correction_method: z
    .string()
    .min(1)
    .nullable()
    .describe(
      'P-value correction method applied to the significance metric, or null when the method reports uncorrected p-values (for example decoupler mlm). Rendered as "Not corrected" when null.',
    ),
  effect_size_label: z
    .string()
    .min(1)
    .describe("Human-readable effect-size label when effect-size metrics are available."),
  significance_label: z
    .string()
    .min(1)
    .describe("Significance metric label when significance metrics are available."),
  effect_size_max: z
    .number()
    .nonnegative()
    .describe("Absolute effect-size bound when effect-size metrics are available."),
  significance_max: z
    .number()
    .nonnegative()
    .describe("Data-driven significance cap (plot y-max) when significance metrics are available."),
});

const EffectOnlyEnrichmentSchema = z.object({
  ...SharedEnrichmentFields,
  has_effect_size: z.literal(true),
  has_significance: z.literal(false),
  correction_method: z
    .null()
    .describe("P-value correction method must be null when significance metrics are unavailable."),
  effect_size_label: z
    .string()
    .min(1)
    .describe("Human-readable effect-size label when effect-size metrics are available."),
  significance_label: z
    .null()
    .describe("Significance metric label must be null when significance metrics are unavailable."),
  effect_size_max: z
    .number()
    .nonnegative()
    .describe("Absolute effect-size bound when effect-size metrics are available."),
  significance_max: z
    .null()
    .describe("Significance cap must be null when significance metrics are unavailable."),
});

const SignificanceOnlyEnrichmentSchema = z.object({
  ...SharedEnrichmentFields,
  has_effect_size: z.literal(false),
  has_significance: z.literal(true),
  correction_method: z
    .string()
    .min(1)
    .nullable()
    .describe(
      'P-value correction method applied to the significance metric, or null when the method reports uncorrected p-values (for example decoupler mlm). Rendered as "Not corrected" when null.',
    ),
  effect_size_label: z
    .null()
    .describe("Effect-size label must be null when effect-size metrics are unavailable."),
  significance_label: z
    .string()
    .min(1)
    .describe("Significance metric label when significance metrics are available."),
  effect_size_max: z
    .null()
    .describe("Effect-size bound must be null when effect-size metrics are unavailable."),
  significance_max: z
    .number()
    .nonnegative()
    .describe("Data-driven significance cap (plot y-max) when significance metrics are available."),
});

const ScoreOnlyEnrichmentSchema = z.object({
  ...SharedEnrichmentFields,
  has_effect_size: z.literal(false),
  has_significance: z.literal(false),
  correction_method: z
    .null()
    .describe("P-value correction method must be null when significance metrics are unavailable."),
  effect_size_label: z
    .null()
    .describe("Effect-size label must be null when effect-size metrics are unavailable."),
  significance_label: z
    .null()
    .describe("Significance metric label must be null when significance metrics are unavailable."),
  effect_size_max: z
    .null()
    .describe("Effect-size bound must be null when effect-size metrics are unavailable."),
  significance_max: z
    .null()
    .describe("Significance cap must be null when significance metrics are unavailable."),
});

/**
 * Method-agnostic enrichment schema.
 *
 * Entry schema for uns/enrichment/enrichment_registry.json -> contrast_enrichments[].
 *
 * Frontend intent:
 * - Discover available gene set enrichment results.
 * - Configure chart labels and metric bounds without loading gene-set-level arrays.
 * - Infer optional visualizations from available metrics and arrays.
 *
 * This intentionally does not discriminate on `enrichment_method` so collaborators
 * can provide outputs from any enrichment engine (for example decoupler GSEA/ORA/
 * MLM/ULM, fgsea, GSVA, or custom pipelines) without schema updates.
 */
export const EnrichmentSchema = z.union([
  EffectAndSignificanceEnrichmentSchema,
  EffectOnlyEnrichmentSchema,
  SignificanceOnlyEnrichmentSchema,
  ScoreOnlyEnrichmentSchema,
]);

/**
 * Entry schema for uns/enrichment/enrichment_registry.json -> cell_activities[].
 *
 * Each entry describes one per-cell activity matrix in obsm/ (cells x gene sets).
 */
export const CellActivitySchema = z.object({
  obsm_key: z
    .string()
    .min(1)
    .describe("obsm key holding the cells x gene sets activity matrix (for example X_activity_mlm_progeny)."),
  method: z
    .string()
    .min(1)
    .describe("Enrichment algorithm identifier as a free-form string."),
  gene_set_collection: z
    .string()
    .min(1)
    .describe("Gene set collection label used to compute the activities."),
  value_label: z
    .string()
    .min(1)
    .describe("Human-readable label for the activity values (for example MLM t-value)."),
  gene_set_ids: z
    .array(z.string().min(1))
    .describe("Gene set ids in column order of the obsm matrix."),
  n_gene_sets: z
    .number()
    .int()
    .nonnegative()
    .describe("Number of gene sets (columns) in the obsm matrix; equals gene_set_ids length."),
});

/**
 * Top-level schemas for enrichment_registry.json arrays.
 */
export const EnrichmentsArraySchema = z.array(EnrichmentSchema);
export const CellActivitiesArraySchema = z.array(CellActivitySchema);

export type Enrichment = z.infer<typeof EnrichmentSchema>;
export type EnrichmentsArray = z.infer<typeof EnrichmentsArraySchema>;
export type CellActivity = z.infer<typeof CellActivitySchema>;
export type CellActivitiesArray = z.infer<typeof CellActivitiesArraySchema>;
