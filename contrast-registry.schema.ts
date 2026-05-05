import { z } from "zod";

const SharedContrastFields = {
  contrast_id: z
    .string()
    .min(1)
    .describe(
      "Stable contrast identifier and folder key (for example de_001).",
    ),
  group_1: z
    .string()
    .min(1)
    .describe(
      "Target group label for this contrast (shown in UI selectors/titles).",
    ),
  group_2: z
    .string()
    .nullable()
    .describe(
      "Reference group label. Null for one-vs-rest and current multiclass outputs where reference is implicit.",
    ),
  test_type: z
    .string()
    .min(1)
    .describe(
      "Frontend grouping label for comparison semantics (for example one_vs_rest, pairwise, multiclass, or custom labels).",
    ),
  subset_column: z
    .string()
    .nullable()
    .describe(
      "Optional obs column used to restrict the test domain (for example region). Null when no subset restriction was applied.",
    ),
  subset_value: z
    .string()
    .nullable()
    .describe(
      "Optional subset value within subset_column (for example Hippocampus). Null when test was run on the full dataset.",
    ),
  contrast_column: z
    .string()
    .min(1)
    .describe(
      "Obs column that defines the compared groups (for example CellType).",
    ),
  feature_type: z
    .string()
    .min(1)
    .describe(
      "Feature modality label for UI display/filtering (for example Gene Expression).",
    ),
  de_method: z
    .string()
    .min(1)
    .describe("DE algorithm identifier as a free-form string."),
  has_effect_size: z
    .boolean()
    .describe("Whether effect-size metrics are available for this contrast."),
  has_significance: z
    .boolean()
    .describe("Whether significance metrics are available for this contrast."),
  has_pct_expr_target: z
    .boolean()
    .describe(
      "Whether target-group percent-expression values are available for this contrast.",
    ),
  has_pct_expr_ref: z
    .boolean()
    .describe(
      "Whether reference/rest percent-expression values are available for this contrast.",
    ),
  top_10_gene_ids: z
    .array(z.string().min(1))
    .describe("Top 10 feature IDs ranked by score in descending order."),
};

const EffectAndSignificanceContrastSchema = z.object({
  ...SharedContrastFields,
  has_effect_size: z.literal(true),
  has_significance: z.literal(true),
  correction_method: z
    .string()
    .min(1)
    .describe(
      "P-value correction method when significance metrics are available.",
    ),
  effect_size_label: z
    .string()
    .min(1)
    .describe(
      "Human-readable effect-size label when effect-size metrics are available.",
    ),
  significance_label: z
    .string()
    .min(1)
    .describe(
      "Significance metric label when significance metrics are available.",
    ),
  effect_size_max: z
    .number()
    .nonnegative()
    .describe(
      "Absolute effect-size bound when effect-size metrics are available.",
    ),
  significance_max: z
    .number()
    .nonnegative()
    .describe("Significance bound when p-value metrics are available."),
});

const EffectOnlyContrastSchema = z.object({
  ...SharedContrastFields,
  has_effect_size: z.literal(true),
  has_significance: z.literal(false),
  correction_method: z
    .null()
    .describe(
      "P-value correction method must be null when significance metrics are unavailable.",
    ),
  effect_size_label: z
    .string()
    .min(1)
    .describe(
      "Human-readable effect-size label when effect-size metrics are available.",
    ),
  significance_label: z
    .null()
    .describe(
      "Significance metric label must be null when significance metrics are unavailable.",
    ),
  effect_size_max: z
    .number()
    .nonnegative()
    .describe(
      "Absolute effect-size bound when effect-size metrics are available.",
    ),
  significance_max: z
    .null()
    .describe(
      "Significance bound must be null when p-value metrics are unavailable.",
    ),
});

const SignificanceOnlyContrastSchema = z.object({
  ...SharedContrastFields,
  has_effect_size: z.literal(false),
  has_significance: z.literal(true),
  correction_method: z
    .string()
    .min(1)
    .describe(
      "P-value correction method when significance metrics are available.",
    ),
  effect_size_label: z
    .null()
    .describe(
      "Effect-size label must be null when effect-size metrics are unavailable.",
    ),
  significance_label: z
    .string()
    .min(1)
    .describe(
      "Significance metric label when significance metrics are available.",
    ),
  effect_size_max: z
    .null()
    .describe(
      "Effect-size bound must be null when effect-size metrics are unavailable.",
    ),
  significance_max: z
    .number()
    .nonnegative()
    .describe("Significance bound when p-value metrics are available."),
});

const ScoreOnlyContrastSchema = z.object({
  ...SharedContrastFields,
  has_effect_size: z.literal(false),
  has_significance: z.literal(false),
  correction_method: z
    .null()
    .describe(
      "P-value correction method must be null when significance metrics are unavailable.",
    ),
  effect_size_label: z
    .null()
    .describe(
      "Effect-size label must be null when effect-size metrics are unavailable.",
    ),
  significance_label: z
    .null()
    .describe(
      "Significance metric label must be null when significance metrics are unavailable.",
    ),
  effect_size_max: z
    .null()
    .describe(
      "Effect-size bound must be null when effect-size metrics are unavailable.",
    ),
  significance_max: z
    .null()
    .describe(
      "Significance bound must be null when p-value metrics are unavailable.",
    ),
});

/**
 * Method-agnostic contrast schema.
 *
 * Contrast entry schema for uns/de/contrast_registry.json -> contrasts[].
 *
 * Frontend intent:
 * - Discover available DE comparisons.
 * - Configure chart labels and metric bounds without loading gene-level arrays.
 * - Infer optional visualizations from available metrics and arrays.
 *
 * This intentionally does not discriminate on `de_method` so collaborators can
 * provide outputs from any DE engine (for example Scanpy, DESeq2, edgeR, MAST,
 * or custom pipelines) without schema updates.
 */
export const ContrastSchema = z.union([
  EffectAndSignificanceContrastSchema,
  EffectOnlyContrastSchema,
  SignificanceOnlyContrastSchema,
  ScoreOnlyContrastSchema,
]);

/**
 * Top-level schema for contrast_registry.json -> contrasts.
 */
export const ContrastsArraySchema = z.array(ContrastSchema);

export type Contrast = z.infer<typeof ContrastSchema>;
export type ContrastsArray = z.infer<typeof ContrastsArraySchema>;
