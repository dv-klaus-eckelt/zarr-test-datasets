import { z } from "zod";

/**
 * Contrast entry schema for uns/de/contrast_registry.json -> contrasts[].
 *
 * Frontend intent:
 * - Discover available DE comparisons.
 * - Configure chart labels and axis bounds without loading gene-level arrays.
 * - Decide which chart types are valid for a given DE method.
 */
const TestTypeSchema = z.enum(["one_vs_rest", "pairwise", "multiclass"]);

const CommonContrastFieldsSchema = z.object({
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
  test_type: TestTypeSchema.describe(
    "Frontend category of comparison semantics: one_vs_rest, pairwise, or multiclass.",
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
  x_axis_label: z
    .string()
    .min(1)
    .describe(
      'Human-readable x-axis label emitted by the generator (currently "log2(Fold Change)" for all methods).',
    ),
  n_features: z
    .number()
    .int()
    .nonnegative()
    .describe(
      "Total number of ranked features available in the corresponding de_xxx folder arrays.",
    ),
  top_10_genes: z
    .array(z.string().min(1))
    .describe("Top 10 gene symbols ranked by score in descending order."),
});

/**
 * Non-logreg contrasts (wilcoxon and t-test) expose p-value based metrics.
 * These entries include y-axis labels and numeric chart bounds.
 */
const NonLogregContrastSchema = CommonContrastFieldsSchema.extend({
  de_method: z
    .enum(["wilcoxon", "t-test"])
    .describe("DE algorithm used to produce this contrast."),
  correction_method: z
    .string()
    .min(1)
    .describe(
      "P-value adjustment method used for significance (for example benjamini-hochberg or bonferroni).",
    ),
  y_axis_label: z
    .string()
    .min(1)
    .describe(
      "Y-axis label for significance scale, usually -log10(FDR) or -log10(p).",
    ),
  lfc_max: z
    .number()
    .nonnegative()
    .describe(
      "Maximum absolute effect-size bound used to keep volcano x-axis symmetric and stable across contrasts.",
    ),
  logp_max: z
    .number()
    .nonnegative()
    .describe(
      "Maximum -log10(p) bound used to stabilize volcano y-axis height.",
    ),
  available_plots: z
    .array(z.enum(["volcano", "dotplot"]))
    .nonempty()
    .describe(
      "Plot types enabled for this contrast; non-logreg supports volcano and dotplot.",
    ),
});

/**
 * Multiclass logreg contrasts do not expose p-value metrics in registry metadata.
 * Therefore correction_method, y_axis_label, lfc_max, and logp_max are null.
 * Note: per-gene arrays still exist on disk with NaN placeholders for unsupported metrics.
 */
const LogregContrastSchema = CommonContrastFieldsSchema.extend({
  de_method: z
    .literal("logreg")
    .describe("Multiclass logistic-regression output emitted per class."),
  correction_method: z
    .null()
    .describe(
      "Null because logreg output here does not include p-value correction metadata.",
    ),
  y_axis_label: z
    .null()
    .describe(
      "Null because volcano-style significance y-axis is not used for multiclass logreg output.",
    ),
  lfc_max: z
    .null()
    .describe(
      "Null because fold-change style x-axis bounds are not defined for this logreg representation.",
    ),
  logp_max: z
    .null()
    .describe(
      "Null because p-value significance bounds are not produced for this logreg representation.",
    ),
  available_plots: z
    .array(z.enum(["bar_chart", "dotplot"]))
    .nonempty()
    .describe(
      "Plot types enabled for this contrast; logreg supports bar_chart and dotplot.",
    ),
});

/**
 * Discriminated union guarantees method-specific nullability at type level.
 */
export const ContrastSchema = z.discriminatedUnion("de_method", [
  NonLogregContrastSchema,
  LogregContrastSchema,
]);

/**
 * Top-level schema for contrast_registry.json -> contrasts.
 */
export const ContrastsArraySchema = z.array(ContrastSchema);

export type Contrast = z.infer<typeof ContrastSchema>;
export type ContrastsArray = z.infer<typeof ContrastsArraySchema>;
