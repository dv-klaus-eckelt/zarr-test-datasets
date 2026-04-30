# /// script
# requires-python = ">=3.12"
# dependencies = [
#     "anndata>=0.12.11",
#     "scanpy>=1.12.1",
#     "zarr>=3.1.6",
# ]
# ///

from pathlib import Path

import anndata as ad
import numpy as np
import scanpy as sc


INPUT_PATH = Path("habib17.h5ad")
OUTPUT_PATH = Path("test-data/habib17-differential-expression-test-data.zarr")
GROUPBY_COLUMN = "CellType"


def _sample_expression_values(x: object, max_items: int = 10000) -> np.ndarray:
    if hasattr(x, "tocoo"):
        values = np.asarray(x.data)
    else:
        values = np.asarray(x).ravel()

    if values.size == 0:
        return values

    if values.size > max_items:
        step = max(1, values.size // max_items)
        values = values[::step]

    return values


def _needs_preprocessing(adata: ad.AnnData) -> tuple[bool, str]:
    if "log1p" in adata.uns:
        return False, "Detected adata.uns['log1p']; matrix appears already log-transformed."

    values = _sample_expression_values(adata.X)
    if values.size == 0:
        return True, "Empty matrix sample; applying preprocessing by default."

    tol = 1e-6
    non_integer_fraction = float(np.mean(np.abs(values - np.round(values)) > tol))
    max_value = float(np.max(values))

    # Heuristic: mostly non-integers with compressed range usually indicates logged data.
    if non_integer_fraction > 0.2 and max_value < 50:
        return (
            False,
            (
                "Expression values look already transformed "
                f"(non-integer fraction={non_integer_fraction:.3f}, max={max_value:.3f})."
            ),
        )

    return (
        True,
        (
            "Expression values look like raw counts "
            f"(non-integer fraction={non_integer_fraction:.3f}, max={max_value:.3f})."
        ),
    )


def main() -> None:
    if not INPUT_PATH.exists():
        raise FileNotFoundError(f"Input file not found: {INPUT_PATH}")

    adata = ad.read_h5ad(INPUT_PATH)

    if GROUPBY_COLUMN not in adata.obs.columns:
        available = ", ".join(adata.obs.columns.astype(str).tolist())
        raise KeyError(
            f"Column '{GROUPBY_COLUMN}' not found in obs. Available columns: {available}"
        )

    # Make sure group labels are categorical for rank_genes_groups.
    adata.obs[GROUPBY_COLUMN] = adata.obs[GROUPBY_COLUMN].astype("category")

    should_preprocess, reason = _needs_preprocessing(adata)
    print(reason)

    if should_preprocess:
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        print("Applied preprocessing: normalize_total + log1p")
    else:
        print("Skipped preprocessing")

    sc.tl.rank_genes_groups(
        adata,
        groupby=GROUPBY_COLUMN,
        method="wilcoxon",
        use_raw=False,
    )

    # Align with the zarr v3 output style used in other dataset creation scripts.
    ad.settings.zarr_write_format = 3
    ad.settings.write_csr_csc_indices_with_min_possible_dtype = True
    ad.settings.auto_shard_zarr_v3 = True

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    adata.write_zarr(OUTPUT_PATH)

    print(f"Wrote DE results to {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
