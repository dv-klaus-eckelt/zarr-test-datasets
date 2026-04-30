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
import scanpy as sc


INPUT_PATH = Path("habib17.h5ad")
OUTPUT_PATH = Path("test-data/habib17-differential-expression-test-data.zarr")
GROUPBY_COLUMN = "CellType"


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

    # Standard preprocessing for differential expression.
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)

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
