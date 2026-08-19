"""Create a zarr v3 store whose index columns use AnnData's `nullable-string-array` encoding.

This is a deliberate A/B partner for `single-cell-v3-test-data.zarr`: same cells, same
genes, same seed, same values. The only difference is the pandas dtype of the index and
the string columns.

- `object` dtype  -> AnnData writes a plain **array** with `encoding-type: string-array`
- `string` dtype  -> AnnData writes a **group** with `encoding-type: nullable-string-array`
                     holding `values` + `mask` children

Newer producers (anndata >= 0.11 with `allow_write_nullable_strings`, and anything built on
pandas nullable/Arrow-backed strings) emit the second form. Readers that assume the index is
a flat array cannot open it.

Run:
    python create_nullable_string_test_data.py
"""

import shutil

import anndata as ad
import numpy as np
import pandas as pd
import zarr
from scipy.sparse import csc_matrix

# Same seed as create-single-cell-data.ipynb so the two stores hold identical values.
np.random.seed(1)

ad.settings.zarr_write_format = 3
ad.settings.write_csr_csc_indices_with_min_possible_dtype = True
ad.settings.auto_shard_zarr_v3 = True
# Opt in to pd.arrays.StringArray output. Without this AnnData refuses to write the
# nullable form at all, which is exactly why the encoding is easy to miss in testing.
ad.settings.allow_write_nullable_strings = True

# The dtype that triggers the nullable encoding. Swap to "string[pyarrow]" for the
# Arrow-backed variant; both produce the same on-disk `values` + `mask` group.
STRING_DTYPE = "string"

gene_map = {
    "ACTB": "ENSG00000075624",  # Housekeeping
    "CD3E": "ENSG00000198851",  # T-cell marker
    "CD79A": "ENSG00000105369",  # B-cell marker
    "GAPDH": "ENSG00000111640",  # Housekeeping
    "GNLY": "ENSG00000115523",  # NK-cell marker
    "LYZ": "ENSG00000090382",  # Monocyte marker
}

gene_symbols = list(gene_map.keys())
ensembl_ids = list(gene_map.values())

n_cells = 100
n_genes = len(gene_symbols)

# Observation metadata (cells). Index is nullable string, the columns stay plain object
# so they still round-trip through the normal categorical path — the index is the only
# thing under test.
obs = pd.DataFrame(index=pd.Index(pd.array([f"Cell_{i:03d}" for i in range(n_cells)], dtype=STRING_DTYPE)))
obs["cell_type"] = np.random.choice(["T-cell", "B-cell", "Monocyte"], size=n_cells).astype(object)
obs["treatment"] = np.random.choice(["treated", "untreated"], size=n_cells).astype(object)

# Variable metadata (genes). Index and the two string columns are nullable strings.
var = pd.DataFrame(index=pd.Index(pd.array(gene_symbols, dtype=STRING_DTYPE)))
var["symbol"] = pd.array(gene_symbols, dtype=STRING_DTYPE)
var["ensembl_gene_id"] = pd.array(ensembl_ids, dtype=STRING_DTYPE)
# Single category on purpose: all codes equal fill_value, so zarr writes no chunk file
# and the chunk request 404s. That is legal zarr and a reader must treat it as fill_value.
var["feature_type"] = np.array(["Random Counts"] * n_genes, dtype=object)

counts_dense = np.random.poisson(lam=2.0, size=(n_cells, n_genes)).astype(np.float32)
adata = ad.AnnData(X=csc_matrix(counts_dense), obs=obs, var=var)

adata.obsm["X_umap"] = np.random.normal(size=(n_cells, 2))
adata.obsm["tsne"] = np.random.normal(size=(n_cells, 2))

# The app locates the identifier columns through these `_index` names.
adata.obs.index.name = "cell_id"
adata.var.index.name = "gene_symbol"

path = "test-data/single-cell-nullable-string-v3-test-data.zarr"
shutil.rmtree(path, ignore_errors=True)

adata.write_zarr(path, chunks=(n_cells, n_genes))
zarr.consolidate_metadata(path)

print(f"Wrote {path}")
