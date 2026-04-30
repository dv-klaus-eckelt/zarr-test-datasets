# Single Cell & Spatial Test Data

Related project: github.com/BIMSBBioinfo/dummy-spatialdata

## Create new testdata

To work with different versions of [anndata](https://pypi.org/project/anndata/) and [zarr](https://pypi.org/project/zarr/) we use [juv](https://github.com/manzt/juv) to have notebooks with their individual dependencies.

```
# Create the notebook with python 3.12
juv init --python=3.12 create-single-cell-data.ipynb

# add dependencies
juv add create-single-cell-data.ipynb anndata scanpy zarr spatialdata

# timestamp dependencies
juv stamp create-single-cell-data.ipynb
```

Run:

```
juv run create-single-cell-data.ipynb
```
