# Single Cell & Spatial Test Data

Related project: github.com/BIMSBBioinfo/dummy-spatialdata

## Create new test data

To work with different versions of [anndata](https://pypi.org/project/anndata/) and [zarr](https://pypi.org/project/zarr/) we use [juv](https://github.com/manzt/juv) to have notebooks with their individual dependencies.

```sh
# Create the notebook with python 3.12
juv init --python=3.12 create-single-cell-data.ipynb

# add dependencies
juv add create-single-cell-data.ipynb anndata scanpy zarr spatialdata

# timestamp dependencies
juv stamp create-single-cell-data.ipynb
```

Run:

```sh
juv run create-single-cell-data.ipynb
```

## Serve test data

Run the local sevrer:

```sh
python -m local_file_server ./test-data
```

Files are served from [localhost:8000](http://localhost:8000/).

Current test data sets in this repository and links for aevidence running locally:

**Single Cell**

- [Single Cell (zarr v2 format)](http://127.0.0.1:8080/app/single_cell?id=sc-zarr-v2&name=Single+Cell+%28zarr+v2+format%29&type=single-cell&url=http%3A%2F%2Flocalhost%3A8000%2Fsingle-cell-v2-test-data.zarr%2F)
- [Single Cell (zarr v3 format)](http://127.0.0.1:8080/app/single_cell?id=sc-zarr-v3&name=Single+Cell+%28zarr+v3+format%29&type=single-cell&url=http%3A%2F%2Flocalhost%3A8000%2Fsingle-cell-v3-test-data.zarr%2F)

**Spatial**

- [Spatial (zarr v2 format)](http://127.0.0.1:8080/app/single_cell?id=spatial-zarr-v2&name=Spatial+%28zarr+v2+format%29&type=spatial&url=http%3A%2F%2Flocalhost%3A8000%2Fspatial-v2-test-data.zarr%2F)
- [Spatial (zarr v3 format)](http://127.0.0.1:8080/app/single_cell?id=spatial-zarr-v3&name=Spatial+%28zarr+v3+format%29&type=spatial&url=http%3A%2F%2Flocalhost%3A8000%2Fspatial-v3-test-data.zarr%2F)
