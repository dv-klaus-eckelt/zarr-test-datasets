# Single Cell & Spatial Test Data

## Serve test data

There are two ways to provide the datasets to Aevidence: via nginx (preferred, through docker) or by spinning up a server in python, see below.

For both, files are served from [localhost:8000](http://localhost:8000/).

Current test data sets in this repository with links to Aevidence:

| **Single Cell** Dataset              | Localhost                                                                                                                                                                                                                                                                            | Public                                                                                                                                                                                                                                                                                                |
| ------------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| zarr v2                              | [Link](http://127.0.0.1:8080/app/single_cell?id=sc-zarr-v2&name=Single+Cell+%28zarr+v2+format%29&type=single-cell&url=http%3A%2F%2Flocalhost%3A8000%2Fsingle-cell-v2-test-data.zarr%2F)                                                                                              | [Link](https://aevidence.dev.app.datavisyn.io/app/single_cell?id=sc-zarr-v2&name=Single+Cell+%28zarr+v2+format%29&type=single-cell&url=http%3A%2F%2Flocalhost%3A8000%2Fsingle-cell-v2-test-data.zarr%2F)                                                                                              |
| zarr v3                              | [Link](http://127.0.0.1:8080/app/single_cell?id=sc-zarr-v3&name=Single+Cell+%28zarr+v3+format%29&type=single-cell&url=http%3A%2F%2Flocalhost%3A8000%2Fsingle-cell-v3-test-data.zarr%2F)                                                                                              | [Link](https://aevidence.dev.app.datavisyn.io/app/single_cell?id=sc-zarr-v3&name=Single+Cell+%28zarr+v3+format%29&type=single-cell&url=http%3A%2F%2Flocalhost%3A8000%2Fsingle-cell-v3-test-data.zarr%2F)                                                                                              |
| zarr v3 with differential expression | [Link](http://127.0.0.1:8080/app/single_cell?id=46652b4c&name=Habib+et+al.+Differential+Expression+Test&primaryAttribute=CellType&storeType=simple&type=single-cell&url=http%3A%2F%2Flocalhost%3A8000%2Fhabib17-differential-expression-test-data-format.zarr%2F&watermark=46652b4c) | [Link](https://aevidence.dev.app.datavisyn.io/app/single_cell?id=46652b4c&name=Habib+et+al.+Differential+Expression+Test&primaryAttribute=CellType&storeType=simple&type=single-cell&url=http%3A%2F%2Flocalhost%3A8000%2Fhabib17-differential-expression-test-data-format.zarr%2F&watermark=46652b4c) |
| zarr v3 with gene set enrichment     | [Link](http://127.0.0.1:8080/app/single_cell?id=7e9d1a3c&name=Habib+et+al.+Gene+Set+Enrichment+Test&primaryAttribute=CellType&storeType=simple&type=single-cell&url=http%3A%2F%2Flocalhost%3A8000%2Fhabib17-enrichment-test-data-format.zarr%2F&watermark=7e9d1a3c)                                                             | [Link](https://aevidence.dev.app.datavisyn.io/app/single_cell?id=7e9d1a3c&name=Habib+et+al.+Gene+Set+Enrichment+Test&primaryAttribute=CellType&storeType=simple&type=single-cell&url=http%3A%2F%2Flocalhost%3A8000%2Fhabib17-enrichment-test-data-format.zarr%2F&watermark=7e9d1a3c)                                                             |

| **Spatial** Dataset | Localhost                                                                                                                                                                        | Public                                                                                                                                                                                            |
| ------------------- | -------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| zarr v2             | [Link](http://127.0.0.1:8080/app/single_cell?id=spatial-zarr-v2&name=Spatial+%28zarr+v2+format%29&type=spatial&url=http%3A%2F%2Flocalhost%3A8000%2Fspatial-v2-test-data.zarr%2F) | [Link](https://aevidence.dev.app.datavisyn.io/app/single_cell?id=spatial-zarr-v2&name=Spatial+%28zarr+v2+format%29&type=spatial&url=http%3A%2F%2Flocalhost%3A8000%2Fspatial-v2-test-data.zarr%2F) |
| zarr v3             | [Link](http://127.0.0.1:8080/app/single_cell?id=spatial-zarr-v3&name=Spatial+%28zarr+v3+format%29&type=spatial&url=http%3A%2F%2Flocalhost%3A8000%2Fspatial-v3-test-data.zarr%2F) | [Link](https://aevidence.dev.app.datavisyn.io/app/single_cell?id=spatial-zarr-v3&name=Spatial+%28zarr+v3+format%29&type=spatial&url=http%3A%2F%2Flocalhost%3A8000%2Fspatial-v3-test-data.zarr%2F) |

ℹ️ If you want to open datasets in a custom Aevidence deployment, use the link generator page at [test-data/list.html](test-data/list.html).

### nginx/Docker

Comes with native support for range requests and concurrent requests:

```
docker compose up
```

### Python server

Run the local server (no .venv/dependencies needed):

```sh
python -m local_file_server ./test-data
```

## Create new test data

Install dependencies and create a virtual environment with [uv](https://docs.astral.sh/uv/):

**Dependencies:**

```
uv sync
source .venv/bin/activate
```

## Notes

Related project: github.com/BIMSBBioinfo/dummy-spatialdata

Using [spatialdata 0.7.3a1](https://github.com/scverse/spatialdata/releases/tag/v0.7.3a1) to support chunks and [ome-zarr 0.16.0](https://github.com/ome/ome-zarr-py/releases/tag/v0.16.0) which introduces sharding.
Pinning dask as there is [no version supported by both spatialdata and ome-zarr](https://github.com/scverse/spatialdata/issues/1059#issuecomment-4352434612).
