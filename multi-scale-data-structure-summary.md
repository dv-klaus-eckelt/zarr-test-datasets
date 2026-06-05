# Multi-Scale Data Structure in SpatialData Zarr Stores

## Zarr Directory Layout

A SpatialData Zarr v2 store has the following structure:

```
store.zarr/
├── .zgroup                          # Root Zarr group metadata (zarr_format: 2)
├── .zattrs                          # SpatialData version & software provenance
├── zmetadata                        # Consolidated metadata (optional, aggregates all .zarray/.zattrs/.zgroup)
│
├── images/                          # Image elements (multiplexed microscopy, H&E, etc.)
│   └── .zgroup                      # Zarr group marker
│
├── shapes/                          # Shape elements (segmentations, ROIs, annotations)
│   └── .zgroup                      # Zarr group marker (if shapes exist)
│
└── tables/                          # Table elements (AnnData-like observation tables)
    └── .zgroup                      # Zarr group marker (if tables exist)
```

## Images: The Multi-Scale Question

### Structure observed in this dataset

Each image is stored as an **independent Zarr group** with a single array at one resolution:

```
images/
└── <element_name>/
    ├── .zgroup                      # zarr_format: 2
    ├── .zattrs                      # multiscales spec + omero channel metadata
    └── 0/                           # Single dataset: one resolution level only
        ├── .zarray                  # Shape, chunks, dtype, compressor
        └── 0/                       # Chunk coordinate: c=0 (first channel axis)
            └── 0/                   # Chunk coordinate: y=0 (first row axis)
                └── 0                # Chunk data file: x=0 (first column axis)
```

**Key observation:** The `datasets` array in `.zattrs` always contains exactly one entry:

### Image arrays: channels, shape, and chunk layout

Images are stored as **3-dimensional arrays** with axes `[C, Y, X]`:

- **C** — channel axis (first axis, interleaved in each chunk)
- **Y** — spatial axis (image rows / height)
- **X** — spatial axis (image columns / width)

Channels are **not** stored in separate directories or arrays. All channels are interleaved along the first dimension of a single Zarr array. For example, an RGB image has shape `[3, H, W]` and a single chunk contains all three color channels for the covered pixel region.

**Chunk file path anatomy.** With `dimension_separator: "/"`, a chunk at grid coordinate
`[c_idx, y_idx, x_idx]` is stored at `<c_idx>/<y_idx>/<x_idx>` under the dataset
directory. For example, the path components:

```
<element_name>/0/0/0/0
                │ │ │ │
                │ │ │ └─ x=0: chunk data file (last axis == file)
                │ │ └─── y=0: row chunk directory
                │ └───── c=0: channel chunk directory
                └─────── dataset path (from multiscales `datasets[].path`)
```

When the chunk size equals the full array shape (i.e. one chunk covers the entire image),
there is only a single chunk file at `0/0/0/0` containing all channels and all spatial
pixels.

**Channel metadata** is stored in `.zattrs` under the `omero` key, not in the array structure:

```json
{
  "omero": {
    "channels": [
      { "label": "r" },
      { "label": "g" },
      { "label": "b" }
    ]
  }
}
```

Each entry in the `channels` array corresponds by position to the first axis of the Zarr
array — channel index 0 maps to `channels[0]`, index 1 to `channels[1]`, and so on.

```json
{
  "multiscales": [
    {
      "datasets": [
        {
          "path": "0",
          "coordinateTransformations": [
            { "scale": [1.0, 1.0, 1.0], "type": "scale" }
          ]
        }
      ]
    }
  ]
}
```

There is never a `1/`, `2/`, etc. subdirectory. Each `.zattrs` declares only a single `path` in `datasets`.

Multi-resolution imagery is achieved _here_ by creating **separate named image elements** — for example, one called `<name>_hires_image` and another called `<name>_lowres_image` — each containing exactly one resolution level and living in independent Zarr groups with no storage-level link between them.

### Standard OME-NGFF multiscale structure (what this dataset does NOT use)

In a true multiscale pyramid, multiple resolution levels are stored under **the same Zarr group**:

```
images/
└── <element_name>/                  # Single multiscale group
    ├── .zgroup
    ├── .zattrs                      # multiscales spec with multiple datasets
    ├── 0/                           # Level 0: full resolution
    │   ├── .zarray
    │   └── ...chunks...
    ├── 1/                           # Level 1: first downsampled level
    │   ├── .zarray
    │   └── ...chunks...
    └── 2/                           # Level 2: further downsampled (if present)
        ├── .zarray
        └── ...chunks...
```

The `.zattrs` lists **multiple datasets**, each with its own path and scale factor relative to the **base resolution** (not the pixel space):

```json
{
  "multiscales": [
    {
      "datasets": [
        {
          "path": "0",
          "coordinateTransformations": [
            { "scale": [1.0, 1.0, 1.0], "type": "scale" }
          ]
        },
        {
          "path": "1",
          "coordinateTransformations": [
            { "scale": [1.0, 2.0, 2.0], "type": "scale" }
          ]
        },
        {
          "path": "2",
          "coordinateTransformations": [
            { "scale": [1.0, 4.0, 4.0], "type": "scale" }
          ]
        }
      ]
    }
  ]
}
```

### Structural comparison

|                                 | This dataset (old convention)                | True OME-NGFF pyramid                        |
| ------------------------------- | -------------------------------------------- | -------------------------------------------- |
| Zarr groups per resolution      | 1 independent group per level                | 1 group for all levels                       |
| `datasets` entries in `.zattrs` | Always 1                                     | ≥2                                           |
| Subdirectories under group      | Only `0/`                                    | `0/`, `1/`, `2/`, ...                        |
| Multi-scale linking             | Naming convention + shared coordinate system | Explicit `datasets` array in one `.zattrs`   |
| Viewer fallback logic           | SpatialData API maps element names to levels | Zarr reader loads appropriate `path` by zoom |

## Shapes

Shape elements store polygon/multipolygon geometries as **Apache Parquet files**:

```
shapes/
└── <element_name>/
    ├── .zgroup
    ├── .zattrs                      # axes + coordinateTransformations + encoding-type
    └── shapes.parquet               # Geometry data (Parquet format)
```

Shape `.zattrs` structure:

```json
{
  "axes": ["x", "y"],
  "coordinateTransformations": [
    { "input":  { "name": "xy" },     "output": { "name": "<coordinate_system>"      }, "type": "identity" },
    { "input":  { "name": "xy" },     "output": { "name": "<cs>_downscaled_hires" }, "scale": [...], "type": "scale" },
    { "input":  { "name": "xy" },     "output": { "name": "<cs>_downscaled_lowres" }, "scale": [...], "type": "scale" }
  ],
  "encoding-type": "ngff:shapes"
}
```

## Tables

Tables follow the AnnData Zarr layout (anndata v0.10+ spec):

```
tables/
└── <table_name>/
    ├── .zgroup
    ├── X/                            # Expression matrix (sparse/dense)
    ├── obs/                          # Observations (if separate from var)
    ├── var/                          # Variables / gene metadata
    ├── obsm/                         # Multi-dimensional observations (e.g., spatial coords, UMAP)
    ├── obsp/                         # Pairwise observations
    ├── varp/                         # Pairwise variables
    ├── layers/                       # Additional data layers
    ├── raw/                          # Raw (pre-filtered) expression data
    └── uns/                          # Unstructured metadata
```

Tables contain coordinate information in `obsm/spatial/` which maps observations to the same coordinate systems used by images and shapes.

## Coordinate Systems: The Linking Mechanism

### Where coordinate systems are defined

There is **no dedicated `coordinate_systems/` directory or array** in the Zarr store. Coordinate systems exist purely as **string labels** in the `output.name` field within each element's `coordinateTransformations` metadata. They are defined in:

| Element type | Metadata location            |
| ------------ | ---------------------------- |
| Images       | `images/<name>/.zattrs`      |
| Shapes       | `shapes/<name>/.zattrs`      |
| Tables       | `tables/<name>/obsm/.zattrs` |

### How coordinate transformations work

Every spatial element declares a chain of coordinate transformations mapping its **intrinsic (pixel/array) coordinates** to a named **output coordinate system**:

```
intrinsic coords → [transform 1] → [transform 2] → ... → output coordinate system (label)
```

The output coordinate system name is an abstract identifier — it has no physical presence in the Zarr store. It is simply a string that multiple elements can reference.

### How elements are aligned

Elements that target the **same output coordinate system name** are understood to describe the same physical space:

```
Element A (image,  pixel coords) → [scale: 10.5 µm/px] → "tissue_section_X"
Element B (image,  pixel coords) → [scale: 105 µm/px] → "tissue_section_X"
Element C (shapes, polygon coords) → [identity]        → "tissue_section_X"
```

Because A, B, and C all resolve to `"tissue_section_X"`, they occupy the same coordinate space. The scale factors in A and B convert pixel coordinates to physical units (micrometers), ensuring the two images overlap correctly despite having different pixel dimensions.

### Intermediate coordinate systems

The `_downscaled_hires` and `_downscaled_lowres` names in transformation chains represent **intermediate coordinate spaces** — they sit between pixel space and the final tissue coordinate system. They exist because images may be stored at a resolution that is already downscaled relative to the raw acquisition resolution, and the chain of transforms accounts for this.

### Transformation chain example (anatomy of one image)

```
pixel coordinates (intrinsic)
       │
       ▼  [identity transform]
  <name>_downscaled_hires   ← intermediate space
       │
       ▼  [scale transform: 1.0, <µm/px>, <µm/px>]
  <tissue_section_name>     ← final coordinate system (shared label)
```

### Transformation chain example (anatomy of one shape)

```
polygon coordinates (intrinsic)
       │
       ├── [identity transform] → <tissue_section_name>           ← directly to tissue space
       ├── [scale transform]    → <tissue_section_name>_downscaled_hires
       └── [scale transform]    → <tissue_section_name>_downscaled_lowres
```

## How Multi-Scale Resolution Is Achieved (Without a True Pyramid)

Since the Zarr storage layer has no multi-scale pyramid, resolution selection is handled at the **SpatialData API level**:

1. Two separate image elements exist per tissue section — one high-resolution, one low-resolution.
2. Both elements share the same output coordinate system name.
3. The SpatialData library infers the relationship: elements targeting the same coordinate system with different pixel dimensions and scale factors form a de facto multi-scale stack.
4. At viewing time, the library selects the appropriate element based on the requested zoom level / field of view.

This is a **convention over configuration** approach — the Zarr store itself has no multi-scale semantics, but the naming scheme, coordinate system labels, and scale factor relationships allow the consuming library to reconstruct the pyramid logically.

## .zattrs Specification Versions

| Field                       | Value in this dataset              | Meaning                                      |
| --------------------------- | ---------------------------------- | -------------------------------------------- |
| `spatialdata_attrs.version` | `"0.2"` (elements), `"0.1"` (root) | SpatialData element/container spec version   |
| `multiscales[].version`     | `"0.4-dev-spatialdata"`            | Pre-release OME-NGFF multiscale spec variant |
| `zarr_format`               | `2`                                | Zarr v2 storage format                       |
