# controlledburn (Python)

Scanline polygon/line/point rasterization with exact coverage
fractions, O(perimeter) in time and memory, producing sparse output.
Python bindings (pybind11) over the pure C++17 core in ../cpp.

## Usage

```python
import shapely
import controlledburn as cb
poly = shapely.box(2.5, 4.5, 6.5, 8.5)
r = cb.burn([shapely.to_wkb(poly)],
            extent=(0, 10, 0, 10),   # (xmin, xmax, ymin, ymax), matches R's extent
            shape=(10, 10))          # (nrow, ncol), numpy-style

r.runs    # interior RLE:       (row, col_start, col_end, id)
r.edges   # boundary fractions: (row, col, fraction, id)
r.lines   # line lengths:       (row, col, length, id)
r.points  # point cells:        (row, col, id)

import pandas as pd
pd.DataFrame(r.edges)  # structured arrays convert directly

# Optional dense consumer (fasterize-style pixel functions):
m = cb.materialize(r, fn="sum", edge_policy="fraction")

# Crop the sparse result to a sub-window and materialize just that tile
# (the counterpart of R's crop_burn() + materialise_chunk()):
tile = r.crop((3, 7, 5, 9)).materialize(fn="sum", edge_policy="fraction")
```

Input can be WKB bytes, shapely geometries (scalar or sequence), a
geopandas GeoSeries/GeoDataFrame, or a mixed sequence with `None`
placeholders. Grid parameters follow the R package's rules when
omitted: `extent` defaults to the geometry bounding box, `shape`
defaults to a 256-cell fit along the longer axis, and `resolution=`
(mutually exclusive with `shape`) computes the shape from a cell size:

```python
r = cb.burn(geoms)                    # bbox extent, 256-cell fit
r = cb.burn(geoms, resolution=0.01)   # cell size in CRS units
```

`materialize()` is the polygon consumer (runs + edges, the fasterize
path). Line, point, and mixed-kind results raise rather than silently
returning an empty array: what a dense line or point raster should
mean is an aggregation question deliberately left unresolved (see the
tracking issue). The sparse tables are the product; dense consumption
of lines and points is a few lines of numpy over them.

Tables use 0-BASED row/col with row 0 at the top and an exclusive
`col_end` (a run covers `[col_start, col_end)`); geometry k gets id = k.
The R package's shim adds 1 to row/col/col_start/id, so its tables are
1-based and inclusive -- the same coverage values, a different index
base. `extent` uses the same ordering as R's
`extent = c(xmin, xmax, ymin, ymax)`, while `shape=(nrow, ncol)` is the
numpy transpose of R's `dimension = c(ncol, nrow)`.

Non-fatal problems (unparseable WKB, GeometryCollection input) are
raised as Python warnings and the offending geometry is skipped.
GeometryCollections must be split into homogeneous groups upstream;
curved types must be linearised.

## Install (development)

From the repository root:

```sh
pip install ./python
```

Requires a C++17 compiler and CMake >= 3.16 (via scikit-build-core).
The build references the core at ../cpp, so it works for git checkouts
and `pip install git+...#subdirectory=python`; publishing an sdist to
PyPI would need the core vendored into this tree first (a sync step
like the R side's tools/sync-core.sh) -- deliberately deferred until a
first release.

## Test

```sh
pip install ./python[test]
pytest python/tests
```
