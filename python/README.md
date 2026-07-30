# controlledburn (Python)

Scanline polygon/line/point rasterization with exact coverage
fractions, O(perimeter) in time and memory, producing sparse output.
Python bindings (pybind11) over the pure C++17 core in ../cpp.

## Usage

```python
import shapely
import controlledburn as cb

r = cb.burn(shapely.box(2.5, 4.5, 6.5, 8.5),
            bounds=(0, 0, 10, 10),   # (xmin, ymin, xmax, ymax), rasterio-style
            shape=(10, 10))          # (nrow, ncol), numpy-style

r.runs    # interior RLE:       (row, col_start, col_end, id)
r.edges   # boundary fractions: (row, col, fraction, id)
r.lines   # line lengths:       (row, col, length, id)
r.points  # point cells:        (row, col, id)

import pandas as pd
pd.DataFrame(r.edges)  # structured arrays convert directly

# Optional dense consumer (fasterize-style pixel functions):
m = cb.materialize(r, fn="sum", edge_policy="fraction")
```

Input can be WKB bytes, shapely geometries (scalar or sequence), a
geopandas GeoSeries/GeoDataFrame, or a mixed sequence with `None`
placeholders. Grid parameters follow the R package's rules when
omitted: `bounds` defaults to the geometry bounding box, `shape`
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

Tables use 1-BASED row/col with row 1 at the top -- identical values to
the R package, so cross-language fixtures compare directly. The
argument conventions are deliberately Python-native and differ from R:
`bounds=(xmin, ymin, xmax, ymax)` vs R's `extent = c(xmin, xmax, ymin,
ymax)`, and `shape=(nrow, ncol)` vs R's `dimension = c(ncol, nrow)`.

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
