# Using controlledburn from Python

controlledburn rasterizes geometry *without allocating a raster*. It walks
the boundary of each shape and emits a compact description of which grid
cells were touched and how much. Nothing is ever materialized unless you
ask for it.

The practical consequence: you can burn 25,000 polygons onto an eight
trillion cell grid in a few seconds, because the cost is proportional to
the perimeter of your geometry, not to the area of your raster.

This guide covers the Python bindings. The C++ core, the R package and
the Python package all share one engine and one output contract, so
values agree exactly across languages.

---

## Contents

1. [Install](#install)
2. [The thirty second version](#the-thirty-second-version)
3. [The four tables](#the-four-tables)
4. [Specifying the grid](#specifying-the-grid)
5. [Getting geometry in](#getting-geometry-in)
6. [Coverage mode and approx mode](#coverage-mode-and-approx-mode)
7. [Polygons](#polygons)
8. [Lines](#lines)
9. [Points](#points)
10. [Materializing to a dense array](#materializing-to-a-dense-array)
11. [Cropping and tiling](#cropping-and-tiling)
12. [Georeferencing the output](#georeferencing-the-output)
13. [pandas, scipy.sparse, xarray](#pandas-scipysparse-xarray)
14. [Scale](#scale)
15. [Sharp edges](#sharp-edges)
16. [R to Python cheatsheet](#r-to-python-cheatsheet)

---

## Install

From a git checkout:

```sh
pip install ./python
```

Needs a C++17 compiler and CMake 3.16 or newer (supplied through
scikit-build-core). The build reads the core from `../cpp`, so installing
straight from GitHub also works:

```sh
pip install "git+https://github.com/hypertidy/controlledburn#subdirectory=python"
```

The only hard runtime dependency is numpy. shapely is used throughout
this guide to make geometry, but controlledburn itself only ever sees
WKB bytes, so anything that can produce WKB will do.

---

## The thirty second version

```python
import shapely
import controlledburn as cb

poly = shapely.box(2.5, 4.5, 6.5, 8.5)

r = cb.burn([shapely.to_wkb(poly)],
            extent=(0, 10, 0, 10),   # (xmin, xmax, ymin, ymax), matches R's extent
            shape=(10, 10))          # (nrow, ncol)

r.runs     # interior cells, run-length encoded
r.edges    # boundary cells with exact coverage fractions
r.lines    # line cells with length-in-cell
r.points   # point cells
```

For this box:

```
r.runs
array([(3, 4, 6, 1), (4, 4, 6, 1), (5, 4, 6, 1)],
      dtype=[('row','<i4'), ('col_start','<i4'), ('col_end','<i4'), ('id','<i4')])

r.edges[:5]
array([(2, 3, 0.25, 1), (2, 4, 0.5, 1), (2, 5, 0.5, 1),
       (2, 6, 0.5, 1), (2, 7, 0.25, 1)],
      dtype=[('row','<i4'), ('col','<i4'), ('fraction','<f4'), ('id','<i4')])
```

Three runs of three fully covered cells each, plus sixteen boundary cells.
Add it up: 9 full cells plus 4 quarters plus 12 halves is 16 cell-areas,
and each cell is 1 x 1, so 16.0 square units. The box is 4 x 4. The
coverage is exact, not sampled.

---

## The four tables

`burn()` returns a `BurnResult` holding four numpy structured arrays.
There is one table per geometry kind, and each table's measure column
means exactly one thing.

| Table | Columns | Measure | Produced by |
|---|---|---|---|
| `runs` | `row, col_start, col_end, id` | none (full coverage) | polygons |
| `edges` | `row, col, fraction, id` | dimensionless, in (0, 1) | polygons |
| `lines` | `row, col, length, id` | CRS units | lines |
| `points` | `row, col, id` | none | points |

Plus `r.extent` and `r.shape`, echoing what you passed in.

Three rules that hold everywhere:

- **Indices are 1-based.** `row` 1 is the first row, `col` 1 is the first
  column. This is a deliberate cross-language contract, not an oversight.
  Subtract 1 to index a numpy array.
- **Row 1 is the top of the grid**, at `ymax`. Same orientation as a
  raster, opposite to a matplotlib default axis.
- **`id` is the 1-based position of the geometry in your input.** Geometry
  `k` (0-based) gets `id = k + 1`.

The separation between tables is the point. A coverage fraction is
dimensionless, a line length is in CRS units, and a point has no measure
at all. Stacking them into one column would silently mix units, so the
package refuses to do it for you.

Mixed input is fine, and populates whichever tables apply:

```python
mixed = cb.burn([shapely.to_wkb(g) for g in [
    shapely.box(2.5, 4.5, 6.5, 8.5),
    shapely.LineString([(0.5, 0.5), (9.5, 7.5)]),
    shapely.Point(1.5, 1.5),
]], extent=(0, 10, 0, 10), shape=(10, 10))

len(mixed.runs), len(mixed.edges), len(mixed.lines), len(mixed.points)
# (3, 16, 16, 1)
```

---

## Specifying the grid

Two arguments, both required, both in Python-native order:

```python
extent = (xmin, xmax, ymin, ymax)   # matches R's extent ordering
shape  = (nrow, ncol)               # numpy ordering
```

`extent` uses the same ordering as the R package's
`extent = c(xmin, xmax, ymin, ymax)`. `shape` follows the numpy idiom
`(nrow, ncol)`, which is the transpose of the R package's
`dimension = c(ncol, nrow)`. The *values* in the output tables are
identical either way.

`shape` is the easiest thing to get backwards. If your result looks
transposed or oddly clipped, check that argument first. Note that ordering
mistakes are silent -- a square grid will not complain at all.

Invalid grids do raise:

```python
cb.burn(wkb, extent=(0, 0, 0, 10), shape=(10, 10))
# ValueError: invalid extent: xmax must be > xmin, ymax must be > ymin

cb.burn(wkb, extent=(0, 10, 0, 10), shape=(0, 10))
# ValueError: ncol and nrow must be positive
```

### Deriving a grid from the geometry

The R package lets you skip both arguments, or give a `resolution`
instead of a dimension. Python does not do this yet (see
[Sharp edges](#sharp-edges)), but it is a few lines:

```python
import numpy as np
import shapely

def resolve_grid(geoms, extent=None, shape=None, resolution=None, size=256):
    """Mirror of the R package's grid parameter resolution."""
    if extent is None:
        g = as_wkb_list(geoms)                 # see the next section
        b = np.atleast_2d(shapely.bounds(shapely.from_wkb(g)))
        extent = (float(np.nanmin(b[:, 0])), float(np.nanmax(b[:, 2])),
                  float(np.nanmin(b[:, 1])), float(np.nanmax(b[:, 3])))
    xmin, xmax, ymin, ymax = (float(v) for v in extent)

    if shape is not None and resolution is not None:
        raise ValueError("specify 'shape' or 'resolution', not both")

    if resolution is not None:
        res = np.broadcast_to(np.asarray(resolution, float), (2,))
        shape = (int(np.ceil((ymax - ymin) / res[1])),
                 int(np.ceil((xmax - xmin) / res[0])))

    if shape is None:
        w, h = xmax - xmin, ymax - ymin
        m = max(w, h)
        shape = (int(np.ceil(size * h / m)), int(np.ceil(size * w / m)))

    return extent, (int(shape[0]), int(shape[1]))
```

```python
resolve_grid(shapely.box(2.5, 4.5, 6.5, 8.5))
# ((2.5, 6.5, 4.5, 8.5), (256, 256))

resolve_grid(shapely.box(2.5, 4.5, 6.5, 8.5), resolution=0.5)
# ((2.5, 6.5, 4.5, 8.5), (8, 8))
```

The default fits at most 256 cells along the longer axis and preserves
aspect ratio, which is the same rule the R package uses.

---

## Getting geometry in

`burn()` wants a sequence of WKB blobs. Both ISO WKB and EWKB are
accepted, in either byte order, and Z or M ordinates are skipped rather
than rejected. `None` entries are skipped (and still consume an id).

**shapely, one at a time:**

```python
cb.burn([shapely.to_wkb(g) for g in geoms], extent=..., shape=...)
```

**shapely, vectorized** (much faster for large collections; `to_wkb` on an
array returns an object array of bytes, which `burn()` accepts directly):

```python
import numpy as np
arr = np.asarray(geoms)
cb.burn(shapely.to_wkb(arr), extent=..., shape=...)
```

**geopandas:**

```python
xmin, ymin, xmax, ymax = gdf.total_bounds
cb.burn(gdf.geometry.to_wkb(), extent=(xmin, xmax, ymin, ymax), shape=(2000, 2000))
```

`gdf.total_bounds` is `(xmin, ymin, xmax, ymax)` (rasterio order), so it
needs reordering to controlledburn's `(xmin, xmax, ymin, ymax)` extent.

**pyogrio, skipping shapely entirely:**

```python
import pyogrio
tbl = pyogrio.read_arrow("countries.gpkg")[1]
wkb = tbl.column("wkb_geometry").to_pylist()
cb.burn(wkb, extent=(-180, 180, -90, 90), shape=(1800, 3600))
```

**Anything else** that can hand you bytes: a WKB column out of DuckDB or
PostGIS, `fiona`, `gdal`, a parquet file of WKB blobs. controlledburn has
no geometry model of its own and no GEOS dependency, so there is nothing
to convert to.

### A coercion helper

If you would rather pass shapely objects around directly, this small
wrapper mirrors the R package's `as_wkb_list()` generic:

```python
def as_wkb_list(x):
    if isinstance(x, (bytes, bytearray, memoryview)):
        return [x]
    if hasattr(x, "geom_type"):                            # one shapely geom
        return [shapely.to_wkb(x)]
    if hasattr(x, "geometry") and hasattr(x, "columns"):   # GeoDataFrame
        x = x.geometry
    if hasattr(x, "to_wkb"):                               # GeoSeries
        return list(x.to_wkb())
    out = []
    for g in x:
        if g is None or isinstance(g, (bytes, bytearray, memoryview)):
            out.append(g)
        elif hasattr(g, "geom_type"):
            out.append(shapely.to_wkb(g))
        else:
            raise TypeError(f"cannot interpret {type(g).__name__} as geometry")
    return out
```

### What is not accepted

- **GeometryCollection** warns and is skipped. Split into homogeneous
  groups first.
- **Curved types** (CircularString, CompoundCurve and friends) must be
  linearised by the caller.
- **Unparseable bytes** warn and are skipped; the remaining geometries
  still burn, and ids are unaffected.

```python
import warnings
with warnings.catch_warnings(record=True) as w:
    warnings.simplefilter("always")
    r = cb.burn([b"\x01\x03\x00", shapely.to_wkb(shapely.box(2, 4, 6, 8))],
                extent=(0, 10, 0, 10), shape=(10, 10))
# UserWarning: geometry 1: ... WKB ...
# the second geometry is still burned, with id 2
```

### What is not validated

controlledburn rasterizes whatever it is given. Self-intersecting rings,
unclosed polygons, repeated vertices, near-degenerate slivers -- all go
through the same fast path with no validity check. If your input has
topological problems that matter for your science you will see them in
the output. If they do not matter, you have saved the cost of checking.
Use `shapely.is_valid` and `shapely.make_valid` when you do want that.

---

## Coverage mode and approx mode

```python
cb.burn(wkb, extent=..., shape=..., mode="coverage")   # default
cb.burn(wkb, extent=..., shape=..., mode="approx")
```

**`coverage`** computes exact analytical coverage fractions for polygon
boundary cells, using the exactextract algorithm. Every partially covered
cell appears in `edges` with its true area fraction.

**`approx`** uses the cell-centre rule, matching fasterize semantics: a
boundary cell is included as a full run cell if and only if the polygon
edge crosses to the left of the cell centre. No `edges` are produced.
It uses a dedicated lightweight sweep that bypasses the coverage walker
entirely, and it is substantially faster.

```python
g = shapely.box(2.5, 4.5, 6.5, 8.5)
for mode in ("coverage", "approx"):
    r = cb.burn([shapely.to_wkb(g)], extent=(0, 10, 0, 10),
                shape=(10, 10), mode=mode)
    print(mode, len(r.runs), "runs,", len(r.edges), "edges")
# coverage 3 runs, 16 edges
# approx   4 runs, 0 edges
```

Lines and points are unaffected by mode.

Choose `coverage` when the fraction is the answer -- area-weighted
zonal statistics, fractional land cover, anything that has to conserve
area. Choose `approx` when you want a mask and you want it now.

---

## Polygons

Interior cells arrive run-length encoded by row, which is why the output
stays small even on enormous grids. Boundary cells arrive individually.

### Recovering total area

```python
def covered_area(r):
    xmin, xmax, ymin, ymax = r.extent
    nrow, ncol = r.shape
    cell = ((xmax - xmin) / ncol) * ((ymax - ymin) / nrow)
    full = (r.runs["col_end"] - r.runs["col_start"] + 1).sum()
    return cell * (full + r.edges["fraction"].sum(dtype="f8"))
```

```python
tri = shapely.Polygon([(13.3, 17.7), (88.1, 22.4), (41.9, 79.2)])
r = cb.burn([shapely.to_wkb(tri)], extent=(0, 100, 0, 100), shape=(41, 37))
covered_area(r), tri.area
# agrees to better than 1 part in 10,000 on an awkward non-square grid
```

### Holes and multipart geometries

Both are handled inside a single geometry, in either ring orientation,
and produce a single `id`:

```python
g = shapely.Polygon([(1, 1), (9, 1), (9, 9), (1, 9)],
                    [[(3, 3), (7, 3), (7, 7), (3, 7)]])
r = cb.burn([shapely.to_wkb(g)], extent=(0, 10, 0, 10), shape=(10, 10))
covered_area(r)   # 48.0, the 64-unit square minus the 16-unit hole
```

### Shared boundaries

Adjacent polygons with a shared edge produce complementary fractions that
sum to exactly 1.0 in every boundary cell. No gaps, no double counting:

```python
left, right = shapely.box(0, 0, 5, 10), shapely.box(5, 0, 10, 10)
r = cb.burn(shapely.to_wkb(np.array([left, right])),
            extent=(0, 10, 0, 10), shape=(20, 20))
m = cb.materialize(r, values=[1.0, 1.0], fn="sum",
                   edge_policy="fraction", background=0.0)
np.unique(np.round(m, 6))
# array([1.])
```

Every cell in the grid gets exactly 1.0. This is the property that makes
coverage mode usable for partitioning a surface between polygons.

### Geometry outside the grid

Handled correctly, including the case where a polygon entirely encloses
the grid:

```python
g = shapely.box(-100, -100, 100, 100)
r = cb.burn([shapely.to_wkb(g)], extent=(0, 10, 0, 10), shape=(5, 5))
len(r.edges), covered_area(r)
# (0, 100.0)  -- every cell fully covered, no boundary cells at all
```

---

## Lines

Lines produce one record per touched cell, with the **absolute length of
the line inside that cell**, in CRS units. Not a fraction, not a
normalised weight.

```python
line = shapely.LineString([(0.5, 0.5), (9.5, 7.5)])
r = cb.burn([shapely.to_wkb(line)], extent=(0, 10, 0, 10), shape=(10, 10))

len(r.lines), r.lines["length"].sum(dtype="f8"), line.length
# (16, 11.4018..., 11.4018...)
```

Length is conserved: summing the table gives back the length of the input
line. That is the invariant to test against when you are checking your own
pipeline.

Because the measure is in CRS units, geometry in degrees gives you length
in degrees. If you want metres, project first. controlledburn does no
geodesy and knows nothing about the CRS.

---

## Points

Points produce one record per touched cell and no measure column, because
a point is either in a cell or it is not.

```python
pts = [shapely.Point(0.5, 9.5), shapely.Point(9.5, 0.5), shapely.Point(15, 5)]
r = cb.burn([shapely.to_wkb(p) for p in pts],
            extent=(0, 10, 0, 10), shape=(10, 10))

r.points
# array([(1, 1, 1), (10, 10, 2)],
#       dtype=[('row','<i4'), ('col','<i4'), ('id','<i4')])
```

The third point falls outside the grid and is dropped silently. Note the
orientation: `(0.5, 9.5)` is near the top-left corner and lands in row 1,
column 1.

Repeated points in the same cell produce repeated records, one per input
geometry, so counting is just a group-by on `(row, col)`.

---

## Materializing to a dense array

Everything so far has been sparse. When you actually want pixels:

```python
m = cb.materialize(r)
```

### The built-in materializer

`cb.materialize()` is the fasterize reconciliation layer. It handles
**polygons only** -- `runs` and `edges`.

```python
cb.materialize(
    result,
    shape=None,              # defaults to result.shape
    values=None,             # burn value per id; default is the id itself
    fn="last",               # first, last, sum, min, max, count, any
    edge_policy="threshold", # threshold | fraction
    threshold=0.5,
    background=np.nan,
)
```

`values[k - 1]` is the value burned for geometry `id = k`. The default
burns the id itself, which is how you get a zone raster.

`fn` decides what happens when several geometries hit the same pixel.
`edge_policy` decides what happens on boundary cells:

- `"fraction"` writes `value * fraction`. Exact and area-conserving.
- `"threshold"` includes the cell whole if `fraction >= threshold`.
  With `threshold=0.5` this approximates the fasterize cell-centre rule
  but is not identical to it. If you want true fasterize parity, use
  `mode="approx"` at burn time instead of thresholding afterwards.

```python
geoms = [shapely.box(2, 4, 6, 8), shapely.box(4, 4, 8, 8)]
r = cb.burn([shapely.to_wkb(g) for g in geoms],
            extent=(0, 10, 0, 10), shape=(10, 10))

m = cb.materialize(r, values=[1.0, 1.0], fn="sum")
np.isfinite(m).sum()          # 24 cells in the union
(m[np.isfinite(m)] == 2).sum() # 8 cells in the overlap
```

Area conservation with `edge_policy="fraction"`:

```python
g = shapely.box(2.5, 4.5, 6.5, 8.5)
r = cb.burn([shapely.to_wkb(g)], extent=(0, 10, 0, 10), shape=(10, 10))
np.nansum(cb.materialize(r, values=[1.0], fn="sum", edge_policy="fraction"))
# 16.0 -- the polygon area, in cell units
```

### Materializing lines and points

`cb.materialize()` currently **ignores `lines` and `points` and returns an
all-background array** for line or point input. It does not warn. This is
a known gap (see [Sharp edges](#sharp-edges)). Until it is fixed, use
numpy directly:

```python
def materialize_all(r, shape=None, background=0.0):
    """Dense sum of all four tables. Mirrors R's materialise_chunk()."""
    nrow, ncol = shape or r.shape
    out = np.full((nrow, ncol), float(background), dtype="f8")

    for run in r.runs:
        out[run["row"] - 1, run["col_start"] - 1:run["col_end"]] += 1.0
    if len(r.edges):
        np.add.at(out, (r.edges["row"] - 1, r.edges["col"] - 1),
                  r.edges["fraction"].astype("f8"))
    if len(r.lines):
        np.add.at(out, (r.lines["row"] - 1, r.lines["col"] - 1),
                  r.lines["length"].astype("f8"))
    if len(r.points):
        np.add.at(out, (r.points["row"] - 1, r.points["col"] - 1), 1.0)
    return out
```

```python
line = shapely.LineString([(0.5, 0.5), (9.5, 7.5)])
r = cb.burn([shapely.to_wkb(line)], extent=(0, 10, 0, 10), shape=(10, 10))

cb.materialize(r, background=0.0).sum()   # 0.0   <- silently empty
materialize_all(r).sum()                  # 11.4018, the line length
```

**Do not run `materialize_all` on mixed-kind input** unless you mean it.
Summing a dimensionless fraction, a length in metres and a point count
into one array mixes three different units. Filter by `id` to one
geometry kind first, or burn each kind separately.

---

## Cropping and tiling

Burn once, then cut as many windows out of the result as you like.
`BurnResult.crop()` filters and clips the four sparse tables to a
sub-window, re-basing the row/col indices to 1 and snapping the window
outward to whole cells. Nothing is materialized, so cropping a burn that
covers a trillion-cell grid is just structured-array filtering. It is the
Python counterpart of the R package's `crop_burn()`.

```python
r = cb.burn([shapely.to_wkb(shapely.box(2.5, 4.5, 6.5, 8.5))],
            extent=(0, 10, 0, 10), shape=(10, 10))

tile = r.crop((3, 7, 5, 9))   # (xmin, xmax, ymin, ymax), same order as extent
tile.shape, tile.extent
# ((4, 4), (3.0, 7.0, 5.0, 9.0))
```

The target window uses the same `(xmin, xmax, ymin, ymax)` ordering as
`extent`, matching the R package's `crop_burn()`.

Because `crop()` returns a `BurnResult` and `materialize()` is also a
method, the tile workflow reads as a single chain -- and nothing dense is
allocated until the final step:

```python
window = (3, 7, 5, 9)
dense = r.crop(window).materialize(fn="sum", edge_policy="fraction")
```

Cropping then materializing gives exactly the same pixels as slicing the
full dense array over that window, without ever allocating the full array:

```python
full = r.materialize(fn="sum", edge_policy="fraction")
np.array_equal(np.nan_to_num(dense), np.nan_to_num(full[1:5, 3:7]))
# True
```

That equivalence is what makes `crop()` the tiling primitive: burn one
huge extent, then loop over windows writing tiles.

```python
def tiles(r, tile_rows, tile_cols):
    """Yield (row_offset, col_offset, dense_tile) over a grid of windows."""
    xmin, xmax, ymin, ymax = r.extent
    nrow, ncol = r.shape
    dx, dy = (xmax - xmin) / ncol, (ymax - ymin) / nrow
    for r0 in range(0, nrow, tile_rows):
        for c0 in range(0, ncol, tile_cols):
            win = (xmin + c0 * dx,
                   xmin + min(c0 + tile_cols, ncol) * dx,
                   ymax - min(r0 + tile_rows, nrow) * dy,
                   ymax - r0 * dy)
            yield r0, c0, r.crop(win).materialize(fn="sum",
                                                  edge_policy="fraction")
```

A window that does not overlap the grid warns and returns empty tables
with `shape == (0, 0)`.

---

## Georeferencing the output

The tables carry grid indices, not coordinates. Two helpers cover almost
every case.

### Affine transform

```python
def transform(r):
    """GDAL/rasterio-style affine (a, b, c, d, e, f)."""
    xmin, xmax, ymin, ymax = r.extent
    nrow, ncol = r.shape
    return ((xmax - xmin) / ncol, 0.0, xmin,
            0.0, -(ymax - ymin) / nrow, ymax)
```

```python
transform(r)   # (1.0, 0.0, 0.0, 0.0, -1.0, 10.0) for extent (0,10,0,10), shape (10,10)
```

Feeding rasterio:

```python
import rasterio
from rasterio.transform import Affine

m = materialize_all(r)
with rasterio.open("out.tif", "w", driver="GTiff",
                   height=r.shape[0], width=r.shape[1], count=1,
                   dtype="float32", crs="EPSG:4326",
                   transform=Affine(*transform(r))) as dst:
    dst.write(m.astype("float32"), 1)
```

The `-` on the y term and the `ymax` origin are what encode "row 1 is the
top". Because controlledburn already emits top-down rows, a materialized
array drops into a raster with no flipping.

### Cell centres

```python
def cell_centers(row, col, extent, shape):
    """1-based (row, col) to (x, y) cell centre coordinates."""
    xmin, xmax, ymin, ymax = extent
    nrow, ncol = shape
    x = xmin + (np.asarray(col) - 0.5) * (xmax - xmin) / ncol
    y = ymax - (np.asarray(row) - 0.5) * (ymax - ymin) / nrow
    return x, y
```

```python
r = cb.burn([shapely.to_wkb(p) for p in
             [shapely.Point(0.5, 9.5), shapely.Point(9.5, 0.5)]],
            extent=(0, 10, 0, 10), shape=(10, 10))
cell_centers(r.points["row"], r.points["col"], r.extent, r.shape)
# (array([0.5, 9.5]), array([9.5, 0.5]))
```

This is the Python equivalent of what the `vaster` package gives R users:
primitive cell-to-xy arithmetic on the same `(row, col)` schema.

---

## pandas, scipy.sparse, xarray

### pandas

Structured arrays convert with no copy step:

```python
import pandas as pd
pd.DataFrame(r.edges)      # columns: row, col, fraction, id
pd.DataFrame(r.lines)      # columns: row, col, length, id
```

Expanding runs to one row per cell, when you want a long table:

```python
runs = r.runs
n = runs["col_end"] - runs["col_start"] + 1
long = pd.DataFrame({
    "row": np.repeat(runs["row"], n),
    "col": np.concatenate([np.arange(a, b + 1)
                           for a, b in zip(runs["col_start"], runs["col_end"])]),
    "id":  np.repeat(runs["id"], n),
})
```

Area-weighted zonal statistics against an existing raster, which is the
job coverage mode exists for:

```python
edges = pd.DataFrame(r.edges)
edges["value"] = raster[edges["row"] - 1, edges["col"] - 1]
weighted = (edges.groupby("id")
                 .apply(lambda d: np.average(d["value"], weights=d["fraction"])))
```

### scipy.sparse

The output is already sparse; this just changes the container:

```python
from scipy.sparse import coo_array

def to_coo(r, kind="polygon"):
    nrow, ncol = r.shape
    rows, cols, vals = [], [], []
    if kind == "polygon":
        runs = r.runs
        if len(runs):
            n = runs["col_end"] - runs["col_start"] + 1
            rows.append(np.repeat(runs["row"], n) - 1)
            cols.append(np.concatenate([np.arange(a, b + 1) for a, b in
                        zip(runs["col_start"], runs["col_end"])]) - 1)
            vals.append(np.ones(n.sum()))
        if len(r.edges):
            rows.append(r.edges["row"] - 1)
            cols.append(r.edges["col"] - 1)
            vals.append(r.edges["fraction"].astype("f8"))
    elif kind == "line":
        rows.append(r.lines["row"] - 1); cols.append(r.lines["col"] - 1)
        vals.append(r.lines["length"].astype("f8"))
    elif kind == "point":
        rows.append(r.points["row"] - 1); cols.append(r.points["col"] - 1)
        vals.append(np.ones(len(r.points)))
    if not rows:
        return coo_array((nrow, ncol))
    return coo_array((np.concatenate(vals),
                      (np.concatenate(rows), np.concatenate(cols))),
                     shape=(nrow, ncol))
```

```python
polys = shapely.buffer(shapely.points(rng.uniform(0, 1000, 500),
                                      rng.uniform(0, 1000, 500)), 5)
r = cb.burn(shapely.to_wkb(polys), extent=(0, 1000, 0, 1000),
            shape=(20000, 20000), mode="approx")
s = to_coo(r)
s.nnz            # 15,544,875 of 400,000,000 cells -- 3.9 percent occupied
```

### xarray

```python
import xarray as xr

nrow, ncol = r.shape
xmin, xmax, ymin, ymax = r.extent
da = xr.DataArray(
    materialize_all(r),
    dims=("y", "x"),
    coords={
        "y": ymax - (np.arange(nrow) + 0.5) * (ymax - ymin) / nrow,
        "x": xmin + (np.arange(ncol) + 0.5) * (xmax - xmin) / ncol,
    },
)
```

Note `y` is descending, which is what you want for a north-up raster and
what rioxarray expects.

---

## Scale

The whole point of the sparse output is that grid size costs almost
nothing. 2,000 buffered points burned in `approx` mode:

| Grid | Cells | Time | Runs emitted |
|---|---|---|---|
| 1,000 x 1,000 | 1M | 0.01s | 42,751 |
| 10,000 x 10,000 | 100M | 0.09s | 428,163 |
| 100,000 x 100,000 | 10B | 1.5s | 4,282,455 |

Ten billion cells in a second and a half. A dense float64 array at that
size would need 80 terabytes.

The real-world case from the R package: Antarctic rock outcrop polygons
(25,954 geometries) against a REMA 2m DEM grid of 2.7 million by 2.9
million pixels, roughly 8 trillion cells. Eight seconds, 438 MB of sparse
output, 99.8 percent sparsity.

Cost scales with perimeter and with the number of rows each geometry
spans, not with area. Doubling the resolution roughly doubles the work;
it does not quadruple it, the way a dense rasterizer does.

The GIL is released around the burn, so `concurrent.futures.ThreadPoolExecutor`
gives real parallelism across chunks of geometry.

---

## Sharp edges

Things that will bite you, in rough order of likelihood.

**1. `shape` is `(nrow, ncol)`.** Numpy order, not `(width, height)`.
Silent when the grid is square.

**2. Indices are 1-based.** `r.runs["row"] - 1` to index numpy. Forgetting
this shifts everything up and left by one cell, which looks almost right.

**3. Row 1 is the top.** If your picture is upside down, you flipped
something that was already correct.

**4. `cb.materialize()` silently ignores lines and points.** It handles
`runs` and `edges` only. Burning a LineString and calling `materialize()`
returns an all-background array with no warning. Use `materialize_all()`
above. This is the single most likely way to get a wrong answer quietly.

**5. `extent` and `shape` are mandatory.** No geometry-derived default, no
`resolution=`. Use `resolve_grid()` above.

**6. `burn()` needs WKB, not shapely objects.** Passing a shapely geometry
gives `TypeError: a bytes-like object is required, not 'Polygon'`. Passing
a single bare `bytes` object instead of a list gives the genuinely
confusing `TypeError: a bytes-like object is required, not 'int'`, because
iterating bytes yields integers. Use `as_wkb_list()` above.

**7. `repr(BurnResult)` prints every table in full.** In a notebook with a
million runs this is unpleasant. Print a summary instead:

```python
def summary(r):
    nrow, ncol = r.shape
    interior = int((r.runs["col_end"] - r.runs["col_start"] + 1).sum())
    touched = interior + len(r.edges) + len(r.lines) + len(r.points)
    ids = np.unique(np.concatenate([r.runs["id"], r.edges["id"],
                                    r.lines["id"], r.points["id"]]))
    lines = [f"<BurnResult> {ncol} x {nrow} grid, {len(ids)} geometries",
             f"  runs:   {len(r.runs)} ({interior} interior cells)",
             f"  edges:  {len(r.edges)} polygon boundary cells"]
    if len(r.lines):
        lines.append(f"  lines:  {len(r.lines)} cells")
    if len(r.points):
        lines.append(f"  points: {len(r.points)} cells")
    lines.append(f"  sparsity: {100 * (1 - touched / (nrow * ncol)):.1f}% empty")
    return "\n".join(lines)
```

```
<BurnResult> 10 x 10 grid, 3 geometries
  runs:   3 (9 interior cells)
  edges:  16 polygon boundary cells
  lines:  16 cells
  points: 1 cells
  sparsity: 58.0% empty
```

**8. Measures are in CRS units.** Degrees in, square degrees and degrees
out. controlledburn does no geodesy.

**9. `edge_policy="threshold"` is not fasterize parity.** For true
cell-for-cell fasterize agreement use `mode="approx"` at burn time.

**10. Independent versions.** The Python package (`0.3.0`) is versioned
separately from the R package (`0.2.0`); they share the same C++ core, so
numeric output matches regardless of the package version numbers. See
`python/CHANGELOG.md` for the Python history.

---

## R to Python cheatsheet

| Concept | R | Python |
|---|---|---|
| Burn | `burn(x, extent, dimension)` | `cb.burn(wkb, extent, shape)` |
| Extent | `c(xmin, xmax, ymin, ymax)` | `(xmin, xmax, ymin, ymax)` |
| Dimensions | `c(ncol, nrow)` | `(nrow, ncol)` |
| Mode | `mode = "coverage" \| "approx"` | `mode="coverage" \| "approx"` |
| Interior | `r$runs` | `r.runs` |
| Boundary | `r$edges` (`$fraction`) | `r.edges` (`["fraction"]`) |
| Lines | `r$lines` (`$length`) | `r.lines` (`["length"]`) |
| Points | `r$points` | `r.points` |
| Grid echo | `r$extent`, `r$dimension` | `r.extent`, `r.shape` |
| Dense | `materialise_chunk(r)` | `cb.materialize(r)` or `r.materialize()` (polygons only) |
| Crop / tile | `crop_burn(r, c(xmin, xmax, ymin, ymax))` | `r.crop((xmin, xmax, ymin, ymax))` |
| Summary | `print(r)` | `summary(r)` recipe above |
| From bbox | `burn(x)` | `resolve_grid()` recipe above |
| By resolution | `burn(x, resolution = 0.5)` | `resolve_grid()` recipe above |
| Input types | geos, sfc, wk_wkb, blob, raw list | WKB bytes |
| Cell to xy | `vaster` package | `cell_centers()` recipe above |

Indices, orientation and all numeric values are identical between the two.
The shared fixtures in `fixtures/` are read by the C++, R and Python test
suites precisely to keep it that way.

---

## See also

- `cpp/README.md` and `doc/cpp-core.md` for the C++ core.
- `vignette("architecture")` in the R package for how the engine works.
- `inst/docs-design/` for the design records, including the unified
  polygon/line/point rasterization design.
- [exactextract](https://github.com/isciences/exactextract) by Daniel
  Baston, the source of the coverage fraction algorithm.
- [fasterize](https://github.com/ecohealthalliance/fasterize) by Noam
  Ross, the origin of the scanline approach.
