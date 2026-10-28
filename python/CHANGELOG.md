# Changelog — controlledburn (Python)

The Python package is versioned independently of the R package. Both
share the same C++ core; coverage fractions, lengths and areas match
across languages, though the index base differs (Python 0-based, R
1-based -- see 0.4.0).

## Unreleased

### Fixes from the shared C++ core

- **Corner passes and near-coincident edges (coverage mode, #25).** A
  cell was overcounted where an edge passes exactly through a cell
  corner, or where two edges nearly coincide (a sliver, a ring doubling
  back along itself), by up to a whole cell. Cells crossed by several
  traversals now get their coverage by Green's theorem, which does not
  depend on the order of the traversal ends. Covered in
  `tests/test_regressions.py` and the shared `corner_pass` fixture.
- **Cell-centre ties (approx mode, #25).** The side of a cell centre an
  edge passes is decided exactly, so a centre two edges both miss by
  1e-13 no longer burns.
- **Holes removed too much area (coverage mode).** A cell crossed only by
  a hole's boundary, inside the exterior, got the hole's negative fraction
  alone and was dropped; it now keeps 1 minus the hole's fraction. Covered
  in `tests/test_regressions.py`, including a cell-by-cell comparison with
  shapely's exact intersection.
- **A cell with edges only along its walls was filled (coverage mode).**
  With lattice-aligned input, a cell whose polygon edges lie along both
  of its walls, crossing the row centre in opposite directions, has zero
  coverage and zero net winding; it got no boundary record and was filled
  as part of the interior run, one whole cell too many. Covered in
  `tests/test_regressions.py` and the shared `wall_notch` fixture.
These landed in the core through the R package (#16, #19) and apply to
the Python package unchanged; they are now covered by Python tests too
(`tests/test_regressions.py`).

- **Lattice snapping (coverage mode).** A vertex a few ulp off a cell edge
  (0.3 on a 0.1 grid) no longer sends the cell walker off the grid until
  memory runs out.
- **Row-centre crossings on a column wall (coverage mode).** Such a
  crossing now registers its winding exactly once, so the row keeps its
  interior runs from that column on.
- **Vertex on a row centre line (approx mode).** The crossing direction
  now comes from the edge's endpoints, so a vertex an ulp off the
  recomputed centre line no longer fills the row to the right edge of the
  grid.
- **Non-finite extents are rejected.** An `inf` bound, or a span that
  overflows double precision, used to reach the walker with an infinite
  cell size and crash the interpreter. `burn()` now raises `ValueError`,
  and so does the core itself.

### Input handling

- The default extent is now derived from the WKB by the core's envelope
  scan (new `_core.bbox_wkb`) on every input path, not only for Arrow
  input. **shapely is no longer needed to omit `extent`.** All-empty
  input raises "no finite coordinates" instead of a misleading
  invalid-extent error with a NumPy all-NaN warning.
- Sliced Arrow arrays (a non-zero offset, e.g. `pa.array(...).slice(...)`)
  are now supported; they used to raise.
- WKB passed as a buffer with items wider than a byte (e.g.
  `memoryview.cast("Q")`) is read in full; only the item count was used
  before, truncating the WKB. Non-contiguous buffers are rejected rather
  than read as if contiguous.
- `extent` and `shape` of the wrong length, and `ncol`/`nrow` beyond
  int32, raise a clear `ValueError` instead of an unpacking error or a
  pybind11 signature error.
- `_core.materialize()` checks that the column arrays of each table have
  equal lengths.

### Documentation

- `burn()` documented geometry ids as `k + 1`; they are `k` since 0.4.0.

## 0.4.0

### Breaking changes

- **0-based indices with an exclusive `col_end`.** Tables now use 0-based
  `row`/`col` with row 0 at the top, and `runs` cover the half-open range
  `[col_start, col_end)` -- so `col_end - col_start` is the run length and
  `buf[row, col_start:col_end]` is the slice. Geometry k gets `id = k`
  (also 0-based). Index a dense array directly; the old `- 1` corrections
  are gone. This moves the `+1` offset out of the shared C++ core (now
  0-based, matching controlledburn-rs and every non-R consumer) and into
  the R package's cpp11 shim, so R output is unchanged. `crop()`,
  `materialize()` (now `values[k]`), `cell_centers()` and the guide's
  recipe snippets are updated to match.

## 0.3.0

### Breaking changes

- **`extent` replaces `bounds`.** The grid-extent argument to `burn()`
  and `resolve_grid()`, and the `BurnResult.bounds` attribute, are
  renamed to `extent`. There is no compatibility alias — `bounds=` and
  `r.bounds` no longer exist.
- **Extent ordering now matches R.** `extent` is
  `(xmin, xmax, ymin, ymax)`, the same ordering as the R package's
  `extent = c(xmin, xmax, ymin, ymax)` (previously the rasterio-style
  `(xmin, ymin, xmax, ymax)`). `shape=(nrow, ncol)` is unchanged — it
  remains the numpy transpose of R's `dimension = c(ncol, nrow)`.

### What's new

- **`BurnResult.crop(target)`** — crop a sparse result to a sub-window,
  the Python counterpart of the R package's `crop_burn()`. Filters and
  clips the `runs`/`edges`/`lines`/`points` tables, re-bases row/col
  indices to 1, snaps the window outward to whole cells, and updates
  `extent`/`shape`. No dense allocation — just structured-array
  filtering. A non-overlapping window warns and returns empty tables
  with `shape == (0, 0)`.
- **`BurnResult.materialize(...)`** — method form of the module-level
  `materialize()`, so the crop/materialize tile workflow reads as a
  chain: `r.crop(window).materialize(fn="sum", edge_policy="fraction")`.
- **Arrow-native input.** `burn()` accepts anything implementing the
  Arrow C data interface (`__arrow_c_array__` / `__arrow_c_stream__`)
  with `(large_)binary` or `geoarrow.wkb` geometry — pyarrow, nanoarrow,
  geoarrow-pyarrow, duckdb, polars. It reads WKB straight from the Arrow
  values buffer through a new zero-copy core entry
  (`_core.burn_wkb_arrow`): no per-geometry Python objects and no
  shapely. Nulls consume an id; chunked streams are concatenated with an
  id offset. Requires `nanoarrow` (optional dependency
  `controlledburn[arrow]`). When `extent` is omitted it is derived from
  the geometry by a WKB-envelope scan in the C++ core
  (`bbox_wkb` / `_core.bbox_wkb_arrow`), so the Arrow path needs no
  shapely at all. `as_wkb_list()` also gains an Arrow fallback branch.

### Core

- **`bbox_wkb()`** added to the shared C++ core: the axis-aligned
  envelope of a set of WKB blobs (skipping nulls, unparseable blobs, and
  non-finite coordinates), so bindings can derive a default grid extent
  without a geometry library.

## 0.2.0

- First tagged Python surface: `burn()` accepting WKB bytes, shapely
  geometries, and geopandas input; the four sparse output tables
  (`runs`, `edges`, `lines`, `points`); `mode="coverage"` and
  `mode="approx"`; `resolve_grid()` and `as_wkb_list()` conveniences;
  a summary `__repr__`; and the polygon-only `materialize()`.
