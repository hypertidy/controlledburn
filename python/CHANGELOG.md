# Changelog — controlledburn (Python)

The Python package is versioned independently of the R package. Both
share the same C++ core, so numeric output matches across languages.

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
