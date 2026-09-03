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

## 0.2.0

- First tagged Python surface: `burn()` accepting WKB bytes, shapely
  geometries, and geopandas input; the four sparse output tables
  (`runs`, `edges`, `lines`, `points`); `mode="coverage"` and
  `mode="approx"`; `resolve_grid()` and `as_wkb_list()` conveniences;
  a summary `__repr__`; and the polygon-only `materialize()`.
