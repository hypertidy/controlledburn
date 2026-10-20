# Sparse rasterization of polygon, line, and point geometry

Rasterizes geometry into a sparse output schema, with one type-pure
table per geometry kind. Polygons produce run-length-encoded interior
cells (`$runs`) and weighted boundary cells (`$edges`, with a coverage
`fraction` in `[0, 1]`). Lines produce one record per touched cell with
the absolute `length` of line within the cell, in CRS units (`$lines`).
Points produce one record per touched cell, with no measure column
(`$points`). See the unified geometry rasterization design doc in
`inst/docs-design/`.

## Usage

``` r
burn(
  x,
  extent = NULL,
  dimension = NULL,
  resolution = NULL,
  mode = c("coverage", "approx")
)

burn_scanline(x, extent = NULL, dimension = NULL, resolution = NULL)
```

## Arguments

- x:

  geometry input, one of:

  - a `geos_geometry` vector (from geos)

  - an `sfc` geometry column (from sf)

  - a `wk_wkb` vector (from wk)

  - a list of raw vectors containing WKB

- extent:

  numeric vector `c(xmin, xmax, ymin, ymax)` defining the raster extent.
  If `NULL` (default), derived from the bounding box of `x`.

- dimension:

  integer vector `c(ncol, nrow)` defining the grid dimensions. If `NULL`
  (default), fitted to the extent with at most 256 cells along the
  longer axis, preserving aspect ratio. Mutually exclusive with
  `resolution`.

- resolution:

  numeric, cell size (scalar for square cells, or `c(dx, dy)`). If
  supplied, `dimension` is computed from `extent / resolution`. Mutually
  exclusive with `dimension`.

- mode:

  character, one of `"coverage"` (default) or `"approx"`. `"coverage"`
  computes exact coverage fractions for polygon boundary cells (the
  `$edges` table). `"approx"` uses the cell-center rule (fasterize
  semantics): a boundary cell is included as a full run cell iff the
  polygon edge crosses to the left of the cell center. No `$edges` are
  produced for polygons. Faster (skips traversal math), but approximate.
  Lines and points are unaffected by mode.

## Value

A list with class `"controlledburn"` containing:

- `runs`:

  data.frame with columns `row`, `col_start`, `col_end`, `id` —
  run-length encoded polygon interior cells (always full coverage).
  Empty for line/point input.

- `edges`:

  data.frame with columns `row`, `col`, `fraction`, `id` — polygon
  boundary cells with partial coverage; `fraction` is in (0, 1). Empty
  for line/point input, or when `mode = "approx"`.

- `lines`:

  data.frame with columns `row`, `col`, `length`, `id` — line cells;
  `length` is the absolute length of the line within the cell, in CRS
  units. Empty for polygon/point input.

- `points`:

  data.frame with columns `row`, `col`, `id` — point cells; no measure
  column (a point is either in a cell or it isn't). Empty for
  polygon/line input.

- `extent`:

  the raster extent (`c(xmin, xmax, ymin, ymax)`)

- `dimension`:

  the grid dimensions (`c(ncol, nrow)`)

Row and column indices are 1-based. Row 1 is the top (ymax) row. The
`id` column is a 1-based index into the input geometry vector.

The measure is computed in the CRS the geometry and grid are specified
in — *not* the geodesic measure on the Earth. Geometry in degrees
produces area in square degrees and length in degrees. For geodesic
measures, project to an appropriate CRS first.

## Details

Mixed geometry-kind input (e.g. a vector containing both polygons and
lines) is permitted and produces tables populated for whichever kinds
are present. Each table's measure column means exactly one thing —
polygon `$edges$fraction` is dimensionless, `$lines$length` is in CRS
units, points have no measure — so combining them in a single
materialized matrix would silently mix units: run separate burns for
each kind before
[`materialize_chunk()`](https://hypertidy.github.io/controlledburn/reference/materialise_chunk.md)
(which sums whatever it is given) or
[`extract_burn()`](https://hypertidy.github.io/controlledburn/reference/extract_burn.md)
with `fun = "sum"` (which refuses mixed kinds).

GeometryCollection input warns and is skipped. Curved types
(CircularString, CompoundCurve, etc.) must be linearised by the caller
before calling.

Memory usage is O(perimeter) rather than O(bounding-box area), so no
`tile_size` parameter is required (in contrast to
[`burn_sparse()`](https://hypertidy.github.io/controlledburn/reference/burn_sparse.md),
which is polygon-only and tile-bounded).

## Examples

``` r
if (requireNamespace("geos", quietly = TRUE)) {
  library(geos)

  # Polygon: produces $runs (interior) and $edges (fraction in [0, 1])
  poly <- as_geos_geometry("POLYGON ((0.5 0.5, 2.5 0.5, 2.5 2.5, 0.5 2.5, 0.5 0.5))")
  r <- burn(poly, extent = c(0, 3, 0, 3), dimension = c(3, 3))

  # Approx mode: fasterize-style cell-center rule, runs only
  r <- burn(poly, extent = c(0, 3, 0, 3), dimension = c(3, 3), mode = "approx")

  # Line: produces $lines with absolute length in CRS units
  line <- as_geos_geometry("LINESTRING (0.5 0.5, 2.5 2.5)")
  r <- burn(line, extent = c(0, 3, 0, 3), dimension = c(3, 3))

  # Point: produces $points with row/col only (no measure column)
  pt <- as_geos_geometry("POINT (1.5 1.5)")
  r <- burn(pt, extent = c(0, 3, 0, 3), dimension = c(3, 3))

  # Defaults: extent from bbox, 256-cell fitted grid
  r <- burn(poly)

  # By resolution
  r <- burn(poly, resolution = 0.01)
}
```
