# Slice a burn into tiles

Cuts a
[`burn()`](https://hypertidy.github.io/controlledburn/reference/burn.md)
result into a regular grid of tiles in one pass over its records,
anchored at the top-left cell of the grid. A record in cell `(row, col)`
goes to the tile that holds it, with indices local to that tile, and a
run is split where it crosses a tile seam. Nothing is recomputed, so
each tile is
[`crop_burn()`](https://hypertidy.github.io/controlledburn/reference/crop_burn.md)
of the whole result to the tile's extent, record for record, and the
materialised tiles put back together are the materialised whole grid.

## Usage

``` r
slice_burn(x, dimension = c(256L, 256L), threads = 1L)
```

## Arguments

- x:

  A `controlledburn` object returned by
  [`burn()`](https://hypertidy.github.io/controlledburn/reference/burn.md).

- dimension:

  integer, the tile size in cells, `c(ncol, nrow)`, or one number for
  square tiles. Tiles on the right and bottom edges are smaller when
  this does not divide the grid.

- threads:

  integer, the number of threads (default `1`; `0` for all available).
  The tiles are the same. Results of fewer than 262,144 records are
  sliced on one thread whatever `threads` is.

## Value

A list of `controlledburn` objects, one per tile holding at least one
record, in row-major tile order (tiles without records are left out).
Each has the usual tables, its own `extent` and `dimension`, and also

- `tile`:

  integer `c(row = , col = )`, the tile's 1-based position in the tile
  grid, row 1 at the top.

- `offset`:

  integer `c(row = , col = )`, the cells before the tile in the whole
  grid: cell `(r, c)` of the tile is cell
  `(r + offset["row"], c + offset["col"])` of `x`.

## Details

Calling
[`crop_burn()`](https://hypertidy.github.io/controlledburn/reference/crop_burn.md)
once per tile gives the same tables but scans every record for every
tile; `slice_burn()` costs about one pass over the records whatever the
number of tiles.

## Examples

``` r
if (requireNamespace("geos", quietly = TRUE)) {
  poly <- geos::as_geos_geometry("POLYGON ((1 1, 9 1, 9 9, 1 9, 1 1))")
  r <- burn(poly, extent = c(0, 10, 0, 10), dimension = c(100L, 100L))
  tiles <- slice_burn(r, c(32L, 32L))
  length(tiles)
  tiles[[1]]$tile
  tiles[[1]]$extent
  ## each tile is the crop of the whole burn to the tile
  identical(materialise_chunk(tiles[[1]]),
            materialise_chunk(r, target = tiles[[1]]$extent))
}
#> [1] TRUE
```
