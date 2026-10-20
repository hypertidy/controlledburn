# Coarse overview of a burn

Aggregate a
[`burn()`](https://hypertidy.github.io/controlledburn/reference/burn.md)
result onto a grid of square tiles of `block` cells and return the mean
coverage per tile (in `[0, 1]` for non-overlapping polygons). Tiles cut
off by the grid boundary are normalised by their actual cell count. Cost
is proportional to the number of records plus the number of tiles: a run
that spans several tiles contributes partial lengths to its two end
tiles, and the fully covered tiles between them are filled with a
difference array, so nothing is expanded to cells whatever the block
size.

## Usage

``` r
overview_burn(x, tiles = 256L, block = NULL)
```

## Arguments

- x:

  A `controlledburn` object returned by
  [`burn()`](https://hypertidy.github.io/controlledburn/reference/burn.md).

- tiles:

  Target number of tiles along the longer grid axis; ignored if `block`
  is given.

- block:

  Tile size in cells. `block = 1` gives the exact coverage, the same as
  [`materialise_chunk()`](https://hypertidy.github.io/controlledburn/reference/materialise_chunk.md)
  for polygons, with `NA` for 0.

## Value

A numeric matrix (tile rows by tile columns) with attributes `extent`
(`c(xmin, xmax, ymin, ymax)` of the tile grid, which extends past the
burn grid on the right and bottom when `block` does not divide its
dimension) and `block`. Untouched tiles are `NA`.

## Details

The result is a matrix in raster order (row 1 is the top), ready for
`ximage::ximage(ov, extent = attr(ov, "extent"))`.

## Examples

``` r
if (requireNamespace("geos", quietly = TRUE)) {
  poly <- geos::as_geos_geometry("POLYGON ((1 1, 9 1, 9 9, 1 9, 1 1))")
  r <- burn(poly, extent = c(0, 10, 0, 10), dimension = c(1000L, 1000L))
  ov <- overview_burn(r, tiles = 20)
  dim(ov)
  attr(ov, "extent")
}
#> [1]  0 10  0 10
```
