# Summarise a burn without materialising it

Compute a compact description of a
[`burn()`](https://hypertidy.github.io/controlledburn/reference/burn.md)
result in time proportional to the number of records in its sparse
tables, never to the number of grid cells, so it is cheap even for a
national grid of billions of cells.

## Usage

``` r
# S3 method for class 'controlledburn'
summary(object, tiles = 256L, ...)

# S3 method for class 'summary.controlledburn'
print(x, ascii = FALSE, ...)
```

## Arguments

- object:

  A `controlledburn` object returned by
  [`burn()`](https://hypertidy.github.io/controlledburn/reference/burn.md).

- tiles:

  Target number of tiles along the longer grid axis for the overview,
  passed to
  [`overview_burn()`](https://hypertidy.github.io/controlledburn/reference/overview_burn.md).

- ...:

  Ignored.

- x:

  A `summary.controlledburn` object.

- ascii:

  Also print the overview as text.

## Value

A list of class `summary.controlledburn`.

## Details

The pieces are:

- `overview`:

  a coarse matrix of mean coverage per tile, see
  [`overview_burn()`](https://hypertidy.github.io/controlledburn/reference/overview_burn.md).
  When the grid is small enough that one tile is one cell, this is the
  exact coverage.

- `row_marginal`, `col_marginal`:

  exact covered-cell totals per grid row and per grid column, at full
  resolution.

- `content`:

  the extent and 1-based row/col range of touched cells, or `NULL` for
  an empty burn.

- `area`:

  per-geometry totals in cells: `interior` (run cells), `boundary`
  (summed edge fractions), `lines` (summed length, in extent units),
  `points` (count).

- `fractions`:

  quantiles of the edge coverage fractions.

- `records`:

  row counts of the sparse tables.

Coverage here means: run cells 1 each, edge cells their fraction, line
and point cells 1 each (presence).

The number of geometries is taken as the largest id present, so
geometries after the last one that touched the grid are not counted.

## Examples

``` r
if (requireNamespace("geos", quietly = TRUE)) {
  poly <- geos::as_geos_geometry(c("POLYGON ((1 1, 9 1, 9 9, 1 9, 1 1))",
                                   "POLYGON ((2 2, 4 3, 3 5, 2 2))"))
  r <- burn(poly, extent = c(0, 10, 0, 10), dimension = c(100L, 100L))
  s <- summary(r)
  s
  s$area
}
#>   id interior boundary lines points
#> 1  1     6400        0     0      0
#> 2  2      216       34     0      0
```
