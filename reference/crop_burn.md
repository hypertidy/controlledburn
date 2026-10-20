# Crop a controlledburn result to a target extent

Filters and clips the sparse tables (runs, edges, lines, points) to a
sub-window defined by a target extent. Returns a new `controlledburn`
object with offset row/col indices and snapped extent. No dense
allocation — just data frame filtering.

## Usage

``` r
crop_burn(x, target)
```

## Arguments

- x:

  A `controlledburn` object returned by
  [`burn()`](https://hypertidy.github.io/controlledburn/reference/burn.md).

- target:

  numeric vector `c(xmin, xmax, ymin, ymax)` defining the target window.
  Snapped outward to cell boundaries.

## Value

A `controlledburn` object covering only the target window, with row/col
indices re-based to 1.

## Examples

``` r
if (requireNamespace("geos", quietly = TRUE)) {
  library(geos)
  poly <- as_geos_geometry("POLYGON ((1 1, 9 1, 9 9, 1 9, 1 1))")
  r <- burn(poly, extent = c(0, 10, 0, 10), dimension = c(100L, 100L))
  sub <- crop_burn(r, target = c(2, 5, 3, 7))
  sub
}
#> <controlledburn> 30 x 40 grid, 1 geometry
#>   runs:   40 (1200 interior cells)
#>   edges:  0 polygon boundary cells
#>   sparsity: 0.0% empty
```
