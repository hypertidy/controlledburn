# Extract values from a controlledburn result at point locations

Queries the sparse tables directly; no dense matrix is created. For
points at cell centres this reproduces
[`materialize_chunk()`](https://hypertidy.github.io/controlledburn/reference/materialise_chunk.md)
exactly:
`extract_burn(x, xy, fun) == materialize_chunk(x, fun)[cbind(row, col)]`

## Usage

``` r
extract_burn(x, xy, fun = c("id", "sum"))
```

## Arguments

- x:

  A `controlledburn` object returned by
  [`burn()`](https://hypertidy.github.io/controlledburn/reference/burn.md).

- xy:

  A two-column matrix or data.frame of x, y coordinates in the CRS of
  the burn.

- fun:

  character, `"id"` (last geometry written wins, `NA` background) or
  `"sum"` (accumulated coverage/length/count, `0` background). Same
  semantics as
  [`materialize_chunk()`](https://hypertidy.github.io/controlledburn/reference/materialise_chunk.md).

## Value

A numeric vector of length `nrow(xy)`. Points outside the grid extent
are `NA` for both funs.

## Details

Interior runs are matched by row equality plus column interval
containment; edges, lines, and points by exact (row, col) match. Cost is
O(points + hits), independent of grid dimension.

## Examples

``` r
if (requireNamespace("geos", quietly = TRUE)) {
  library(geos)
  poly <- as_geos_geometry("POLYGON ((0 0, 10 0, 10 10, 0 10, 0 0))")
  r <- burn(poly, extent = c(0, 10, 0, 10), dimension = c(10L, 10L))
  # value at three cell centres (fun = "id": geometry id, else NA)
  extract_burn(r, cbind(c(1.5, 5.5, 8.5), c(1.5, 5.5, 8.5)), fun = "id")
}
#> [1] 1 1 1
```
