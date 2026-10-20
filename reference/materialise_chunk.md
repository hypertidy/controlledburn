# Materialise a controlledburn sparse result into a dense matrix

Converts the sparse run/edge/line tables from
[`burn()`](https://hypertidy.github.io/controlledburn/reference/burn.md)
into a dense matrix. With `target`, only the cells inside a sub-window
of the burn grid are materialised: the sparse tables are first filtered
and re-indexed by
[`crop_burn()`](https://hypertidy.github.io/controlledburn/reference/crop_burn.md)
and only the sub-window is allocated, so a chunk of a very large burn
grid can be realised without ever allocating the full grid.

## Usage

``` r
materialise_chunk(x, fun = c("sum", "id"), target = NULL)

materialize_chunk(x, fun = c("sum", "id"), target = NULL)
```

## Arguments

- x:

  A `controlledburn` object returned by
  [`burn()`](https://hypertidy.github.io/controlledburn/reference/burn.md).

- fun:

  character, the pixel function:

  `"sum"`

  :   (default) Polygon runs contribute 1, edges contribute their
      coverage fraction, lines contribute length, points contribute 1.
      Values accumulate additively across overlapping geometries.
      Background is 0. Mixed geometry kinds are summed as-is, although
      the measures have different units.

  `"id"`

  :   Write the geometry id (`$runs$id`, `$edges$id`, etc.) into each
      cell. When geometries overlap, the last id wins. Background is
      `NA`.

- target:

  optional numeric `c(xmin, xmax, ymin, ymax)`; when given, materialise
  only this sub-window of the burn grid.

## Value

A numeric matrix with dimensions `nrow x ncol` of the burn grid, or of
the snapped `target` window.

## Details

The sub-window is snapped outward to the burn grid's cell boundaries
(see
[`crop_burn()`](https://hypertidy.github.io/controlledburn/reference/crop_burn.md)),
so the returned matrix covers at least `target`. Use
`crop_burn(x, target)` to recover the snapped grid's `extent` and
`dimension` (for plotting with
[`rasterImage()`](https://rdrr.io/r/graphics/rasterImage.html), say).
The materialised chunk is identical to the corresponding sub-matrix of
the full materialisation.

## Examples

``` r
if (requireNamespace("geos", quietly = TRUE)) {
  poly <- geos::as_geos_geometry("POLYGON ((1 1, 9 1, 9 9, 1 9, 1 1))")
  r <- burn(poly, extent = c(0, 10, 0, 10), dimension = c(100L, 100L))
  full <- materialise_chunk(r)
  chunk <- materialise_chunk(r, target = c(2, 5, 3, 7))
  dim(chunk)
  crop_burn(r, c(2, 5, 3, 7))$extent
}
#> [1] 2 5 3 7
```
