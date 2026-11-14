# Spatial predicates of polygons against grid cells

Computes the DE-9IM intersection matrix of each polygon against each
grid cell its boundary touches, from a finished
[`burn()`](https://hypertidy.github.io/controlledburn/reference/burn.md)
and the same geometry. With it, every cell's relation to every polygon
is known: a cell the boundary does not reach is either wholly inside the
polygon (it is in one of the polygon's runs: matrix `"212FF1FF2"`,
covered) or wholly outside (no record: `"FF2FF1212"`, disjoint).

## Usage

``` r
relate_cells(x, burn)
```

## Arguments

- x:

  the geometry that was burned (any input
  [`burn()`](https://hypertidy.github.io/controlledburn/reference/burn.md)
  accepts).

- burn:

  the `controlledburn` result of `burn(x, ...)`, in either mode and with
  any options. Only its runs, grid and `id_offset` are used. A result
  cropped with
  [`crop_burn()`](https://hypertidy.github.io/controlledburn/reference/crop_burn.md)
  or sliced with
  [`slice_burn()`](https://hypertidy.github.io/controlledburn/reference/slice_burn.md)
  gives the relations on its window's grid.

## Value

A data.frame with one row per (polygon, cell) the polygon's boundary
touches, sorted by `id`, `row`, `col`:

- `row`, `col`:

  the cell, 1-based, row 1 at the top.

- `id`:

  the geometry, as in the burn.

- `matrix`:

  the DE-9IM matrix as nine characters, the polygon's interior, boundary
  and exterior (rows) against the cell's (columns).

- `intersects`:

  they share a point.

- `interiors_meet`:

  the polygon covers part of the cell's area.

- `covers`:

  no point of the cell is outside the polygon.

- `touches`:

  they share a point but no area.

- `covered_by`:

  no point of the polygon is outside the cell.

## Details

The cells returned include those the boundary only touches, along a wall
or at a corner, which have no record in the burn. So an exact "all
touched" set of a polygon is its runs plus the relations with
`intersects`, and the cells it covers are its runs plus the relations
with `covers`. Contacts are decided with exact predicates on the input
coordinates, not from the burn's coverage fractions: the burn treats
coverage within 1e-6 of 0 or 1 as 0 or 1, so a cell holding a sliver of
a polygon may have no record, or be in a run, and is still reported here
as touched or partly covered.

Polygons are assumed valid (rings that do not cross, holes inside their
shell, the parts of a multipolygon meeting only at points). Lines,
points and empty geometries give no relations.

## Examples

``` r
if (requireNamespace("geos", quietly = TRUE)) {
  ## a square on cell walls: its 4 cells are covered, and the 12 around
  ## it touch it without sharing any area
  sq <- geos::as_geos_geometry("POLYGON ((1 1, 3 1, 3 3, 1 3, 1 1))")
  r <- burn(sq, extent = c(0, 4, 0, 4), dimension = c(4L, 4L))
  rel <- relate_cells(sq, r)
  table(rel$matrix)
  sum(rel$touches)

  ## exact all-touched cells of a polygon: its runs plus the cells
  ## its boundary meets
  poly <- geos::as_geos_geometry("POLYGON ((0.5 0.5, 2.5 0.7, 2 3.2, 0.5 0.5))")
  r <- burn(poly, extent = c(0, 4, 0, 4), dimension = c(4L, 4L))
  rel <- relate_cells(poly, r)
  rel[rel$intersects, c("row", "col")]
}
#>    row col
#> 1    1   2
#> 2    1   3
#> 3    2   2
#> 4    2   3
#> 5    3   1
#> 6    3   2
#> 7    3   3
#> 8    4   1
#> 9    4   2
#> 10   4   3
```
