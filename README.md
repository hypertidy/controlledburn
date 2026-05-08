
<!-- README.md is generated from README.Rmd. Please edit that file -->

# controlledburn

<!-- badges: start -->

[![R-CMD-check](https://github.com/hypertidy/controlledburn/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/hypertidy/controlledburn/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Rasterize geometry without materializing any pixel values.
controlledburn produces sparse tables for polygon, line, and point input
— one type-pure table per geometry kind.

For polygons: run-length-encoded interior cells and boundary cells with
exact partial coverage (`fraction` in \[0, 1\]). For lines: per-cell
absolute length in CRS units. For points: per-cell records with no
measure column.

The scanline algorithm is O(perimeter) in time and memory; no dense
matrix is allocated.

## Installation

``` r
remotes::install_github("hypertidy/controlledburn")
```

## Usage

``` r
library(controlledburn)
library(geos)

poly <- as_geos_geometry("POLYGON ((1 1, 9 1, 9 9, 1 9, 1 1))")

# Sparse output — no dense matrix allocated
r <- burn_scanline(poly, extent = c(0, 10, 0, 10), dimension = c(20L, 20L))
r
#> <controlledburn> 20 x 20 grid, 1 geometry
#>   runs:   74 (256 interior cells)
#>   edges:  0 polygon boundary cells
#>   sparsity: 36.0% empty

# Materialise only when you need it
mat <- materialise_chunk(r, c(0, 10, 0, 10))
```

### Default grid parameters

With no extent or dimension, controlledburn derives both from the
geometry:

``` r
r <- burn_scanline(poly)
# extent from wk::wk_bbox(), 256 cells on the long axis
```

Or specify resolution:

``` r
r <- burn_scanline(poly, resolution = 0.5)
```

### Lines and points

Lines produce a `$lines` table with absolute length within each cell (in
CRS units, not a fraction):

``` r
line <- as_geos_geometry("LINESTRING (0 5, 10 5)")
r <- burn_scanline(line, extent = c(0, 10, 0, 10), dimension = c(20L, 20L))
r$lines
#>    row col length id
#> 1   11   1    0.5  1
#> 2   11   2    0.5  1
#> 3   11   3    0.5  1
#> 4   11   4    0.5  1
#> 5   11   5    0.5  1
#> 6   11   6    0.5  1
#> 7   11   7    0.5  1
#> 8   11   8    0.5  1
#> 9   11   9    0.5  1
#> 10  11  10    0.5  1
#> 11  11  11    0.5  1
#> 12  11  12    0.5  1
#> 13  11  13    0.5  1
#> 14  11  14    0.5  1
#> 15  11  15    0.5  1
#> 16  11  16    0.5  1
#> 17  11  17    0.5  1
#> 18  11  18    0.5  1
#> 19  11  19    0.5  1
#> 20  11  20    0.5  1
```

Points produce a `$points` table with one record per cell hit (no
measure column — a point is either in a cell or it isn’t):

``` r
pts <- as_geos_geometry(c("POINT (2 3)", "POINT (7 8)"))
r <- burn_scanline(pts, extent = c(0, 10, 0, 10), dimension = c(20L, 20L))
r$points
#>   row col id
#> 1  15   5  1
#> 2   5  15  2
```

### Geometry input

`burn_scanline()` accepts `geos_geometry`, `sfc` (sf), `wk::wkb()`,
`blob`, or a list of raw WKB vectors. The raw-WKB path is compatible
with `vapour::vapour_read_geometry()` and `gdalraster::GDALVector`
output. For `terra::vect()` input, round-trip via
`geos::as_geos_geometry()`.

### Shared boundary complementarity

Adjacent polygons with shared edges produce complementary coverage
fractions that sum to exactly 1.0 in every boundary cell — no gaps, no
overlaps:

``` r
left  <- as_geos_geometry("POLYGON ((0 0, 5 0, 5 10, 0 10, 0 0))")
right <- as_geos_geometry("POLYGON ((5 0, 10 0, 10 10, 5 10, 5 0))")

r <- burn_scanline(c(left, right), extent = c(0,10,0,10), dimension = c(20L,20L))

# Coverage sums to 1.0 in every touched cell
mat1 <- materialise_chunk(r, id = 1)
mat2 <- materialise_chunk(r, id = 2)
max(mat1 + mat2)
#> [1] 1
#> [1] 1
```

## Output format

`burn_scanline()` returns a list with class `"controlledburn"`:

- **`runs`**: `data.frame(row, col_start, col_end, id)` — polygon
  interior cells (full coverage), run-length encoded by row.
- **`edges`**: `data.frame(row, col, fraction, id)` — polygon boundary
  cells with partial coverage; `fraction` is in (0, 1).
- **`lines`**: `data.frame(row, col, length, id)` — line cells; `length`
  is the absolute length of the line within the cell, in CRS units.
- **`points`**: `data.frame(row, col, id)` — point cells; no measure
  column (a point is either in a cell or it isn’t).
- **`extent`**: `c(xmin, xmax, ymin, ymax)`
- **`dimension`**: `c(ncol, nrow)`

Tables are populated for whichever geometry kinds are in the input. Each
table’s measure column means exactly one thing — `$edges$fraction` is
dimensionless, `$lines$length` is in CRS units, points have no measure.
This separation is deliberate: the three measures are different
mathematical objects and combining them in one column would silently mix
units.

`burn_sparse()` is polygon-only and returns just `$runs` and `$edges`.
Line and point input there is rejected with an error pointing at
`burn_scanline()`.

This is the natural output of scanline rasterization — no dense matrix
is allocated until `materialise_chunk()` is called.

## Performance

Scanline algorithm scales with perimeter, not area. The comparison table
is polygon-only — `burn_sparse()` is the older bbox-bounded exactextract
path and does not accept line or point input:

| Shape       | Resolution | Scanline | Dense (burn_sparse) | Speedup |
|-------------|------------|----------|---------------------|---------|
| Star        | 3200×3200  | 13 ms    | 225 ms              | 17×     |
| Jagged      | 3200×3200  | 15 ms    | 136 ms              | 9×      |
| NC counties | 2000×800   | 29 ms    | 61 ms               | 2×      |

Memory for real-world grids (CGAZ at 32K×16K, ~500M cells): ~50 MB
sparse vs ~2 GB dense.

## Real-world example: global administrative boundaries

The CGAZ (Geo Boundaries) ADM0 dataset is roughly
<!-- TODO: a few hundred MB
of WKB across ~250 polygons of varying complexity --> a useful stress
test for both shape complexity and scale.

``` r
v <- new(gdalraster::GDALVector, "/vsicurl/https://github.com/mdsumner/geoboundaries/releases/download/latest/geoBoundariesCGAZ_ADM0.parquet")
v$returnGeomAs ## WKB
#> [1] "WKB"
gcol <- v$getGeometryColumn()
v$setIgnoredFields( setdiff(v$getFieldNames(), gcol))
wkbgeom <- wk::wkb(v$fetch(-1)[[gcol]])
v$close()

system.time(burn_scanline(wkbgeom))
#>    user  system elapsed 
#>   0.513   0.014   0.527

system.time(r1 <- burn_scanline(wkbgeom, dimension = c(8192, 4096)))
#>    user  system elapsed 
#>   0.926   0.016   0.943
str(r1)
#> List of 6
#>  $ runs     :'data.frame':   81149 obs. of  4 variables:
#>   ..$ row      : int [1:81149] 1067 1068 1069 1070 1071 1072 1073 1074 1075 1076 ...
#>   ..$ col_start: int [1:81149] 5712 5706 5706 5704 5704 5703 5702 5702 5701 5699 ...
#>   ..$ col_end  : int [1:81149] 5712 5712 5715 5717 5719 5719 5719 5719 5718 5718 ...
#>   ..$ id       : int [1:81149] 1 1 1 1 1 1 1 1 1 1 ...
#>  $ edges    :'data.frame':   469461 obs. of  4 variables:
#>   ..$ row     : int [1:469461] 1065 1066 1066 1066 1066 1066 1066 1066 1067 1067 ...
#>   ..$ col     : int [1:469461] 5712 5707 5708 5709 5710 5711 5712 5713 5705 5706 ...
#>   ..$ fraction: num [1:469461] 0.0126 0.1343 0.0792 0.296 0.1872 ...
#>   ..$ id      : int [1:469461] 1 1 1 1 1 1 1 1 1 1 ...
#>  $ lines    :'data.frame':   0 obs. of  4 variables:
#>   ..$ row   : int(0) 
#>   ..$ col   : int(0) 
#>   ..$ length: num(0) 
#>   ..$ id    : int(0) 
#>  $ points   :'data.frame':   0 obs. of  3 variables:
#>   ..$ row: int(0) 
#>   ..$ col: int(0) 
#>   ..$ id : int(0) 
#>  $ extent   : num [1:4] -180 180 -90 83.6
#>  $ dimension: int [1:2] 8192 4096
#>  - attr(*, "class")= chr "controlledburn"
```

<!-- TODO: short observation about what the timing demonstrates —
e.g. "X seconds for the full world at 8K×4K, scaling linearly with
perimeter as the resolution increases." -->

``` r
system.time(r1 <- burn_scanline(wkbgeom, dimension = c(8192, 4096) * 20))
#>    user  system elapsed 
#>  16.975   0.828  17.804
pryr::object_size(r1)
#> 278.57 MB
tibble::as_tibble(r1$runs)
#> # A tibble: 2,402,331 × 4
#>      row col_start col_end    id
#>    <int>     <int>   <int> <int>
#>  1 21300    114229  114229     1
#>  2 21301    114228  114231     1
#>  3 21302    114227  114232     1
#>  4 21303    114226  114232     1
#>  5 21304    114225  114232     1
#>  6 21305    114224  114232     1
#>  7 21306    114210  114215     1
#>  8 21306    114221  114234     1
#>  9 21307    114209  114235     1
#> 10 21308    114209  114236     1
#> # ℹ 2,402,321 more rows
tibble::as_tibble(r1$edges)
#> # A tibble: 12,006,186 × 4
#>      row    col fraction    id
#>    <int>  <int>    <dbl> <int>
#>  1 21299 114228 0.49016      1
#>  2 21299 114229 0.52762      1
#>  3 21299 114230 0.029348     1
#>  4 21300 114227 0.14377      1
#>  5 21300 114228 0.91647      1
#>  6 21300 114230 0.95528      1
#>  7 21300 114231 0.61698      1
#>  8 21300 114232 0.36007      1
#>  9 21301 114226 0.44189      1
#> 10 21301 114227 0.95691      1
#> # ℹ 12,006,176 more rows
r1[c("extent", "dimension")]
#> $extent
#> [1] -180.00000  180.00000  -90.00000   83.63339
#> 
#> $dimension
#> [1] 163840  81920
```

## History

controlledburn was derived from
[fasterize](https://github.com/ecohealthalliance/fasterize) by Noam Ross
(EcoHealth Alliance), removing Armadillo and raster package
dependencies. See [NEWS](NEWS.md) for the version history; the design
record lives in `inst/docs-design/`.

## See also

- [vaster](https://github.com/hypertidy/vaster) — primitive grid cell ↔
  xy operations; consumes the `(row, col, …)` schema this package emits.

- [silicate](https://github.com/hypertidy/silicate) — the
  primitives-first geometry stance this package follows (segments and
  vertices as first-class objects).

- [polymer2](https://github.com/hypertidy/polymer2) — sparse geometry
  overlay; consumes the `(row, col, fraction, id)` schema this package
  emits for polygons.

- [exactextractr](https://github.com/isciences/exactextractr) — raster
  extraction with polygon weights; the source of the exactextract C++
  algorithm vendored here.

- <!-- TODO: any other neighbours worth flagging for line/point
    consumers — `lazysf`, `wk`, `geos`, `gdalraster`? -->

## Code of Conduct

Please note that the controlledburn project is released with a
[Contributor Code of
Conduct](https://contributor-covenant.org/version/2/1/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
