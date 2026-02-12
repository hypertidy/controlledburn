
<!-- README.md is generated from README.Rmd. Please edit that file -->

# controlledburn

<!-- badges: start -->

[![R-CMD-check](https://github.com/hypertidy/controlledburn/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/hypertidy/controlledburn/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

Rasterize polygons without materializing any pixel values.
controlledburn computes exact coverage fractions for polygon-grid
intersections and returns results as sparse tables: run-length encoded
interior cells and individually weighted boundary cells.

The scanline algorithm is O(perimeter) in both time and memory.

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
#>   runs:  74 (256 interior cells)
#>   edges: 0 boundary cells
#>   sparsity: 36.0% empty
#> <controlledburn> 20 x 20 grid, 1 geometry
#>   runs:  N (interior cells)
#>   edges: N boundary cells

# Materialise only when you need it
mat <- materialise_chunk(r)
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

### Geometry input

Accepts `geos_geometry`, `sfc`, `wk::wkb()`, `blob`, or raw WKB list
(compatible with vapour/gdalraster):

Note that `geos::as_geos_geometry()` provides interchange for
`terra::vect()`.

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

`burn_scanline()` and `burn_sparse()` return a list with class
`"controlledburn"` containing:

- **`runs`**: `data.frame(row, col_start, col_end, id)` — interior cells
  with coverage = 1.0, run-length encoded by row.
- **`edges`**: `data.frame(row, col, weight, id)` — boundary cells with
  exact partial coverage in (0, 1).
- **`extent`**: `c(xmin, xmax, ymin, ymax)`
- **`dimension`**: `c(ncol, nrow)`

This is the natural output of scanline rasterization — no dense matrix
is allocated until `materialise_chunk()` is called.

## Performance

Scanline algorithm scales with perimeter, not area:

| Shape       | Resolution | Scanline | Dense (burn_sparse) | Speedup |
|-------------|------------|----------|---------------------|---------|
| Star        | 3200×3200  | 13 ms    | 225 ms              | 17×     |
| Jagged      | 3200×3200  | 15 ms    | 136 ms              | 9×      |
| NC counties | 2000×800   | 29 ms    | 61 ms               | 2×      |

Memory for real-world grids (CGAZ at 32K×16K, ~500M cells): ~50 MB
sparse vs ~2 GB dense.

## Fast

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
#>   0.524   0.015   0.538

system.time(r1 <- burn_scanline(wkbgeom, dimension = c(8192, 4096)))
#>    user  system elapsed 
#>   0.953   0.000   0.953
str(r1)
#> List of 4
#>  $ runs     :'data.frame':   81149 obs. of  4 variables:
#>   ..$ row      : int [1:81149] 1067 1068 1069 1070 1071 1072 1073 1074 1075 1076 ...
#>   ..$ col_start: int [1:81149] 5712 5706 5706 5704 5704 5703 5702 5702 5701 5699 ...
#>   ..$ col_end  : int [1:81149] 5712 5712 5715 5717 5719 5719 5719 5719 5718 5718 ...
#>   ..$ id       : int [1:81149] 1 1 1 1 1 1 1 1 1 1 ...
#>  $ edges    :'data.frame':   469461 obs. of  4 variables:
#>   ..$ row   : int [1:469461] 1065 1066 1066 1066 1066 1066 1066 1066 1067 1067 ...
#>   ..$ col   : int [1:469461] 5712 5707 5708 5709 5710 5711 5712 5713 5705 5706 ...
#>   ..$ weight: num [1:469461] 0.0126 0.1343 0.0792 0.296 0.1872 ...
#>   ..$ id    : int [1:469461] 1 1 1 1 1 1 1 1 1 1 ...
#>  $ extent   : num [1:4] -180 180 -90 83.6
#>  $ dimension: int [1:2] 8192 4096
#>  - attr(*, "class")= chr "controlledburn"
```

Really fast.

``` r
system.time(r1 <- burn_scanline(wkbgeom, dimension = c(8192, 4096) * 20))
#>    user  system elapsed 
#>  17.005   0.625  17.629
pryr::object_size(r1)
#> 278.56 MB
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
#>      row    col   weight    id
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
(EcoHealth Alliance), removing Armadillo and raster package dependencies
to return bare scanline indexes.

Version 0.1.0 is a complete rewrite using the exactextract algorithm
(Daniel Baston, vendored from
[exactextractr](https://github.com/isciences/exactextractr)) for exact
coverage fractions via a new O(perimeter) scanline sweep.

See also: [vaster](https://github.com/hypertidy/vaster) for grid cell
indexing, [exactextractr](https://github.com/isciences/exactextractr)
for raster extraction with polygon weights.

## Code of Conduct

Please note that the controlledburn project is released with a
[Contributor Code of
Conduct](https://contributor-covenant.org/version/2/1/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
