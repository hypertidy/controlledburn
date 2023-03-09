
<!-- README.md is generated from README.Rmd. Please edit that file -->

# controlledburn

<!-- badges: start -->

[![R-CMD-check](https://github.com/hypertidy/controlledburn/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/hypertidy/controlledburn/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of controlledburn is to rasterize without materializing any
pixel values.

A very fast rasterization algorithm for polygons in
[fasterize](https://github.com/ecohealthalliance/fasterize) does the
following:

- restructures polygons to edges
- considers only non-horizontal polygon edges and indexes them in y,x
  order (slope always down)
- determines all affected raster rows and scans these for edge *start*
  and *end* column
- burns polygon value or identity into every pixel between those starts
  and ends

controlledburn avoids that last step. With these rasterization functions
the return value is not a set of pixels or a side-effect of a
materialized raster file but simply those row and column start and end
indexes.

This is an expression of the *cell abstraction* wish item
[fasterize/issues/11](https://github.com/ecohealthalliance/fasterize/issues/11).

## TODO

- [ ] move to cpp11
- [ ] rasterize lines
  [fasterize/issues/30](https://github.com/ecohealthalliance/fasterize/issues/30)
- [ ] formats for import (wk, geos, grd, rct, triangles etc.)
- [ ] streaming with wkb/xy unpack with wk
- [ ] provide output options (see next section)
- [ ] port back into fasterize, with options for efficiently writing out
  to a tiled and sparse GeoTIFF
- [x] points is too easy, see vaster::cell_from_xy
- [x] name the package
- [x] copy logic from fasterize, and remove Armadillo array handling
- [x] remove use of raster objects, in favour of input extent and
  dimension
- [x] remove all trace of the raster package
- [x] implement return of the ‘yline, xpix’ and polygon ID info to user
  (see below)
- [x] make return of ylin,xpix structure efficient (CollectorList.h ftw)

Lines is working but still has some problems. For polygons it’s pretty
good, see [fasterize
\#6](https://github.com/ecohealthalliance/fasterize/issues/6) for a
remaining issue.

## Outputs

Internal function `burn_polygon()` has arguments `sf`, `extent`,
`dimension`. The first is a sf polygons data frame, extent is
`c(xmin, xmax, ymin, ymax)` and dimension is `c(ncol, nrow)`. In
fasterize you need an actual non-materialized *RasterLayer* object, but
all that was really used for was the six numbers extent, dimension and
as the shell for the in-memory output.

The output of `burn_polygon()` is a list of four-element indexes
`start,end,row,poly_id` - these are zero-based atm because they reflect
the underlying C++ code. Examples shown here flatten this to a 4-column
matrix (and add 1).

These options are still in play for what the interface/s could do:

- record presence of polygon (this is all we have currently) OR ID of
  polygon
- tools to format this meaningfully, and plot lazily (see example for
  quick plot)
- tools to materialize as actual raster data

I also found this real world example for which this is relevant,
discussed in PROJ for very fast lookup for large non-materialized
(highly compressed) grids by Thomas Knudsen:

<https://github.com/OSGeo/PROJ/issues/1461#issuecomment-491501992>

## Installation

You can install the development version of controlledburn like so:

``` r
remotes::install_github("hypertidy/controlledburn")
```

## Example

This is a basic example, this is fast, and shows that it works.

``` r
pols <- silicate::inlandwaters
library(vaster)
#> 
#> Attaching package: 'vaster'
#> The following object is masked from 'package:stats':
#> 
#>     ts
## define a raster (xmin, xmax, ymin, ymax), (ncol, nrow)
ext <- unlist(lapply(silicate::sc_vertex(pols), range))
dm <- c(50, 40)
r <- controlledburn:::burn_polygon(pols, extent = ext,
                               dimension = dm)

## our index is triplets of start,end,line where the polygon edge was detected - 
## this essentially an rle by scanline of start,end polygon coverage
index <- matrix(unlist(r, use.names = F), ncol = 4L, byrow = TRUE) + 1 ## plus one because 0-index internally

## plot just the start and ends of each scanline detected
xy0 <- vaster::xy_from_cell(dm, ext, vaster::cell_from_row_col(dm, index[,c(3, 3)], index[,1:2]))
plot(xy0)
plot(silicate::PATH0(pols), add = TRUE)
```

<img src="man/figures/README-example-1.png" width="100%" />

``` r

## expand out to every cell
  
cr <- do.call(rbind, apply(index, 1,\(.x) cbind(seq(.x[1], .x[2]), .x[3], .x[4])))
xy <- vaster::xy_from_cell(dm, ext, vaster::cell_from_row_col(dm, cr[,2], cr[,1]))
plot(xy, pch = 19, cex = .3, col = palr::d_pal(cr[,3]))
plot(silicate::PATH0(pols), add = TRUE)
```

<img src="man/figures/README-example-2.png" width="100%" />

It scales to very large tasks, with small output.

``` r
dm <- c(500000, 400000)
system.time(r <- controlledburn:::burn_polygon(pols, extent = ext,
                               dimension = dm))
#>    user  system elapsed 
#>   0.913   0.068   0.982
length(r)
#> [1] 989153

## consider a prod(dm) raster of type double (or even bool) obviously would compress again but why churn
pryr::object_size(r)  
#> 71.22 MB
```

The following is inefficient, but shows that we get the right result.

``` r
dm <- c(500, 400)
system.time(r <- controlledburn:::burn_polygon(pols, extent = ext,
                               dimension = dm))
#>    user  system elapsed 
#>   0.002   0.001   0.003

index <- matrix(unlist(r, use.names = F), ncol = 4L, byrow = TRUE) + 1 ## plus one because 0-index internally

## now go inefficient, this is every column,row index, then converted to cell, converted to xy
cr <- do.call(rbind, apply(index, 1, \(.x) cbind(seq(.x[1], .x[2]), .x[3], .x[4])))
xy <- vaster::xy_from_cell(dm, ext, vaster::cell_from_row_col(dm, cr[,2], cr[,1]))
plot(xy, pch = ".", col = palr::d_pal(cr[,3]))
```

<img src="man/figures/README-bad-1.png" width="100%" />

``` r

rr <- terra::rast(cbind(xy, 0), type = "xyz")

rr[terra::cellFromXY(rr, xy)] <- 1
terra::plot(rr ,col = "firebrick", asp = NA)
plot(silicate::SC0(pols), add = TRUE)
```

<img src="man/figures/README-bad-2.png" width="100%" />

## Code of Conduct

Please note that the controlledburn project is released with a
[Contributor Code of
Conduct](https://contributor-covenant.org/version/2/1/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
