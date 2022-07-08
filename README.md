
<!-- README.md is generated from README.Rmd. Please edit that file -->

# laserize

<!-- badges: start -->

[![R-CMD-check](https://github.com/hypertidy/laserize/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/hypertidy/laserize/actions/workflows/R-CMD-check.yaml)
<!-- badges: end -->

The goal of laserize is to rasterize without materializing any pixel
values.

WARNING: does not do anything yet, it’s just me learning C++. Do not
use, VERY VERY WIP

This is an expression of my WISHUE “cell abstraction”
[fasterize/issues/11](https://github.com/ecohealthalliance/fasterize/issues/11).

-   [x] copy logic from fasterize, and remove Armadillo array handling
-   [x] remove use of raster objects, in favour of input extent and
    dimension
-   [x] remove all trace of the raster package
-   [ ] implement return of the ‘yline, xpix’ and polygon ID info to
    user (see below)
-   [ ] move to cpp11
-   [ ] rasterize lines and points
    [fasterize/issues/30](https://github.com/ecohealthalliance/fasterize/issues/30)
-   [ ] consider formats other than sf (wk, geos, grd, rct, triangles
    etc.)

## Outputs

-   two options, record presence of polygon OR ID of polygon
-   a row-indexed (*yline*) set of edge instances (start, end *xpix*)
    along scanlines with the two options
-   tools to format this meaningfully, and plot lazily
-   tools to materialize as actual raster data

For the record, I wanted this facility before I read this issue - but
here’s a real world example, discussed in PROJ for very fast lookup for
large non-materialized (highly compressed) grids by Thomas Knudsen:

<https://github.com/OSGeo/PROJ/issues/1461#issuecomment-491501992>

## Installation

You can install the development version of laserize like so:

``` r
# FILL THIS IN! HOW CAN PEOPLE INSTALL YOUR DEV PACKAGE?
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(laserize)
## basic example code
```

## Code of Conduct

Please note that the laserize project is released with a [Contributor
Code of
Conduct](https://contributor-covenant.org/version/2/1/CODE_OF_CONDUCT.html).
By contributing to this project, you agree to abide by its terms.
