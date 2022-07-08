
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
-   [ ] remove use of raster objects, in favour of input extent and
    dimension
-   [ ] remove all trace of the raster package
-   [ ] move to cpp11
-   [ ] ensure return value/s of intermediate steps (see below) rather
    than materialized rasters
-   [ ] rasterize lines and points
    [fasterize/issues/30](https://github.com/ecohealthalliance/fasterize/issues/30)
-   [ ] consider formats other than sf (wk, geos, grd, rct, triangles
    etc.)

## Outputs

-   two options, record presence of polygon OR ID of polygon
-   a row-indexed set of edge instances (start, end) along scanlines
    with the two options
-   tools to format this meaningfully, and plot lazily
-   tools to materialize as actual raster data

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
