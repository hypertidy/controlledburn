# Deprecated: use burn() instead

`burn_sparse()` has been removed. Use
[`burn()`](https://hypertidy.github.io/controlledburn/reference/burn.md)
instead, which produces identical output for polygon input and also
supports lines and points.

## Usage

``` r
burn_sparse(...)
```

## Arguments

- ...:

  ignored

## Value

This function always errors; it never returns a value. Call
[`burn()`](https://hypertidy.github.io/controlledburn/reference/burn.md)
instead.

## Examples

``` r
# burn_sparse() is deprecated and errors; use burn() instead
try(burn_sparse())
#> Warning: burn_sparse() has been removed. Use burn() instead (identical output for polygons, plus line and point support).
#> Error in burn_sparse() : 
#>   burn_sparse() is no longer available. Use burn().
```
