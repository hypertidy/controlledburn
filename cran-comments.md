## Resubmission

This is version 0.3.1. It replaces 0.3.0, which is still pending in the
new-submissions queue. After submitting 0.3.0 we found two bugs that we
did not want to reach CRAN:

* `burn()` with a non-finite extent (e.g. `Inf`) reached the C++ core with
  an infinite cell size and crashed R ("C stack overflow"). The extent is
  now validated in R and in the core.
* `materialise_chunk(fun = "sum")` undercounted cells where polygon
  interiors overlap.

Both have regression tests. No other user-facing changes beyond
documentation corrections.

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

## Test environments

* local: Ubuntu 24.04, R release
* GitHub Actions: macOS (release), Windows (release),
  Ubuntu (devel, release, oldrel-1)
* win-builder: R-devel
* R-hub v2: linux, windows, macos, plus the clang-asan and valgrind
  containers (the package compiles vendored C++)

## Notes

* The package vendors a GEOS-free subset of Daniel Baston's exactextract
  C++ library (credited in Authors@R with `cph` roles and cited in the
  Description).
* There are no reverse dependencies (new submission).
