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
