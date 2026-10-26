## Update

This is version 0.3.2. It fixes the one issue CRAN's checks show for 0.3.1.

## The fix

0.3.1 shows a WARN on the two GCC 16 r-devel flavors
(r-devel-linux-x86_64-debian-gcc and r-devel-linux-x86_64-fedora-gcc):

    exactextract/coordinate.h:41:31: warning: array subscript -1 is
    outside array bounds of 'exactextract::Coordinate [576460752303423487]'
    [-Warray-bounds=]

inlined from `analytical_covered_fraction()` via `walk_ring()` and
`process_polygon()`. The other eight flavors are OK, including
r-devel-linux-x86_64-fedora-clang, so it is specific to GCC 16's
-Warray-bounds analysis.

The code tested ring closedness on a copied vector; GCC 16 cannot prove the
copy is non-empty and so reads `poly.back()` as a possible index -1. The
input is already known to have at least two points, so closedness is now
tested on the input vector before copying. No behaviour change; the existing
regression tests still pass.

Verified on the r-hub `gcc16` container (g++ 16.2.1, the same toolchain as
the CRAN fedora-gcc flavor): the 0.3.1 sources reproduce the warning and
0.3.2 compiles with no warnings. (See tools/check-gcc16.sh.)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Test environments

* local: Ubuntu 24.04, R release
* GitHub Actions: macOS (release), Windows (release),
  Ubuntu (devel, release, oldrel-1)
* win-builder: R-devel
* R-hub v2: gcc16 (reproduces and confirms the fix above), plus linux,
  windows, macos and the clang-asan / valgrind containers (the package
  compiles vendored C++)

## Notes

* The package vendors a GEOS-free subset of Daniel Baston's exactextract
  C++ library (credited in Authors@R with `cph` roles and cited in the
  Description).
* There are no reverse dependencies.
