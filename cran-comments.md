## controlledburn 0.3.2

This is version 0.3.2. It fixes the issue CRAN's checks show for 0.3.1.

Verified on the r-hub `gcc16` container (g++ 16.2.1, the same toolchain as
the CRAN fedora-gcc flavor): the 0.3.1 sources reproduce the warning and
0.3.2 compiles with no warnings. (See tools/check-gcc16.sh.)

Thanks! 

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

