# controlledburn: Sparse Geometry Rasterization with Exact Coverage Fractions

Rasterize polygon, line, and point geometry without materializing any
pixel values. For polygons it computes exact coverage fractions for
polygon-grid intersections, returning sparse tables: run-length encoded
interior cells and individually weighted boundary cells. The scanline
algorithm is O(perimeter) in both time and memory, using the scan line
fill attributed to Wylie et al. (1967)
[doi:10.1145/1465611.1465619](https://doi.org/10.1145/1465611.1465619) .
Exact coverage fractions use the 'exactextract' algorithm by Daniel
Baston, vendored from the 'exactextractr' package
[doi:10.32614/CRAN.package.exactextractr](https://doi.org/10.32614/CRAN.package.exactextractr)
. Derived from the 'fasterize' package by Noam Ross.

## See also

Useful links:

- <https://github.com/hypertidy/controlledburn>

- Report bugs at <https://github.com/hypertidy/controlledburn/issues>

## Author

**Maintainer**: Michael Sumner <mdsumner@gmail.com>
([ORCID](https://orcid.org/0000-0002-2471-7511))

Authors:

- Noam Ross ([ORCID](https://orcid.org/0000-0002-2136-0000))

Other contributors:

- EcoHealth Alliance (Original fasterize scanline logic) \[copyright
  holder\]

- Jim Hester (Wrote CollectorList.h for early version) \[contributor\]

- Daniel Baston (exactextract C++ algorithm (vendored)) \[copyright
  holder\]

- ISciences, LLC (exactextract C++ algorithm (vendored)) \[copyright
  holder\]
