# controlledburn 0.1.0

Complete rewrite of controlledburn using the exactextract algorithm (Daniel
Baston, vendored from exactextractr) for exact polygon-grid coverage fractions.

## What's new

* `burn_scanline()`: O(perimeter) scanline sweep with winding-number interior
  classification and exact boundary coverage fractions. No dense matrix
  allocation — output is sparse runs + edges tables.

* `burn_sparse()`: Reference implementation using the exactextract dense
  algorithm, compressed to the same sparse output format.

* `materialise_chunk()`: Opt-in expansion to dense matrix or vector, with
  per-polygon-id filtering.

* Default grid parameters: extent derived from geometry bbox via `wk::wk_bbox()`,
  dimension auto-fitted to 256 cells on the long axis preserving aspect ratio,
  or specified as `resolution`.

* Geometry input via `wk::wkb()`, `geos_geometry`, `sfc`, `blob`, or raw WKB
  list (compatible with vapour/gdalraster output).

* Moved from Rcpp to cpp11, using libgeos for GEOS access.

## Internals

* Scanline algorithm: lightweight walk using `Box::crossing()` directly (no
 `Cell` class allocation), winding-count interior classification, analytical
  single-traversal coverage via `perimeter_distance()`.

* Validated against burn_sparse across 52 test cases: simple shapes, NC
  counties, shared boundaries, edge cases (grid-aligned edges, slivers,
  extent clipping, degenerate shapes).

* O(perimeter) scaling confirmed by benchmark: 17× faster than dense at
  3200×3200 resolution for complex shapes. Memory: sparse output ~50 MB vs
  ~2 GB dense for real-world 32K×16K grids.

## History

* Derived from fasterize (Noam Ross, EcoHealth Alliance) — original scanline
  edge logic for polygon rasterization.

* Previous controlledburn (0.0.2) returned binary in/out run-length indexes
  without coverage fractions. That code is archived at
  `archive0-2026-02-13` branch.

# controlledburn 0.0.2

* Basic function working, return a list of triplets of zero-index
  start,end,line index.

* Raster package objects and use removed.

* Converted from fasterize edge scanline logic, to not materialize any
  raster data. Armadillo usage removed.
