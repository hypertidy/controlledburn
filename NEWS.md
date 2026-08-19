# controlledburn 0.2.0

## Breaking changes

* **`burn()` replaces `burn_scanline()`** as the main entry point.
  `burn_scanline()` remains as a deprecated wrapper.

* **`burn_sparse()` removed.** The function is now a deprecated stub
  that errors with a message directing users to `burn()`. All
  GEOS-dependent code has been removed; the package no longer depends
  on `libgeos`.

* Polygon boundary cells: `$edges$weight` renamed to `$edges$fraction`
  (done in a previous dev cycle, retained here for completeness).

## What's new

* `extract_burn()` to extract values at points

* **`mode = c("coverage", "approx")`** parameter on `burn()`.
  `"coverage"` (default) computes exact coverage fractions for polygon
  boundary cells via the analytical traversal engine. `"approx"` uses
  the cell-centre rule (fasterize semantics): boundary cells are
  included as full run cells iff the cell centre is inside the polygon.
  No `$edges` produced in approx mode. Lines and points are unaffected
  by mode.

* **`crop_burn(x, target)`** extracts a sub-window from a
  `controlledburn` result. Filters and clips runs/edges/lines/points
  to the target extent (snapped outward to cell boundaries), re-bases
  row/col indices. Pure R data frame filtering — no dense allocation.

* **Lightweight approx sweep.** Approx mode uses a dedicated
  edge-row intersection sweep (~120 lines of C++) that bypasses the
  exactextract walker entirely. On CGAZ (218 countries, 10.1M
  vertices), this is 1.6–14.5× faster than the walker-based path and
  crosses over fasterize at ~134 million cells. At 537M cells,
  controlledburn is 6.8× faster than fasterize; at 8.6 billion cells
  (where fasterize runs out of memory), the sweep completes in 2.4
  seconds.

* **Cell-for-cell fasterize parity.** On NC counties at 2000×800, the
  lightweight sweep produces identical output to fasterize (zero
  discrepant cells). Boundary conventions match fasterize:
  left-inclusive, top-inclusive, half-open interval for horizontal
  edges at cell-centre y-coordinates.

* **C++ core** (`cpp/`): pure C++17 library with zero external
  dependencies, built and tested independently via CMake/CTest.
  `tools/sync-core.sh` derives R package sources from the canonical
  core.

* **Python bindings** (`python/`): pybind11 + scikit-build-core.
  `burn()` accepts WKB bytes and supports `mode="coverage"` and
  `mode="approx"`. 33 pytest cases pass.

* **Shared parity fixtures** (`fixtures/`): CSV files with WKT, WKB
  hex, grid specs, and expected results. Read by C++, R, and Python
  test suites for cross-language consistency.

* **CI** for all three surfaces: C++ ctest, Python pytest, R CMD check.

## Dependencies

* **Dropped**: `libgeos` (LinkingTo and Imports). The vendored
  exactextract subset is trimmed to 9 GEOS-free analytical geometry
  files.
* **Retained**: `cpp11` (LinkingTo), `wk` (Imports).
* Test dependencies (`sf`, `fasterize`) are optional and gated behind
  `skip_if_not_installed()`. NC county test fixture is bundled as
  serialised WKB — no `sf` or `vapour` needed at test time.

## Internals

* Unified geometry rasterization: `burn()` handles polygon, line, and
  point input through a single entry point with type-pure output
  tables (`$runs`, `$edges`, `$lines`, `$points`).

* `process_polygon_approx()`: the lightweight sweep. For each polygon
  edge, computes x-intercepts at each row's y_mid, accumulates winding
  per row, sweeps left-to-right to emit runs. Top-inclusive half-open
  interval `(ya, yb]` matches the walker's crossing convention.

* `process_polygon()`: the full walker path for coverage mode. Uses
  the exactextract traversal engine for exact analytical coverage
  fractions.

* Edge zoo (`tests/testthat/test-edge-zoo.R`) pins canonical
  rasterizer edge cases by category. Each test pins a *convention*,
  not just a numerical expectation.

* `vignette("architecture")` documents the full development story:
  fasterize origins, exactextract integration, GEOS replacement,
  dual-mode engine, and scaling characteristics.

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
