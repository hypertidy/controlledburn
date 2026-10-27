# Changelog

## controlledburn 0.3.2

- Fixed a GCC 16 `-Warray-bounds` warning (“array subscript -1 is
  outside array bounds”) in the vendored exactextract code, reported by
  CRAN’s r-devel gcc flavors. Ring closedness in
  `analytical_covered_fraction()` is now tested on the input coordinate
  vector, which the compiler can prove is non-empty, rather than on a
  copy it cannot. No change in behaviour or output.

- Removed the quotes around the exactextract algorithm name in the
  package Description, per the CRAN maintainers’ request (the
  `'exactextractr'` package reference is retained).

## controlledburn 0.3.1

CRAN release: 2026-10-24

- Fixed: `materialise_chunk(fun = "sum")` assigned 1 to interior run
  cells instead of adding 1, so where polygon interiors overlapped the
  sum undercounted (and disagreed with `extract_burn(fun = "sum")`,
  which documents that it reproduces
  [`materialise_chunk()`](https://hypertidy.github.io/controlledburn/reference/materialise_chunk.md)).
  Interiors now accumulate like edges, lines and points.

- Fixed: a non-finite extent (`Inf`, or a span that overflows double
  precision) reached the core with an infinite cell size and crashed R
  with a C stack overflow.
  [`burn()`](https://hypertidy.github.io/controlledburn/reference/burn.md)
  now errors, and the core’s own grid validation rejects non-finite
  bounds for every binding.

- [`?burn`](https://hypertidy.github.io/controlledburn/reference/burn.md)
  no longer says that
  [`materialise_chunk()`](https://hypertidy.github.io/controlledburn/reference/materialise_chunk.md)
  refuses mixed geometry kinds or refers to an `id =` filter it doesn’t
  have. It sums whatever it is given; `extract_burn(fun = "sum")` is the
  one that refuses mixed kinds.

- The row-centre crossing fix below
  ([\#19](https://github.com/hypertidy/controlledburn/issues/19)) is now
  also part of the shared parity fixtures (`fixtures/`), so the C++, R
  and Python suites all check it.

## controlledburn 0.3.0

- Fixed: in coverage mode, a vertex lying a few ulp off a cell edge
  (e.g. a vertex at 0.3 on a grid with resolution 0.1, where the edge is
  computed as 3 \* 0.1 = 0.30000000000000004) sent the cell walker off
  the grid along the extrapolated edge until memory was exhausted; where
  it did return, a degenerate corner exit made `left_hand_area()`
  over-count the cell. Vertices within 1e-9 of a cell (relative to
  resolution) of an edge are now snapped onto it before the walk; a
  crossing computed beyond a segment’s end is clamped to the cell; and a
  non-converging ring walk throws a clear error instead of `bad_alloc`.
  The same input crashes exactextractr, so this is inherited from the
  exactextract traversal. Fix lives in `cpp/src/scanline.cpp` and is
  synced to `src/`.

- Fixed: in coverage mode, an edge crossing a row’s centre line exactly
  on a column wall registered no winding delta in either adjacent cell,
  and the row lost its interior runs from that column on. Pre-existing
  (the Falkland Islands from Natural Earth 110m lose two cells at 0.5
  degrees), but the lattice snapping above puts many more crossings
  exactly on walls and made it far more visible. The crossing test is
  now half-open at the row centre in coverage mode as it already was in
  approx mode. Natural Earth 110m at 0.01 degree is now 0.035 square
  degrees short of the planar area, versus 0.50 before; no geometry gets
  worse.

- Fixed: in approx mode, a vertex sitting on a row’s centre line (as a
  user computes it, one ulp from the same line recomputed from the
  extent) could flip the sign of that edge’s intercept, after which the
  row filled from the polygon to the right edge of the grid. The
  crossing direction now comes from the edge’s own endpoints. Found by
  the property-based suite in controlledburn-rs.

- [`materialise_chunk()`](https://hypertidy.github.io/controlledburn/reference/materialise_chunk.md)
  gains a `target` argument: materialise only a sub-window of the burn
  grid, via
  [`crop_burn()`](https://hypertidy.github.io/controlledburn/reference/crop_burn.md),
  without allocating the full grid. The chunk is identical to the
  corresponding sub-matrix of the full materialisation.

- Internal: the vendored C++ core now emits 0-based indices with an
  exclusive `col_end` and `id = k` (aligning with controlledburn-rs and
  every non-R consumer). The cpp11 shim (`src/scanline_shim.cpp`) adds 1
  to row/col/col_start/id to restore the package’s 1-based, inclusive
  contract, so R output is unchanged.

## controlledburn 0.2.0

### Breaking changes

- **[`burn()`](https://hypertidy.github.io/controlledburn/reference/burn.md)
  replaces
  [`burn_scanline()`](https://hypertidy.github.io/controlledburn/reference/burn.md)**
  as the main entry point.
  [`burn_scanline()`](https://hypertidy.github.io/controlledburn/reference/burn.md)
  remains as a deprecated wrapper.

- **[`burn_sparse()`](https://hypertidy.github.io/controlledburn/reference/burn_sparse.md)
  removed.** The function is now a deprecated stub that errors with a
  message directing users to
  [`burn()`](https://hypertidy.github.io/controlledburn/reference/burn.md).
  All GEOS-dependent code has been removed; the package no longer
  depends on `libgeos`.

- Polygon boundary cells: `$edges$weight` renamed to `$edges$fraction`
  (done in a previous dev cycle, retained here for completeness).

### What’s new

- [`extract_burn()`](https://hypertidy.github.io/controlledburn/reference/extract_burn.md)
  to extract values at points

- **`mode = c("coverage", "approx")`** parameter on
  [`burn()`](https://hypertidy.github.io/controlledburn/reference/burn.md).
  `"coverage"` (default) computes exact coverage fractions for polygon
  boundary cells via the analytical traversal engine. `"approx"` uses
  the cell-centre rule (fasterize semantics): boundary cells are
  included as full run cells iff the cell centre is inside the polygon.
  No `$edges` produced in approx mode. Lines and points are unaffected
  by mode.

- **`crop_burn(x, target)`** extracts a sub-window from a
  `controlledburn` result. Filters and clips runs/edges/lines/points to
  the target extent (snapped outward to cell boundaries), re-bases
  row/col indices. Pure R data frame filtering — no dense allocation.

- **Lightweight approx sweep.** Approx mode uses a dedicated edge-row
  intersection sweep (~120 lines of C++) that bypasses the exactextract
  walker entirely. On CGAZ (218 countries, 10.1M vertices), this is
  1.6–14.5× faster than the walker-based path and crosses over fasterize
  at ~134 million cells. At 537M cells, controlledburn is 6.8× faster
  than fasterize; at 8.6 billion cells (where fasterize runs out of
  memory), the sweep completes in 2.4 seconds.

- **Cell-for-cell fasterize parity.** On NC counties at 2000×800, the
  lightweight sweep produces identical output to fasterize (zero
  discrepant cells). Boundary conventions match fasterize:
  left-inclusive, top-inclusive, half-open interval for horizontal edges
  at cell-centre y-coordinates.

- **C++ core** (`cpp/`): pure C++17 library with zero external
  dependencies, built and tested independently via CMake/CTest.
  `tools/sync-core.sh` derives R package sources from the canonical
  core.

- **Python bindings** (`python/`): pybind11 + scikit-build-core.
  [`burn()`](https://hypertidy.github.io/controlledburn/reference/burn.md)
  accepts WKB bytes and supports `mode="coverage"` and `mode="approx"`.
  33 pytest cases pass.

- **Shared parity fixtures** (`fixtures/`): CSV files with WKT, WKB hex,
  grid specs, and expected results. Read by C++, R, and Python test
  suites for cross-language consistency.

- **CI** for all three surfaces: C++ ctest, Python pytest, R CMD check.

### Dependencies

- **Dropped**: `libgeos` (LinkingTo and Imports). The vendored
  exactextract subset is trimmed to 9 GEOS-free analytical geometry
  files.
- **Retained**: `cpp11` (LinkingTo), `wk` (Imports).
- Test dependencies (`sf`, `fasterize`) are optional and gated behind
  `skip_if_not_installed()`. NC county test fixture is bundled as
  serialised WKB — no `sf` or `vapour` needed at test time.

### Internals

- Unified geometry rasterization:
  [`burn()`](https://hypertidy.github.io/controlledburn/reference/burn.md)
  handles polygon, line, and point input through a single entry point
  with type-pure output tables (`$runs`, `$edges`, `$lines`, `$points`).

- `process_polygon_approx()`: the lightweight sweep. For each polygon
  edge, computes x-intercepts at each row’s y_mid, accumulates winding
  per row, sweeps left-to-right to emit runs. Top-inclusive half-open
  interval `(ya, yb]` matches the walker’s crossing convention.

- `process_polygon()`: the full walker path for coverage mode. Uses the
  exactextract traversal engine for exact analytical coverage fractions.

- Edge zoo (`tests/testthat/test-edge-zoo.R`) pins canonical rasterizer
  edge cases by category. Each test pins a *convention*, not just a
  numerical expectation.

- [`vignette("architecture")`](https://hypertidy.github.io/controlledburn/articles/architecture.md)
  documents the full development story: fasterize origins, exactextract
  integration, GEOS replacement, dual-mode engine, and scaling
  characteristics.

## controlledburn 0.1.0

Complete rewrite of controlledburn using the exactextract algorithm
(Daniel Baston, vendored from exactextractr) for exact polygon-grid
coverage fractions.

### What’s new

- [`burn_scanline()`](https://hypertidy.github.io/controlledburn/reference/burn.md):
  O(perimeter) scanline sweep with winding-number interior
  classification and exact boundary coverage fractions. No dense matrix
  allocation — output is sparse runs + edges tables.

- [`burn_sparse()`](https://hypertidy.github.io/controlledburn/reference/burn_sparse.md):
  Reference implementation using the exactextract dense algorithm,
  compressed to the same sparse output format.

- [`materialise_chunk()`](https://hypertidy.github.io/controlledburn/reference/materialise_chunk.md):
  Opt-in expansion to dense matrix or vector, with per-polygon-id
  filtering.

- Default grid parameters: extent derived from geometry bbox via
  [`wk::wk_bbox()`](https://paleolimbot.github.io/wk/reference/wk_bbox.html),
  dimension auto-fitted to 256 cells on the long axis preserving aspect
  ratio, or specified as `resolution`.

- Geometry input via
  [`wk::wkb()`](https://paleolimbot.github.io/wk/reference/wkb.html),
  `geos_geometry`, `sfc`, `blob`, or raw WKB list (compatible with
  vapour/gdalraster output).

- Moved from Rcpp to cpp11, using libgeos for GEOS access.

### Internals

- Scanline algorithm: lightweight walk using `Box::crossing()` directly
  (no `Cell` class allocation), winding-count interior classification,
  analytical single-traversal coverage via `perimeter_distance()`.

- Validated against burn_sparse across 52 test cases: simple shapes, NC
  counties, shared boundaries, edge cases (grid-aligned edges, slivers,
  extent clipping, degenerate shapes).

- O(perimeter) scaling confirmed by benchmark: 17× faster than dense at
  3200×3200 resolution for complex shapes. Memory: sparse output ~50 MB
  vs ~2 GB dense for real-world 32K×16K grids.

### History

- Derived from fasterize (Noam Ross, EcoHealth Alliance) — original
  scanline edge logic for polygon rasterization.

- Previous controlledburn (0.0.2) returned binary in/out run-length
  indexes without coverage fractions. That code is archived at
  `archive0-2026-02-13` branch.

## controlledburn 0.0.2

- Basic function working, return a list of triplets of zero-index
  start,end,line index.

- Raster package objects and use removed.

- Converted from fasterize edge scanline logic, to not materialize any
  raster data. Armadillo usage removed.
