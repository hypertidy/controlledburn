# controlledburn (development version)

## Aligned with the controlledburn Rust crate 0.3.1

The Rust crate (hypertidy/controlledburn-rs, on crates.io) is now the
canonical implementation, and this package follows it. The R package
still builds from the C++ core, which now has the crate's options and
operations; `cpp/tests/test_canon_ops.cpp` checks each one record for
record against the crate's own answers (`fixtures/canon-ops.txt`, written
by `tools/canon-ops` from crate 0.3.1).

* **Changed defaults (records, not cells).** `burn()` now gives one record
  per cell for a multi-part geometry and merges adjacent runs, as the
  crate does by default:
  * `merge_parts = TRUE`: the parts of a MULTIPOLYGON were burned one by
    one, so a cell two parts cross had two edges (two halves either side
    of a shared border came out as `0.5 + 0.5`, not a full cell), and an
    overlap stacked an edge on a run. The parts' fractions are now summed
    per cell, a sum reaching 1 becomes a run, and an edge in a cell
    another part covers is dropped. Overlapping parts are capped at a full
    cell, not unioned. MULTILINESTRING parts are summed per cell the same
    way.
  * `coalesce_runs = TRUE`: a fully covered boundary cell was its own
    one-cell run; adjacent runs on a row are now one run.

  Materialised values are unchanged except where multipolygon parts meet
  or overlap. `merge_parts = FALSE, coalesce_runs = FALSE` gives the
  previous records exactly.

* `burn()` gains `id_offset` (ids are `id_offset` plus the position, so
  burns of slices of a collection bind together), `threads` (burn on
  several threads, with the same result record for record) and `target`
  (burn only the cells of the grid overlapping a window: the result is
  `crop_burn(burn(...), target)` but costs about the geometry near the
  window, so a large grid can be burned tile by tile). The options used
  are kept in the result's `$options`.

* New `slice_burn()`: cut a burn into a regular grid of tiles in one pass
  over its records (each tile equals `crop_burn()` to that tile), with
  `threads`.

* New `relate_cells()`: the DE-9IM matrix of each polygon against each
  grid cell its boundary touches, from a finished burn and the same
  geometry, with `intersects`, `covers`, `touches`, `interiors_meet` and
  `covered_by` columns. Every other cell is covered (in a run) or
  disjoint (no record), so together they give every cell's relation; the
  exact all-touched set of a polygon is its runs plus the relations that
  intersect. Contacts are decided with exact predicates, so cells the
  boundary only touches along a wall or at a corner, and slivers within
  the burn's 1e-6 coverage tolerance, are reported. Tested against
  `geos::geos_relate()` on every cell.

* `crop_burn()` and `materialise_chunk()` now run in C++ (the core's
  `crop()`); `materialise_chunk()` no longer loops over records in R.

* The burn now notes, as warnings, rings of fewer than 4 coordinates once
  closed (an exterior skips its part, a hole only itself), lines of one
  coordinate, and geometries with non-finite coordinates, which are
  skipped; unclosed rings are closed. Before, these were skipped silently
  or processed as given.

* C++ core: `BurnOptions` (with `parity()` for the earlier records),
  `burn_par()`/`burn_wkb_par()`, `burn_window()`/`burn_wkb_window()`,
  `coalesce_runs()`, `merge_cells()`, `crop()`, `slice_tiles()`/
  `slice_tiles_par()` (`ops.hpp`) and `relate_cells()` (`relate.hpp`), all
  ports of the crate's. A `BurnMode` passed where options are expected
  gives the crate's defaults, so C++ callers of `burn(geoms, grid, mode)`
  get coalesced runs and merged parts; use `BurnOptions::parity(mode)` for
  the previous records. The package links with `-lpthread` on Unix.

* Python: built from crate 0.3.1 (was 0.2.4), with its defaults;
  `burn()` gains `merge_parts`, `coalesce_runs` and `id_offset`.

* The Python package (`python/`) now builds from the controlledburn Rust
  crate (PyO3 + maturin) instead of this C++ core, with an unchanged API
  and identical output on the canonical tables. The R package still
  builds from the C++ core.

* The core no longer uses fused multiply-add. Every core source file
  includes `no_fma.h` first, which turns off floating-point contraction by
  pragma (clang, GCC, MSVC), so the R, Python and CMake builds all get it
  without a compiler flag (`-ffp-contract=off` in `Makevars` is a
  non-portable flag for R CMD check). Without it, clang on arm64 macOS and
  GCC on FMA-capable targets (aarch64, `-march=x86-64-v3`) moved run
  boundaries on knife-edge inputs, so results differed by platform. The
  Rust crate never fuses. Found by the canonical-table check
  (`cpp/tests/test_canon_tables.cpp`), which now passes with GCC and clang
  on FMA targets. All core translation units must agree: turning
  contraction off in the scanline file alone gave more differences, not
  fewer.

* Fixed (coverage mode): a cell was overcounted where a polygon edge passes
  exactly through a cell corner, or where two edges nearly coincide (#25,
  found by the controlledburn-rs property suite). An edge through a corner
  gave a crossing an ulp off it and a sliver traversal along the next
  cell's wall; the cell walker now exits through the corner itself. A cell
  crossed by several traversals got its coverage from exactextract's chain
  chaser, which depends on the order of the traversal ends around the cell;
  rounding can swap two nearly coincident ends and the chaser then returned
  the complement, a whole cell too much. Such cells, and a single traversal
  lying along the cell walls, now use Green's theorem: each traversal adds
  its own term, plus the cell width times the covered length of its left
  wall, read from the row. A zero-width spike out to a corner vertex and
  back no longer leaves a 1e-16 coverage record that blanks a covered cell.
  The issue's ring covered 6.6016 instead of 6.1578; a 4e-13 sliver covered
  up to 11 whole cells. A self-intersecting bowtie with its closing edge on
  a column wall now covers the exact 32 cells, not 26.

* Fixed (approx mode): which side of a cell centre an edge passes is now an
  exact orientation test, not a comparison of the rounded intercept with
  the centre. Two edges of a ring doubling back along a line, both 1e-13
  right of a centre, burned that cell when one intercept rounded onto it.

# controlledburn 0.3.2

* Fixed a GCC 16 `-Warray-bounds` warning ("array subscript -1 is outside
  array bounds") in the vendored exactextract code, reported by CRAN's
  r-devel gcc flavors. Ring closedness in `analytical_covered_fraction()` is
  now tested on the input coordinate vector, which the compiler can prove is
  non-empty, rather than on a copy it cannot. No change in behaviour or output.

* Removed the quotes around the exactextract algorithm name in the package
  Description, per the CRAN maintainers' request (the `'exactextractr'`
  package reference is retained).

# controlledburn 0.3.1

* Fixed: `materialise_chunk(fun = "sum")` assigned 1 to interior run cells
  instead of adding 1, so where polygon interiors overlapped the sum
  undercounted (and disagreed with `extract_burn(fun = "sum")`, which
  documents that it reproduces `materialise_chunk()`). Interiors now
  accumulate like edges, lines and points.

* Fixed: a non-finite extent (`Inf`, or a span that overflows double
  precision) reached the core with an infinite cell size and crashed R
  with a C stack overflow. `burn()` now errors, and the core's own grid
  validation rejects non-finite bounds for every binding.

* `?burn` no longer says that `materialise_chunk()` refuses mixed geometry
  kinds or refers to an `id =` filter it doesn't have. It sums whatever it
  is given; `extract_burn(fun = "sum")` is the one that refuses mixed
  kinds.

* The row-centre crossing fix below (#19) is now also part of the shared
  parity fixtures (`fixtures/`), so the C++, R and Python suites all check
  it.

# controlledburn 0.3.0

* Fixed: in coverage mode, a vertex lying a few ulp off a cell edge (e.g. a
  vertex at 0.3 on a grid with resolution 0.1, where the edge is computed as
  3 * 0.1 = 0.30000000000000004) sent the cell walker off the grid along the
  extrapolated edge until memory was exhausted; where it did return, a
  degenerate corner exit made `left_hand_area()` over-count the cell.
  Vertices within 1e-9 of a cell (relative to resolution) of an edge are now
  snapped onto it before the walk; a crossing computed beyond a segment's end
  is clamped to the cell; and a non-converging ring walk throws a clear error
  instead of `bad_alloc`. The same input crashes exactextractr, so this is
  inherited from the exactextract traversal. Fix lives in
  `cpp/src/scanline.cpp` and is synced to `src/`.

* Fixed: in coverage mode, an edge crossing a row's centre line exactly on
  a column wall registered no winding delta in either adjacent cell, and the
  row lost its interior runs from that column on. Pre-existing (the Falkland
  Islands from Natural Earth 110m lose two cells at 0.5 degrees), but the
  lattice snapping above puts many more crossings exactly on walls and made
  it far more visible. The crossing test is now half-open at the row centre
  in coverage mode as it already was in approx mode. Natural Earth 110m at
  0.01 degree is now 0.035 square degrees short of the planar area, versus
  0.50 before; no geometry gets worse.

* Fixed: in approx mode, a vertex sitting on a row's centre line (as a user
  computes it, one ulp from the same line recomputed from the extent) could
  flip the sign of that edge's intercept, after which the row filled from
  the polygon to the right edge of the grid. The crossing direction now
  comes from the edge's own endpoints. Found by the property-based suite in
  controlledburn-rs.

* `materialise_chunk()` gains a `target` argument: materialise only a
sub-window of the burn grid, via `crop_burn()`, without allocating the full
grid. The chunk is identical to the corresponding sub-matrix of the full
materialisation.

* Internal: the vendored C++ core now emits 0-based indices with an
  exclusive `col_end` and `id = k` (aligning with controlledburn-rs and
  every non-R consumer). The cpp11 shim (`src/scanline_shim.cpp`) adds 1
  to row/col/col_start/id to restore the package's 1-based, inclusive
  contract, so R output is unchanged.

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
