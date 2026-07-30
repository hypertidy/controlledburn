# controlledburn core: status and next steps

## Where it stands

Branches `cpp-core` and `approx-mode` (on top of `cpp-core`) reconcile
controlledburn and fasterize around a single engine with two burn modes.

### C++ core (`cpp/`)

Pure C++17 library, zero external dependencies.

- **Scanline engine** (`src/scanline.cpp`): walker, ring processing,
  winding sweep. O(perimeter) time and memory.
- **WKB reader** (`src/wkb.cpp`): ~230-line zero-dependency reader
  (ISO + EWKB, both byte orders, Z/M skipped).
- **exactextract subset** (`src/ee/`): GEOS-free analytical coverage
  math vendored from exactextract.
- **`materialize.hpp`**: optional consumer — fasterize pixel-function
  semantics (Sum, First, Last, Min, Max, Count, Any) over the sparse
  output. EdgePolicy controls boundary cell treatment when
  materializing.

No GEOS, no Armadillo, no R, no materialized raster in the core.

#### `BurnMode` enum (`burn.hpp`)

Two modes, selectable per burn call:

- **`Coverage`** (default): exact coverage fractions via analytical
  traversal. Boundary cells appear in `$edges` with a fraction in
  (0, 1). Total burned coverage equals exact polygon area.
- **`Approx`**: cell-center rule (fasterize semantics). Skips the
  exactextract traversal math entirely. For each boundary cell,
  interpolates the polygon edge's x-crossing at the row's y-midpoint
  and classifies the cell as inside iff the crossing is left of (or at)
  the cell center x. Inside cells become runs; outside cells are
  dropped. No `$edges` for polygons. Left-inclusive boundary convention.

Lines and points are unaffected by mode.

### R package

- **`burn()`** is the main entry point (renamed from `burn_scanline()`).
  Signature: `burn(x, extent, dimension, resolution, mode)` where
  `mode = c("coverage", "approx")`.
- **`burn_scanline()`** remains as a deprecated wrapper calling
  `burn(..., mode = "coverage")`.
- **`burn_sparse()`**: older polygon-only path using GEOS + full
  exactextract dense intermediate. The only remaining GEOS user.
  Deprecation candidate.
- `cpp_scanline_burn` is a thin cpp11 shim over the core
  (`src/scanline_shim.cpp`). `tools/sync-core.sh` derives the build
  copies from canonical `cpp/`: public headers to
  `inst/include/controlledburn/` (downstream `LinkingTo: controlledburn`
  works), core sources to `src/core_*.cpp` compiled against the
  existing `src/exactextract/` objects.
- 196 tests pass (0 failures, 0 warnings, 1 expected skip).

### Python bindings (`python/`)

pybind11 + scikit-build-core. `burn()` takes WKB bytes (shapely
`to_wkb` works directly), returns numpy structured arrays;
`materialize()` returns dense 2D arrays. Python-native argument
conventions (`bounds` rasterio-style, `shape` numpy-style); identical
1-based table values to R. `mode` parameter supports both
`"coverage"` and `"approx"`. 33 pytest cases pass (including 6
approx mode tests).

Install: `pip install git+https://github.com/hypertidy/controlledburn#subdirectory=python`

## CI

Three GitHub Actions workflows, all triggering on `main`, `master`,
and `cpp-core` branches:

- **C++ core tests** (`.github/workflows/cpp-core-test.yaml`): cmake +
  ctest on ubuntu-latest and macos-latest.
- **Python tests** (`.github/workflows/python-test.yaml`): pytest on
  ubuntu-latest and macos-latest, Python 3.11 and 3.12.
- **R CMD check** (`.github/workflows/R-CMD-check-check.yaml`):
  macOS-latest (release) and ubuntu-latest (release).

## Parity fixtures

Shared test geometries and expected results live in `fixtures/`:

- **`geometries.csv`**: WKT, WKB hex, and grid spec per test case.
- **`expected.csv`**: expected aggregate results (covered area, edge
  emptiness, line length, point count, tolerance).

Fixture-reading tests exist for all three surfaces:
- C++ (`cpp/tests/test_parity_fixtures.cpp`)
- R (`tests/testthat/test-parity-fixtures.R`)
- Python (`python/tests/test_parity_fixtures.py`)

## Performance

### Approx mode: two code paths

Approx mode has a dedicated lightweight sweep
(`process_polygon_approx`) that bypasses the exactextract walker
entirely. For each polygon edge, it computes x-intercepts at each
row's y_mid, accumulates winding per row, and sweeps left-to-right
to emit runs. No CellRecords, no traversal coordinates, no
BoundaryCellRecords — just edge math and winding. ~120 lines of C++.

Coverage mode still uses the full walker for exact coverage fractions.

### Scaling: controlledburn vs fasterize

Benchmarked on CGAZ (218 country polygons, 10.1M vertices) at
increasing grid resolutions. controlledburn produces sparse output
(runs only in approx mode); fasterize allocates a dense raster.

| Grid | Cells | cb approx | fasterize | Winner |
|------|-------|-----------|-----------|--------|
| 256×128 | 33K | 1.0s | 0.5s | fasterize |
| 1024×512 | 524K | 1.0s | 0.3s | fasterize |
| 4096×2048 | 8.4M | 1.1s | 0.4s | fasterize |
| **16384×8192** | **134M** | **1.2s** | **2.8s** | **cb** |
| 32768×16384 | 537M | 1.4s | 9.7s | cb |
| 65536×32768 | 2.1B | 1.7s | OOM | cb only |
| 131072×65536 | 8.6B | 2.4s | OOM | cb only |

**Crossover at ~134 million cells.** Below that, fasterize's simpler
per-cell cost wins (2–3×). Above it, controlledburn's O(perimeter)
memory and the lightweight sweep win — at 537M cells, cb is 6.8×
faster; fasterize OOMs above ~2B cells while cb reaches 8.6B in
2.4 seconds.

The ~1.0s floor at small grids is the cost of 10M edge–row
intersections. The sweep scales primarily with polygon perimeter
(total edge count × rows touched), barely increasing with grid
resolution (1.0s → 2.4s across a 260× increase in cell count).

### Fasterize parity (NC counties, approx mode)

At 2000×800 on North Carolina's 100 counties: 14 discrepant cells
out of 841,558 (0.002%), all at polygon boundaries where an edge
grazes a cell center at floating-point precision. Consistent with
the boundary handling uncertainties documented in
[hypertidy/fasterize#6](https://github.com/hypertidy/fasterize/issues/6)
and [OSGeo/gdal#14615](https://github.com/OSGeo/gdal/issues/14615).
On constructed geometries (aligned rectangles, offset rectangles,
triangles, holes), approx mode is cell-for-cell identical to fasterize.

## Validation

- **Core**: 17 tests — 12 original invariant tests (area conservation
  across rectangles, triangles, holes, multipolygons, beyond-extent,
  lines, points, WKB round trips, materialize, degenerate input) plus
  5 Approx mode tests (aligned rect, offset rect, beyond-extent, hole,
  lines unchanged).
- **Parity fixtures**: 10 cross-language cases driven by shared CSV.
- **R**: 218 tests pass, including 22 approx-mode tests (output
  contract, geometry cases, 5 fasterize parity comparisons).
  NC county fixture bundled as WKB (no sf/vapour needed at test time).
  Malformed WKB warns-and-skips rather than hard GEOS error.
- **Python**: 33 pytest cases including 6 approx-mode tests;
  cross-language parity with R confirmed.

## Next steps

1. **Port target/snap/clamp into `materialize.hpp`** using the
   resurrected `test-materialise.R` as the spec (chunked windowed reads
   for Python too).
4. **Decide `burn_sparse` deprecation**: the only remaining GEOS user.
   Removing it drops libgeos entirely and shrinks vendored exactextract
   to the nine GEOS-free files.
5. **Python sdist/PyPI**: vendor the core into `python/` (sync step
   like the R side) when a first release is wanted.
6. **fasterize convergence**: point fasterize at the core (engine +
   pixel-function consumer both exist now); approx mode provides the
   center-rule sweep.
7. Later: Rust binding candidate (same WKB-in, tables-out contract).
