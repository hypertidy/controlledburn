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
1-based table values to R. 17 pytest cases pass.

Install: `pip install git+https://github.com/hypertidy/controlledburn#subdirectory=python`

**Note:** `mode` parameter not yet wired through Python bindings.

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

## Validation

- **Core**: 17 tests — 12 original invariant tests (area conservation
  across rectangles, triangles, holes, multipolygons, beyond-extent,
  lines, points, WKB round trips, materialize, degenerate input) plus
  5 Approx mode tests (aligned rect, offset rect, beyond-extent, hole,
  lines unchanged).
- **Parity fixtures**: 10 cross-language cases driven by shared CSV.
- **R**: 196 tests pass. Shim parity with the pre-shim engine
  confirmed. Malformed WKB warns-and-skips rather than hard GEOS error.
- **Python**: 17 pytest cases; cross-language parity with R confirmed.

## Known limitations

- **Approx mode horizontal-edge boundary rows**: when a polygon
  boundary is exactly at a cell center's y-coordinate (horizontal edge
  at y_mid), those rows don't generate winding deltas. This means
  polygons whose extent aligns with cell centers may miss boundary
  rows. The coverage mode is unaffected. Refinement candidate.

## Next steps

1. **Wire approx mode through Python bindings**: add `mode` kwarg to
   the Python `burn()`.
2. **Refine approx mode boundary handling**: investigate the horizontal-
   edge-at-y_mid case and compare against fasterize output for parity.
3. **Port target/snap/clamp into `materialize.hpp`** using the
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
