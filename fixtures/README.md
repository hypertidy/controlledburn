# Parity fixtures

Shared test geometries and expected results for cross-language parity testing.
Read by C++ ctest, R testthat, and Python pytest.

## Files

- **geometries.csv** -- WKT geometries with grid specs and labels.
- **expected.csv** -- expected aggregate results per test case (total area,
  run count, edge count, etc.).

## Contract

Each row in `geometries.csv` defines a test case with:

| Column     | Description                                              |
|------------|----------------------------------------------------------|
| `case`     | unique label                                             |
| `wkt`      | well-known text (one geometry per row)                   |
| `xmin`     | grid extent xmin                                         |
| `ymin`     | grid extent ymin                                         |
| `xmax`     | grid extent xmax                                         |
| `ymax`     | grid extent ymax                                         |
| `ncol`     | grid columns                                             |
| `nrow`     | grid rows                                                |

Each row in `expected.csv` defines the expected results:

| Column           | Description                                        |
|------------------|----------------------------------------------------|
| `case`           | matches `geometries.csv`                           |
| `covered_area`   | total burned coverage in CRS units²                |
| `n_runs`         | number of full-cell run records (0 = "don't check")|
| `n_edges`        | number of edge records (0 = "don't check")         |
| `edges_empty`    | 1 if edges must be empty, 0 otherwise              |
| `line_length`    | total line length (NA for non-line cases)          |
| `n_points`       | number of point records (NA for non-point cases)   |
| `tol_rel`        | relative tolerance for area check                  |

## Canonical tables

- **canon-tables.txt** -- record-level tables written by the Rust crate,
  hypertidy/controlledburn-rs (`tests/canon.rs`), which is the canonical
  implementation of the burn algorithm. 240 generated cases on
  non-dyadic grids, with vertices steered onto cell corners, walls and
  centre lines (where ties are decided), holes, multipolygons, lines and
  points, and every record the Rust crate emits for them in both modes
  with run coalescing off. `cpp/tests/test_canon_tables.cpp` requires the
  C++ core to emit the same records in the same order, with fractions
  and lengths equal as bits.

  The file is not edited here. When the Rust crate changes the algorithm
  on purpose, it regenerates the file (`CANON_WRITE=1 cargo test --release
  --test canon`), the new file is copied here, and the change is ported to
  `cpp/` until the test passes. R and Python see the same core through
  `tools/sync-core.sh`, so they need no copy of this check.
