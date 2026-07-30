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
