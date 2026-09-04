# test_parity_fixtures.py -- cross-language parity tests driven by shared
# fixtures in fixtures/. Identical fixtures are read by C++ and R test
# suites so any regression is caught across all three surfaces.

import csv
import pathlib

import numpy as np
import pytest
import shapely

import controlledburn as cb

FIXTURES_DIR = pathlib.Path(__file__).resolve().parents[2] / "fixtures"


def load_fixtures():
    geoms = {}
    with open(FIXTURES_DIR / "geometries.csv") as f:
        for row in csv.DictReader(f):
            geoms[row["case"]] = row

    expected = {}
    with open(FIXTURES_DIR / "expected.csv") as f:
        for row in csv.DictReader(f):
            expected[row["case"]] = row

    fixtures = []
    for case, g in geoms.items():
        e = expected[case]
        fixtures.append((case, g, e))
    return fixtures


def covered_area(r, extent, shape):
    # extent is (xmin, xmax, ymin, ymax), matching R's extent ordering
    nrow, ncol = shape
    cell = ((extent[1] - extent[0]) / ncol) * ((extent[3] - extent[2]) / nrow)
    full = cell * (r.runs["col_end"] - r.runs["col_start"]).sum()
    frac = cell * r.edges["fraction"].sum(dtype=np.float64)
    return full + frac


@pytest.fixture(params=load_fixtures(), ids=lambda x: x[0])
def fixture(request):
    return request.param


def test_parity(fixture):
    case, g, e = fixture
    geom = shapely.from_wkt(g["wkt"])
    wkb = shapely.to_wkb(geom)

    extent = (float(g["xmin"]), float(g["xmax"]),
              float(g["ymin"]), float(g["ymax"]))
    shape = (int(g["nrow"]), int(g["ncol"]))

    r = cb.burn([wkb], extent=extent, shape=shape)

    if e["covered_area"] != "NA":
        expected_area = float(e["covered_area"])
        tol = float(e["tol_rel"])
        actual = covered_area(r, extent, shape)
        assert actual == pytest.approx(expected_area, rel=tol), \
            f"{case}: area {actual} != {expected_area}"

    if e["edges_empty"] != "NA" and int(e["edges_empty"]) == 1:
        assert len(r.edges) == 0, f"{case}: expected no edges"

    if e["line_length"] != "NA":
        expected_len = float(e["line_length"])
        tol = float(e["tol_rel"])
        actual = r.lines["length"].sum(dtype=np.float64)
        assert actual == pytest.approx(expected_len, rel=tol), \
            f"{case}: line length {actual} != {expected_len}"

    if e["n_points"] != "NA":
        assert len(r.points) == int(e["n_points"]), \
            f"{case}: point count {len(r.points)} != {e['n_points']}"


# The checks above (areas, counts, line length) pass under any consistent
# index base. The tests below pin the *0-based, half-open, id = k* contract
# on known fixtures, so a regression -- e.g. reintroducing the R shim's +1
# in the core -- is caught as an exact table mismatch, not silently.

_GEOMS = {g["case"]: g for g in
          csv.DictReader(open(FIXTURES_DIR / "geometries.csv"))}


def _burn_fixture(case, mode="coverage"):
    g = _GEOMS[case]
    wkb = shapely.to_wkb(shapely.from_wkt(g["wkt"]))
    extent = (float(g["xmin"]), float(g["xmax"]),
              float(g["ymin"]), float(g["ymax"]))
    shape = (int(g["nrow"]), int(g["ncol"]))
    return cb.burn([wkb], extent=extent, shape=shape, mode=mode)


def test_runs_table_zero_based_coverage():
    # beyond_extent: whole 5x5 grid covered -> one clean run per row,
    # columns [0, ncol) exclusive. Pins row 0 = top, col_end exclusive,
    # id = 0 for the single geometry.
    r = _burn_fixture("beyond_extent")
    expected = np.array(
        [(0, 0, 5, 0), (1, 0, 5, 0), (2, 0, 5, 0), (3, 0, 5, 0), (4, 0, 5, 0)],
        dtype=r.runs.dtype,
    )
    assert np.array_equal(r.runs, expected)
    assert len(r.edges) == 0


def test_runs_table_zero_based_approx():
    # aligned_rect box x[2,6] y[4,8] on a 10x10 grid, approx mode -> clean
    # per-row runs: rows 2..5 (0-based, row 0 at top), cols [2, 6) = 4 wide.
    r = _burn_fixture("aligned_rect", mode="approx")
    expected = np.array(
        [(2, 2, 6, 0), (3, 2, 6, 0), (4, 2, 6, 0), (5, 2, 6, 0)],
        dtype=r.runs.dtype,
    )
    assert np.array_equal(r.runs, expected)


def test_points_table_zero_based():
    # points_mixed is a single MULTIPOINT -> all points carry id 0; the
    # off-grid point is dropped. (0, 0) is top-left, (9, 9) bottom-right.
    r = _burn_fixture("points_mixed")
    expected = np.array([(0, 0, 0), (9, 9, 0)], dtype=r.points.dtype)
    assert np.array_equal(r.points, expected)
