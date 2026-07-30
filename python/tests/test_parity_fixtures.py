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


def covered_area(r, bounds, shape):
    nrow, ncol = shape
    cell = ((bounds[2] - bounds[0]) / ncol) * ((bounds[3] - bounds[1]) / nrow)
    full = cell * (r.runs["col_end"] - r.runs["col_start"] + 1).sum()
    frac = cell * r.edges["fraction"].sum(dtype=np.float64)
    return full + frac


@pytest.fixture(params=load_fixtures(), ids=lambda x: x[0])
def fixture(request):
    return request.param


def test_parity(fixture):
    case, g, e = fixture
    geom = shapely.from_wkt(g["wkt"])
    wkb = shapely.to_wkb(geom)

    bounds = (float(g["xmin"]), float(g["ymin"]),
              float(g["xmax"]), float(g["ymax"]))
    shape = (int(g["nrow"]), int(g["ncol"]))

    r = cb.burn([wkb], bounds=bounds, shape=shape)

    if e["covered_area"] != "NA":
        expected_area = float(e["covered_area"])
        tol = float(e["tol_rel"])
        actual = covered_area(r, bounds, shape)
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
