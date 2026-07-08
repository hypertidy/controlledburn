# test_burn.py -- invariant tests for the Python bindings, mirroring
# cpp/tests/test_burn.cpp, plus shapely WKB integration.
#
# The core invariant throughout: total burned coverage (full run cells +
# edge fractions, times cell area) equals the exact polygon area.
#
# Copyright (c) 2025 Michael Sumner
# Licensed under Apache License 2.0

import numpy as np
import pytest
import shapely

import controlledburn as cb

G10 = dict(bounds=(0, 0, 10, 10), shape=(10, 10))


def covered_area(r, bounds, shape):
    nrow, ncol = shape
    cell = ((bounds[2] - bounds[0]) / ncol) * ((bounds[3] - bounds[1]) / nrow)
    full = cell * (r.runs["col_end"] - r.runs["col_start"] + 1).sum()
    frac = cell * r.edges["fraction"].sum(dtype=np.float64)
    return full + frac


def wkb(geom):
    return shapely.to_wkb(geom)


def test_aligned_rectangle():
    g = shapely.box(2, 4, 6, 8)
    r = cb.burn([wkb(g)], **G10)
    assert len(r.edges) == 0
    assert len(r.runs) > 0
    assert covered_area(r, G10["bounds"], G10["shape"]) == pytest.approx(16.0)
    assert (r.runs["id"] == 1).all()


def test_offset_rectangle():
    g = shapely.box(2.5, 4.5, 6.5, 8.5)
    r = cb.burn([wkb(g)], **G10)
    assert covered_area(r, G10["bounds"], G10["shape"]) == pytest.approx(16.0)
    quarters = np.isclose(r.edges["fraction"], 0.25).sum()
    halves = np.isclose(r.edges["fraction"], 0.5).sum()
    assert quarters == 4
    assert halves == 12


def test_triangle_awkward_grid():
    tri = shapely.Polygon([(13.3, 17.7), (88.1, 22.4), (41.9, 79.2)])
    bounds, shape = (0, 0, 100, 100), (41, 37)
    r = cb.burn([wkb(tri)], bounds=bounds, shape=shape)
    assert covered_area(r, bounds, shape) == pytest.approx(tri.area, rel=1e-4)


def test_hole_both_orientations():
    outer = [(1, 1), (9, 1), (9, 9), (1, 9)]
    hole_ccw = [(3, 3), (7, 3), (7, 7), (3, 7)]
    for hole in (hole_ccw, hole_ccw[::-1]):
        g = shapely.Polygon(outer, [hole])
        r = cb.burn([wkb(g)], **G10)
        assert covered_area(r, G10["bounds"], G10["shape"]) == pytest.approx(48.0)


def test_multipolygon_disjoint():
    g = shapely.MultiPolygon([shapely.box(1, 1, 3, 3), shapely.box(6, 6, 9, 9)])
    r = cb.burn([wkb(g)], **G10)
    assert covered_area(r, G10["bounds"], G10["shape"]) == pytest.approx(13.0)
    assert (r.runs["id"] == 1).all()  # one geometry -> one id


def test_beyond_extent():
    # Exercises the padding-column winding case: edges entirely outside
    # the grid must still contribute winding to grid rows.
    g = shapely.box(-100, -100, 100, 100)
    r = cb.burn([wkb(g)], bounds=(0, 0, 10, 10), shape=(5, 5))
    assert len(r.edges) == 0
    assert covered_area(r, (0, 0, 10, 10), (5, 5)) == pytest.approx(100.0)


def test_line_length_conserved():
    g = shapely.LineString([(0.5, 0.5), (9.5, 7.5)])
    r = cb.burn([wkb(g)], **G10)
    assert len(r.runs) == 0 and len(r.edges) == 0
    assert r.lines["length"].sum(dtype=np.float64) == pytest.approx(
        g.length, rel=1e-4
    )


def test_points_binning():
    pts = [shapely.Point(0.5, 9.5), shapely.Point(9.5, 0.5), shapely.Point(15, 5)]
    r = cb.burn([wkb(p) for p in pts], **G10)
    assert len(r.points) == 2  # out-of-extent point dropped
    # 1-based, row 1 at top: matches the R package contract exactly
    assert tuple(r.points[0]) == (1, 1, 1)
    assert tuple(r.points[1]) == (10, 10, 2)


def test_ids_and_multiple_geometries():
    geoms = [shapely.box(2, 4, 6, 8), shapely.box(2.5, 4.5, 6.5, 8.5)]
    r = cb.burn([wkb(g) for g in geoms], **G10)
    assert set(np.unique(r.runs["id"])) <= {1, 2}
    assert set(np.unique(r.edges["id"])) == {2}  # only the offset box has edges


def test_geometrycollection_warns_and_skips():
    g = shapely.GeometryCollection([shapely.Point(0.5, 0.5)])
    with pytest.warns(UserWarning, match="GeometryCollection"):
        r = cb.burn([wkb(g)], **G10)
    assert len(r.runs) == 0 and len(r.points) == 0


def test_bad_wkb_warns_and_skips():
    with pytest.warns(UserWarning, match="WKB"):
        r = cb.burn([b"\x01\x03\x00", wkb(shapely.box(2, 4, 6, 8))], **G10)
    # second geometry still burned with id 2
    assert (r.runs["id"] == 2).all()
    assert covered_area(r, G10["bounds"], G10["shape"]) == pytest.approx(16.0)


def test_none_entries_skipped():
    r = cb.burn([None, wkb(shapely.box(2, 4, 6, 8))], **G10)
    assert (r.runs["id"] == 2).all()


def test_materialize_sum_and_overlap():
    geoms = [shapely.box(2, 4, 6, 8), shapely.box(4, 4, 8, 8)]
    r = cb.burn([wkb(g) for g in geoms], **G10)
    m = cb.materialize(r, values=[1.0, 1.0], fn="sum")
    assert m.shape == (10, 10)
    touched = np.isfinite(m)
    assert touched.sum() == 24  # union
    assert (m[touched] == 2.0).sum() == 8  # overlap
    assert np.nansum(m) == pytest.approx(32.0)


def test_materialize_fraction_conserves_area():
    g = shapely.box(2.5, 4.5, 6.5, 8.5)
    r = cb.burn([wkb(g)], **G10)
    m = cb.materialize(r, values=[1.0], fn="sum", edge_policy="fraction")
    assert np.nansum(m) == pytest.approx(16.0)  # area / cell_area


def test_materialize_ids_default_and_background():
    g = shapely.box(2, 4, 6, 8)
    r = cb.burn([wkb(g)], **G10)
    m = cb.materialize(r, background=0.0)
    assert m.sum() == pytest.approx(16.0)  # id 1 burned into 16 cells
    assert m.max() == 1.0


def test_pandas_roundtrip():
    pd = pytest.importorskip("pandas")
    g = shapely.box(2.5, 4.5, 6.5, 8.5)
    r = cb.burn([wkb(g)], **G10)
    df = pd.DataFrame(r.edges)
    assert list(df.columns) == ["row", "col", "fraction", "id"]
    assert len(df) == len(r.edges)


def test_invalid_grid_raises():
    with pytest.raises(Exception, match="extent|positive"):
        cb.burn([wkb(shapely.box(0, 0, 1, 1))], bounds=(0, 0, 0, 10), shape=(10, 10))
