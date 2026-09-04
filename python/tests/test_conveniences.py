# test_conveniences.py -- tests for the ergonomic layer added on top of
# the core bindings: input coercion (as_wkb_list), grid parameter
# resolution (resolve_grid, optional extent/shape/resolution on burn),
# the summary __repr__, and materialize() handling of line/point/mixed
# input.
#
# Copyright (c) 2025 Michael Sumner
# Licensed under Apache License 2.0

import numpy as np
import pytest
import shapely

import controlledburn as cb

G10 = dict(extent=(0, 10, 0, 10), shape=(10, 10))
BOX = shapely.box(2.5, 4.5, 6.5, 8.5)
LINE = shapely.LineString([(0.5, 0.5), (9.5, 7.5)])


def covered_area(r):
    xmin, xmax, ymin, ymax = r.extent
    nrow, ncol = r.shape
    cell = ((xmax - xmin) / ncol) * ((ymax - ymin) / nrow)
    full = (r.runs["col_end"] - r.runs["col_start"]).sum()
    return cell * (full + r.edges["fraction"].sum(dtype=np.float64))


# ---- input coercion --------------------------------------------------

def test_burn_accepts_single_shapely_geometry():
    r = cb.burn(BOX, **G10)
    assert covered_area(r) == pytest.approx(16.0)


def test_burn_accepts_list_of_shapely_geometries():
    r = cb.burn([BOX, shapely.box(0, 0, 1, 1)], **G10)
    assert set(np.unique(np.concatenate([r.runs["id"], r.edges["id"]]))) == {0, 1}


def test_burn_accepts_bare_wkb_bytes():
    r = cb.burn(shapely.to_wkb(BOX), **G10)
    assert covered_area(r) == pytest.approx(16.0)


def test_burn_accepts_shapely_array():
    arr = shapely.to_wkb(np.array([BOX, BOX]))
    r = cb.burn(arr, **G10)
    assert covered_area(r) == pytest.approx(32.0)


def test_burn_accepts_mixed_wkb_shapely_none():
    r = cb.burn([None, shapely.to_wkb(BOX), shapely.box(0, 0, 1, 1)], **G10)
    ids = set(np.unique(np.concatenate([r.runs["id"], r.edges["id"]])))
    assert ids == {1, 2}  # None consumed id 0


def test_burn_accepts_geopandas_geoseries():
    gpd = pytest.importorskip("geopandas")
    gs = gpd.GeoSeries([BOX])
    r = cb.burn(gs, **G10)
    assert covered_area(r) == pytest.approx(16.0)


def test_as_wkb_list_rejects_nonsense():
    with pytest.raises(TypeError, match="cannot interpret"):
        cb.as_wkb_list([42])


# ---- grid parameter resolution ---------------------------------------

def test_extent_default_from_geometry():
    r = cb.burn(BOX, shape=(10, 10))
    assert r.extent == (2.5, 6.5, 4.5, 8.5)


def test_shape_default_256_fit():
    # taller than wide: 256 on the long (y) axis, aspect preserved
    r = cb.burn(shapely.box(0, 0, 5, 10))
    assert r.shape == (256, 128)
    assert r.extent == (0.0, 5.0, 0.0, 10.0)


def test_resolution_scalar():
    r = cb.burn(BOX, extent=(0, 10, 0, 10), resolution=0.5)
    assert r.shape == (20, 20)


def test_resolution_anisotropic():
    r = cb.burn(BOX, extent=(0, 10, 0, 10), resolution=(0.5, 1.0))
    assert r.shape == (10, 20)  # (nrow from dy, ncol from dx)


def test_resolution_and_shape_conflict():
    with pytest.raises(ValueError, match="not both"):
        cb.burn(BOX, extent=(0, 10, 0, 10), shape=(10, 10), resolution=0.5)


def test_resolution_matches_r_ceiling_rule():
    # R: dimension <- ceiling(c(dx_extent, dy_extent) / resolution)
    extent, shape = cb.resolve_grid(BOX, extent=(0, 10, 0, 10), resolution=3.0)
    assert shape == (4, 4)  # ceil(10/3)


def test_invalid_extent_still_raise():
    with pytest.raises(ValueError, match="extent"):
        cb.burn(shapely.to_wkb(BOX), extent=(0, 0, 0, 10), shape=(10, 10))


def test_invalid_shape_still_raises():
    with pytest.raises(ValueError, match="positive"):
        cb.burn(shapely.to_wkb(BOX), extent=(0, 10, 0, 10), shape=(0, 10))


# ---- repr -------------------------------------------------------------

def test_repr_is_summary_not_dump():
    r = cb.burn(BOX, **G10)
    s = repr(r)
    assert "<BurnResult>" in s
    assert "10 x 10 grid" in s
    assert "1 geometry" in s
    assert "16 polygon boundary cells" in s
    assert "sparsity" in s
    assert "dtype" not in s          # no structured-array dump
    assert len(s) < 500


def test_repr_mixed_kinds_lists_all_tables():
    r = cb.burn([BOX, LINE, shapely.Point(1.5, 1.5)], **G10)
    s = repr(r)
    assert "3 geometries" in s
    assert "lines:" in s and "points:" in s


# ---- materialize: polygon-only contract ------------------------------
# Dense line/point semantics are deliberately unresolved (see the
# tracking issue): materialize() must fail LOUDLY on non-polygon input
# rather than silently returning an all-background array, which was the
# pre-0.2.0 behaviour.

def test_materialize_lines_raises_not_silent():
    r = cb.burn(LINE, **G10)
    with pytest.raises(NotImplementedError, match="line"):
        cb.materialize(r)


def test_materialize_points_raises_not_silent():
    r = cb.burn(shapely.Point(1.5, 1.5), **G10)
    with pytest.raises(NotImplementedError, match="point"):
        cb.materialize(r)


def test_materialize_mixed_kinds_raises():
    r = cb.burn([BOX, LINE], **G10)
    with pytest.raises(ValueError, match="mixed geometry kinds"):
        cb.materialize(r)


def test_line_tables_remain_the_product():
    # the sparse table is the supported dense path: a few lines of numpy
    r = cb.burn(LINE, **G10)
    out = np.zeros(r.shape)
    np.add.at(out, (r.lines["row"] - 1, r.lines["col"] - 1),
              r.lines["length"].astype("f8"))
    assert out.sum() == pytest.approx(LINE.length, rel=1e-4)


def test_materialize_polygon_path_unchanged():
    # polygon-only input still goes through the C++ core
    r = cb.burn(BOX, **G10)
    m = cb.materialize(r, values=[1.0], fn="sum", edge_policy="fraction",
                       background=0.0)
    assert m.sum() == pytest.approx(16.0)


def test_materialize_empty_result():
    # an empty burn materializes to all-background via the polygon path
    r = cb.burn([shapely.to_wkb(shapely.Point(50, 50))], **G10)  # off-grid
    m = cb.materialize(r, background=0.0)
    assert m.sum() == 0.0


def test_materialize_bad_fn():
    r = cb.burn(BOX, **G10)
    with pytest.raises(ValueError, match="fn must be"):
        cb.materialize(r, fn="mean")
