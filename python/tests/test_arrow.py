# test_arrow.py -- Arrow C data interface input path: burn() consuming
# WKB directly from Arrow (large_)binary / geoarrow.wkb buffers via
# nanoarrow, with no shapely on the burn path. Parity is checked against
# the WKB-list path (which uses the same C++ core).
#
# Copyright (c) 2025 Michael Sumner
# Licensed under Apache License 2.0

import numpy as np
import pytest
import shapely

import controlledburn as cb

pa = pytest.importorskip("pyarrow")
pytest.importorskip("nanoarrow")

G10 = dict(extent=(0, 10, 0, 10), shape=(10, 10))
BOX = shapely.box(2.5, 4.5, 6.5, 8.5)
TWO = [shapely.box(2, 2, 6, 6), shapely.box(4, 4, 8, 8)]


def wkbs(geoms):
    return [shapely.to_wkb(g) for g in geoms]


def _same(a, b):
    for k in ("runs", "edges", "lines", "points"):
        assert np.array_equal(getattr(a, k), getattr(b, k)), k


# ---- parity with the WKB-list path ----------------------------------

def test_arrow_binary_matches_wkb_list():
    ws = wkbs(TWO)
    _same(cb.burn(ws, **G10), cb.burn(pa.array(ws, type=pa.binary()), **G10))


def test_arrow_large_binary_matches():
    ws = wkbs(TWO)
    _same(cb.burn(ws, **G10), cb.burn(pa.array(ws, type=pa.large_binary()), **G10))


def test_arrow_approx_mode_matches():
    ws = wkbs([BOX])
    ref = cb.burn(ws, mode="approx", **G10)
    got = cb.burn(pa.array(ws, type=pa.binary()), mode="approx", **G10)
    assert np.array_equal(ref.runs, got.runs)
    assert len(got.edges) == 0


def test_arrow_lines_match():
    line = shapely.to_wkb(shapely.LineString([(0.5, 0.5), (9.5, 7.5)]))
    ref = cb.burn([line], **G10)
    got = cb.burn(pa.array([line], type=pa.binary()), **G10)
    assert np.array_equal(ref.lines, got.lines)


# ---- nulls and chunking ---------------------------------------------

def test_arrow_nulls_consume_an_id():
    ws = wkbs(TWO)
    r = cb.burn(pa.array([ws[0], None, ws[1]], type=pa.binary()), **G10)
    ids = {int(i) for i in np.unique(np.concatenate([r.runs["id"], r.edges["id"]]))}
    assert ids == {0, 2}          # the None consumed id 1


def test_arrow_chunked_stream_offsets_ids():
    ws = wkbs(TWO)
    ref = cb.burn(ws, **G10)
    ca = pa.chunked_array([[ws[0]], [ws[1]]], type=pa.binary())
    r = cb.burn(ca, **G10)
    ids = {int(i) for i in np.unique(np.concatenate([r.runs["id"], r.edges["id"]]))}
    assert ids == {0, 1}
    assert np.array_equal(r.runs, ref.runs)


# ---- coercion + validation ------------------------------------------

def test_as_wkb_list_accepts_arrow():
    ws = wkbs([shapely.box(2, 2, 6, 6)])
    lst = cb.as_wkb_list(pa.array([ws[0], None], type=pa.binary()))
    assert lst[1] is None
    assert bytes(lst[0]) == ws[0]


def test_arrow_extent_derived_from_geometry():
    # box(2,2,6,6) + box(4,4,8,8) -> envelope x[2,8], y[2,8]
    r = cb.burn(pa.array(wkbs(TWO), type=pa.binary()))   # no extent
    assert r.extent == (2.0, 8.0, 2.0, 8.0)              # (xmin, xmax, ymin, ymax)


def test_arrow_default_grid_matches_shapely_path():
    ws = wkbs(TWO)
    r_arrow = cb.burn(pa.array(ws, type=pa.binary()))    # extent+shape derived, no shapely
    r_list = cb.burn(ws)                                 # extent via shapely, same rule
    assert r_arrow.extent == r_list.extent
    assert r_arrow.shape == r_list.shape
    assert np.array_equal(r_arrow.runs, r_list.runs)


def test_arrow_all_null_extent_raises():
    with pytest.raises(ValueError, match="no finite coordinates"):
        cb.burn(pa.array([None, None], type=pa.binary()))


def test_arrow_native_encoding_rejected():
    # a non-binary Arrow array is not WKB -> a clear error, not garbage
    arr = pa.array([{"x": 1.0, "y": 2.0}],
                   type=pa.struct([("x", pa.float64()), ("y", pa.float64())]))
    with pytest.raises((NotImplementedError, ValueError)):
        cb.burn(arr, **G10)


def test_arrow_values_are_zero_copy_view():
    ws = wkbs(TWO)
    data, offsets, valid, n = next(cb._arrow_chunks(pa.array(ws, type=pa.binary())))
    assert data.dtype == np.uint8
    assert offsets.dtype == np.int64 and offsets.shape[0] == n + 1
    assert valid is None


# ---- real GeoArrow sources (skipped when the libs are absent) --------

def test_geopandas_to_arrow_wkb_matches():
    gpd = pytest.importorskip("geopandas")
    arr = gpd.GeoSeries(TWO).to_arrow(geometry_encoding="WKB")
    _same(cb.burn(wkbs(TWO), **G10), cb.burn(arr, **G10))


def test_geopandas_to_arrow_derives_extent():
    gpd = pytest.importorskip("geopandas")
    r = cb.burn(gpd.GeoSeries(TWO).to_arrow(geometry_encoding="WKB"))
    assert r.extent == (2.0, 8.0, 2.0, 8.0)


def test_geopandas_native_geoarrow_rejected():
    gpd = pytest.importorskip("geopandas")
    arr = gpd.GeoSeries(TWO).to_arrow(geometry_encoding="geoarrow")
    with pytest.raises((NotImplementedError, ValueError)):
        cb.burn(arr, **G10)


def test_geoarrow_pyarrow_wkb_extension_matches():
    ga = pytest.importorskip("geoarrow.pyarrow")
    gawkb = ga.wkb().wrap_array(pa.array(wkbs(TWO), type=pa.binary()))
    _same(cb.burn(wkbs(TWO), **G10), cb.burn(gawkb, **G10))
