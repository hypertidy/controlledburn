# test_regressions.py -- Python-side coverage of core fixes made on the R
# side (lattice snapping, #16; row-centre crossings, #19) plus input
# validation at the binding boundary. The engine is shared, so these are
# the same cases as tests/testthat/test-lattice-snap.R and
# cpp/tests/test_burn.cpp, checked through the Python surface.

import struct
import sys
import warnings

import numpy as np
import pytest

import controlledburn as cb

shapely = pytest.importorskip("shapely")


def covered_area(r):
    xmin, xmax, ymin, ymax = r.extent
    nrow, ncol = r.shape
    cell = (xmax - xmin) / ncol * (ymax - ymin) / nrow
    full = int((r.runs["col_end"] - r.runs["col_start"]).sum())
    return cell * (full + float(r.edges["fraction"].sum(dtype="f8")))


def tri_area(a, b, c):
    return abs((b[0] - a[0]) * (c[1] - a[1]) - (c[0] - a[0]) * (b[1] - a[1])) / 2


# ---- lattice snap (#16) --------------------------------------------------

@pytest.mark.parametrize("wkt,pts", [
    ("POLYGON ((0.3 0.3, 4.7 0.3, 2.5 6.2, 0.3 0.3))",
     ((0.3, 0.3), (4.7, 0.3), (2.5, 6.2))),
    ("POLYGON ((4.7 0.3, 2.5 6.2, 0.3 0.3, 4.7 0.3))",
     ((0.3, 0.3), (4.7, 0.3), (2.5, 6.2))),
    ("POLYGON ((0.3 0.3, 4.7 0.3, 2.5 6.25, 0.3 0.3))",
     ((0.3, 0.3), (4.7, 0.3), (2.5, 6.25))),
])
def test_vertex_ulp_off_cell_edge(wkt, pts):
    # 0.3 parsed vs 3 * 0.1 computed differ by ~5e-17; before snapping the
    # coverage walk left the grid and ran until bad_alloc.
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        r = cb.burn(shapely.from_wkt(wkt), extent=(0, 10, 0, 10), shape=(100, 100))
    assert len(r.runs) > 0
    assert covered_area(r) == pytest.approx(tri_area(*pts), rel=1e-6)


# ---- row-centre crossing on a column wall (#19, coverage) ----------------

FALKLANDS = ("POLYGON ((0.8 1.15, 2.0 1.75, 2.85 1.5, 3.45 1.9, 4.25 1.45, "
             "3.95 1.1, 2.6 0.8, 2.15 1.15, 1.3 0.7, 0.8 1.15))")


@pytest.mark.parametrize("n", [20, 100, 200, 1000])
def test_wall_crossing_keeps_winding(n):
    # At 20 columns the first edge crosses y = 1.25 at x = 1.0, exactly on a
    # column wall; with a strict crossing test the row lost interior cells.
    g = shapely.from_wkt(FALKLANDS)
    r = cb.burn(g, extent=(0, 10, 0, 10), shape=(n, n))
    assert covered_area(r) == pytest.approx(g.area, rel=1e-6)
    assert g.area == pytest.approx(2.12875)


def test_wall_crossing_materialized_fraction_conserves_area():
    g = shapely.from_wkt(FALKLANDS)
    r = cb.burn(g, extent=(0, 10, 0, 10), shape=(20, 20))
    m = r.materialize(fn="sum", values=[1.0], edge_policy="fraction",
                      background=0.0)
    assert m.sum() * 0.25 == pytest.approx(2.12875, rel=1e-6)


# ---- approx: vertex on a row centre line (#19, approx) -------------------

def test_approx_centre_line_vertex_does_not_fill_row():
    # The middle vertex sits on row 24's centre line as a user computes it,
    # one ulp below the line recomputed from the extent. The sign test used
    # to call the downward edge upward and the row filled to the grid edge.
    g = shapely.from_wkt(
        "POLYGON ((-179.89654756364214 -88.83333333333334, "
        "-180.94076559960791 -92.3680897750401, "
        "-179.6 -88.518152273978, "
        "-179.89654756364214 -88.83333333333334))")
    r = cb.burn(g, extent=(-180, -175, -90, -80.66666666666667),
                shape=(28, 15), mode="approx")
    assert len(r.runs) == 2
    assert (r.runs["col_end"] <= 2).all()


# ---- grid validation ----------------------------------------------------

SQUARE = "POLYGON ((1 1, 3 1, 3 3, 1 3, 1 1))"


@pytest.mark.parametrize("extent", [
    (0, np.inf, 0, 4), (-np.inf, 4, 0, 4), (0, 4, 0, np.inf),
    (0, np.nan, 0, 4), (np.nan, 4, 0, 4),
])
def test_non_finite_extent_raises(extent):
    # inf used to reach the core with an infinite cell size and segfault
    with pytest.raises(ValueError, match="extent"):
        cb.burn(shapely.from_wkt(SQUARE), extent=extent, shape=(4, 4))


def test_overflowing_extent_span_raises_from_core():
    # finite bounds whose span overflows: caught by the core's validation
    with pytest.raises(ValueError, match="overflow"):
        cb.burn(shapely.from_wkt(SQUARE), extent=(-1e308, 1e308, 0, 4),
                shape=(4, 4))


def test_core_rejects_non_finite_extent_directly():
    # the private binding validates too, so no caller can crash the core
    wkb = [shapely.to_wkb(shapely.from_wkt(SQUARE))]
    with pytest.raises(ValueError, match="finite"):
        cb._core.burn_wkb(wkb, 0, 0, np.inf, 4, 4, 4)


def test_extent_wrong_length_raises():
    with pytest.raises(ValueError, match="xmin, xmax, ymin, ymax"):
        cb.burn(shapely.from_wkt(SQUARE), extent=(0, 4, 0), shape=(4, 4))


def test_shape_wrong_length_raises():
    with pytest.raises(ValueError, match="nrow, ncol"):
        cb.burn(shapely.from_wkt(SQUARE), extent=(0, 4, 0, 4), shape=(4,))


def test_shape_beyond_int32_raises_clearly():
    with pytest.raises(ValueError, match="int32"):
        cb.burn(shapely.from_wkt(SQUARE), extent=(0, 4, 0, 4),
                shape=(4, 2**31))


# ---- extent derivation without shapely ----------------------------------

def _wkb_polygon(coords):
    # little-endian ISO WKB polygon, one ring
    out = struct.pack("<BII", 1, 3, 1) + struct.pack("<I", len(coords))
    for x, y in coords:
        out += struct.pack("<dd", x, y)
    return out


def test_extent_from_wkb_needs_no_shapely(monkeypatch):
    wkb = _wkb_polygon([(1, 2), (5, 2), (5, 7), (1, 7), (1, 2)])
    monkeypatch.setitem(sys.modules, "shapely", None)   # import would fail
    r = cb.burn([wkb, None])
    assert r.extent == (1.0, 5.0, 2.0, 7.0)
    assert covered_area(r) == pytest.approx(20.0, rel=1e-6)


def test_extent_from_all_empty_geometry_raises_clearly():
    with pytest.raises(ValueError, match="no finite coordinates"):
        cb.burn([shapely.from_wkt("POLYGON EMPTY"), None])


def test_extent_matches_shapely_bounds():
    gs = shapely.from_wkt(["POLYGON ((-3 1, 2 1, 2 9, -3 1))",
                           "LINESTRING (0 -4, 6 0)", "POINT (7 2)"])
    (xmin, xmax, ymin, ymax), _ = cb.resolve_grid(gs)
    b = shapely.total_bounds(gs)
    assert (xmin, ymin, xmax, ymax) == tuple(float(v) for v in b)


# ---- buffer handling at the binding boundary -----------------------------

def test_wkb_from_memoryview_with_wide_items():
    # a contiguous buffer of any item size carries its full byte length
    # (previously only the item count was passed, truncating the WKB)
    wkb = shapely.to_wkb(shapely.from_wkt(SQUARE))
    pad = wkb + b"\x00" * (-len(wkb) % 8)            # trailing bytes are ignored
    r8 = cb.burn([memoryview(wkb)], extent=(0, 4, 0, 4), shape=(4, 4))
    with warnings.catch_warnings():
        warnings.simplefilter("error")
        r64 = cb.burn([memoryview(pad).cast("Q")], extent=(0, 4, 0, 4),
                      shape=(4, 4))
    assert covered_area(r8) == covered_area(r64) == pytest.approx(4.0)


def test_non_contiguous_buffer_rejected():
    # a strided view used to be read as if contiguous
    wkb = shapely.to_wkb(shapely.from_wkt(SQUARE))
    strided = memoryview(bytes(b for c in wkb for b in (c, 0)))[::2]
    assert bytes(strided) == wkb
    with pytest.raises(ValueError, match="contiguous"):
        cb.burn([strided], extent=(0, 4, 0, 4), shape=(4, 4))


def test_materialize_core_rejects_ragged_columns():
    i = np.zeros(2, dtype=np.int32)
    with pytest.raises(ValueError, match="differ in length"):
        cb._core.materialize(i, i[:1], i, i, i[:0], i[:0],
                             np.zeros(0, np.float32), i[:0],
                             4, 4, [], "last", "threshold", 0.5, np.nan)
