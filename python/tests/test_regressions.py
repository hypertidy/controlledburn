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


# ---- holes: a hole's boundary cells keep the exterior's cover -------------

HOLE_EXT = [(0.5, 0.5), (9.5, 0.5), (9.5, 9.5), (0.5, 9.5)]


@pytest.mark.parametrize("hole", [
    [(4.2, 4.2), (5.8, 4.2), (5.8, 5.8), (4.2, 5.8)],   # across four cells
    [(4.2, 4.2), (4.8, 4.2), (4.8, 4.8), (4.2, 4.8)],   # within one cell
    [(4.2, 4.6), (5.8, 4.6), (5.0, 5.4)],               # off the row centres
    [(4.2, 4.2), (4.2, 5.8), (5.8, 5.8), (5.8, 4.2)],   # the other winding
])
def test_hole_removes_exactly_its_area(hole):
    # A cell crossed only by the hole's boundary used to get the hole's
    # negative fraction alone and be dropped (77 cells instead of 78.44
    # for the first hole).
    g = shapely.Polygon(HOLE_EXT, [hole])
    r = cb.burn(g, extent=(0, 10, 0, 10), shape=(10, 10))
    assert covered_area(r) == pytest.approx(g.area, rel=1e-6)


def test_hole_cells_match_exact_intersection():
    g = shapely.Polygon([(0, 0), (10, 0), (10, 10), (0, 10)],
                        [[(1.3, 1.3), (3.7, 1.3), (3.7, 3.7), (1.3, 3.7)],
                         [(6.1, 6.6), (8.9, 6.6), (7.5, 8.8)]])
    r = cb.burn(g, extent=(0, 10, 0, 10), shape=(10, 10))
    m = r.materialize(fn="sum", values=[1.0], edge_policy="fraction", background=0.0)
    for row in range(10):
        for col in range(10):
            cell = shapely.box(col, 9 - row, col + 1, 10 - row)
            assert m[row, col] == pytest.approx(g.intersection(cell).area, abs=1e-6), (row, col)


def test_many_holes_cover_their_area():
    holes = []
    for i in range(20):
        for j in range(20):
            cx, cy, rr = 2.5 + 5 * i + 0.013, 2.5 + 5 * j + 0.029, 1.7
            holes.append([(cx - rr, cy), (cx, cy - rr), (cx + rr, cy), (cx, cy + rr)])
    g = shapely.Polygon([(0, 0), (100, 0), (100, 100), (0, 100)], holes)
    r = cb.burn(g, extent=(0, 100, 0, 100), shape=(333, 333))
    assert covered_area(r) == pytest.approx(g.area, rel=1e-6)


def test_wall_notch_cell_is_a_boundary_cell():
    # Lattice-aligned notch from a raster-derived waterbody (DEA Waterbodies
    # uid q9g5b3eyr_v1). Cell [4, 5] x [2, 3] has edges along both walls,
    # crossing the row centre in opposite directions: zero coverage, zero
    # net winding. It got no boundary record and was filled as interior.
    g = shapely.Polygon([(4, 2), (2, 7), (6, 2), (5, 2), (5, 3), (3, 5.5), (4, 2.5)])
    r = cb.burn(g, extent=(0, 10, 0, 10), shape=(10, 10))
    assert covered_area(r) == pytest.approx(g.area, rel=1e-6)
    m = r.materialize(fn="sum", values=[1.0], edge_policy="fraction", background=0.0)
    for row in range(10):
        for col in range(10):
            cell = shapely.box(col, 9 - row, col + 1, 10 - row)
            assert m[row, col] == pytest.approx(g.intersection(cell).area, abs=1e-6), (row, col)


def test_zero_width_spike_cell_stays_covered():
    # Found by the controlledburn-rs property suite: the ring runs in
    # through the top-left corner of cell [4.8, 5.1] x [1.2, 1.5] to its
    # centre and straight back out, enclosing nothing. The cell is inside.
    g = shapely.Polygon([
        (2.2955346720684804, -0.28831957426678434), (4.2871984006926445, -0.6968703145658899),
        (6.183801073287196, -1.4711665488619694), (5.24135385306532, 1.4474746273177665),
        (5.275730794077829, 1.5426958381891471), (5.824891418164306, 2.019485421020327),
        (5.151353853065319, 1.696766993828962), (3.5999999999999996, 2.6999999999999997),
        (4.95, 1.3499999999999996), (4.8, 1.4999999999999996),
    ])
    r = cb.burn(g, extent=(0, 6, 0, 5.3999999999999995), shape=(18, 20))
    m = r.materialize(fn="sum", values=[1.0], edge_policy="fraction", background=0.0)
    assert m[13, 16] == pytest.approx(1.0)
    assert covered_area(r) == pytest.approx(3.026275470389206, rel=1e-6)


CORNER_CLIP_CASES = [
    # controlledburn#25: cell (12, 18) came out 1 instead of 0.0943
    ("POLYGON ((13.49307281327689 0.280951025273552, 13.774064548591133 0.1953615710514231, "
     "14.150447247739818 0.19536157105142304, 14.454947247739819 0.2945346261241548, "
     "14.571255898165475 0.455, 14.454947247739819 0.6154653738758452, "
     "17.207761779641103 4.932666216410308, 12.6 1.4000000000000004, "
     "13.649999999999999 0.34999999999999964, 6.0614012054372814 0.45500000000000046, "
     "13.49307281327689 0.280951025273552))",
     (0, 20.299999999999997, 0, 9.1), (13, 29)),
    # controlledburn-rs#6: cell (4, 4) came out 0.739 instead of 0.455
    ("POLYGON ((2.9404541633227357 6.235726447802992, 7.140490483291759 -0.6880522184656641, "
     "12.79922599805937 6.649999999999999, 13.222799970968294 7.999704574445243, "
     "12.79922599805937 8.351528741555857, 8.3 10.399999999999999, 12.05 6.649999999999999, "
     "2.9404541633227357 6.235726447802992))",
     (-1.7, 35.8, 2.9, 17.9), (6, 15)),
    # controlledburn-rs#6: cell (3, 0) came out 0.839 instead of 0.339
    ("POLYGON ((0.22857142857142856 0.657142857142857, 0.3 0.7285714285714284, "
     "0.5385327833222915 0.5669218770762636, 0.3826046772369731 0.8968780211084503, "
     "0.44285714285714284 0.9150607669103901, 0.6141351860119378 1.2110801232422874, "
     "0.3 0.9407750526246758, 0.37142857142857144 0.7999999999999999, "
     "0.22857142857142856 0.657142857142857))",
     (0.3, 1.0142857142857142, 0.3, 1.157142857142857), (6, 5)),
    # a zero-width spike out of a cell and back: exit and next entry an ulp
    # apart sent the area chase round the whole cell (rs property suite)
    ("POLYGON ((10.449999999999999 1.3499999999999996, 10.1 0.29999999999999982, "
     "10.800000000000001 2.3999999999999995, 13.447616999963875 1.037844845134956, "
     "15.411473092828203 2.2720230998577735, 11.15 2.4078014931373346, "
     "10.800000000000001 2.4917171551219539, 10.229283250618847 2.4917171551219539, "
     "3.0359419646823742 4.3208089369106828, 9.9390703882563258 2.2720230998577735, "
     "10.449999999999999 1.3499999999999996))",
     (0.3, 14.3, 0.3, 7.999999999999999), (11, 20)),
    ("POLYGON ((1000.1 -2000.5130338951737, 1000.1819562805505 -2000.5770261925379, "
     "1000.2234724980999 -2000.6086258946541, 1000.2686486160532 -2000.5978003193729, "
     "1000.3031415361594 -2000.5500000000002, 1000.7 -2000.1000000000001, 1000.85 -1999.95, "
     "1000.25 -2000.5500000000002, 1000.2234724980999 -2000.3413741053462, "
     "999.56959736624071 -1999.7508032852861, 999.25351910076768 -1999.9946853427875, "
     "1000.1 -2000.5130338951737))",
     (1000.1, 1002.8000000000001, -2000.7, -1996.2), (15, 9)),
]


@pytest.mark.parametrize("wkt,extent,shape", CORNER_CLIP_CASES)
def test_corner_clip_cells_match_exact_intersection(wkt, extent, shape):
    # An edge clipping a cell corner is recorded as a traversal with entry
    # and exit an ulp or two apart; compared exactly it counted as movement
    # and the cell gained a whole cell of coverage. Every start vertex.
    xmin, xmax, ymin, ymax = extent
    nrow, ncol = shape
    dx, dy = (xmax - xmin) / ncol, (ymax - ymin) / nrow
    ring = list(shapely.from_wkt(wkt).exterior.coords)[:-1]
    for k in range(len(ring)):
        g = shapely.Polygon(ring[k:] + ring[:k])
        r = cb.burn(g, extent=extent, shape=shape)
        m = r.materialize(fn="sum", values=[1.0], edge_policy="fraction", background=0.0)
        for row in range(nrow):
            for col in range(ncol):
                cell = shapely.box(xmin + col * dx, ymax - (row + 1) * dy, xmin + (col + 1) * dx, ymax - row * dy)
                want = g.intersection(cell).area / (dx * dy)
                assert m[row, col] == pytest.approx(want, abs=1e-5), (k, row, col)
