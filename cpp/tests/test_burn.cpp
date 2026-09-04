// test_burn.cpp -- invariant tests for the controlledburn core
//
// The primary invariant throughout: total burned coverage (full cells
// from runs + fractions from edges, times cell area) equals the exact
// polygon area. This is what "exact coverage fractions" promises.
//
// Copyright (c) 2025 Michael Sumner
// Licensed under Apache License 2.0

#include "controlledburn/controlledburn.hpp"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <vector>

using namespace controlledburn;

static int failures = 0;

#define CHECK(cond) do { \
    if (!(cond)) { \
        std::printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); \
        failures++; \
    } \
} while (0)

#define CHECK_NEAR(a, b, tol) do { \
    double aa = (a), bb = (b); \
    if (std::fabs(aa - bb) > (tol)) { \
        std::printf("FAIL %s:%d: %s (%.12g) != %s (%.12g), diff %.3g\n", \
                    __FILE__, __LINE__, #a, aa, #b, bb, std::fabs(aa - bb)); \
        failures++; \
    } \
} while (0)

// Total covered area implied by the sparse output, in CRS units.
static double covered_area(const BurnResult& r, const GridSpec& gs) {
    double cell = gs.dx() * gs.dy();
    double total = 0.0;
    for (const auto& run : r.runs) {
        total += cell * (run.col_end - run.col_start + 1);
    }
    for (const auto& e : r.edges) {
        total += cell * e.fraction;
    }
    return total;
}

static Geometry make_polygon(std::vector<CoordSeq> rings) {
    Geometry g;
    g.kind = GeomKind::Polygon;
    g.polygons.push_back({std::move(rings)});
    return g;
}

// ---- Tests ----

// Axis-aligned rectangle exactly on cell boundaries: pure runs, no edges.
static void test_aligned_rectangle() {
    GridSpec gs{0, 0, 10, 10, 10, 10};
    // Rectangle covering cells cols 3..6, rows (y 4..8 -> grid rows 3..6)
    Geometry g = make_polygon({{{2, 4}, {6, 4}, {6, 8}, {2, 8}, {2, 4}}});

    BurnResult r = burn({g}, gs);
    CHECK(r.edges.empty());
    CHECK(!r.runs.empty());
    CHECK_NEAR(covered_area(r, gs), 16.0, 1e-9);
    for (const auto& run : r.runs) CHECK(run.id == 1);
}

// Half-cell offset rectangle: edge cells with fraction 0.5 / 0.25.
static void test_offset_rectangle() {
    GridSpec gs{0, 0, 10, 10, 10, 10};
    Geometry g = make_polygon({{{2.5, 4.5}, {6.5, 4.5}, {6.5, 8.5}, {2.5, 8.5}, {2.5, 4.5}}});

    BurnResult r = burn({g}, gs);
    CHECK(!r.edges.empty());
    CHECK_NEAR(covered_area(r, gs), 16.0, 1e-6);

    // Corner cells are quarter-covered, side cells half-covered.
    int quarters = 0, halves = 0;
    for (const auto& e : r.edges) {
        if (std::fabs(e.fraction - 0.25f) < 1e-6f) quarters++;
        if (std::fabs(e.fraction - 0.5f) < 1e-6f) halves++;
    }
    CHECK(quarters == 4);
    CHECK(halves == 12);
}

// Triangle: coverage sum equals exact area.
static void test_triangle() {
    GridSpec gs{0, 0, 100, 100, 37, 41}; // deliberately awkward cell sizes
    Geometry g = make_polygon({{{13.3, 17.7}, {88.1, 22.4}, {41.9, 79.2}, {13.3, 17.7}}});
    double area = std::fabs(signed_area(g.polygons[0].rings[0]));

    BurnResult r = burn({g}, gs);
    CHECK_NEAR(covered_area(r, gs), area, area * 1e-4);
}

// Polygon with a hole: coverage = outer area minus hole area.
// Hole orientation should not matter (normalised internally).
static void test_hole() {
    GridSpec gs{0, 0, 10, 10, 20, 20};
    CoordSeq outer = {{1, 1}, {9, 1}, {9, 9}, {1, 9}, {1, 1}};          // CCW, area 64
    CoordSeq hole_ccw = {{3, 3}, {7, 3}, {7, 7}, {3, 7}, {3, 3}};       // area 16
    CoordSeq hole_cw = {{3, 3}, {3, 7}, {7, 7}, {7, 3}, {3, 3}};

    for (const CoordSeq& hole : {hole_ccw, hole_cw}) {
        Geometry g = make_polygon({outer, hole});
        BurnResult r = burn({g}, gs);
        CHECK_NEAR(covered_area(r, gs), 64.0 - 16.0, 1e-6);
    }
}

// Disjoint multipolygon: components must not bleed winding into each other.
static void test_multipolygon() {
    GridSpec gs{0, 0, 10, 10, 10, 10};
    Geometry g;
    g.kind = GeomKind::MultiPolygon;
    g.polygons.push_back({{{{1, 1}, {3, 1}, {3, 3}, {1, 3}, {1, 1}}}});   // area 4
    g.polygons.push_back({{{{6, 6}, {9, 6}, {9, 9}, {6, 9}, {6, 6}}}});   // area 9

    BurnResult r = burn({g}, gs);
    CHECK_NEAR(covered_area(r, gs), 13.0, 1e-6);
    for (const auto& run : r.runs) CHECK(run.id == 1); // one geometry, one id
}

// Polygon larger than the grid: every cell fully covered.
// Exercises the padding-column winding fix (edges entirely outside the
// grid must still contribute winding deltas to grid rows).
static void test_beyond_extent() {
    GridSpec gs{0, 0, 10, 10, 5, 5};
    Geometry g = make_polygon({{{-100, -100}, {100, -100}, {100, 100}, {-100, 100}, {-100, -100}}});

    BurnResult r = burn({g}, gs);
    CHECK(r.edges.empty());
    CHECK_NEAR(covered_area(r, gs), 100.0, 1e-9); // whole grid
}

// Polygon straddling one grid edge: partial rows, correct clipped area.
static void test_straddle() {
    GridSpec gs{0, 0, 10, 10, 10, 10};
    Geometry g = make_polygon({{{-5, 2.5}, {4.5, 2.5}, {4.5, 6.5}, {-5, 6.5}, {-5, 2.5}}});
    // Clipped to grid: x in [0, 4.5], y in [2.5, 6.5] -> area 18

    BurnResult r = burn({g}, gs);
    CHECK_NEAR(covered_area(r, gs), 18.0, 1e-6);
}

// Diagonal line: per-cell lengths sum to the full line length.
static void test_line_length() {
    GridSpec gs{0, 0, 10, 10, 10, 10};
    Geometry g;
    g.kind = GeomKind::LineString;
    g.lines.push_back({{0.5, 0.5}, {9.5, 7.5}});
    double len = std::hypot(9.0, 7.0);

    BurnResult r = burn({g}, gs);
    CHECK(r.runs.empty());
    CHECK(r.edges.empty());
    double total = 0.0;
    for (const auto& l : r.lines) total += l.length;
    CHECK_NEAR(total, len, 1e-4);
}

// Points: correct cell binning, out-of-extent dropped.
static void test_points() {
    GridSpec gs{0, 0, 10, 10, 10, 10};
    Geometry g;
    g.kind = GeomKind::MultiPoint;
    g.points.push_back({{0.5, 9.5}});   // top-left cell: row 1, col 1
    g.points.push_back({{9.5, 0.5}});   // bottom-right cell: row 10, col 10
    g.points.push_back({{15.0, 5.0}});  // outside: dropped

    BurnResult r = burn({g}, gs);
    CHECK(r.points.size() == 2);
    CHECK(r.points[0].row == 1 && r.points[0].col == 1);
    CHECK(r.points[1].row == 10 && r.points[1].col == 10);
}

// Hand-encoded little-endian WKB POLYGON round trip through burn_wkb.
static void test_wkb() {
    // POLYGON ((2 4, 6 4, 6 8, 2 8, 2 4))
    std::vector<uint8_t> wkb;
    auto push_u32 = [&](uint32_t v) {
        for (int i = 0; i < 4; i++) wkb.push_back(static_cast<uint8_t>(v >> (8 * i)));
    };
    auto push_f64 = [&](double d) {
        uint64_t v;
        std::memcpy(&v, &d, 8);
        for (int i = 0; i < 8; i++) wkb.push_back(static_cast<uint8_t>(v >> (8 * i)));
    };
    wkb.push_back(1);    // little endian
    push_u32(3);         // Polygon
    push_u32(1);         // 1 ring
    push_u32(5);         // 5 coords
    double coords[5][2] = {{2, 4}, {6, 4}, {6, 8}, {2, 8}, {2, 4}};
    for (auto& c : coords) { push_f64(c[0]); push_f64(c[1]); }

    GridSpec gs{0, 0, 10, 10, 10, 10};
    BurnResult r = burn_wkb({{wkb.data(), wkb.size()}}, gs);
    CHECK(r.notes.empty());
    CHECK_NEAR(covered_area(r, gs), 16.0, 1e-9);

    // Malformed input is a note, not a crash.
    std::vector<uint8_t> bad = {1, 3, 0};
    BurnResult rb = burn_wkb({{bad.data(), bad.size()}}, gs);
    CHECK(rb.notes.size() == 1);
    CHECK(rb.runs.empty());
}

// bbox_wkb: union envelope over WKB blobs, skipping nulls and unparseable.
static void test_bbox_wkb() {
    auto poly_wkb = [](double x0, double y0, double x1, double y1) {
        std::vector<uint8_t> wkb;
        auto push_u32 = [&](uint32_t v) {
            for (int i = 0; i < 4; i++)
                wkb.push_back(static_cast<uint8_t>(v >> (8 * i)));
        };
        auto push_f64 = [&](double d) {
            uint64_t v;
            std::memcpy(&v, &d, 8);
            for (int i = 0; i < 8; i++)
                wkb.push_back(static_cast<uint8_t>(v >> (8 * i)));
        };
        wkb.push_back(1);    // little endian
        push_u32(3);         // Polygon
        push_u32(1);         // 1 ring
        push_u32(5);         // 5 coords
        double c[5][2] = {{x0, y0}, {x1, y0}, {x1, y1}, {x0, y1}, {x0, y0}};
        for (auto& p : c) { push_f64(p[0]); push_f64(p[1]); }
        return wkb;
    };
    std::vector<uint8_t> a = poly_wkb(2, 4, 6, 8);
    std::vector<uint8_t> b = poly_wkb(-1, 0, 3, 5);
    std::vector<uint8_t> bad = {1, 3, 0};

    std::vector<WKBSpan> spans = {
        {a.data(), a.size()},
        {nullptr, 0},               // null span skipped
        {b.data(), b.size()},
        {bad.data(), bad.size()},   // unparseable skipped
    };
    BBox bb = bbox_wkb(spans);
    CHECK(bb.valid);
    CHECK_NEAR(bb.xmin, -1.0, 1e-12);
    CHECK_NEAR(bb.ymin, 0.0, 1e-12);
    CHECK_NEAR(bb.xmax, 6.0, 1e-12);
    CHECK_NEAR(bb.ymax, 8.0, 1e-12);

    // No parseable geometry -> invalid envelope.
    std::vector<WKBSpan> none = {{nullptr, 0}, {bad.data(), bad.size()}};
    CHECK(!bbox_wkb(none).valid);
}

// Materialize: fasterize-style consumption of the sparse output.
static void test_materialize() {
    GridSpec gs{0, 0, 10, 10, 10, 10};
    Geometry a = make_polygon({{{2, 4}, {6, 4}, {6, 8}, {2, 8}, {2, 4}}});   // 16 cells
    Geometry b = make_polygon({{{4, 4}, {8, 4}, {8, 8}, {4, 8}, {4, 4}}});   // 16 cells, overlaps 8

    BurnResult r = burn({a, b}, gs);

    std::vector<double> buf(100, std::nan(""));
    MaterializeOptions opts;
    opts.fn = PixelFn::Sum;
    materialize(r, buf.data(), 10, 10, {1.0, 1.0}, opts);

    double total = 0.0;
    int touched = 0, twos = 0;
    for (double v : buf) {
        if (!std::isnan(v)) { total += v; touched++; if (v == 2.0) twos++; }
    }
    CHECK(touched == 24);        // union of the two rectangles
    CHECK(twos == 8);            // overlap cells counted twice
    CHECK_NEAR(total, 32.0, 1e-9);

    // Fraction policy conserves exact area for a half-offset polygon.
    Geometry c = make_polygon({{{2.5, 4.5}, {6.5, 4.5}, {6.5, 8.5}, {2.5, 8.5}, {2.5, 4.5}}});
    BurnResult rc = burn({c}, gs);
    std::vector<double> buf2(100, std::nan(""));
    MaterializeOptions frac;
    frac.fn = PixelFn::Sum;
    frac.edge_policy = EdgePolicy::Fraction;
    materialize(rc, buf2.data(), 10, 10, {1.0}, frac);
    double total2 = 0.0;
    for (double v : buf2) if (!std::isnan(v)) total2 += v;
    CHECK_NEAR(total2, 16.0, 1e-6); // area / cell_area
}

// ---- Approx mode tests ----

// Aligned rectangle: approx mode should produce runs only, same total area.
static void test_approx_aligned_rectangle() {
    GridSpec gs{0, 0, 10, 10, 10, 10};
    Geometry g = make_polygon({{{2, 4}, {6, 4}, {6, 8}, {2, 8}, {2, 4}}});

    BurnResult r = burn({g}, gs, BurnMode::Approx);
    CHECK(r.edges.empty());
    CHECK(!r.runs.empty());
    CHECK_NEAR(covered_area(r, gs), 16.0, 1e-9);
}

// Offset rectangle: approx mode gives runs only (no edges), area is
// approximate -- cell centers inside get full cells, so total area
// differs from exact but is deterministic.
static void test_approx_offset_rectangle() {
    GridSpec gs{0, 0, 10, 10, 10, 10};
    // 2.5..6.5 x 4.5..8.5: cell centers at 0.5, 1.5, ..., 9.5
    // x: centers 2.5, 3.5, 4.5, 5.5 are inside (cols 3-6); 6.5 is on boundary
    // y: centers 4.5, 5.5, 6.5, 7.5 are inside (rows 3-6); 8.5 is on boundary
    // Expect 4*4=16 cells if boundary cells with center ON edge are excluded,
    // but the exact count depends on the center-rule decision for edge-on-center.
    Geometry g = make_polygon({{{2.5, 4.5}, {6.5, 4.5}, {6.5, 8.5}, {2.5, 8.5}, {2.5, 4.5}}});

    BurnResult r = burn({g}, gs, BurnMode::Approx);
    CHECK(r.edges.empty());  // no edges in approx mode
    CHECK(!r.runs.empty());
    // Area should be a whole number of cells (each run cell = 1 cell area)
    double area = covered_area(r, gs);
    double cell_area = gs.dx() * gs.dy();
    double n_cells = area / cell_area;
    CHECK_NEAR(n_cells, std::round(n_cells), 1e-9);
}

// Beyond-extent polygon: approx mode should still fill the entire grid.
static void test_approx_beyond_extent() {
    GridSpec gs{0, 0, 10, 10, 5, 5};
    Geometry g = make_polygon({{{-100, -100}, {100, -100}, {100, 100}, {-100, 100}, {-100, -100}}});

    BurnResult r = burn({g}, gs, BurnMode::Approx);
    CHECK(r.edges.empty());
    CHECK_NEAR(covered_area(r, gs), 100.0, 1e-9);
}

// Hole: approx mode should still respect holes.
static void test_approx_hole() {
    GridSpec gs{0, 0, 10, 10, 20, 20};
    CoordSeq outer = {{1, 1}, {9, 1}, {9, 9}, {1, 9}, {1, 1}};
    CoordSeq hole = {{3, 3}, {7, 3}, {7, 7}, {3, 7}, {3, 3}};
    Geometry g = make_polygon({outer, hole});

    BurnResult r = burn({g}, gs, BurnMode::Approx);
    CHECK(r.edges.empty());
    // Approx area should be close to 48 (64 - 16), within a few cells
    double area = covered_area(r, gs);
    double cell_area = gs.dx() * gs.dy();
    CHECK(area > 40.0);
    CHECK(area < 56.0);
    // Should be a whole number of cells
    double n_cells = area / cell_area;
    CHECK_NEAR(n_cells, std::round(n_cells), 1e-9);
}

// Approx mode does not affect lines or points.
static void test_approx_line_unchanged() {
    GridSpec gs{0, 0, 10, 10, 10, 10};
    Geometry g;
    g.kind = GeomKind::LineString;
    g.lines.push_back({{0.5, 0.5}, {9.5, 7.5}});
    double len = std::hypot(9.0, 7.0);

    BurnResult r = burn({g}, gs, BurnMode::Approx);
    double total = 0.0;
    for (const auto& l : r.lines) total += l.length;
    CHECK_NEAR(total, len, 1e-4);
}

// Degenerate and empty inputs must not crash.
static void test_degenerate() {
    GridSpec gs{0, 0, 10, 10, 10, 10};

    Geometry empty;
    empty.kind = GeomKind::Polygon;
    BurnResult r1 = burn({empty}, gs);
    CHECK(r1.runs.empty() && r1.edges.empty());

    Geometry tiny = make_polygon({{{1, 1}, {2, 2}, {1, 1}}}); // < 4 coords
    BurnResult r2 = burn({tiny}, gs);
    CHECK(r2.runs.empty() && r2.edges.empty());

    Geometry sub = make_polygon({{{3.2, 3.2}, {3.4, 3.2}, {3.4, 3.4}, {3.2, 3.4}, {3.2, 3.2}}});
    BurnResult r3 = burn({sub}, gs); // polygon inside one cell
    CHECK(r3.runs.empty());
    CHECK(r3.edges.size() == 1);
    CHECK_NEAR(covered_area(r3, gs), 0.04, 1e-6);
}

int main() {
    test_aligned_rectangle();
    test_offset_rectangle();
    test_triangle();
    test_hole();
    test_multipolygon();
    test_beyond_extent();
    test_straddle();
    test_line_length();
    test_points();
    test_wkb();
    test_bbox_wkb();
    test_materialize();
    test_degenerate();

    test_approx_aligned_rectangle();
    test_approx_offset_rectangle();
    test_approx_beyond_extent();
    test_approx_hole();
    test_approx_line_unchanged();

    if (failures == 0) {
        std::printf("all tests passed\n");
        return 0;
    }
    std::printf("%d failure(s)\n", failures);
    return 1;
}
