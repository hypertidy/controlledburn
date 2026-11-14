// test_canon_ops.cpp -- the core's ports of the Rust crate's operations,
// record for record against the crate.
//
// fixtures/canon-ops.txt is written by tools/canon-ops (controlledburn
// crate 0.3.2) from the inputs of fixtures/canon-tables.txt: burns with
// the default options (coalesced runs, merged parts) and others, burns
// onto windows, slice_tiles, merge_cells, crop and relate_cells, in both
// modes. This test recomputes each with the C++ core and requires the
// same records in the same order, fractions and lengths compared as bits.
// It also checks that the threaded burn and slicing equal the serial ones.
//
// Copyright (c) 2026 Michael Sumner
// Licensed under Apache License 2.0

#include "controlledburn/controlledburn.hpp"

#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

using namespace controlledburn;

namespace {

struct Op {
    std::string header;  // the op line, without "op "
    std::vector<std::string> records;
};

struct Case {
    std::string name;
    GridSpec grid{};
    std::vector<std::vector<uint8_t>> wkb;
    std::vector<Op> ops;
};

double from_bits(const std::string& h) {
    uint64_t u = std::stoull(h, nullptr, 16);
    double d;
    std::memcpy(&d, &u, sizeof d);
    return d;
}

std::string bits(double d) {
    uint64_t u;
    std::memcpy(&u, &d, sizeof u);
    char buf[24];
    std::snprintf(buf, sizeof buf, "%016llx", static_cast<unsigned long long>(u));
    return buf;
}

std::vector<uint8_t> unhex(const std::string& h) {
    std::vector<uint8_t> out(h.size() / 2);
    for (size_t i = 0; i < out.size(); i++)
        out[i] = static_cast<uint8_t>(std::stoul(h.substr(2 * i, 2), nullptr, 16));
    return out;
}

std::string f32_hex(float f) {
    uint32_t u;
    std::memcpy(&u, &f, sizeof u);
    char buf[16];
    std::snprintf(buf, sizeof buf, "%08x", u);
    return buf;
}

std::string grid_line(const GridSpec& g) {
    return bits(g.xmin) + " " + bits(g.ymin) + " " + bits(g.xmax) + " " + bits(g.ymax) + " " +
           std::to_string(g.ncol) + " " + std::to_string(g.nrow);
}

void records(std::vector<std::string>& out, const BurnResult& r, bool notes) {
    for (const auto& x : r.runs)
        out.push_back("run " + std::to_string(x.row) + " " + std::to_string(x.col_start) + " " +
                      std::to_string(x.col_end) + " " + std::to_string(x.id));
    for (const auto& x : r.edges)
        out.push_back("edge " + std::to_string(x.row) + " " + std::to_string(x.col) + " " +
                      f32_hex(x.fraction) + " " + std::to_string(x.id));
    for (const auto& x : r.lines)
        out.push_back("line " + std::to_string(x.row) + " " + std::to_string(x.col) + " " +
                      f32_hex(x.length) + " " + std::to_string(x.id));
    for (const auto& x : r.points)
        out.push_back("point " + std::to_string(x.row) + " " + std::to_string(x.col) + " " +
                      std::to_string(x.id));
    if (notes)
        for (const auto& n : r.notes) out.push_back("note " + std::to_string(n.geom_index));
}

std::vector<Case> read_cases(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        std::printf("cannot open %s\n", path.c_str());
        return {};
    }
    std::vector<Case> cases;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        std::string tag;
        ss >> tag;
        if (tag == "case") {
            cases.emplace_back();
            ss >> cases.back().name;
        } else if (tag == "grid" && cases.back().ops.empty()) {
            std::string a, b, c, d;
            int nc, nr;
            ss >> a >> b >> c >> d >> nc >> nr;
            cases.back().grid = GridSpec{from_bits(a), from_bits(b), from_bits(c), from_bits(d), nc, nr};
        } else if (tag == "wkb") {
            std::string h;
            ss >> h;
            cases.back().wkb.push_back(unhex(h));
        } else if (tag == "op") {
            cases.back().ops.push_back(Op{line.substr(3), {}});
        } else if (tag == "end") {
        } else {
            cases.back().ops.back().records.push_back(line);
        }
    }
    return cases;
}

BurnMode mode_of(const std::string& m) { return m == "approx" ? BurnMode::Approx : BurnMode::Coverage; }

// The C++ core's answer for one op.
std::vector<std::string> run_op(const Case& c, const std::vector<Geometry>& geoms,
                                const std::vector<WKBSpan>& spans, const std::string& header) {
    std::istringstream ss(header);
    std::string kind, m;
    ss >> kind >> m;
    BurnMode mode = mode_of(m);
    const GridSpec& g = c.grid;
    BurnOptions opts(mode);
    std::vector<std::string> got;
    if (kind == "burn") {
        std::string variant;
        ss >> variant;
        BurnOptions o = opts;
        if (variant == "merge_only") o = BurnOptions::parity(mode).with_merge_parts(true);
        if (variant == "coalesce_only")
            o = BurnOptions::parity(mode).with_coalesce_runs(true).with_id_offset(7);
        records(got, burn_wkb(spans, g, o), true);
    } else if (kind == "window") {
        uint32_t r0, r1, c0, c1;
        ss >> r0 >> r1 >> c0 >> c1;
        records(got, burn_window(geoms, g, r0, r1, c0, c1, opts), true);
        // the WKB entry point gives the same
        std::vector<std::string> w;
        records(w, burn_wkb_window(spans, g, r0, r1, c0, c1, opts), true);
        if (w != got) got.push_back("burn_wkb_window differs from burn_window");
    } else if (kind == "tiles") {
        uint32_t tc, tr;
        ss >> tc >> tr;
        BurnResult whole = burn(geoms, g, opts);
        for (const auto& t : slice_tiles(whole, g, tc, tr)) {
            got.push_back("tile " + std::to_string(t.tile_row) + " " + std::to_string(t.tile_col) + " " +
                          std::to_string(t.row_off) + " " + std::to_string(t.col_off) + " " +
                          grid_line(t.grid));
            records(got, t.result, false);
        }
    } else if (kind == "merge_cells") {
        BurnResult r = burn(geoms, g, BurnOptions::parity(mode));
        merge_cells(r);
        records(got, r, false);
    } else if (kind == "crop") {
        std::string a, b, cc, d;
        ss >> a >> b >> cc >> d;
        BurnResult whole = burn(geoms, g, opts);
        BurnResult sub;
        GridSpec sg{};
        if (crop(whole, g, from_bits(a), from_bits(b), from_bits(cc), from_bits(d), sub, sg)) {
            got.push_back("subgrid " + grid_line(sg));
            records(got, sub, false);
        }
    } else if (kind == "relate") {
        BurnResult whole = burn(geoms, g, opts);
        for (const auto& rel : relate_cells(whole, geoms, g, opts)) {
            got.push_back("rel " + std::to_string(rel.row) + " " + std::to_string(rel.col) + " " +
                          std::to_string(rel.id) + " " + rel.matrix.str());
        }
    } else {
        got.push_back("unknown op " + kind);
    }
    return got;
}

int failures = 0;

void check(bool ok, const char* what) {
    if (!ok) {
        std::printf("FAIL: %s\n", what);
        failures++;
    }
}

bool same(const BurnResult& a, const BurnResult& b) {
    std::vector<std::string> x, y;
    records(x, a, true);
    records(y, b, true);
    return x == y;
}

// A ring of k vertices around (cx, cy), radius wobbled.
CoordSeq ring(double cx, double cy, double r, int k, double phase) {
    CoordSeq v;
    for (int i = 0; i < k; i++) {
        double a = 2.0 * M_PI * i / k;
        double rr = r * (0.6 + 0.4 * std::fabs(std::sin(phase + 1.7 * i)));
        v.push_back({cx + rr * std::cos(a), cy + rr * std::sin(a)});
    }
    v.push_back(v[0]);
    return v;
}

// WKB of `g` in either byte order, with `extra` (0-2) zero ordinates after
// x and y: ISO codes (1000, 3000), or EWKB flags and an SRID when `ewkb`.
struct WkbWriter {
    bool big;
    size_t extra;
    bool ewkb;
    std::vector<uint8_t> b;

    void u32(uint32_t v) {
        for (int i = 0; i < 4; i++)
            b.push_back(static_cast<uint8_t>(big ? v >> (8 * (3 - i)) : v >> (8 * i)));
    }
    void f64(double d) {
        uint64_t u;
        std::memcpy(&u, &d, 8);
        for (int i = 0; i < 8; i++)
            b.push_back(static_cast<uint8_t>(big ? u >> (8 * (7 - i)) : u >> (8 * i)));
    }
    void head(uint32_t base, bool top) {
        b.push_back(big ? 0 : 1);
        if (ewkb) {
            uint32_t flags = (extra >= 1 ? 0x80000000u : 0) | (extra >= 2 ? 0x40000000u : 0);
            u32(base | flags | (top ? 0x20000000u : 0));
            if (top) u32(4326);
        } else {
            u32(base + (extra == 0 ? 0 : extra == 1 ? 1000 : 3000));
        }
    }
    void pt(const Coord& c) {
        f64(c.x);
        f64(c.y);
        for (size_t i = 0; i < extra; i++) f64(0.0);
    }
    void seq(const CoordSeq& s) {
        u32(static_cast<uint32_t>(s.size()));
        for (const auto& c : s) pt(c);
    }
    void poly(const Polygon& p) {
        u32(static_cast<uint32_t>(p.rings.size()));
        for (const auto& r : p.rings) seq(r);
    }
    std::vector<uint8_t> write(const Geometry& g) {
        b.clear();
        switch (g.kind) {
        case GeomKind::Point: head(1, true); pt(g.points[0][0]); break;
        case GeomKind::LineString: head(2, true); seq(g.lines[0]); break;
        case GeomKind::Polygon: head(3, true); poly(g.polygons[0]); break;
        case GeomKind::MultiPoint:
            head(4, true);
            u32(static_cast<uint32_t>(g.points.size()));
            for (const auto& p : g.points) { head(1, false); pt(p[0]); }
            break;
        case GeomKind::MultiLineString:
            head(5, true);
            u32(static_cast<uint32_t>(g.lines.size()));
            for (const auto& l : g.lines) { head(2, false); seq(l); }
            break;
        case GeomKind::MultiPolygon:
            head(6, true);
            u32(static_cast<uint32_t>(g.polygons.size()));
            for (const auto& p : g.polygons) { head(3, false); poly(p); }
            break;
        }
        return b;
    }
};

// burn_wkb_par equals burn_wkb, with the blobs of `geoms` in both byte
// orders, ISO and EWKB, with Z and M, and with empty and unparseable
// blobs among them; and each blob's work read from its bytes equals the
// work of the parsed geometry.
void check_wkb_threads(std::vector<Geometry> geoms, const GridSpec& g) {
    // kinds the polygon-and-line set lacks: a lone point, a NaN point (POINT
    // EMPTY), an empty polygon, and a multi-line whose last coordinate alone
    // reaches higher (38 coordinates, so the sample's step of 3 misses none
    // but the last)
    for (int i = 0; i < 40; i++) {
        double x = 5.0 + 9.0 * i, y = 20.0 + 4.0 * i;
        Geometry p; p.kind = GeomKind::Point; p.points.push_back({{x, y}});
        Geometry e; e.kind = GeomKind::Point; e.points.push_back({{std::nan(""), y}});
        Geometry ep; ep.kind = GeomKind::Polygon; ep.polygons.push_back(Polygon{});
        Geometry ml; ml.kind = GeomKind::MultiLineString;
        ml.lines.push_back(ring(x, y, 2.0, 8, i));
        ml.lines.push_back(CoordSeq{});
        CoordSeq open;
        for (int k = 0; k < 38; k++) open.push_back({x + 0.2 * k, k == 37 ? y + 5.0 : y});
        ml.lines.push_back(open);
        for (auto* q : {&p, &e, &ep, &ml}) geoms.insert(geoms.begin() + 61 * i, *q);
    }
    std::vector<std::vector<uint8_t>> blobs;
    for (size_t k = 0; k < geoms.size(); k++) {
        if (k % 97 == 13) blobs.push_back({});                  // empty: skipped
        else if (k % 89 == 7) blobs.push_back({1, 3, 0, 0, 0, 7});  // truncated: a note
        else {
            WkbWriter w{k % 2 == 1, k % 3, k % 5 == 2, {}};
            blobs.push_back(w.write(geoms[k]));
        }
    }
    std::vector<WKBSpan> spans;
    for (const auto& b : blobs) spans.push_back(WKBSpan{b.data(), b.size()});
    spans[spans.size() / 2] = WKBSpan{nullptr, 0};             // null: skipped

    size_t total = 0, bad_work = 0;
    for (const auto& s : spans) {
        size_t w = detail::burn_work_wkb(s, g);
        total += w;
        size_t want = 1;
        if (s.data && s.size) {
            try { want = detail::burn_work(parse_wkb(s.data, s.size), g); } catch (const WKBParseError&) {}
        }
        if (w != want) bad_work++;
    }
    check(bad_work == 0, "work read from WKB equals work of the parsed geometry");
    check(total > 8 * 8192, "WKB burn is large enough for threads");

    for (BurnMode mode : {BurnMode::Coverage, BurnMode::Approx}) {
        for (bool merge : {true, false}) {
            BurnOptions o = BurnOptions(mode).with_merge_parts(merge).with_id_offset(3);
            BurnResult serial = burn_wkb(spans, g, o);
            check(!serial.notes.empty(), "the truncated blobs give notes");
            for (unsigned t : {2u, 5u, 0u}) {
                check(same(burn_wkb_par(spans, g, o, t), serial), "burn_wkb_par equals burn_wkb");
            }
        }
    }
    // errors as burn_wkb's: a bad grid, a negative offset, an id beyond int32
    GridSpec bad{0.0, 0.0, 10.0, 10.0, 0, 10};
    std::vector<std::pair<GridSpec, BurnOptions>> errs = {
        {bad, BurnOptions()},
        {g, BurnOptions().with_id_offset(-1)},
        {g, BurnOptions().with_id_offset(INT32_MAX)},
    };
    for (const auto& eo : errs) {
        std::string want, got;
        try { burn_wkb(spans, eo.first, eo.second); } catch (const std::exception& e) { want = e.what(); }
        try { burn_wkb_par(spans, eo.first, eo.second, 4); } catch (const std::exception& e) { got = e.what(); }
        check(!want.empty() && got == want, "burn_wkb_par errors as burn_wkb");
    }
    std::printf("wkb threads: %zu blobs, work %zu\n", spans.size(), total);
}

// The threaded burn and slicing equal the serial ones (inputs large
// enough that the threads are used).
void check_threads() {
    GridSpec g{0.0, 0.0, 370.0, 230.0, 1530, 910};
    std::vector<Geometry> geoms;
    for (int i = 0; i < 3000; i++) {
        double x = std::fmod(i * 37.31, 370.0), y = std::fmod(i * 17.77, 230.0);
        Geometry geo;
        switch (i % 4) {
        case 0:
            geo.kind = GeomKind::Polygon;
            geo.polygons.push_back(Polygon{{ring(x, y, 3.0, 40, i)}});
            break;
        case 1:
            geo.kind = GeomKind::MultiPolygon;
            geo.polygons.push_back(Polygon{{ring(x, y, 3.0, 30, i), ring(x, y, 0.9, 12, i)}});
            geo.polygons.push_back(Polygon{{ring(x + 2.4, y, 2.0, 20, i + 1)}});
            break;
        case 2:
            geo.kind = GeomKind::LineString;
            geo.lines.push_back(ring(x, y, 2.5, 25, i));
            break;
        default:
            geo.kind = GeomKind::MultiPoint;
            for (int k = 0; k < 7; k++) geo.points.push_back({{x + 0.7 * k, y - 0.4 * k}});
        }
        geoms.push_back(geo);
    }
    for (BurnMode mode : {BurnMode::Coverage, BurnMode::Approx}) {
        for (bool merge : {true, false}) {
            BurnOptions o = BurnOptions(mode).with_merge_parts(merge).with_id_offset(3);
            BurnResult serial = burn(geoms, g, o);
            for (unsigned t : {2u, 5u, 0u}) {
                check(same(burn_par(geoms, g, o, t), serial), "burn_par equals burn");
            }
        }
    }
    check_wkb_threads(geoms, g);
    BurnResult r = burn(geoms, g, BurnOptions());
    size_t n = r.runs.size() + r.edges.size() + r.lines.size() + r.points.size();
    // fine tiles over a dense burn, so slice_tiles_par takes its threaded path
    GridSpec fine{0.0, 0.0, 370.0, 230.0, 6120, 3640};
    BurnResult rf = burn(geoms, fine, BurnOptions());
    size_t nf = rf.runs.size() + rf.edges.size() + rf.lines.size() + rf.points.size();
    check(nf > 262144, "fine burn is large enough for threaded slicing");
    for (unsigned t : {2u, 7u, 0u}) {
        auto a = slice_tiles(rf, fine, 97, 61);
        auto b = slice_tiles_par(rf, fine, 97, 61, t);
        bool ok = a.size() == b.size();
        for (size_t i = 0; ok && i < a.size(); i++) {
            ok = a[i].tile_row == b[i].tile_row && a[i].tile_col == b[i].tile_col &&
                 same(a[i].result, b[i].result) && grid_line(a[i].grid) == grid_line(b[i].grid);
        }
        check(ok, "slice_tiles_par equals slice_tiles");
    }
    std::printf("threads: %zu geometries, %zu and %zu records\n", geoms.size(), n, nf);
}

} // namespace

int main(int argc, char** argv) {
    std::string dir = argc > 1 ? argv[1] : "fixtures";
    auto cases = read_cases(dir + "/canon-ops.txt");
    if (cases.size() < 100) {
        std::printf("FAIL: expected the canonical cases, read %zu\n", cases.size());
        return 1;
    }
    size_t nops = 0, nrec = 0;
    for (const auto& c : cases) {
        std::vector<WKBSpan> spans;
        std::vector<Geometry> geoms;
        for (const auto& b : c.wkb) {
            spans.push_back(WKBSpan{b.data(), b.size()});
            geoms.push_back(parse_wkb(b.data(), b.size()));
        }
        for (const auto& op : c.ops) {
            auto got = run_op(c, geoms, spans, op.header);
            nops++;
            nrec += op.records.size();
            if (got != op.records) {
                size_t i = 0;
                while (i < got.size() && i < op.records.size() && got[i] == op.records[i]) i++;
                std::printf("FAIL %s [%s]: %zu records, canon %zu; first difference at %zu: '%s' vs canon '%s'\n",
                            c.name.c_str(), op.header.c_str(), got.size(), op.records.size(), i,
                            i < got.size() ? got[i].c_str() : "(none)",
                            i < op.records.size() ? op.records[i].c_str() : "(none)");
                failures++;
            }
        }
    }
    std::printf("%zu cases, %zu ops, %zu canonical records\n", cases.size(), nops, nrec);
    check_threads();
    std::printf("%d failures\n", failures);
    return failures == 0 ? 0 : 1;
}
