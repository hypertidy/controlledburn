// Experiment 01: GEOS GridIntersectionFractions (exactextract in GEOS >= 3.14,
// what GDAL's `gdal raster zonal-stats --pixels fractional` calls) against the
// controlledburn C++ core, cell by cell, on two corpora:
//
//   - burnscars' oracle cases (fixtures/oracle in controlledburn-rs): 50
//     polygons with exactextract 0.3.0 coverage and GDAL's cell-centre
//     rasterize as external references;
//   - controlledburn's canon-tables cases (fixtures/canon-tables.txt): 240
//     generated cases on non-dyadic grids with vertices steered onto cell
//     corners, walls and centre lines (polygonal geometries only here).
//
// Then a timing run on one large polygon: GEOS fills a dense float buffer the
// size of the geometry's window; controlledburn returns sparse runs and edges.
//
// Build: see run.sh. Nothing here is part of the R package.

#include <controlledburn/controlledburn.hpp>
#include <geos_c.h>

#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <map>
#include <set>
#include <sstream>
#include <string>
#include <vector>

using namespace controlledburn;

static std::vector<uint8_t> unhex(const std::string& h) {
    std::vector<uint8_t> b(h.size() / 2);
    for (size_t i = 0; i < b.size(); ++i) b[i] = (uint8_t)std::stoul(h.substr(2 * i, 2), nullptr, 16);
    return b;
}
static double hexbits(const std::string& h) {
    uint64_t u = std::stoull(h, nullptr, 16);
    double d;
    std::memcpy(&d, &u, 8);
    return d;
}
static uint32_t wkb_type(const std::vector<uint8_t>& w) {
    uint32_t t;
    std::memcpy(&t, w.data() + 1, 4);
    if (w[0] == 0) t = __builtin_bswap32(t);
    return (t & 0xffff) % 1000;
}

struct Case {
    std::string name;
    GridSpec grid;
    std::vector<std::vector<uint8_t>> wkb;
};

// Dense fractions of one geometry from controlledburn (row 0 at ymax).
static std::vector<double> cb_dense(const std::vector<uint8_t>& wkb, const GridSpec& g, BurnMode mode) {
    std::vector<WKBSpan> spans{{wkb.data(), wkb.size()}};
    BurnResult r = burn_wkb(spans, g, BurnOptions(mode));
    std::vector<double> d((size_t)g.ncol * g.nrow, 0.0);
    for (auto& run : r.runs)
        for (int c = run.col_start; c < run.col_end; ++c) d[(size_t)run.row * g.ncol + c] += 1.0;
    for (auto& e : r.edges) d[(size_t)e.row * g.ncol + e.col] += e.fraction;
    return d;
}

static std::string g_current;

static std::vector<double> geos_dense(GEOSContextHandle_t h, const std::vector<uint8_t>& wkb, const GridSpec& g,
                                      bool& ok) {
    if (std::getenv("CB_TRACE")) { std::fprintf(stderr, "geos: %s\n", g_current.c_str()); std::fflush(stderr); }
    GEOSWKBReader* rd = GEOSWKBReader_create_r(h);
    GEOSGeometry* geom = GEOSWKBReader_read_r(h, rd, wkb.data(), wkb.size());
    GEOSWKBReader_destroy_r(h, rd);
    std::vector<float> buf((size_t)g.ncol * g.nrow, 0.0f);
    ok = geom && GEOSGridIntersectionFractions_r(h, geom, g.xmin, g.ymin, g.xmax, g.ymax, g.ncol, g.nrow, buf.data());
    if (geom) GEOSGeom_destroy_r(h, geom);
    return std::vector<double>(buf.begin(), buf.end());
}

struct Tally {
    long cases = 0, cells = 0, geos_fail = 0;
    long diff_1e6 = 0, diff_1e3 = 0, zero_vs_full = 0, full_vs_zero = 0;
    double max_abs = 0;
    std::string worst;
};

static void compare(const char* label, Tally& t, const std::string& name, const std::vector<double>& a,
                    const std::vector<double>& b, int ncol, FILE* detail) {
    for (size_t i = 0; i < a.size(); ++i) {
        double d = std::fabs(a[i] - b[i]);
        if (d > 1e-6) {
            ++t.diff_1e6;
            if (detail)
                std::fprintf(detail, "%s,%s,%zu,%zu,%.9g,%.9g\n", label, name.c_str(), i / ncol, i % ncol, a[i], b[i]);
        }
        if (d > 1e-3) ++t.diff_1e3;
        if (a[i] > 1 - 1e-6 && b[i] < 1e-6) ++t.full_vs_zero;
        if (a[i] < 1e-6 && b[i] > 1 - 1e-6) ++t.zero_vs_full;
        if (d > t.max_abs) { t.max_abs = d; t.worst = name; }
    }
    t.cells += (long)a.size();
}

static void report(const char* what, const Tally& t) {
    std::printf("%-44s cases %4ld  cells %8ld  |diff|>1e-6 %5ld  >1e-3 %5ld  full/0 %3ld  0/full %3ld  max %.3g (%s)%s\n",
                what, t.cases, t.cells, t.diff_1e6, t.diff_1e3, t.full_vs_zero, t.zero_vs_full, t.max_abs,
                t.worst.c_str(), t.geos_fail ? (" GEOS failures " + std::to_string(t.geos_fail)).c_str() : "");
}

static std::vector<std::string> split(const std::string& s, char c) {
    std::vector<std::string> out;
    std::stringstream ss(s);
    std::string x;
    while (std::getline(ss, x, c)) out.push_back(x);
    return out;
}

int main(int argc, char** argv) {
    if (argc < 3) {
        std::fprintf(stderr, "usage: compare <oracle_dir> <canon-tables.txt> [detail.csv]\n");
        return 2;
    }
    std::string odir = argv[1];
    FILE* detail = argc > 3 ? std::fopen(argv[3], "w") : nullptr;
    if (detail) std::fprintf(detail, "comparison,case,row,col,a,b\n");
    GEOSContextHandle_t h = GEOS_init_r();
    std::printf("GEOS %s\n\n", GEOSversion());

    // ---- oracle corpus -------------------------------------------------------
    std::vector<Case> oracle;
    {
        std::ifstream f(odir + "/cases.csv");
        std::string line;
        std::getline(f, line);
        while (std::getline(f, line)) {
            auto v = split(line, ',');
            Case c;
            c.name = v[0];
            c.wkb.push_back(unhex(v[1]));
            c.grid = GridSpec{std::stod(v[2]), std::stod(v[3]), std::stod(v[4]), std::stod(v[5]), std::stoi(v[6]),
                              std::stoi(v[7])};
            oracle.push_back(c);
        }
    }
    std::map<std::string, std::map<std::pair<int, int>, double>> ee;
    std::map<std::string, std::set<std::pair<int, int>>> gdal_approx;
    {
        std::ifstream f(odir + "/coverage.csv");
        std::string line;
        std::getline(f, line);
        while (std::getline(f, line)) {
            auto v = split(line, ',');
            ee[v[0]][{std::stoi(v[1]), std::stoi(v[2])}] = std::stod(v[3]);
        }
        std::ifstream g(odir + "/approx.csv");
        std::getline(g, line);
        while (std::getline(g, line)) {
            auto v = split(line, ',');
            gdal_approx[v[0]].insert({std::stoi(v[1]), std::stoi(v[2])});
        }
    }

    Tally t_cb_geos, t_geos_ee, t_cb_ee, t_approx;
    for (auto& c : oracle) {
        bool ok;
        g_current = c.name;
        auto cb = cb_dense(c.wkb[0], c.grid, BurnMode::Coverage);
        auto gs = geos_dense(h, c.wkb[0], c.grid, ok);
        std::vector<double> eed(cb.size(), 0.0);
        for (auto& kv : ee[c.name]) eed[(size_t)kv.first.first * c.grid.ncol + kv.first.second] = kv.second;
        if (!ok) ++t_cb_geos.geos_fail;
        compare("cb_vs_geos", t_cb_geos, c.name, cb, gs, c.grid.ncol, detail);
        compare("geos_vs_ee030", t_geos_ee, c.name, gs, eed, c.grid.ncol, detail);
        compare("cb_vs_ee030", t_cb_ee, c.name, cb, eed, c.grid.ncol, detail);
        ++t_cb_geos.cases, ++t_geos_ee.cases, ++t_cb_ee.cases;

        auto ap = cb_dense(c.wkb[0], c.grid, BurnMode::Approx);
        std::vector<double> gd(ap.size(), 0.0);
        for (auto& rc : gdal_approx[c.name]) gd[(size_t)rc.first * c.grid.ncol + rc.second] = 1.0;
        compare("cb_approx_vs_gdal", t_approx, c.name, ap, gd, c.grid.ncol, detail);
        ++t_approx.cases;
    }
    std::printf("oracle corpus (burnscars, 50 cases)\n");
    report("  coverage: controlledburn vs GEOS", t_cb_geos);
    report("  coverage: GEOS vs exactextract 0.3.0", t_geos_ee);
    report("  coverage: controlledburn vs exactextract 0.3.0", t_cb_ee);
    report("  approx:   controlledburn vs GDAL rasterize", t_approx);

    // ---- canon corpus --------------------------------------------------------
    std::vector<Case> canon;
    {
        std::ifstream f(argv[2]);
        std::string line;
        while (std::getline(f, line)) {
            std::stringstream ss(line);
            std::string tag;
            ss >> tag;
            if (tag == "case") {
                canon.emplace_back();
                ss >> canon.back().name;
            } else if (tag == "grid") {
                std::string a, b, c2, d;
                int nc, nr;
                ss >> a >> b >> c2 >> d >> nc >> nr;
                canon.back().grid = GridSpec{hexbits(a), hexbits(b), hexbits(c2), hexbits(d), nc, nr};
            } else if (tag == "wkb") {
                std::string hx;
                ss >> hx;
                canon.back().wkb.push_back(unhex(hx));
            }
        }
    }
    Tally t_canon;
    long n_geoms = 0;
    for (auto& c : canon)
        for (size_t k = 0; k < c.wkb.size(); ++k) {
            uint32_t ty = wkb_type(c.wkb[k]);
            if (ty != 3 && ty != 6) continue;
            bool ok;
            g_current = c.name + "#" + std::to_string(k);
            auto cb = cb_dense(c.wkb[k], c.grid, BurnMode::Coverage);
            auto gs = geos_dense(h, c.wkb[k], c.grid, ok);
            if (!ok) { ++t_canon.geos_fail; continue; }
            compare("canon_cb_vs_geos", t_canon, c.name + "#" + std::to_string(k), cb, gs, c.grid.ncol, detail);
            ++n_geoms;
        }
    t_canon.cases = n_geoms;
    std::printf("\ncanon corpus (controlledburn, polygonal geometries; 'cases' counts geometries)\n");
    report("  coverage: controlledburn vs GEOS", t_canon);

    // ---- timing --------------------------------------------------------------
    // A wavy star of N vertices on an n x n grid: perimeter-heavy, mostly interior.
    std::printf("\ntiming: one polygon, star of N vertices filling an n x n grid (best of 3)\n");
    std::printf("%8s %7s %12s %12s %10s %14s %14s\n", "N", "n", "cb ms", "GEOS ms", "ratio", "cb records", "GEOS buf MB");
    for (auto [N, n] : std::vector<std::pair<int, int>>{{1000, 1000}, {10000, 4000}, {100000, 10000}, {100000, 20000}}) {
        std::vector<uint8_t> w(1 + 4 + 4 + 4 + (size_t)(N + 1) * 16);
        size_t p = 0;
        w[p++] = 1;
        uint32_t t = 3, nr = 1, np = N + 1;
        std::memcpy(&w[p], &t, 4), p += 4;
        std::memcpy(&w[p], &nr, 4), p += 4;
        std::memcpy(&w[p], &np, 4), p += 4;
        for (int i = 0; i <= N; ++i) {
            double a = 2 * M_PI * (i % N) / N;
            double r = 0.45 * n * (1 + 0.08 * std::sin(37 * a) + 0.03 * std::sin(311 * a));
            double xy[2] = {n / 2.0 + 0.3 + r * std::cos(a), n / 2.0 + 0.7 + r * std::sin(a)};
            std::memcpy(&w[p], xy, 16), p += 16;
        }
        GridSpec g{0, 0, (double)n, (double)n, n, n};
        double best_cb = 1e30, best_g = 1e30;
        size_t recs = 0;
        for (int rep = 0; rep < 3; ++rep) {
            auto t0 = std::chrono::steady_clock::now();
            std::vector<WKBSpan> spans{{w.data(), w.size()}};
            BurnResult r = burn_wkb(spans, g, BurnOptions());
            auto t1 = std::chrono::steady_clock::now();
            recs = r.runs.size() + r.edges.size();
            best_cb = std::min(best_cb, std::chrono::duration<double, std::milli>(t1 - t0).count());
        }
        for (int rep = 0; rep < 3; ++rep) {
            auto t0 = std::chrono::steady_clock::now();
            bool ok;
            auto gs = geos_dense(h, w, g, ok);
            auto t1 = std::chrono::steady_clock::now();
            best_g = std::min(best_g, std::chrono::duration<double, std::milli>(t1 - t0).count());
        }
        std::printf("%8d %7d %12.1f %12.1f %10.1f %14zu %14.0f\n", N, n, best_cb, best_g, best_g / best_cb, recs,
                    4.0 * n * n / 1e6);
    }

    if (detail) std::fclose(detail);
    GEOS_finish_r(h);
    return 0;
}
