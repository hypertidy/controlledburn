// Experiment 05: GDALRasterizeGeometries (cell-centre rule, and ALL_TOUCHED)
// against controlledburn's approx and coverage burns, on the same grid and
// geometries. Both sides end in the same dense Byte buffer, so the
// controlledburn timings include materialising its runs and edges into it.
// Cell agreement is reported alongside the times. GDAL's time is the
// GDALRasterizeGeometries call into a MEM dataset (zero-initialised when
// created, outside the timing); controlledburn's is burn plus materialise into
// a zeroed buffer (zeroing outside the timing).
#include <controlledburn/controlledburn.hpp>
#include <gdal.h>
#include <gdal_alg.h>
#include <ogr_api.h>
#include <cpl_conv.h>
#include <chrono>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <random>
#include <vector>

using namespace controlledburn;
using clk = std::chrono::steady_clock;
static double ms(clk::time_point a, clk::time_point b) { return std::chrono::duration<double, std::milli>(b - a).count(); }

static std::vector<uint8_t> ring_wkb(const std::vector<std::pair<double, double>>& pts) {
    std::vector<uint8_t> w(13 + pts.size() * 16);
    w[0] = 1;
    uint32_t t = 3, nr = 1, np = (uint32_t)pts.size();
    std::memcpy(&w[1], &t, 4), std::memcpy(&w[5], &nr, 4), std::memcpy(&w[9], &np, 4);
    for (size_t i = 0; i < pts.size(); ++i) std::memcpy(&w[13 + 16 * i], &pts[i], 16);
    return w;
}

struct Corpus { const char* name; int n; std::vector<std::vector<uint8_t>> wkb; };

static Corpus star(int N, int n) {
    std::vector<std::pair<double, double>> p;
    for (int i = 0; i <= N; ++i) {
        double a = 2 * M_PI * (i % N) / N, r = 0.45 * n * (1 + 0.08 * std::sin(37 * a) + 0.03 * std::sin(311 * a));
        p.push_back({n / 2.0 + 0.3 + r * std::cos(a), n / 2.0 + 0.7 + r * std::sin(a)});
    }
    return {"one star, 100k vertices", n, {ring_wkb(p)}};
}
static Corpus blobs(int count, int n) {
    std::mt19937_64 rng(42);
    std::uniform_real_distribution<double> u(0, 1);
    Corpus c{"20k random 24-gons, r 2-30 cells", n, {}};
    for (int k = 0; k < count; ++k) {
        double cx = u(rng) * n, cy = u(rng) * n, r = 2 + 28 * u(rng);
        std::vector<std::pair<double, double>> p;
        for (int i = 0; i <= 24; ++i) {
            double a = 2 * M_PI * (i % 24) / 24, rr = r * (0.7 + 0.3 * std::sin(5 * a + k));
            p.push_back({cx + rr * std::cos(a), cy + rr * std::sin(a)});
        }
        c.wkb.push_back(ring_wkb(p));
    }
    return c;
}

int main() {
    GDALAllRegister();
    std::printf("GDAL %s\n", GDALVersionInfo("RELEASE_NAME"));
    std::printf("%-34s %-34s %10s %12s %12s\n", "corpus", "method", "ms", "cells set", "vs GDAL");
    for (Corpus c : {star(100000, 10000), blobs(20000, 10000)}) {
        int n = c.n;
        GridSpec g{0, 0, (double)n, (double)n, n, n};
        std::vector<OGRGeometryH> geoms;
        for (auto& w : c.wkb) {
            OGRGeometryH h = nullptr;
            OGR_G_CreateFromWkb(w.data(), nullptr, &h, (int)w.size());
            geoms.push_back(h);
        }
        std::vector<WKBSpan> spans;
        for (auto& w : c.wkb) spans.push_back({w.data(), w.size()});

        auto gdal_burn = [&](bool all_touched, std::vector<uint8_t>& out) {
            out.assign((size_t)n * n, 0);
            GDALDatasetH ds = GDALCreate(GDALGetDriverByName("MEM"), "", n, n, 1, GDT_Byte, nullptr);
            double gt[6] = {0, 1, 0, (double)n, 0, -1};
            GDALSetGeoTransform(ds, gt);
            int band = 1;
            std::vector<double> burn(geoms.size(), 1.0);
            const char* o[] = {all_touched ? "ALL_TOUCHED=TRUE" : nullptr, nullptr};
            auto t0 = clk::now();
            GDALRasterizeGeometries(ds, 1, &band, (int)geoms.size(), geoms.data(), nullptr, nullptr, burn.data(),
                                    (char**)o, nullptr, nullptr);
            auto t1 = clk::now();
            GDALRasterIO(GDALGetRasterBand(ds, 1), GF_Read, 0, 0, n, n, out.data(), n, n, GDT_Byte, 0, 0);
            GDALClose(ds);
            return ms(t0, t1);
        };
        auto cb_burn = [&](BurnMode mode, unsigned threads, std::vector<uint8_t>& out) {
            out.assign((size_t)n * n, 0);
            auto t0 = clk::now();
            BurnResult r = burn_wkb_par(spans, g, BurnOptions(mode), threads);
            for (auto& run : r.runs) std::memset(&out[(size_t)run.row * n + run.col_start], 1, run.col_end - run.col_start);
            for (auto& e : r.edges)
                if (e.fraction > 0) out[(size_t)e.row * n + e.col] = 1;
            auto t1 = clk::now();
            return ms(t0, t1);
        };
        auto count = [](const std::vector<uint8_t>& a) { long s = 0; for (auto v : a) s += v; return s; };
        auto differ = [](const std::vector<uint8_t>& a, const std::vector<uint8_t>& b) {
            long d = 0; for (size_t i = 0; i < a.size(); ++i) d += a[i] != b[i]; return d; };

        std::vector<uint8_t> gc, ga, out;
        auto best = [](auto f) { double b = 1e30; for (int i = 0; i < 3; ++i) b = std::min(b, f()); return b; };
        double t_gc = best([&] { return gdal_burn(false, gc); });
        std::printf("%-34s %-34s %10.1f %12ld %12s\n", c.name, "GDAL, cell centre", t_gc, count(gc), "-");
        for (unsigned th : {1u, 8u}) {
            double t = best([&] { return cb_burn(BurnMode::Approx, th, out); });
            char m[64]; std::snprintf(m, sizeof m, "controlledburn approx, %u thread%s", th, th > 1 ? "s" : "");
            std::printf("%-34s %-34s %10.1f %12ld %10ld %s\n", "", m, t, count(out), differ(out, gc), "cells differ");
        }
        double t_ga = best([&] { return gdal_burn(true, ga); });
        std::printf("%-34s %-34s %10.1f %12ld %12s\n", "", "GDAL, ALL_TOUCHED", t_ga, count(ga), "-");
        for (unsigned th : {1u, 8u}) {
            double t = best([&] { return cb_burn(BurnMode::Coverage, th, out); });
            char m[64]; std::snprintf(m, sizeof m, "controlledburn fraction > 0, %u thr", th);
            std::printf("%-34s %-34s %10.1f %12ld %10ld %s\n", "", m, t, count(out), differ(out, ga), "cells differ");
        }
        // Classify the ALL_TOUCHED disagreements by the exact area of the cell
        // inside the union of the geometries that reach it (OGR/GEOS).
        long g_only = 0, cb_only = 0, g_only_zero = 0, cb_only_pos = 0, checked = 0;
        double max_g_only_area = 0;
        for (size_t i = 0; i < out.size(); ++i) {
            if (out[i] == ga[i]) continue;
            bool gdal_set = ga[i] != 0;
            (gdal_set ? g_only : cb_only)++;
            if (checked >= 400) continue;
            ++checked;
            int row = (int)(i / n), col = (int)(i % n);
            double x0 = col, y1 = n - row, x1 = col + 1, y0 = n - row - 1;
            char wkt[256];
            std::snprintf(wkt, sizeof wkt, "POLYGON ((%d %d, %d %d, %d %d, %d %d, %d %d))", (int)x0, (int)y0, (int)x1,
                          (int)y0, (int)x1, (int)y1, (int)x0, (int)y1, (int)x0, (int)y0);
            OGRGeometryH cell = nullptr;
            char* p = wkt;
            OGR_G_CreateFromWkt(&p, nullptr, &cell);
            double area = 0;
            for (auto h : geoms) {
                OGREnvelope e;
                OGR_G_GetEnvelope(h, &e);
                if (e.MaxX < x0 || e.MinX > x1 || e.MaxY < y0 || e.MinY > y1) continue;
                OGRGeometryH ix = OGR_G_Intersection(h, cell);
                if (ix) { area = std::max(area, OGR_G_Area(ix)); OGR_G_DestroyGeometry(ix); }
            }
            OGR_G_DestroyGeometry(cell);
            if (gdal_set) { if (area == 0) ++g_only_zero; max_g_only_area = std::max(max_g_only_area, area); }
            else if (area > 0) ++cb_only_pos;
        }
        std::printf("%-34s ALL_TOUCHED disagreements: %ld set only by GDAL, %ld only by controlledburn; of the first %ld checked,\n"
                    "%-34s %ld GDAL-only cells have zero intersection area (largest GDAL-only area %.3g), %ld cb-only cells have area > 0\n",
                    "", g_only, cb_only, checked, "", g_only_zero, max_g_only_area, cb_only_pos);
        for (auto h : geoms) OGR_G_DestroyGeometry(h);
    }
    return 0;
}
