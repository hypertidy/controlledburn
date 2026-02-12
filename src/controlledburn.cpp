// controlledburn.cpp — cpp11 bridge between R and vendored exactextract algorithm
//
// Accepts WKB geometry bytes and grid parameters from R, computes exact
// coverage fractions using the exactextract algorithm, and returns results
// in sparse two-table format (runs + edges).

#include <cpp11.hpp>
#include <cpp11/list.hpp>
#include <cpp11/integers.hpp>
#include <cpp11/doubles.hpp>
#include <cpp11/raws.hpp>
#include <cpp11/strings.hpp>

#include <cstdarg>

#include "libgeos.h"
#include "dense_to_sparse.h"

#include "exactextract/raster_cell_intersection.h"
#include "exactextract/grid.h"
#include "exactextract/box.h"
#include "exactextract/geos_utils.h"

// Suppress GEOS notice/error messages (CRAN requirement: no stdout/stderr)
static void geos_notice_handler(const char* /*fmt*/, ...) {}
static void geos_error_handler(const char* fmt, ...) {
    // Capture the error message for R
    char buf[1024];
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(buf, sizeof(buf), fmt, ap);
    va_end(ap);
    cpp11::stop("GEOS error: %s", buf);
}

// RAII wrapper for GEOS context handle
class GEOSContextGuard {
public:
    GEOSContextGuard() {
        ctx_ = GEOS_init_r();
        if (!ctx_) cpp11::stop("Failed to create GEOS context");
        GEOSContext_setNoticeHandler_r(ctx_, geos_notice_handler);
        GEOSContext_setErrorHandler_r(ctx_, geos_error_handler);
    }
    ~GEOSContextGuard() {
        if (ctx_) GEOS_finish_r(ctx_);
    }
    GEOSContextHandle_t get() const { return ctx_; }

    // non-copyable
    GEOSContextGuard(const GEOSContextGuard&) = delete;
    GEOSContextGuard& operator=(const GEOSContextGuard&) = delete;
private:
    GEOSContextHandle_t ctx_;
};

// RAII wrapper for GEOSGeometry*
class GEOSGeomGuard {
public:
    GEOSGeomGuard(GEOSContextHandle_t ctx, GEOSGeometry* g) : ctx_(ctx), g_(g) {}
    ~GEOSGeomGuard() { if (g_) GEOSGeom_destroy_r(ctx_, g_); }
    GEOSGeometry* get() const { return g_; }
    bool valid() const { return g_ != nullptr; }

    GEOSGeomGuard(const GEOSGeomGuard&) = delete;
    GEOSGeomGuard& operator=(const GEOSGeomGuard&) = delete;
private:
    GEOSContextHandle_t ctx_;
    GEOSGeometry* g_;
};

// Compute exact coverage fractions for a set of WKB geometries on a grid.
//
// Returns a list with two elements:
//   $runs: data.frame(row, col_start, col_end, id) — interior cells (weight ≈ 1)
//   $edges: data.frame(row, col, weight, id) — boundary cells (0 < weight < 1)
//
// Grid coordinates are 1-based, row-major (row 1 is the top row).
[[cpp11::register]]
cpp11::writable::list cpp_burn_sparse(
    cpp11::list wkb_list,
    double xmin, double ymin, double xmax, double ymax,
    int ncol, int nrow
) {
    if (ncol <= 0 || nrow <= 0) {
        cpp11::stop("ncol and nrow must be positive");
    }
    if (xmax <= xmin || ymax <= ymin) {
        cpp11::stop("Invalid extent: xmax must be > xmin, ymax must be > ymin");
    }

    double dx = (xmax - xmin) / ncol;
    double dy = (ymax - ymin) / nrow;

    // Set up the full raster grid
    exactextract::Grid<exactextract::bounded_extent> full_grid(
        exactextract::Box(xmin, ymin, xmax, ymax), dx, dy);

    GEOSContextGuard geos_guard;
    GEOSContextHandle_t ctx = geos_guard.get();

    // Create a WKB reader (reused across all geometries)
    GEOSWKBReader* wkb_reader = GEOSWKBReader_create_r(ctx);
    if (!wkb_reader) cpp11::stop("Failed to create WKB reader");

    int n_geoms = wkb_list.size();

    // Accumulate results across all geometries
    std::vector<GridRun> all_runs;
    std::vector<GridEdge> all_edges;

    for (int k = 0; k < n_geoms; k++) {
        cpp11::raws wkb_raw(wkb_list[k]);
        int wkb_size = wkb_raw.size();
        if (wkb_size == 0) continue;

        const unsigned char* wkb_data = reinterpret_cast<const unsigned char*>(RAW(wkb_raw));

        // Deserialize WKB to GEOSGeometry
        GEOSGeomGuard geom(ctx,
            GEOSWKBReader_read_r(ctx, wkb_reader, wkb_data, static_cast<size_t>(wkb_size)));

        if (!geom.valid()) {
            cpp11::warning("Failed to parse WKB for geometry %d, skipping", k + 1);
            continue;
        }

        if (GEOSisEmpty_r(ctx, geom.get())) {
            continue;
        }

        // Compute coverage fractions (subgrid only)
        // raster_cell_intersection returns Raster<float> by value
        try {
            auto rci = exactextract::raster_cell_intersection(full_grid, ctx, geom.get());

            if (rci.grid().empty()) continue;

            // Compute offset of subgrid within full raster grid
            const auto& sub_grid = rci.grid();
            size_t row_off = static_cast<size_t>(
                std::round((full_grid.ymax() - sub_grid.ymax()) / dy));
            size_t col_off = static_cast<size_t>(
                std::round((sub_grid.xmin() - full_grid.xmin()) / dx));

            size_t sub_rows = sub_grid.rows();
            size_t sub_cols = sub_grid.cols();

            // Build a contiguous float buffer for dense_to_sparse
            std::vector<float> buf(sub_rows * sub_cols);
            for (size_t i = 0; i < sub_rows; i++) {
                for (size_t j = 0; j < sub_cols; j++) {
                    buf[i * sub_cols + j] = rci(i, j);
                }
            }

            int id = k + 1; // 1-based geometry id
            SparseResult sp = dense_to_sparse(
                buf.data(), sub_rows, sub_cols, row_off, col_off, id);

            // Append to accumulated results
            all_runs.insert(all_runs.end(), sp.runs.begin(), sp.runs.end());
            all_edges.insert(all_edges.end(), sp.edges.begin(), sp.edges.end());

        } catch (const std::exception& e) {
            cpp11::warning("Error processing geometry %d: %s, skipping", k + 1, e.what());
            continue;
        }
    }

    // Clean up WKB reader
    GEOSWKBReader_destroy_r(ctx, wkb_reader);

    // Build R data.frames
    // runs table
    size_t n_runs = all_runs.size();
    cpp11::writable::integers runs_row(n_runs);
    cpp11::writable::integers runs_col_start(n_runs);
    cpp11::writable::integers runs_col_end(n_runs);
    cpp11::writable::integers runs_id(n_runs);

    for (size_t i = 0; i < n_runs; i++) {
        runs_row[i] = all_runs[i].row;
        runs_col_start[i] = all_runs[i].col_start;
        runs_col_end[i] = all_runs[i].col_end;
        runs_id[i] = all_runs[i].id;
    }

    cpp11::writable::list runs_df(4);
    runs_df[0] = static_cast<SEXP>(runs_row);
    runs_df[1] = static_cast<SEXP>(runs_col_start);
    runs_df[2] = static_cast<SEXP>(runs_col_end);
    runs_df[3] = static_cast<SEXP>(runs_id);
    runs_df.attr("names") = cpp11::writable::strings({"row", "col_start", "col_end", "id"});
    runs_df.attr("class") = "data.frame";
    runs_df.attr("row.names") = cpp11::writable::integers({NA_INTEGER, -static_cast<int>(n_runs)});

    // edges table
    size_t n_edges = all_edges.size();
    cpp11::writable::integers edges_row(n_edges);
    cpp11::writable::integers edges_col(n_edges);
    cpp11::writable::doubles edges_weight(n_edges);
    cpp11::writable::integers edges_id(n_edges);

    for (size_t i = 0; i < n_edges; i++) {
        edges_row[i] = all_edges[i].row;
        edges_col[i] = all_edges[i].col;
        edges_weight[i] = static_cast<double>(all_edges[i].weight);
        edges_id[i] = all_edges[i].id;
    }

    cpp11::writable::list edges_df(4);
    edges_df[0] = static_cast<SEXP>(edges_row);
    edges_df[1] = static_cast<SEXP>(edges_col);
    edges_df[2] = static_cast<SEXP>(edges_weight);
    edges_df[3] = static_cast<SEXP>(edges_id);
    edges_df.attr("names") = cpp11::writable::strings({"row", "col", "weight", "id"});
    edges_df.attr("class") = "data.frame";
    edges_df.attr("row.names") = cpp11::writable::integers({NA_INTEGER, -static_cast<int>(n_edges)});

    // Return list with both tables
    cpp11::writable::list result(2);
    result[0] = static_cast<SEXP>(runs_df);
    result[1] = static_cast<SEXP>(edges_df);
    result.attr("names") = cpp11::writable::strings({"runs", "edges"});

    return result;
}
