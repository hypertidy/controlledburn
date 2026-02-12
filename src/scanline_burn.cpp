// scanline_burn.cpp — scanline polygon rasterization with exact coverage fractions
//
// Item 6: edge case fixes:
// - Padding column winding: edges outside the grid (padding columns) now
//   contribute winding deltas to grid rows, fixing beyond-extent polygons.
// - Sweep start condition: prev_col > -2 allows runs after padding cells.
// - Analytical single-traversal coverage via perimeter_distance (Item 3).
// - No Cell class allocation (Item 2).
// - MULTIPOLYGON components processed independently (each gets own sweep).
//
// Copyright (c) 2025 Michael Sumner
// Licensed under Apache License 2.0

#include <cpp11.hpp>
#include <cpp11/list.hpp>
#include <cpp11/integers.hpp>
#include <cpp11/doubles.hpp>
#include <cpp11/raws.hpp>
#include <cpp11/strings.hpp>

#include <cstdarg>
#include <map>
#include <algorithm>
#include <cmath>

#include "libgeos.h"
#include "dense_to_sparse.h"
#include "analytical_coverage.h"

#include "exactextract/grid.h"
#include "exactextract/box.h"
#include "exactextract/geos_utils.h"
#include "exactextract/coordinate.h"
#include "exactextract/side.h"
#include "exactextract/crossing.h"
#include "exactextract/traversal_areas.h"
#include "exactextract/measures.h"

// ---- GEOS context management ----

static void sl_geos_notice_handler(const char* /*fmt*/, ...) {}
static void sl_geos_error_handler(const char* fmt, ...) {
    char buf[1024];
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(buf, sizeof(buf), fmt, ap);
    va_end(ap);
    cpp11::stop("GEOS error: %s", buf);
}

class SLGEOSContextGuard {
public:
    SLGEOSContextGuard() {
        ctx_ = GEOS_init_r();
        if (!ctx_) cpp11::stop("Failed to create GEOS context");
        GEOSContext_setNoticeHandler_r(ctx_, sl_geos_notice_handler);
        GEOSContext_setErrorHandler_r(ctx_, sl_geos_error_handler);
    }
    ~SLGEOSContextGuard() { if (ctx_) GEOS_finish_r(ctx_); }
    GEOSContextHandle_t get() const { return ctx_; }
    SLGEOSContextGuard(const SLGEOSContextGuard&) = delete;
    SLGEOSContextGuard& operator=(const SLGEOSContextGuard&) = delete;
private:
    GEOSContextHandle_t ctx_;
};

class SLGEOSGeomGuard {
public:
    SLGEOSGeomGuard(GEOSContextHandle_t ctx, GEOSGeometry* g) : ctx_(ctx), g_(g) {}
    ~SLGEOSGeomGuard() { if (g_) GEOSGeom_destroy_r(ctx_, g_); }
    GEOSGeometry* get() const { return g_; }
    bool valid() const { return g_ != nullptr; }
    SLGEOSGeomGuard(const SLGEOSGeomGuard&) = delete;
    SLGEOSGeomGuard& operator=(const SLGEOSGeomGuard&) = delete;
private:
    GEOSContextHandle_t ctx_;
    GEOSGeometry* g_;
};

// ---- Lightweight traversal tracking ----

using exactextract::Coordinate;
using exactextract::Side;
using exactextract::Box;
using exactextract::Crossing;

struct LightTraversal {
    std::vector<Coordinate> coords;   // entry → intermediates → exit
    Side entry_side = Side::NONE;
    Side exit_side = Side::NONE;

    bool traversed() const {
        return entry_side != Side::NONE && exit_side != Side::NONE;
    }

    bool is_closed_ring() const {
        return coords.size() >= 3 &&
               coords.front() == coords.back();
    }

    bool multiple_unique_coordinates() const {
        for (size_t i = 1; i < coords.size(); i++) {
            if (coords[0] != coords[i]) return true;
        }
        return false;
    }
};

// Per-cell traversal data, keyed by (row, col) in the infinite_extent grid
struct CellRecord {
    Box box;
    std::vector<LightTraversal> traversals;

    CellRecord() : box(0, 0, 0, 0) {}
    explicit CellRecord(const Box& b) : box(b) {}
};

// Point classification relative to a box
enum class Location { INSIDE, BOUNDARY, OUTSIDE };

static inline Location point_location(const Box& box, const Coordinate& c) {
    if (box.strictly_contains(c)) return Location::INSIDE;
    if (box.contains(c)) return Location::BOUNDARY;
    return Location::OUTSIDE;
}

// ---- Per-cell boundary data for the winding sweep ----

struct BoundaryCellRecord {
    int col;               // 0-based column in full grid
    float coverage;        // accumulated coverage fraction (signed)
    int winding_delta;     // accumulated winding contribution
};

// ---- Scanline algorithm ----

static void walk_ring(
    std::vector<Coordinate> coords,
    bool is_ccw,
    bool is_exterior,
    const exactextract::Grid<exactextract::infinite_extent>& grid,
    double /* dy (unused, kept for signature compat) */,
    std::vector<std::vector<BoundaryCellRecord>>& row_data,
    size_t sub_rows,
    size_t sub_cols,
    size_t row_off,
    size_t col_off
) {
    if (coords.size() < 4) return;

    // Normalise to CCW for correct coverage fraction semantics
    if (!is_ccw) {
        std::reverse(coords.begin(), coords.end());
    }

    float coverage_factor = is_exterior ? 1.0f : -1.0f;
    int winding_factor = is_exterior ? 1 : -1;

    // Per-cell traversal data (keyed by grid row/col in infinite_extent grid)
    using CellKey = std::pair<size_t, size_t>;
    std::map<CellKey, CellRecord> cells;

    auto get_or_create = [&](size_t r, size_t c) -> CellRecord& {
        auto key = std::make_pair(r, c);
        auto it = cells.find(key);
        if (it == cells.end()) {
            Box box = grid_cell(grid, r, c);
            auto result = cells.emplace(key, CellRecord(box));
            return result.first->second;
        }
        return it->second;
    };

    // ---- Lightweight walk (replaces Cell-based walk) ----
    //
    // This replicates the core loop from raster_cell_intersection.cpp::process_line
    // using Box::crossing() directly instead of Cell::take().

    size_t pos = 0;
    size_t row = grid.get_row(coords.front().y);
    size_t col = grid.get_column(coords.front().x);

    // Storage for interpolated exit coordinate (persists across cell transitions)
    Coordinate exit_coord_buf(0, 0);
    const Coordinate* last_exit = nullptr;

    while (pos < coords.size()) {
        CellRecord& cr = get_or_create(row, col);
        const Box& box = cr.box;

        // Start a new traversal for this cell visit
        LightTraversal trav;

        while (pos < coords.size()) {
            const Coordinate* next = last_exit ? last_exit : &coords[pos];
            const Coordinate* prev_original = pos > 0 ? &coords[pos - 1] : nullptr;

            if (trav.coords.empty()) {
                // First coordinate for this traversal — enter the cell
                trav.entry_side = box.side(*next);
                trav.coords.push_back(*next);
                if (last_exit) { last_exit = nullptr; } else { pos++; }
                continue;
            }

            Location loc = point_location(box, *next);

            if (loc != Location::OUTSIDE) {
                // Inside or on boundary — add to traversal
                trav.coords.push_back(*next);
                if (last_exit) { last_exit = nullptr; } else { pos++; }
            } else {
                // Outside — compute exit crossing using Box::crossing()
                // Use prev_original for robustness (same as Cell::take)
                Crossing x = prev_original ?
                    box.crossing(*prev_original, *next) :
                    box.crossing(trav.coords.back(), *next);

                trav.coords.push_back(x.coord());
                trav.exit_side = x.side();

                // If exit coord differs from the target coord, use it as
                // the entry point for the next cell
                if (x.coord() != *next) {
                    exit_coord_buf = x.coord();
                    last_exit = &exit_coord_buf;
                }
                break;
            }
        }

        // Force exit if stuck on boundary (same as Cell::force_exit)
        if (trav.exit_side == Side::NONE && !trav.coords.empty()) {
            const Coordinate& last = trav.coords.back();
            if (point_location(box, last) == Location::BOUNDARY) {
                trav.exit_side = box.side(last);
            }
        }

        bool exited = (trav.exit_side != Side::NONE);
        bool incomplete = exited && (trav.entry_side == Side::NONE);

        // Handle incomplete initial traversal: polygon started inside this cell,
        // ring hasn't closed yet. Push coords to end for reprocessing.
        if (incomplete) {
            for (const auto& c : trav.coords) {
                coords.push_back(c);
            }
        }

        // Store the traversal (even incomplete ones, they'll be filtered later)
        cr.traversals.push_back(std::move(trav));

        // Move to next cell based on exit side
        if (exited) {
            Side exit_s = cr.traversals.back().exit_side;
            switch (exit_s) {
                case Side::TOP:    row--; break;
                case Side::BOTTOM: row++; break;
                case Side::LEFT:   col--; break;
                case Side::RIGHT:  col++; break;
                default: break;
            }
        }
    }

    // ---- Compute coverage fractions and winding ----

    for (auto& kv : cells) {
        size_t r = kv.first.first;
        size_t c = kv.first.second;
        CellRecord& cr = kv.second;

        // Map from infinite_extent grid coords to subgrid coords.
        // Skip padding ROWS — they don't affect any grid row's winding.
        if (r < 1) continue;
        size_t sub_r = r - 1;
        if (sub_r >= sub_rows) continue;

        // Determine column mapping.
        // Padding COLUMNS still carry winding deltas for their grid row
        // (e.g. a polygon edge entirely outside the grid but crossing rows).
        bool in_grid_cols;
        int full_col;
        if (c < 1) {
            // Left padding column — virtual column before grid
            full_col = static_cast<int>(col_off) - 1;
            in_grid_cols = false;
        } else {
            size_t sub_c = c - 1;
            if (sub_c >= sub_cols) {
                // Right padding column — virtual column after grid
                full_col = static_cast<int>(col_off + sub_cols);
                in_grid_cols = false;
            } else {
                full_col = static_cast<int>(col_off + sub_c);
                in_grid_cols = true;
            }
        }

        // Filter to valid traversals (proper enter + exit, or closed ring)
        std::vector<LightTraversal*> valid;
        for (auto& t : cr.traversals) {
            if (t.traversed() && t.multiple_unique_coordinates()) {
                valid.push_back(&t);
            } else if (t.entry_side == Side::NONE && t.is_closed_ring()) {
                // Closed ring entirely within this cell
                valid.push_back(&t);
            }
        }

        if (valid.empty()) continue;

        // ---- Coverage fraction (only for in-grid cells) ----
        float frac = 0.0f;

        if (in_grid_cols) {
            if (valid.size() == 1 && valid[0]->entry_side == Side::NONE
                && valid[0]->is_closed_ring()) {
                frac = static_cast<float>(
                    controlledburn::closed_ring_covered_fraction(cr.box, valid[0]->coords));
            } else if (valid.size() == 1) {
                frac = static_cast<float>(
                    controlledburn::analytical_covered_fraction(
                        cr.box, valid[0]->coords,
                        valid[0]->entry_side, valid[0]->exit_side));
            } else {
                std::vector<const std::vector<Coordinate>*> coord_lists;
                for (auto* t : valid) {
                    coord_lists.push_back(&t->coords);
                }
                double cell_area = cr.box.area();
                if (cell_area > 0.0) {
                    frac = static_cast<float>(
                        exactextract::left_hand_area(cr.box, coord_lists) / cell_area);
                }
            }
        }

        // ---- Store coverage (if nonzero) and winding deltas ----
        //
        // IMPORTANT: winding deltas must be stored even when coverage is zero.
        // A traversal along a cell wall (e.g. vertical edge on a grid line)
        // has zero area but still crosses the row center, contributing to the
        // winding count that classifies interior cells.

        auto& row_vec = row_data[sub_r];

        // Helper: find or create the BoundaryCellRecord for this column
        auto find_or_create = [&]() -> BoundaryCellRecord& {
            for (auto& rec : row_vec) {
                if (rec.col == full_col) return rec;
            }
            row_vec.push_back({full_col, 0.0f, 0});
            return row_vec.back();
        };

        if (frac != 0.0f) {
            BoundaryCellRecord& rec = find_or_create();
            rec.coverage += coverage_factor * frac;
        }

        // Winding deltas from traversals
        for (auto* t : valid) {
            if (!t->traversed()) continue; // closed rings don't contribute winding
            if (t->coords.size() < 2) continue;

            double entry_y = t->coords.front().y;
            double exit_y  = t->coords.back().y;
            double y_mid = (cr.box.ymin + cr.box.ymax) / 2.0;

            bool crosses = (entry_y > y_mid && exit_y < y_mid) ||
                           (entry_y < y_mid && exit_y > y_mid);
            if (!crosses) continue;

            int delta = (entry_y > y_mid) ? -1 : +1; // downward = -1, upward = +1
            delta *= winding_factor;

            BoundaryCellRecord& rec = find_or_create();
            rec.winding_delta += delta;
        }
    }
}


// ---- Process geometry: per-POLYGON processing with padding-aware sweep ----
//
// For MULTIPOLYGON/GEOMETRYCOLLECTION, each POLYGON component is processed
// independently with its own subgrid, row_data, and winding sweep. This
// prevents winding from one disjoint component bleeding into another's
// boundary cells (which would incorrectly promote partial coverage to 1.0).

static void process_geometry(
    GEOSContextHandle_t ctx,
    const GEOSGeometry* g,
    const exactextract::Grid<exactextract::bounded_extent>& full_grid,
    double dx, double dy,
    int poly_id,
    std::vector<GridRun>& all_runs,
    std::vector<GridEdge>& all_edges
) {
    using namespace exactextract;

    int type = GEOSGeomTypeId_r(ctx, g);

    if (type == GEOS_GEOMETRYCOLLECTION || type == GEOS_MULTIPOLYGON) {
        for (int i = 0; i < GEOSGetNumGeometries_r(ctx, g); i++) {
            process_geometry(ctx, GEOSGetGeometryN_r(ctx, g, i),
                             full_grid, dx, dy, poly_id, all_runs, all_edges);
        }
        return;
    }

    if (type != GEOS_POLYGON) return;

    auto component_boxes = geos_get_component_boxes(ctx, g);
    Box region = Box::make_empty();
    for (const auto& box : component_boxes) {
        if (!box.intersects(full_grid.extent())) continue;
        Box isect = full_grid.extent().intersection(box);
        if (region.empty()) {
            region = isect;
        } else if (!region.contains(isect)) {
            region = region.expand_to_include(isect);
        }
    }
    if (region.empty()) return;

    auto subgrid_bounded = full_grid.shrink_to_fit(region);
    auto subgrid = make_infinite(subgrid_bounded);
    if (subgrid.empty()) return;

    size_t sub_rows = subgrid.rows() - 2;
    size_t sub_cols = subgrid.cols() - 2;

    size_t row_off = static_cast<size_t>(
        std::round((full_grid.ymax() - subgrid_bounded.ymax()) / dy));
    size_t col_off = static_cast<size_t>(
        std::round((subgrid_bounded.xmin() - full_grid.xmin()) / dx));

    std::vector<std::vector<BoundaryCellRecord>> row_data(sub_rows);

    // Exterior ring
    {
        const GEOSGeometry* ring = GEOSGetExteriorRing_r(ctx, g);
        const GEOSCoordSequence* seq = GEOSGeom_getCoordSeq_r(ctx, ring);
        auto coords = read(ctx, seq);
        bool is_ccw = geos_is_ccw(ctx, seq);

        walk_ring(coords, is_ccw, true, subgrid, dy,
                  row_data, sub_rows, sub_cols, row_off, col_off);
    }

    // Holes
    int n_holes = GEOSGetNumInteriorRings_r(ctx, g);
    for (int h = 0; h < n_holes; h++) {
        const GEOSGeometry* ring = GEOSGetInteriorRingN_r(ctx, g, h);
        const GEOSCoordSequence* seq = GEOSGeom_getCoordSeq_r(ctx, ring);
        auto coords = read(ctx, seq);
        bool is_ccw = geos_is_ccw(ctx, seq);

        walk_ring(coords, is_ccw, false, subgrid, dy,
                  row_data, sub_rows, sub_cols, row_off, col_off);
    }

    // ---- Winding sweep: emit runs and edges per row ----

    float tol = 1e-6f;

    for (size_t sr = 0; sr < sub_rows; sr++) {
        auto& row_vec = row_data[sr];
        if (row_vec.empty()) continue;

        std::sort(row_vec.begin(), row_vec.end(),
            [](const BoundaryCellRecord& a, const BoundaryCellRecord& b) {
                return a.col < b.col;
            });

        // Merge same-column entries
        std::vector<BoundaryCellRecord> merged;
        for (auto& rec : row_vec) {
            if (!merged.empty() && merged.back().col == rec.col) {
                merged.back().coverage += rec.coverage;
                merged.back().winding_delta += rec.winding_delta;
            } else {
                merged.push_back(rec);
            }
        }

        int winding = 0;
        int prev_col = -2;
        int full_row = static_cast<int>(row_off + sr) + 1;

        for (auto& mc : merged) {
            // Emit interior run between previous boundary cell and this one.
            // prev_col > -2 means at least one cell has been seen (including padding).
            if (winding != 0 && prev_col > -2 && mc.col > prev_col + 1) {
                all_runs.push_back({
                    full_row,
                    prev_col + 1 + 1,
                    mc.col - 1 + 1,
                    poly_id
                });
            }

            float w = mc.coverage;
            if (w > tol && w < (1.0f - tol)) {
                all_edges.push_back({full_row, mc.col + 1, w, poly_id});
            } else if (w >= (1.0f - tol)) {
                all_runs.push_back({full_row, mc.col + 1, mc.col + 1, poly_id});
            }

            winding += mc.winding_delta;
            prev_col = mc.col;
        }
    }
}


// ---- cpp11 entry point ----

[[cpp11::register]]
cpp11::writable::list cpp_scanline_burn(
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

    exactextract::Grid<exactextract::bounded_extent> full_grid(
        exactextract::Box(xmin, ymin, xmax, ymax), dx, dy);

    SLGEOSContextGuard geos_guard;
    GEOSContextHandle_t ctx = geos_guard.get();

    GEOSWKBReader* wkb_reader = GEOSWKBReader_create_r(ctx);
    if (!wkb_reader) cpp11::stop("Failed to create WKB reader");

    int n_geoms = wkb_list.size();

    std::vector<GridRun> all_runs;
    std::vector<GridEdge> all_edges;

    for (int k = 0; k < n_geoms; k++) {
        cpp11::raws wkb_raw(wkb_list[k]);
        int wkb_size = wkb_raw.size();
        if (wkb_size == 0) continue;

        const unsigned char* wkb_data = reinterpret_cast<const unsigned char*>(RAW(wkb_raw));

        SLGEOSGeomGuard geom(ctx,
            GEOSWKBReader_read_r(ctx, wkb_reader, wkb_data, static_cast<size_t>(wkb_size)));

        if (!geom.valid()) {
            cpp11::warning("Failed to parse WKB for geometry %d, skipping", k + 1);
            continue;
        }

        if (GEOSisEmpty_r(ctx, geom.get())) continue;

        try {
            int poly_id = k + 1;
            process_geometry(ctx, geom.get(), full_grid, dx, dy,
                             poly_id, all_runs, all_edges);
        } catch (const std::exception& e) {
            cpp11::warning("Error processing geometry %d: %s, skipping", k + 1, e.what());
            continue;
        }
    }

    GEOSWKBReader_destroy_r(ctx, wkb_reader);

    // Build R data.frames
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

    cpp11::writable::list result(2);
    result[0] = static_cast<SEXP>(runs_df);
    result[1] = static_cast<SEXP>(edges_df);
    result.attr("names") = cpp11::writable::strings({"runs", "edges"});

    return result;
}
