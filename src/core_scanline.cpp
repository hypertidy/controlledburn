// scanline.cpp -- scanline polygon/line/point rasterization with exact
// coverage fractions. Pure C++17: no GEOS, no R, no materialized raster.
//
// The walker (walk_polyline), ring processing (walk_ring), and winding
// sweep are carried over verbatim from the R package's scanline_burn.cpp;
// only the geometry access layer changed (plain coordinate vectors in
// place of GEOS ring accessors, shoelace orientation in place of
// geos_is_ccw, coordinate min/max in place of GEOS envelopes).
//
// Edge cases handled (see the R package history for provenance):
// - Padding column winding: edges outside the grid (padding columns)
//   contribute winding deltas to grid rows, fixing beyond-extent polygons.
// - Sweep start condition: prev_col > -2 allows runs after padding cells.
// - Analytical single-traversal coverage via perimeter_distance.
// - MULTIPOLYGON components processed independently (each gets own sweep).
//
// Copyright (c) 2025 Michael Sumner
// Licensed under Apache License 2.0

#include "controlledburn/burn.hpp"
#include "controlledburn/geometry.hpp"
#include "controlledburn/output.hpp"
#include "controlledburn/wkb.hpp"

#include <map>
#include <algorithm>
#include <cmath>
#include <utility>

#include "exactextract/grid.h"
#include "exactextract/box.h"
#include "exactextract/coordinate.h"
#include "exactextract/side.h"
#include "exactextract/crossing.h"
#include "exactextract/traversal_areas.h"
#include "exactextract/measures.h"
#include "analytical_coverage.h"

namespace controlledburn {

double signed_area(const CoordSeq& ring) {
    // Standard shoelace: positive = counter-clockwise. Note this is the
    // OPPOSITE sign convention to polygon_signed_area in
    // ee/analytical_coverage.h, whose sign is internal to the coverage
    // math (it mirrors exactextract::area).
    if (ring.size() < 3) return 0.0;
    double sum = 0.0;
    double x0 = ring[0].x;
    for (size_t i = 1; i + 1 < ring.size(); i++) {
        double x = ring[i].x - x0;
        double y_next = ring[i + 1].y;
        double y_prev = ring[i - 1].y;
        sum += x * (y_next - y_prev);
    }
    return sum / 2.0;
}

namespace {

// ---- Lightweight traversal tracking ----

using exactextract::Coordinate;
using exactextract::Side;
using exactextract::Box;
using exactextract::Crossing;

struct LightTraversal {
  std::vector<Coordinate> coords;   // entry -> intermediates -> exit
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

  // Approx mode only: winding delta from traversals whose y-midpoint
  // crossing is to the LEFT of the cell center x. Used to compute the
  // winding number at the cell center for center-rule classification.
  int winding_delta_left;
};

// ---- Scanline algorithm ----

// Type alias for the walker's per-visit cell map.
// Keyed by (row, col) in the infinite_extent grid.
using CellKey = std::pair<size_t, size_t>;
using CellMap = std::map<CellKey, CellRecord>;

// Walk a polyline through grid cells, recording per-cell traversals.
//
// Geometry-agnostic core extracted from walk_ring. Steps along the polyline
// using Box::crossing() to compute cell-edge intersections, accumulating one
// or more LightTraversal records per cell visited. Returns a map keyed by
// (row, col) in the supplied infinite_extent grid.
//
// `coords` is taken by value because the walker may push to it (see the
// `incomplete` reprocess branch -- the case where the polyline begins
// strictly inside a cell rather than on a boundary).
//
// `closed` selects how to handle that incomplete-initial-traversal case:
//   - true  (polygon ring): push coords to the end so the start cell gets
//           re-entered with a proper entry side on the second pass. This
//           stitches the start back to the closing edge of the ring.
//   - false (open polyline): leave the partial traversal in place. The
//           consumer accepts traversals with entry_side == NONE as valid
//           (the line genuinely starts inside the cell with no prior edge).
//
// Replicates the core loop from raster_cell_intersection.cpp::process_line.
static CellMap walk_polyline(
    std::vector<Coordinate> coords,
    const exactextract::Grid<exactextract::infinite_extent>& grid,
    bool closed
) {
  CellMap cells;
  if (coords.empty()) return cells;

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
        // First coordinate for this traversal -- enter the cell
        trav.entry_side = box.side(*next);
        trav.coords.push_back(*next);
        if (last_exit) { last_exit = nullptr; } else { pos++; }
        continue;
      }

      Location loc = point_location(box, *next);

      if (loc != Location::OUTSIDE) {
        // Inside or on boundary -- add to traversal
        trav.coords.push_back(*next);
        if (last_exit) { last_exit = nullptr; } else { pos++; }
      } else {
        // Outside -- compute exit crossing using Box::crossing()
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

    // Handle incomplete initial traversal: walk started inside this cell,
    // boundary hasn't closed yet. For closed rings, push coords to end for
    // reprocessing -- this stitches the start back to a proper entry on
    // the second pass. For open polylines, leave the partial traversal as
    // is; the consumer accepts entry_side==NONE as valid for line input.
    if (incomplete && closed) {
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

  return cells;
}

// Walk a polygon ring and emit boundary-cell records (coverage + winding).
//
// Thin wrapper around walk_polyline + polygon-specific post-processing:
//   1. Reject malformed ring (< 4 coords for a closed ring)
//   2. Normalise to CCW for correct coverage fraction semantics
//   3. Walk the polyline to populate the cells map
//   4. Compute coverage fractions and winding deltas, write into row_data
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
    size_t col_off,
    BurnMode mode = BurnMode::Coverage
) {
  if (coords.size() < 4) return;

  // Normalise to CCW for correct coverage fraction semantics
  if (!is_ccw) {
    std::reverse(coords.begin(), coords.end());
  }

  float coverage_factor = is_exterior ? 1.0f : -1.0f;
  int winding_factor = is_exterior ? 1 : -1;

  CellMap cells = walk_polyline(std::move(coords), grid, /*closed=*/true);

  // ---- Compute coverage fractions and winding ----

  for (auto& kv : cells) {
    size_t r = kv.first.first;
    size_t c = kv.first.second;
    CellRecord& cr = kv.second;

    // Map from infinite_extent grid coords to subgrid coords.
    // Skip padding ROWS -- they don't affect any grid row's winding.
    if (r < 1) continue;
    size_t sub_r = r - 1;
    if (sub_r >= sub_rows) continue;

    // Determine column mapping.
    // Padding COLUMNS still carry winding deltas for their grid row
    // (e.g. a polygon edge entirely outside the grid but crossing rows).
    bool in_grid_cols;
    int full_col;
    if (c < 1) {
      // Left padding column -- virtual column before grid
      full_col = static_cast<int>(col_off) - 1;
      in_grid_cols = false;
    } else {
      size_t sub_c = c - 1;
      if (sub_c >= sub_cols) {
        // Right padding column -- virtual column after grid
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

    // ---- Coverage fraction (only for in-grid cells, Coverage mode) ----
    float frac = 0.0f;

    if (in_grid_cols && mode == BurnMode::Coverage) {
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
      row_vec.push_back({full_col, 0.0f, 0, 0});
      return row_vec.back();
    };

    if (frac != 0.0f) {
      BoundaryCellRecord& rec = find_or_create();
      rec.coverage += coverage_factor * frac;
    }

    // Winding deltas from traversals
    double x_mid = (cr.box.xmin + cr.box.xmax) / 2.0;

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

      // Approx mode: classify whether this crossing is left of cell center.
      // Interpolate x at y_mid along the traversal path.
      if (mode == BurnMode::Approx) {
        double x_at_mid = x_mid; // fallback
        const auto& tc = t->coords;
        for (size_t i = 1; i < tc.size(); i++) {
          double y0 = tc[i-1].y, y1 = tc[i].y;
          if ((y0 <= y_mid && y1 >= y_mid) || (y0 >= y_mid && y1 <= y_mid)) {
            if (y1 != y0) {
              double frac_along = (y_mid - y0) / (y1 - y0);
              x_at_mid = tc[i-1].x + frac_along * (tc[i].x - tc[i-1].x);
            } else {
              x_at_mid = (tc[i-1].x + tc[i].x) / 2.0;
            }
            break;
          }
        }
        // Use <= for left-inclusive convention: an edge crossing exactly
        // at the cell center x is considered "left of center", matching
        // the standard rasterization rule (left/bottom inclusive,
        // right/top exclusive).
        if (x_at_mid <= x_mid) {
          rec.winding_delta_left += delta;
        }
      }
    }
  }
}


// ---- Coordinate conversion at the geometry boundary ----
//
// The walker owns and may mutate its coordinate vector (the incomplete
// initial traversal reprocess branch pushes to it), so a copy at the
// boundary is inherent to the design, not overhead added by the
// conversion.

static std::vector<Coordinate> to_coordinates(const CoordSeq& seq) {
  std::vector<Coordinate> out;
  out.reserve(seq.size());
  for (const auto& c : seq) out.emplace_back(c.x, c.y);
  return out;
}

static Box ring_envelope(const CoordSeq& ring) {
  Box b = Box::make_empty();
  if (ring.empty()) return b;
  double xmin = ring[0].x, xmax = ring[0].x;
  double ymin = ring[0].y, ymax = ring[0].y;
  for (const auto& c : ring) {
    if (c.x < xmin) xmin = c.x;
    if (c.x > xmax) xmax = c.x;
    if (c.y < ymin) ymin = c.y;
    if (c.y > ymax) ymax = c.y;
  }
  return Box(xmin, ymin, xmax, ymax);
}

// ---- Process polygon: per-POLYGON processing with padding-aware sweep ----
//
// Single-polygon path. For MULTIPOLYGON each component is processed
// independently with its own subgrid, row_data, and winding sweep. This
// prevents winding from one disjoint component bleeding into another's
// boundary cells (which would incorrectly promote partial coverage to 1.0).

static void process_polygon(
    const Polygon& poly,
    const exactextract::Grid<exactextract::bounded_extent>& full_grid,
    double dx, double dy,
    int poly_id,
    std::vector<GridRun>& all_runs,
    std::vector<GridEdge>& all_edges,
    BurnMode mode = BurnMode::Coverage
) {
  using namespace exactextract;

  if (poly.rings.empty() || poly.rings[0].size() < 4) return;

  // Region of interest: exterior-ring envelope clipped to the grid.
  // Holes lie inside the exterior, so envelope(exterior) ==
  // envelope(polygon); this replaces geos_get_component_boxes.
  Box env = ring_envelope(poly.rings[0]);
  if (!env.intersects(full_grid.extent())) return;
  Box region = full_grid.extent().intersection(env);
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

  // Exterior ring, then holes. Orientation from shoelace sign (replaces
  // geos_is_ccw); walk_ring normalises to CCW internally.
  for (size_t r = 0; r < poly.rings.size(); r++) {
    const CoordSeq& ring = poly.rings[r];
    if (ring.size() < 4) continue;
    bool is_exterior = (r == 0);
    walk_ring(to_coordinates(ring), is_ccw(ring), is_exterior, subgrid, dy,
              row_data, sub_rows, sub_cols, row_off, col_off, mode);
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
        merged.back().winding_delta_left += rec.winding_delta_left;
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

      if (mode == BurnMode::Coverage) {
        float w = mc.coverage;
        if (w > tol && w < (1.0f - tol)) {
          all_edges.push_back({full_row, mc.col + 1, w, poly_id});
        } else if (w >= (1.0f - tol)) {
          all_runs.push_back({full_row, mc.col + 1, mc.col + 1, poly_id});
        }
      } else {
        // Approx: cell center is inside iff winding at center != 0.
        // Winding at center = winding (from left) + delta_left (crossings
        // left of center x within this cell).
        int winding_at_center = winding + mc.winding_delta_left;
        if (winding_at_center != 0) {
          all_runs.push_back({full_row, mc.col + 1, mc.col + 1, poly_id});
        }
      }

      winding += mc.winding_delta;
      prev_col = mc.col;
    }
  }
}

// ---- Process line: per-LINESTRING length-in-cell rasterization ----
//
// Walks the line through the grid using walk_polyline (closed=false) and
// emits one (row, col, length, id) record per cell touched. Length is the
// absolute length of the line within the cell, in CRS units -- sum of
// segment lengths across all traversal records for that cell.

static void process_line(
    const CoordSeq& line,
    const exactextract::Grid<exactextract::bounded_extent>& full_grid,
    int line_id,
    std::vector<GridLine>& all_lines
) {
  using namespace exactextract;

  if (line.size() < 2) return;

  // Walk the line directly on the full grid (no subgrid shrinking):
  // horizontal or vertical lines have degenerate bounding boxes that
  // Box::empty() reports as empty, and the walker's cost is O(cells
  // touched) regardless of grid dimensions.
  auto grid = make_infinite(full_grid);
  if (grid.empty()) return;

  size_t n_rows = grid.rows();
  size_t n_cols = grid.cols();

  CellMap cells = walk_polyline(to_coordinates(line), grid, /*closed=*/false);

  for (const auto& kv : cells) {
    size_t r = kv.first.first;
    size_t c = kv.first.second;
    const CellRecord& cr = kv.second;

    // Skip padding ring (top/bottom row, left/right col of infinite_extent).
    if (r < 1 || r >= n_rows - 1) continue;
    if (c < 1 || c >= n_cols - 1) continue;

    double total_length = 0.0;
    for (const auto& t : cr.traversals) {
      const auto& tc = t.coords;
      if (tc.size() < 2) continue;
      for (size_t i = 1; i < tc.size(); i++) {
        double sx = tc[i].x - tc[i-1].x;
        double sy = tc[i].y - tc[i-1].y;
        total_length += std::sqrt(sx * sx + sy * sy);
      }
    }

    if (total_length > 0.0) {
      // Map infinite_extent (r, c) to 1-based bounded grid; padding=1
      // absorbs the +1.
      all_lines.push_back({
        static_cast<int>(r), static_cast<int>(c),
        static_cast<float>(total_length),
        line_id
      });
    }
  }
}

// ---- Process point: compute cell index, emit one record ----
//
// No fractional intersection for the 0-D case: a point is either in a
// cell or it isn't, so points carry no weight column. Points outside
// the grid extent are dropped silently, consistent with lines/polygons.

static void process_point(
    const CoordSeq& pt,
    const exactextract::Grid<exactextract::bounded_extent>& full_grid,
    int point_id,
    std::vector<GridPoint>& all_points
) {
  if (pt.empty()) return;
  double x = pt[0].x;
  double y = pt[0].y;
  if (!std::isfinite(x) || !std::isfinite(y)) return; // POINT EMPTY

  // Boundary inclusion: a point exactly on xmax or ymin is assigned to
  // the boundary row/col by Grid::get_row / get_column's special cases.
  const auto& ext = full_grid.extent();
  if (x < ext.xmin || x > ext.xmax || y < ext.ymin || y > ext.ymax) return;

  size_t r = full_grid.get_row(y);
  size_t c = full_grid.get_column(x);

  all_points.push_back({
    static_cast<int>(r) + 1, // 1-based
    static_cast<int>(c) + 1,
    point_id
  });
}

// ---- Process geometry: dispatch by kind ----

static void process_geometry(
    const Geometry& g,
    const exactextract::Grid<exactextract::bounded_extent>& full_grid,
    double dx, double dy,
    int geom_id,
    BurnResult& out,
    BurnMode mode = BurnMode::Coverage
) {
  // Multi types iterate the same containers as their singular
  // counterparts; each component is processed independently (polygons
  // get their own winding sweep).
  for (const auto& poly : g.polygons) {
    process_polygon(poly, full_grid, dx, dy, geom_id, out.runs, out.edges, mode);
  }
  for (const auto& line : g.lines) {
    process_line(line, full_grid, geom_id, out.lines);
  }
  for (const auto& pt : g.points) {
    process_point(pt, full_grid, geom_id, out.points);
  }
}

} // anonymous namespace

// ---- Public entry points ----

BurnResult burn(const std::vector<Geometry>& geoms, const GridSpec& gs,
                BurnMode mode) {
  gs.validate();

  double dx = gs.dx();
  double dy = gs.dy();

  exactextract::Grid<exactextract::bounded_extent> full_grid(
      exactextract::Box(gs.xmin, gs.ymin, gs.xmax, gs.ymax), dx, dy);

  BurnResult out;
  for (size_t k = 0; k < geoms.size(); k++) {
    const Geometry& g = geoms[k];
    if (g.empty()) continue;
    try {
      process_geometry(g, full_grid, dx, dy, static_cast<int>(k) + 1, out, mode);
    } catch (const std::exception& e) {
      out.notes.push_back({static_cast<int32_t>(k) + 1,
                           std::string("error processing geometry: ") + e.what()});
    }
  }
  return out;
}

BurnResult burn_wkb(const std::vector<WKBSpan>& wkb, const GridSpec& gs,
                    BurnMode mode) {
  gs.validate();

  double dx = gs.dx();
  double dy = gs.dy();

  exactextract::Grid<exactextract::bounded_extent> full_grid(
      exactextract::Box(gs.xmin, gs.ymin, gs.xmax, gs.ymax), dx, dy);

  BurnResult out;
  for (size_t k = 0; k < wkb.size(); k++) {
    if (!wkb[k].data || wkb[k].size == 0) continue;

    Geometry g;
    try {
      g = parse_wkb(wkb[k].data, wkb[k].size);
    } catch (const WKBParseError& e) {
      out.notes.push_back({static_cast<int32_t>(k) + 1,
                           std::string("failed to parse WKB: ") + e.what()});
      continue;
    }
    if (g.empty()) continue;

    try {
      process_geometry(g, full_grid, dx, dy, static_cast<int>(k) + 1, out, mode);
    } catch (const std::exception& e) {
      out.notes.push_back({static_cast<int32_t>(k) + 1,
                           std::string("error processing geometry: ") + e.what()});
    }
  }
  return out;
}

} // namespace controlledburn
