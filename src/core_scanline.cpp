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

#include "no_fma.h"  // first: no fused multiply-add (see header)
#include "controlledburn/burn.hpp"
#include "controlledburn/geometry.hpp"
#include "controlledburn/output.hpp"
#include "controlledburn/wkb.hpp"

#include <atomic>
#include <climits>
#include <exception>
#include <limits>
#include <map>
#include <string>
#include <thread>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <utility>

#include "exactextract/grid.h"
#include "exactextract/box.h"
#include "exactextract/coordinate.h"
#include "exactextract/side.h"
#include "exactextract/crossing.h"
#include "exactextract/traversal_areas.h"
#include "exactextract/measures.h"
#include "analytical_coverage.h"
#include "predicates.h"

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

using detail::orient2d_sign;

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

  // The ring (0 the exterior) the record comes from, so the sweep can
  // tell which rings pass through a cell and which enclose it.
  int ring;
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
// Snap vertices that lie within a tiny tolerance of a cell edge exactly
// onto that edge, using the same arithmetic grid_cell() uses to place the
// edges. A vertex parsed from input (0.3) and an edge computed from the
// grid (3 * 0.1 = 0.30000000000000004) can differ by a few ulp; left
// unsnapped, such a vertex is classified OUTSIDE the cell it visually sits
// on the corner of, Box::crossing() extrapolates past the segment end, and
// the walk leaves the grid and never returns (or, if it does, exits a cell
// through the wrong side at a corner and left_hand_area() over-counts).
// Tolerance is relative to the cell size; only degenerate near-boundary
// vertices move, by less than 1e-9 of a cell.
static void snap_to_lattice(
    std::vector<Coordinate>& coords,
    const exactextract::Grid<exactextract::infinite_extent>& grid
) {
  const double dx = grid.dx(), dy = grid.dy();
  if (!(dx > 0.0) || !(dy > 0.0)) return;
  const double x0 = grid.xmin(), x1 = grid.xmax();
  const double y0 = grid.ymin(), y1 = grid.ymax();
  // infinite_extent grids carry one padding row/col on each side
  const double ncol = static_cast<double>(grid.cols()) - 2.0;
  const double nrow = static_cast<double>(grid.rows()) - 2.0;
  const double tol_x = dx * 1e-9, tol_y = dy * 1e-9;
  for (auto& c : coords) {
    double kx = std::round((c.x - x0) / dx);
    if (kx >= 0.0 && kx <= ncol) {
      double bx = (kx == ncol) ? x1 : x0 + kx * dx;
      if (std::fabs(c.x - bx) <= tol_x) c.x = bx;
    }
    double ky = std::round((y1 - c.y) / dy);
    if (ky >= 0.0 && ky <= nrow) {
      double by = (ky == nrow) ? y0 : y1 - ky * dy;
      if (std::fabs(c.y - by) <= tol_y) c.y = by;
    }
  }
}

// Move c (on the boundary of box) onto the corner of box it lies within
// tol_x, tol_y of, if any. Corners of padding cells are at -/+DBL_MAX and
// are never within tolerance.
static inline Coordinate snap_to_corner(const Coordinate& c, const Box& box,
                                        double tol_x, double tol_y) {
  double x, y;
  if (std::fabs(c.x - box.xmin) <= tol_x)      x = box.xmin;
  else if (std::fabs(c.x - box.xmax) <= tol_x) x = box.xmax;
  else return c;
  if (std::fabs(c.y - box.ymin) <= tol_y)      y = box.ymin;
  else if (std::fabs(c.y - box.ymax) <= tol_y) y = box.ymax;
  else return c;
  return Coordinate(x, y);
}

static CellMap walk_polyline(
    std::vector<Coordinate> coords,
    const exactextract::Grid<exactextract::infinite_extent>& grid,
    bool closed
) {
  CellMap cells;
  if (coords.empty()) return cells;
  snap_to_lattice(coords, grid);
  const size_t n_input = coords.size();
  // Same tolerance as snap_to_lattice, from the grid's cell size (a padding
  // cell's box is unbounded).
  const double tol_x = grid.dx() * 1e-9, tol_y = grid.dy() * 1e-9;

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
        const Coordinate& from = prev_original ? *prev_original : trav.coords.back();
        Crossing x = box.crossing(from, *next);
        Coordinate exit_c = x.coord();
        Side exit_s = x.side();

        // Rounding guard. A vertex can sit a few ulp outside this cell when
        // the cell edge (computed from the grid origin and resolution) and
        // the vertex (parsed from input) round differently, e.g. a vertex at
        // 0.3 against an edge at 3 * 0.1. point_location() then says OUTSIDE
        // and Box::crossing() extrapolates the line past the end of the
        // segment to find an "exit", after which the walk follows the
        // extended line off the grid and never returns. If the crossing is
        // beyond the segment end, clamp the target onto the box boundary and
        // exit there instead.
        {
          double dx_n = next->x - from.x, dy_n = next->y - from.y;
          double dx_x = exit_c.x - from.x, dy_x = exit_c.y - from.y;
          double d_next = dx_n * dx_n + dy_n * dy_n;
          double d_exit = dx_x * dx_x + dy_x * dy_x;
          if (!std::isfinite(exit_c.x) || !std::isfinite(exit_c.y) ||
              d_exit > d_next * (1.0 + 1e-9)) {
            exit_c = Coordinate(std::min(std::max(next->x, box.xmin), box.xmax),
                                std::min(std::max(next->y, box.ymin), box.ymax));
            // prefer the side the target actually overshoots
            if (next->x < box.xmin)      exit_s = Side::LEFT;
            else if (next->x > box.xmax) exit_s = Side::RIGHT;
            else if (next->y < box.ymin) exit_s = Side::BOTTOM;
            else                         exit_s = Side::TOP;
          }
        }

        // Corner guard: an edge through a cell corner gives a crossing an
        // ulp or two off the corner, on the wall either side of it. The next
        // cell then starts a sliver traversal from there to the corner, or
        // from the corner to a point an ulp down the wall, where the true
        // path only touches the corner. Its two ends can be out of order
        // around the cell boundary, which the old chain chaser read as a
        // slab of the cell or all of it (#25). Exit through the corner
        // itself; the sliver becomes a single-point touch, which the
        // coverage pass already ignores.
        exit_c = snap_to_corner(exit_c, box, tol_x, tol_y);

        trav.coords.push_back(exit_c);
        trav.exit_side = exit_s;

        // If exit coord differs from the target coord, use it as
        // the entry point for the next cell
        if (exit_c != *next) {
          exit_coord_buf = exit_c;
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
      // A closed ring is reprocessed at most once per starting cell; if the
      // coordinate list keeps growing the walk is not converging.
      if (coords.size() > 4 * n_input + 64) {
        throw std::runtime_error("ring walk did not converge (vertex on cell boundary?)");
      }
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

// A traversal that counts: a proper entry and exit with movement, or a
// closed ring entirely within the cell.
static inline bool valid_traversal(const LightTraversal& t) {
  return (t.traversed() && t.multiple_unique_coordinates()) ||
         (t.entry_side == Side::NONE && t.is_closed_ring());
}

// Whether every segment of tc runs along a wall of box. Such a traversal
// does not split the cell: whether the cell is covered depends on the rest
// of the ring (green_fraction reads it from the row). The single-traversal
// formula would take the side left of travel, all of the cell or none of
// it. A ring doubling back along a wall puts the two directions in
// different cells as the walk happens to run, the grid cell and the padding
// row below it say, and the grid cell came out fully covered (#25).
static bool along_walls(const Box& box, const std::vector<Coordinate>& tc) {
  for (size_t i = 1; i < tc.size(); i++) {
    const Coordinate& a = tc[i - 1];
    const Coordinate& b = tc[i];
    bool on_x = a.x == b.x && (a.x == box.xmin || a.x == box.xmax);
    bool on_y = a.y == b.y && (a.y == box.ymin || a.y == box.ymax);
    if (!on_x && !on_y) return false;
  }
  return true;
}

// Coverage fraction of a cell crossed by several traversals of one CCW
// ring, by Green's theorem: the covered area is the sum over the traversals
// of the integral of (x - xmax) dy, plus the width times the covered length
// of the left wall, which is minus left_dy, the net height of the ring's
// traversals left of the cell in the same row.
//
// Every traversal adds its own term, so the result does not depend on the
// order of the traversal ends around the cell boundary. The chain chaser
// this replaces (exactextract::left_hand_area) did: two nearly coincident
// traversals (a sliver, or an edge doubling back along itself) can cross a
// wall an ulp apart in the wrong order, and the chaser then returned the
// complement of the sliver, a whole cell too much (#25). Same as
// controlledburn-rs's green_fraction, and the identity of its
// CoverageWindow.
static double green_fraction(const Box& box,
                             const std::vector<const LightTraversal*>& valid,
                             double left_dy) {
  const double w = box.xmax - box.xmin, h = box.ymax - box.ymin;
  if (!(w > 0.0 && h > 0.0)) return 0.0;
  double a = 0.0;
  for (const auto* t : valid) {
    const auto& tc = t->coords;
    for (size_t i = 1; i < tc.size(); i++) {
      double xm = 0.5 * ((tc[i - 1].x - box.xmax) + (tc[i].x - box.xmax));
      a += xm * (tc[i].y - tc[i - 1].y);
    }
  }
  // Magnitude, as in the single-traversal formula: a lobe of a
  // self-intersecting ring that runs clockwise after normalisation is still
  // covered under the non-zero winding rule.
  return std::min(std::fabs((a - w * left_dy) / (w * h)), 1.0);
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
    int ring,
    const exactextract::Grid<exactextract::infinite_extent>& grid,
    double /* dy (unused, kept for signature compat) */,
    std::vector<std::vector<BoundaryCellRecord>>& row_data,
    size_t sub_rows,
    size_t sub_cols,
    size_t row_off,
    size_t col_off,
    size_t keep_lo,
    size_t keep_hi,
    BurnMode mode = BurnMode::Coverage
) {
  if (coords.size() < 4) return;

  // Normalise to CCW for correct coverage fraction semantics
  if (!is_ccw) {
    std::reverse(coords.begin(), coords.end());
  }

  const bool is_exterior = (ring == 0);
  float coverage_factor = is_exterior ? 1.0f : -1.0f;
  int winding_factor = is_exterior ? 1 : -1;

  CellMap cells = walk_polyline(std::move(coords), grid, /*closed=*/true);

  // Net height of the ring's traversals in the cells left of each cell in
  // its row (see green_fraction). The map is ordered by (row, col), so
  // this is a running sum reset at each row.
  std::map<CellKey, double> left_dy;
  {
    size_t cur_row = static_cast<size_t>(-1);
    double acc = 0.0;
    for (const auto& kv : cells) {
      if (kv.first.first != cur_row) { cur_row = kv.first.first; acc = 0.0; }
      left_dy.emplace_hint(left_dy.end(), kv.first, acc);
      for (const auto& t : kv.second.traversals) {
        if (!valid_traversal(t) || !t.traversed()) continue;
        acc += t.coords.back().y - t.coords.front().y;
      }
    }
  }

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
    // Rows are independent from here on: a window sweeps only its rows.
    if (sub_r < keep_lo || sub_r >= keep_hi) continue;

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
    std::vector<const LightTraversal*> valid;
    for (const auto& t : cr.traversals) {
      if (valid_traversal(t)) valid.push_back(&t);
    }

    if (valid.empty()) continue;

    // ---- Coverage fraction (only for in-grid cells, Coverage mode) ----
    float frac = 0.0f;

    if (in_grid_cols && mode == BurnMode::Coverage) {
      if (valid.size() == 1 && valid[0]->entry_side == Side::NONE
            && valid[0]->is_closed_ring()) {
        frac = static_cast<float>(
          controlledburn::closed_ring_covered_fraction(cr.box, valid[0]->coords));
      } else if (valid.size() == 1 && !along_walls(cr.box, valid[0]->coords)) {
        frac = static_cast<float>(
          controlledburn::analytical_covered_fraction(
            cr.box, valid[0]->coords,
            valid[0]->entry_side, valid[0]->exit_side));
      } else {
        // several traversals, or one along the walls
        frac = static_cast<float>(green_fraction(cr.box, valid, left_dy[kv.first]));
      }
    }

    // ---- Store coverage (if nonzero) and winding deltas ----
    //
    // IMPORTANT: winding deltas must be stored even when coverage is zero.
    // A traversal along a cell wall (e.g. vertical edge on a grid line)
    // has zero area but still crosses the row center, contributing to the
    // winding count that classifies interior cells.

    // Records are pushed per visit (one for coverage, one for winding) and
    // merged in the sweep after a stable sort by column, so contributions
    // add in visit order, as in controlledburn-rs. A shared record per
    // column (the earlier find_or_create) lost which ring each part came
    // from, which the hole handling in the sweep needs; it also created a
    // record for crossings that cancel within a visit, which split an
    // interior run and dropped the cell.
    auto& row_vec = row_data[sub_r];

    int visit_delta = 0;
    int visit_delta_left = 0;

    // Winding deltas from traversals
    double x_mid = (cr.box.xmin + cr.box.xmax) / 2.0;

    for (auto* t : valid) {
      if (!t->traversed()) continue; // closed rings don't contribute winding
      if (t->coords.size() < 2) continue;

      double entry_y = t->coords.front().y;
      double exit_y  = t->coords.back().y;
      double y_mid = (cr.box.ymin + cr.box.ymax) / 2.0;

      // Half-open at the row centre in both modes: a point exactly on
      // y_mid counts as above it. A traversal that ends on the centre line
      // (an edge crossing it on a column wall) and the next one that starts
      // there then register the crossing exactly once, in one cell. With a
      // strict test on both sides in coverage mode neither did, and the row
      // lost its winding from that column on (interior runs vanished); this
      // is common once vertices sit on the lattice. Coverage runs are only
      // filled in the gaps between boundary cells, so a delta on a boundary
      // cell cannot double-count against its edge fraction.
      bool crosses = (entry_y >= y_mid && exit_y < y_mid) ||
                     (entry_y < y_mid && exit_y >= y_mid);
      if (!crosses) continue;

      int delta = (entry_y >= y_mid) ? -1 : +1; // downward = -1, upward = +1
      delta *= winding_factor;

      visit_delta += delta;

      // Approx mode: classify whether this crossing is left of cell center.
      // Interpolate x at y_mid along the traversal path. Skip horizontal
      // segments (y0 == y1) — they lie ON y_mid but don't contribute a
      // directional crossing; the real crossing is on an adjacent
      // non-horizontal segment.
      if (mode == BurnMode::Approx) {
        double x_at_mid = x_mid; // fallback
        const auto& tc = t->coords;
        for (size_t i = 1; i < tc.size(); i++) {
          double y0 = tc[i-1].y, y1 = tc[i].y;
          if (y0 == y1) continue; // skip horizontal segments
          if ((y0 <= y_mid && y1 >= y_mid) || (y0 >= y_mid && y1 <= y_mid)) {
            double frac_along = (y_mid - y0) / (y1 - y0);
            x_at_mid = tc[i-1].x + frac_along * (tc[i].x - tc[i-1].x);
            break;
          }
        }
        // Use <= for left-inclusive convention: an edge crossing exactly
        // at the cell center x is considered "left of center", matching
        // the standard rasterization rule (left/bottom inclusive,
        // right/top exclusive).
        if (x_at_mid <= x_mid) {
          visit_delta_left += delta;
        }
      }
    }
    // A cell with zero coverage still needs a record when the ring runs
    // through it, from one point on its boundary to another: otherwise the
    // sweep fills it with the winding from its left. Edges along both walls
    // of a cell (lattice-aligned input) give zero coverage and zero net
    // winding, one crossing on each wall cancelling, and the cell was
    // filled as interior. A traversal that leaves where it entered
    // (touching a corner, or a zero-width spike in and back out) encloses
    // nothing, and its coverage says nothing about the cell: no record, the
    // sweep fills it from the winding. That holds for coverage that is zero
    // only up to rounding, too: a spike out to a corner vertex and back can
    // return an ulp from where it entered, with a coverage of 1e-16, and
    // that record kept a fully covered cell empty (#25). Below the sweep's
    // tolerance (1e-6) such a record can only blank the cell, so it needs a
    // pass-through.
    const bool tiny = std::fabs(frac) < 1e-6f;
    bool passes_through = false;
    if (tiny && in_grid_cols && mode == BurnMode::Coverage) {
      const double tol_x = (cr.box.xmax - cr.box.xmin) * 1e-9;
      const double tol_y = (cr.box.ymax - cr.box.ymin) * 1e-9;
      for (auto* t : valid) {
        const Coordinate& a = t->coords.front();
        const Coordinate& b = t->coords.back();
        if (std::fabs(a.x - b.x) > tol_x || std::fabs(a.y - b.y) > tol_y) {
          passes_through = true;
          break;
        }
      }
    }
    if (!tiny || passes_through) {
      row_vec.push_back({full_col, coverage_factor * frac, 0, 0, ring});
    }
    if (visit_delta != 0 || visit_delta_left != 0) {
      row_vec.push_back({full_col, 0.0f, visit_delta, visit_delta_left, ring});
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

// ---- A window of the grid (burn_window) ----
//
// Rows [r0, r1) and columns [c0, c1) of the full grid, and `near`, its box
// widened by two cells each side (lines are walked only there).
struct Window {
  int r0, r1, c0, c1;
  exactextract::Box near;

  Window(const GridSpec& gs, uint32_t row0, uint32_t row1, uint32_t col0, uint32_t col1)
      : r0(static_cast<int>(row0)), r1(static_cast<int>(row1)),
        c0(static_cast<int>(col0)), c1(static_cast<int>(col1)),
        near(exactextract::Box::make_empty()) {
    double dx = gs.dx(), dy = gs.dy();
    auto x = [&](double c) { return gs.xmin + c * dx; };
    auto y = [&](double r) { return gs.ymax - r * dy; };
    double rr0 = static_cast<double>(row0) - 2.0, rr1 = static_cast<double>(row1) + 2.0;
    double cc0 = static_cast<double>(col0) - 2.0, cc1 = static_cast<double>(col1) + 2.0;
    near = exactextract::Box(x(cc0), y(rr1), x(cc1), y(rr0));
  }

  bool has(int row, int col) const {
    return row >= r0 && row < r1 && col >= c0 && col < c1;
  }
};

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
    const Window* win = nullptr,
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

  // A window keeps the sub-grid the whole burn uses, so the walk, the snap
  // and the fractions are its own, and sweeps only the window's rows; the
  // caller cuts the columns.
  size_t keep_lo = 0, keep_hi = sub_rows;
  if (win) {
    size_t wc0 = static_cast<size_t>(win->c0), wc1 = static_cast<size_t>(win->c1);
    if (wc1 <= col_off || wc0 >= col_off + sub_cols) return;
    size_t wr0 = static_cast<size_t>(win->r0), wr1 = static_cast<size_t>(win->r1);
    keep_lo = std::min(wr0 > row_off ? wr0 - row_off : 0, sub_rows);
    keep_hi = std::min(wr1 > row_off ? wr1 - row_off : 0, sub_rows);
  }
  if (keep_lo >= keep_hi) return;

  std::vector<std::vector<BoundaryCellRecord>> row_data(sub_rows);

  // Exterior ring, then holes. Orientation from shoelace sign (replaces
  // geos_is_ccw); walk_ring normalises to CCW internally.
  for (size_t r = 0; r < poly.rings.size(); r++) {
    const CoordSeq& ring = poly.rings[r];
    if (ring.size() < 4) continue;
    walk_ring(to_coordinates(ring), is_ccw(ring), static_cast<int>(r), subgrid, dy,
              row_data, sub_rows, sub_cols, row_off, col_off, keep_lo, keep_hi, mode);
  }

  // ---- Winding sweep: emit runs and edges per row ----

  float tol = 1e-6f;

  // Winding per ring, so a boundary cell gains the cover of rings that
  // enclose it without passing through it (the exterior around a hole's
  // boundary cell). Allocated once; each row resets only what it touched.
  const size_t n_rings = poly.rings.size();
  std::vector<int> ring_winding(n_rings, 0);
  std::vector<int> touched;
  // Rings passing through each merged cell, flat, with a span per cell;
  // each ring appears once per cell with its summed winding delta.
  std::vector<std::pair<int, int>> through;
  std::vector<std::pair<size_t, size_t>> spans;
  auto sign = [](int r) { return r == 0 ? 1 : -1; };

  for (size_t sr = keep_lo; sr < keep_hi; sr++) {
    auto& row_vec = row_data[sr];
    if (row_vec.empty()) continue;

    // Stable, so same-column contributions merge in visit order.
    std::stable_sort(row_vec.begin(), row_vec.end(),
                     [](const BoundaryCellRecord& a, const BoundaryCellRecord& b) {
                       return a.col < b.col;
                     });

    // Merge same-column entries, keeping the rings through each cell.
    std::vector<BoundaryCellRecord> merged;
    through.clear();
    spans.clear();
    for (auto& rec : row_vec) {
      if (!merged.empty() && merged.back().col == rec.col) {
        merged.back().coverage += rec.coverage;
        merged.back().winding_delta += rec.winding_delta;
        merged.back().winding_delta_left += rec.winding_delta_left;
        auto& span = spans.back();
        bool found = false;
        for (size_t i = span.first; i < span.second; i++) {
          if (through[i].first == rec.ring) {
            through[i].second += rec.winding_delta;
            found = true;
            break;
          }
        }
        if (!found) {
          through.push_back({rec.ring, rec.winding_delta});
          span.second++;
        }
      } else {
        merged.push_back(rec);
        spans.push_back({through.size(), through.size() + 1});
        through.push_back({rec.ring, rec.winding_delta});
      }
    }

    int winding = 0;
    // Sum of sign(r) over rings with non-zero winding: +1 inside the
    // exterior, -1 inside a hole. Kept up to date as windings change, so
    // a cell costs only the rings through it, not all the rings.
    int inside = 0;
    int prev_col = -2;
    int full_row = static_cast<int>(row_off + sr);

    for (size_t mi = 0; mi < merged.size(); mi++) {
      const auto& mc = merged[mi];
      const size_t s0 = spans[mi].first, s1 = spans[mi].second;
      // Emit interior run between previous boundary cell and this one.
      // prev_col > -2 means at least one cell has been seen (including padding).
      if (winding != 0 && prev_col > -2 && mc.col > prev_col + 1) {
        all_runs.push_back({
          full_row,
          prev_col + 1,
          mc.col,
          poly_id
        });
      }

      if (mode == BurnMode::Coverage) {
        // +1 inside the exterior, -1 inside a hole, for rings that do not
        // pass through this cell.
        int enclosing = inside;
        for (size_t i = s0; i < s1; i++) {
          if (ring_winding[through[i].first] != 0) enclosing -= sign(through[i].first);
        }
        float w = enclosing != 0 ? mc.coverage + static_cast<float>(enclosing) : mc.coverage;
        if (w > tol && w < (1.0f - tol)) {
          all_edges.push_back({full_row, mc.col, w, poly_id});
        } else if (w >= (1.0f - tol)) {
          all_runs.push_back({full_row, mc.col, mc.col + 1, poly_id});
        }
      } else {
        // Approx: cell center is inside iff winding at center != 0.
        // Winding at center = winding (from left) + delta_left (crossings
        // left of center x within this cell).
        int winding_at_center = winding + mc.winding_delta_left;
        if (winding_at_center != 0) {
          all_runs.push_back({full_row, mc.col, mc.col + 1, poly_id});
        }
      }

      winding += mc.winding_delta;
      prev_col = mc.col;
      for (size_t i = s0; i < s1; i++) {
        const int r = through[i].first;
        const int before = ring_winding[r];
        const int after = before + through[i].second;
        if (before == 0 && after != 0) {
          inside += sign(r);
          touched.push_back(r);
        } else if (before != 0 && after == 0) {
          inside -= sign(r);
        }
        ring_winding[r] = after;
      }
    }
    for (int r : touched) ring_winding[r] = 0;
    touched.clear();
  }
}

// ---- Lightweight approx sweep ----
//
// orient2d_sign: exact orientation, in predicates.h (shared with relate.cpp).

// Direct edge-row intersection sweep for Approx mode. Bypasses the
// walker entirely: for each polygon edge, compute x-intercepts at
// each row's y_mid, accumulate winding per row, sweep left-to-right
// to emit runs. O(edges * rows_touched) with minimal per-cell
// overhead -- no CellRecords, no traversal coordinates, no
// BoundaryCellRecords.

static void process_polygon_approx(
    const Polygon& poly,
    const GridSpec& gs,
    int poly_id,
    std::vector<GridRun>& all_runs,
    const Window* win = nullptr
) {
  if (poly.rings.empty() || poly.rings[0].size() < 4) return;

  double dx = gs.dx();
  double dy = gs.dy();

  // Clip to grid extent
  Box env = ring_envelope(poly.rings[0]);
  Box grid_box(gs.xmin, gs.ymin, gs.xmax, gs.ymax);
  if (!env.intersects(grid_box)) return;

  // Row range that the polygon can touch (0-based).
  // Row 0 is the TOP row: y_mid = ymax - 0.5*dy
  // A window sweeps only its own rows.
  int lo = win ? win->r0 : 0;
  int hi = win ? win->r1 - 1 : gs.nrow - 1;
  int row_lo = std::max(lo, static_cast<int>(
      std::floor((gs.ymax - env.ymax) / dy)));
  int row_hi = std::min(hi, static_cast<int>(
      std::floor((gs.ymax - env.ymin) / dy)));
  if (row_lo > row_hi) return;

  int n_sweep_rows = row_hi - row_lo + 1;

  struct Intercept {
    double x;
    int delta;
    // Column of the crossing (clamped to the grid) and whether it lies left
    // of or on that column's centre, decided exactly.
    int col;
    bool left;
  };
  std::vector<std::vector<Intercept>> row_intercepts(
      static_cast<size_t>(n_sweep_rows));

  for (size_t ri = 0; ri < poly.rings.size(); ri++) {
    const CoordSeq& ring = poly.rings[ri];
    if (ring.size() < 4) continue;

    bool exterior = (ri == 0);
    int wf = is_ccw(ring) ? 1 : -1;
    if (!exterior) wf = -wf;

    for (size_t i = 0; i < ring.size() - 1; i++) {
      double x0 = ring[i].x, y0 = ring[i].y;
      double x1 = ring[i + 1].x, y1 = ring[i + 1].y;

      if (y0 == y1) continue; // horizontal edge

      // Ensure ya < yb for row range computation
      double ya = y0, yb = y1, xa = x0, xb = x1;
      if (ya > yb) { std::swap(ya, yb); std::swap(xa, xb); }

      // Which rows have y_mid in (ya, yb]? Top-inclusive, bottom-
      // exclusive — matching the walker's half-open crossing convention
      // (entry_y >= y_mid && exit_y < y_mid, which gives ya < y_mid <= yb).
      // y_mid(r) = ymax - (r + 0.5) * dy
      // ya < y_mid  →  r < (ymax - ya)/dy - 0.5
      // y_mid <= yb →  r >= (ymax - yb)/dy - 0.5
      int e_row_lo = std::max(row_lo, static_cast<int>(
          std::ceil((gs.ymax - yb) / dy - 0.5 - 1e-10)));
      int e_row_hi = std::min(row_hi, static_cast<int>(
          std::ceil((gs.ymax - ya) / dy - 0.5 - 1e-10) - 1));
      if (e_row_lo > e_row_hi) continue;

      double inv_dy_edge = 1.0 / (y1 - y0);

      for (int r = e_row_lo; r <= e_row_hi; r++) {
        double y_mid = gs.ymax - (r + 0.5) * dy;
        double t = (y_mid - y0) * inv_dy_edge;
        double x_int = x0 + t * (x1 - x0);

        // Crossing direction from the edge's own endpoints, not from y0
        // against the recomputed y_mid: the row-range test above admits an
        // edge whose top endpoint is on the centre line (within 1e-10), but
        // that endpoint can sit an ulp below y_mid as recomputed from ymax
        // and dy, and (y0 >= y_mid) then tags a downward edge as upward. Its
        // intercept gets the wrong sign, the winding never returns to zero,
        // and the row fills to the grid edge.
        int delta = (y1 < y0) ? -1 : +1;
        delta *= wf;

        // Left of or on the centre (the left-inclusive rule), decided
        // exactly from the edge's endpoints taken bottom to top, not from
        // the rounded intercept: two edges through nearly the same point (a
        // sliver, or a ring doubling back along itself) give intercepts that
        // round either side of a centre they both miss, and the cell burned.
        int col = static_cast<int>(std::floor((x_int - gs.xmin) / dx));
        if (col < 0) col = 0;
        if (col >= gs.ncol) col = gs.ncol - 1;
        double x_mid = gs.xmin + (col + 0.5) * dx;
        bool left = orient2d_sign(xa, ya, xb, yb, x_mid, y_mid) <= 0;

        row_intercepts[static_cast<size_t>(r - row_lo)].push_back(
            {x_int, delta, col, left});
      }
    }
  }

  // Sweep: sort intercepts per row, accumulate winding, emit runs
  for (int r = row_lo; r <= row_hi; r++) {
    auto& intercepts = row_intercepts[static_cast<size_t>(r - row_lo)];
    if (intercepts.empty()) continue;

    // By column, then side of its centre, then x: only the side of each
    // centre matters, and this keeps the order consistent with the exact
    // side test where rounded intercepts tie or cross. Stable, as in
    // controlledburn-rs.
    std::stable_sort(intercepts.begin(), intercepts.end(),
              [](const Intercept& a, const Intercept& b) {
                if (a.col != b.col) return a.col < b.col;
                if (a.left != b.left) return a.left;
                return a.x < b.x;
              });

    int winding = 0;
    int run_start = -1;

    for (const auto& ix : intercepts) {
      int col = ix.col;
      bool left_of_center = ix.left;

      if (left_of_center) {
        int new_winding = winding + ix.delta;
        if (winding != 0 && new_winding == 0) {
          // Leaving polygon before cell center: close run before this cell
          int col_end = col;  // 0-based, exclusive (last col is col-1)
          if (run_start >= 0 && run_start < col_end)
            all_runs.push_back({r, run_start, col_end, poly_id});
          run_start = -1;
        } else if (winding == 0 && new_winding != 0) {
          // Entering polygon: this cell is inside
          run_start = col;
        }
        winding = new_winding;
      } else {
        int new_winding = winding + ix.delta;
        if (winding != 0 && new_winding == 0) {
          // Leaving polygon after cell center: include this cell
          if (run_start < 0) run_start = col;
          int col_end_1 = std::min(col + 1, gs.ncol);  // clamp to grid
          if (run_start < col_end_1)
            all_runs.push_back({r, run_start, col_end_1, poly_id});
          run_start = -1;
        } else if (winding == 0 && new_winding != 0) {
          // Entering polygon after cell center: start after this cell
          run_start = col + 1;
        }
        winding = new_winding;
      }
    }

    // Polygon extends beyond grid right edge
    if (winding != 0 && run_start >= 0 && run_start < gs.ncol) {
      all_runs.push_back({r, run_start, gs.ncol, poly_id});
    }
  }
}


// ---- Process line: per-LINESTRING length-in-cell rasterization ----
//
// Walks the line through the grid using walk_polyline (closed=false) and
// emits one (row, col, length, id) record per cell touched. Length is the
// absolute length of the line within the cell, in CRS units -- sum of
// segment lengths across all traversal records for that cell.

// The parts of `line` to walk for the cells in `near`, as [start, end]
// vertex index pairs: each run of segments whose box meets `near`, started
// at the last vertex before it that lies inside one cell, clear of every
// wall (or at the line's start), and runs that then overlap joined. From
// such a vertex on, a walk started there and the whole walk are the same
// (the walker computes each crossing from the segment's own vertices); the
// vertices a part's walk starts and ends at belong to segments that miss
// `near`, so where it differs from the whole walk is not the window's.
// Same as near_runs in controlledburn-rs.
static std::vector<std::pair<size_t, size_t>> near_runs(
    const CoordSeq& line,
    const exactextract::Grid<exactextract::bounded_extent>& grid,
    const exactextract::Box& b
) {
  auto meets = [&](const Coord& p, const Coord& q) {
    return !(std::max(p.x, q.x) < b.xmin || std::min(p.x, q.x) > b.xmax ||
             std::max(p.y, q.y) < b.ymin || std::min(p.y, q.y) > b.ymax);
  };
  const double dx = grid.dx(), dy = grid.dy();
  // well clear of the walls: of rounding in `xmin + k * dx` and of the
  // walker's snap (1e-9 of a cell)
  auto clear = [](double v, double d) {
    double f = v / d;
    return std::fabs(f - std::round(f)) * d > 1e-6 * d + 1e-13 * std::fabs(v);
  };
  auto inside_one_cell = [&](const Coord& p) {
    return clear(p.x - grid.xmin(), dx) && clear(grid.ymax() - p.y, dy);
  };
  std::vector<std::pair<size_t, size_t>> runs;
  size_t i = 0;
  const size_t n = line.size() - 1;
  while (i < n) {
    if (!meets(line[i], line[i + 1])) {
      i++;
      continue;
    }
    size_t start = i;
    while (start > 0 && !inside_one_cell(line[start])) start--;
    size_t end = i + 1;
    while (end < n && meets(line[end], line[end + 1])) end++;
    if (!runs.empty() && start <= runs.back().second) {
      runs.back().second = end;
    } else {
      runs.push_back({start, end});
    }
    i = end;
  }
  return runs;
}

static void process_line(
    const CoordSeq& line,
    const exactextract::Grid<exactextract::bounded_extent>& full_grid,
    int line_id,
    std::vector<GridLine>& all_lines,
    const exactextract::Box* near = nullptr
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

  std::vector<std::pair<size_t, size_t>> parts;
  if (near) {
    parts = near_runs(line, full_grid, *near);
  } else {
    parts.push_back({0, line.size() - 1});
  }

  // (padded row, col) -> summed length. One walk visits a cell in one
  // record; with several parts (a window's), a cell's lengths are summed
  // over the parts in line order, so cells near the window get the whole
  // walk's lengths, bit for bit.
  std::map<CellKey, double> lengths;
  for (const auto& part : parts) {
    CoordSeq piece(line.begin() + static_cast<std::ptrdiff_t>(part.first),
                   line.begin() + static_cast<std::ptrdiff_t>(part.second) + 1);
    CellMap cells = walk_polyline(to_coordinates(piece), grid, /*closed=*/false);

    for (const auto& kv : cells) {
      size_t r = kv.first.first;
      size_t c = kv.first.second;
      const CellRecord& cr = kv.second;

      // Skip padding ring (top/bottom row, left/right col of infinite_extent).
      if (r < 1 || r >= n_rows - 1) continue;
      if (c < 1 || c >= n_cols - 1) continue;

      double& total_length = lengths[kv.first];
      for (const auto& t : cr.traversals) {
        const auto& tc = t.coords;
        if (tc.size() < 2) continue;
        for (size_t i = 1; i < tc.size(); i++) {
          double sx = tc[i].x - tc[i-1].x;
          double sy = tc[i].y - tc[i-1].y;
          total_length += std::sqrt(sx * sx + sy * sy);
        }
      }
    }
  }

  // Row-major record order (the map's).
  for (const auto& kv : lengths) {
    if (kv.second > 0.0) {
      // Map infinite_extent (r, c) to the 0-based bounded grid; padding=1
      // means subtracting 1 recovers the bounded index.
      all_lines.push_back({
        static_cast<int>(kv.first.first) - 1, static_cast<int>(kv.first.second) - 1,
        static_cast<float>(kv.second),
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
    static_cast<int>(r), // 0-based
    static_cast<int>(c),
    point_id
  });
}

// ---- Process geometry: dispatch by kind ----

using FullGrid = exactextract::Grid<exactextract::bounded_extent>;

// The sweep's tolerance on a cell's coverage (process_polygon).
constexpr float kTol = 1e-6f;

static void note(BurnResult& out, int32_t id, std::string msg) {
  out.notes.push_back({id, std::move(msg)});
}

// Keep the records from run_start, edge_start and line_start on that lie
// in `win`, cutting runs at its columns.
static void clip(BurnResult& out, const Window& win, size_t run_start,
                 size_t edge_start, size_t line_start) {
  size_t w = run_start;
  for (size_t i = run_start; i < out.runs.size(); i++) {
    GridRun r = out.runs[i];
    r.col_start = std::max(r.col_start, win.c0);
    r.col_end = std::min(r.col_end, win.c1);
    if (r.row >= win.r0 && r.row < win.r1 && r.col_start < r.col_end) {
      out.runs[w++] = r;
    }
  }
  out.runs.resize(w);
  w = edge_start;
  for (size_t i = edge_start; i < out.edges.size(); i++) {
    if (win.has(out.edges[i].row, out.edges[i].col)) out.edges[w++] = out.edges[i];
  }
  out.edges.resize(w);
  w = line_start;
  for (size_t i = line_start; i < out.lines.size(); i++) {
    if (win.has(out.lines[i].row, out.lines[i].col)) out.lines[w++] = out.lines[i];
  }
  out.lines.resize(w);
}

// Merge adjacent runs of runs[start..), which are ordered by (row,
// col_start); the runs before start are left alone.
static void coalesce_sorted_from(std::vector<GridRun>& runs, size_t start) {
  size_t w = start;
  for (size_t i = start; i < runs.size(); i++) {
    GridRun cur = runs[i];
    if (w > start) {
      GridRun& last = runs[w - 1];
      if (last.id == cur.id && last.row == cur.row && cur.col_start <= last.col_end) {
        if (cur.col_end > last.col_end) last.col_end = cur.col_end;
        continue;
      }
    }
    runs[w++] = cur;
  }
  runs.resize(w);
}

// One edge per cell for the edges of one geometry, edges[start..): the
// parts' fractions summed per cell (in part order), a cell whose sum
// reaches 1 made a one-cell run.
static void merge_part_edges(BurnResult& out, size_t start) {
  auto& edges = out.edges;
  std::stable_sort(edges.begin() + static_cast<std::ptrdiff_t>(start), edges.end(),
                   [](const GridEdge& a, const GridEdge& b) {
                     return a.row != b.row ? a.row < b.row : a.col < b.col;
                   });
  size_t w = start;
  for (size_t i = start; i < edges.size(); i++) {
    GridEdge cur = edges[i];
    if (w > start) {
      GridEdge& last = edges[w - 1];
      if (last.row == cur.row && last.col == cur.col) {
        last.fraction += cur.fraction;
        continue;
      }
    }
    edges[w++] = cur;
  }
  edges.resize(w);
  w = start;
  for (size_t i = start; i < edges.size(); i++) {
    GridEdge e = edges[i];
    if (e.fraction >= 1.0f - kTol) {
      out.runs.push_back({e.row, e.col, e.col + 1, e.id});
    } else {
      edges[w++] = e;
    }
  }
  edges.resize(w);
}

// Drop the edges of one geometry (edges[edge_start..), sorted by cell)
// that lie in one of its runs (runs[run_start..), sorted and coalesced).
static void drop_edges_in_runs(BurnResult& out, size_t run_start, size_t edge_start) {
  const auto& runs = out.runs;
  auto& edges = out.edges;
  size_t j = run_start;
  size_t w = edge_start;
  for (size_t i = edge_start; i < edges.size(); i++) {
    GridEdge e = edges[i];
    while (j < runs.size() &&
           (runs[j].row < e.row || (runs[j].row == e.row && runs[j].col_end <= e.col))) {
      j++;
    }
    bool covered = j < runs.size() && runs[j].row == e.row && runs[j].col_start <= e.col;
    if (!covered) edges[w++] = e;
  }
  edges.resize(w);
}

// One record per cell for the line parts of one geometry, lines[start..):
// lengths summed per cell, in part order.
static void merge_part_lines(BurnResult& out, size_t start) {
  auto& lines = out.lines;
  std::stable_sort(lines.begin() + static_cast<std::ptrdiff_t>(start), lines.end(),
                   [](const GridLine& a, const GridLine& b) {
                     return a.row != b.row ? a.row < b.row : a.col < b.col;
                   });
  size_t w = start;
  for (size_t i = start; i < lines.size(); i++) {
    GridLine cur = lines[i];
    if (w > start) {
      GridLine& last = lines[w - 1];
      if (last.row == cur.row && last.col == cur.col) {
        last.length += cur.length;
        continue;
      }
    }
    lines[w++] = cur;
  }
  lines.resize(w);
}

// A note for each ring the engines skip as too short: fewer than 4
// coordinates once closed. A short exterior skips its part, holes and all;
// a short hole skips only itself. Empty rings are skipped silently.
static void note_short_rings(const std::vector<Polygon>& polys, int32_t id, BurnResult& out) {
  for (size_t pi = 0; pi < polys.size(); pi++) {
    const auto& rings = polys[pi].rings;
    for (size_t ri = 0; ri < rings.size(); ri++) {
      size_t n = rings[ri].size();
      if (n == 0 || n >= 4) continue;
      if (ri == 0) {
        note(out, id, "skipped polygon part " + std::to_string(pi) + ": exterior ring has " +
                          std::to_string(n) + " coordinates when closed, 4 needed");
        break;
      }
      note(out, id, "skipped hole " + std::to_string(ri) + " of polygon part " +
                        std::to_string(pi) + ": " + std::to_string(n) +
                        " coordinates when closed, 4 needed");
    }
  }
}

static void process_geometry(
    const Geometry& g,
    const FullGrid& full_grid,
    const GridSpec& gs,
    double dx, double dy,
    int geom_id,
    BurnResult& out,
    const BurnOptions& opts,
    const Window* win
) {
  const size_t run_start = out.runs.size();
  const size_t edge_start = out.edges.size();
  const size_t line_start = out.lines.size();
  const auto& polys = g.polygons;
  note_short_rings(polys, geom_id, out);
  for (const auto& poly : polys) {
    if (opts.mode == BurnMode::Approx) {
      process_polygon_approx(poly, gs, geom_id, out.runs, win);
    } else {
      process_polygon(poly, full_grid, dx, dy, geom_id, out.runs, out.edges, win);
    }
  }
  if (win) clip(out, *win, run_start, edge_start, line_start);
  const bool merge = opts.merge_parts && polys.size() > 1;
  if (merge) merge_part_edges(out, edge_start);
  if ((opts.coalesce_runs || merge) && out.runs.size() > run_start) {
    // One polygon's runs come out of the sweep ordered by (row,
    // col_start); a multipolygon's are in that order only per part, so
    // sort them together first, or a later part's run left of an earlier
    // part's on a shared row is swallowed by the merge.
    if (polys.size() > 1) {
      std::stable_sort(out.runs.begin() + static_cast<std::ptrdiff_t>(run_start), out.runs.end(),
                       [](const GridRun& a, const GridRun& b) {
                         return a.row != b.row ? a.row < b.row : a.col_start < b.col_start;
                       });
    }
    coalesce_sorted_from(out.runs, run_start);
  }
  if (merge) drop_edges_in_runs(out, run_start, edge_start);

  for (size_t i = 0; i < g.lines.size(); i++) {
    if (g.lines[i].size() == 1) {
      note(out, geom_id, "skipped line part " + std::to_string(i) + ": 1 coordinate, 2 needed");
    }
    process_line(g.lines[i], full_grid, geom_id, out.lines, win ? &win->near : nullptr);
  }
  if (win) clip(out, *win, out.runs.size(), out.edges.size(), line_start);
  if (opts.merge_parts && g.lines.size() > 1) merge_part_lines(out, line_start);

  const size_t point_start = out.points.size();
  for (const auto& pt : g.points) {
    process_point(pt, full_grid, geom_id, out.points);
  }
  if (win) {
    size_t w = point_start;
    for (size_t i = point_start; i < out.points.size(); i++) {
      if (win->has(out.points[i].row, out.points[i].col)) out.points[w++] = out.points[i];
    }
    out.points.resize(w);
  }
}

// WKB's POINT EMPTY: both coordinates NaN.
static bool is_empty_point(const CoordSeq& p) {
  return p.empty() || (std::isnan(p[0].x) && std::isnan(p[0].y));
}

// As Geometry::is_empty in controlledburn-rs, by kind.
static bool geometry_is_empty(const Geometry& g) {
  switch (g.kind) {
  case GeomKind::Point:
  case GeomKind::MultiPoint:
    for (const auto& p : g.points) if (!is_empty_point(p)) return false;
    return true;
  case GeomKind::LineString:
    return g.lines.empty() || g.lines[0].empty();
  case GeomKind::MultiLineString:
    return g.lines.empty();
  case GeomKind::Polygon:
    return g.polygons.empty() || g.polygons[0].rings.empty();
  case GeomKind::MultiPolygon:
    return g.polygons.empty();
  }
  return g.empty();
}

static bool finite(const Coord& c) { return std::isfinite(c.x) && std::isfinite(c.y); }

// Every coordinate finite, empty points aside.
static bool geometry_is_finite(const Geometry& g) {
  for (const auto& p : g.points)
    for (const auto& c : p)
      if (!finite(c) && !(std::isnan(c.x) && std::isnan(c.y))) return false;
  for (const auto& l : g.lines)
    for (const auto& c : l)
      if (!finite(c)) return false;
  for (const auto& p : g.polygons)
    for (const auto& r : p.rings)
      for (const auto& c : r)
        if (!finite(c)) return false;
  return true;
}

static bool ring_needs_closing(const CoordSeq& r) {
  return r.size() >= 2 && (r.front().x != r.back().x || r.front().y != r.back().y);
}

// Burn one geometry: skip empties, note non-finite ones, close unclosed
// rings, and turn engine errors into notes.
static void burn_one(const Geometry& g0, int32_t id, const FullGrid& full_grid,
                     const GridSpec& gs, BurnResult& out, const BurnOptions& opts,
                     const Window* win) {
  if (geometry_is_empty(g0)) return;
  if (!geometry_is_finite(g0)) {
    note(out, id, "skipped geometry with non-finite coordinates");
    return;
  }
  bool needs = false;
  for (const auto& p : g0.polygons)
    for (const auto& r : p.rings)
      if (ring_needs_closing(r)) needs = true;
  Geometry closed;
  if (needs) {
    closed = g0;
    for (auto& p : closed.polygons)
      for (auto& r : p.rings)
        if (ring_needs_closing(r)) r.push_back(r.front());
  }
  const Geometry& g = needs ? closed : g0;
  try {
    process_geometry(g, full_grid, gs, gs.dx(), gs.dy(), id, out, opts, win);
  } catch (const std::exception& e) {
    note(out, id, std::string("error processing geometry: ") + e.what());
  }
}

// The id of geometry k, checked.
static int32_t geom_id(const BurnOptions& opts, size_t k) {
  if (k > static_cast<size_t>(INT32_MAX) ||
      static_cast<int64_t>(opts.id_offset) + static_cast<int64_t>(k) > INT32_MAX) {
    throw std::invalid_argument("geometry " + std::to_string(k) + " with id_offset " +
                                std::to_string(opts.id_offset) + ": id beyond i32::MAX");
  }
  return opts.id_offset + static_cast<int32_t>(k);
}

static void validate(const GridSpec& gs, const BurnOptions& opts, size_t n) {
  gs.validate();
  if (opts.id_offset < 0) {
    throw std::invalid_argument("id_offset " + std::to_string(opts.id_offset) + " is negative");
  }
  if (n > 0) geom_id(opts, n - 1);
}

static FullGrid make_full_grid(const GridSpec& gs) {
  return FullGrid(exactextract::Box(gs.xmin, gs.ymin, gs.xmax, gs.ymax), gs.dx(), gs.dy());
}

// A parsed WKB blob: its geometry, or the parse error, or nothing (an
// empty blob, skipped silently).
struct Parsed {
  Geometry g;
  std::string error;
  bool skip = false;
  bool failed = false;
};

static Parsed parse_one(const WKBSpan& s) {
  Parsed p;
  if (!s.data || s.size == 0) {
    p.skip = true;
    return p;
  }
  try {
    p.g = parse_wkb(s.data, s.size);
  } catch (const WKBParseError& e) {
    p.failed = true;
    p.error = e.what();
  }
  return p;
}

static void burn_parsed(const Parsed& p, int32_t id, const FullGrid& full_grid,
                        const GridSpec& gs, BurnResult& out, const BurnOptions& opts,
                        const Window* win) {
  if (p.skip) return;
  if (p.failed) {
    note(out, id, "failed to parse WKB: " + p.error);
    return;
  }
  burn_one(p.g, id, full_grid, gs, out, opts, win);
}

// ---- Threads ----

// burn_par uses at most one thread per this much work (a geometry's work
// being its coordinates plus the grid rows it spans).
constexpr size_t kWorkPerThread = size_t(1) << 13;
// Chunks of geometries per thread.
constexpr size_t kChunksPerThread = 16;

static unsigned resolve_threads(unsigned threads) {
  if (threads != 0) return threads;
  unsigned n = std::thread::hardware_concurrency();
  return n == 0 ? 1 : n;
}

// A geometry's estimated burn cost: one, plus its coordinates, plus the
// grid rows its extent spans, the extent taken from at most 16
// coordinates of each part. As `work` in controlledburn-rs.
static size_t work(const Geometry& g, const GridSpec& grid) {
  size_t n = 0;
  double lo = std::numeric_limits<double>::infinity();
  double hi = -std::numeric_limits<double>::infinity();
  auto part = [&](const CoordSeq& c) {
    n += c.size();
    size_t step = (c.size() + 15) / 16 + (c.empty() ? 1 : 0);
    for (size_t i = 0; i < c.size(); i += step) {
      lo = std::min(lo, c[i].y);
      hi = std::max(hi, c[i].y);
    }
    if (!c.empty()) {
      lo = std::min(lo, c.back().y);
      hi = std::max(hi, c.back().y);
    }
  };
  // Points: one part of all the point coordinates, as Rust's points().
  CoordSeq pts;
  for (const auto& p : g.points) pts.insert(pts.end(), p.begin(), p.end());
  part(pts);
  for (const auto& l : g.lines) part(l);
  for (const auto& p : g.polygons) {
    if (!p.rings.empty()) part(p.rings[0]);
    for (size_t r = 1; r < p.rings.size(); r++) n += p.rings[r].size();
  }
  lo = std::max(lo, grid.ymin);
  hi = std::min(hi, grid.ymax);
  size_t rows = 0;
  if (hi > lo) {
    rows = static_cast<size_t>(std::min((hi - lo) / grid.dy(), static_cast<double>(grid.nrow))) + 1;
  }
  return 1 + n + rows;
}

// Burn items [0, n) on several threads: `cost(k)` estimates item k's
// work, `burn_k(k, out)` burns it. Contiguous chunks of about equal work,
// costliest first, taken by the threads in turn, then concatenated in
// order. Returns false (nothing done) when one thread would do.
template <typename Cost, typename BurnK>
static bool run_parallel(size_t n, unsigned threads_in, Cost cost_of, BurnK burn_k,
                         BurnResult& out) {
  size_t threads = resolve_threads(threads_in);
  std::vector<size_t> cost(n);
  size_t total = 0;
  for (size_t k = 0; k < n; k++) {
    cost[k] = cost_of(k);
    total += cost[k];
  }
  threads = std::min(threads, total / kWorkPerThread);
  if (threads < 2 || n < 2) return false;

  size_t target = (total + threads * kChunksPerThread - 1) / (threads * kChunksPerThread);
  std::vector<std::pair<size_t, size_t>> chunks;
  size_t start = 0, acc = 0;
  for (size_t k = 0; k < n; k++) {
    acc += cost[k];
    if (acc >= target) {
      chunks.push_back({start, k + 1});
      start = k + 1;
      acc = 0;
    }
  }
  if (start < n) chunks.push_back({start, n});
  std::vector<size_t> chunk_cost(chunks.size(), 0);
  for (size_t c = 0; c < chunks.size(); c++)
    for (size_t k = chunks[c].first; k < chunks[c].second; k++) chunk_cost[c] += cost[k];
  std::vector<size_t> order(chunks.size());
  for (size_t i = 0; i < order.size(); i++) order[i] = i;
  std::stable_sort(order.begin(), order.end(),
                   [&](size_t a, size_t b) { return chunk_cost[a] > chunk_cost[b]; });

  std::vector<BurnResult> parts(chunks.size());
  std::atomic<size_t> next(0);
  std::vector<std::exception_ptr> errors(std::min(threads, chunks.size()));
  std::vector<std::thread> pool;
  for (size_t t = 0; t < errors.size(); t++) {
    pool.emplace_back([&, t]() {
      try {
        for (;;) {
          size_t i = next.fetch_add(1, std::memory_order_relaxed);
          if (i >= order.size()) break;
          size_t c = order[i];
          for (size_t k = chunks[c].first; k < chunks[c].second; k++) burn_k(k, parts[c]);
        }
      } catch (...) {
        errors[t] = std::current_exception();
      }
    });
  }
  for (auto& th : pool) th.join();
  for (auto& e : errors)
    if (e) std::rethrow_exception(e);

  size_t nr = 0, ne = 0, nl = 0, np = 0;
  for (const auto& p : parts) {
    nr += p.runs.size(); ne += p.edges.size(); nl += p.lines.size(); np += p.points.size();
  }
  out.runs.reserve(nr); out.edges.reserve(ne); out.lines.reserve(nl); out.points.reserve(np);
  for (auto& p : parts) {
    out.runs.insert(out.runs.end(), p.runs.begin(), p.runs.end());
    out.edges.insert(out.edges.end(), p.edges.begin(), p.edges.end());
    out.lines.insert(out.lines.end(), p.lines.begin(), p.lines.end());
    out.points.insert(out.points.end(), p.points.begin(), p.points.end());
    for (auto& nt : p.notes) out.notes.push_back(std::move(nt));
    p = BurnResult();
  }
  return true;
}

static void check_window(const GridSpec& gs, uint32_t row0, uint32_t row1,
                         uint32_t col0, uint32_t col1) {
  if (row0 >= row1 || col0 >= col1 || row1 > static_cast<uint32_t>(gs.nrow) ||
      col1 > static_cast<uint32_t>(gs.ncol)) {
    throw std::invalid_argument(
        "window rows " + std::to_string(row0) + ".." + std::to_string(row1) + ", cols " +
        std::to_string(col0) + ".." + std::to_string(col1) +
        " must be non-empty and within " + std::to_string(gs.nrow) + " rows, " +
        std::to_string(gs.ncol) + " cols");
  }
}

} // anonymous namespace

// ---- Public entry points ----

BurnResult burn(const std::vector<Geometry>& geoms, const GridSpec& gs,
                const BurnOptions& opts) {
  validate(gs, opts, geoms.size());
  FullGrid full_grid = make_full_grid(gs);
  BurnResult out;
  for (size_t k = 0; k < geoms.size(); k++) {
    burn_one(geoms[k], geom_id(opts, k), full_grid, gs, out, opts, nullptr);
  }
  return out;
}

BurnResult burn_wkb(const std::vector<WKBSpan>& wkb, const GridSpec& gs,
                    const BurnOptions& opts) {
  gs.validate();
  if (opts.id_offset < 0) {
    throw std::invalid_argument("id_offset " + std::to_string(opts.id_offset) + " is negative");
  }
  FullGrid full_grid = make_full_grid(gs);
  BurnResult out;
  for (size_t k = 0; k < wkb.size(); k++) {
    if (!wkb[k].data || wkb[k].size == 0) continue;
    int32_t id = geom_id(opts, k);
    burn_parsed(parse_one(wkb[k]), id, full_grid, gs, out, opts, nullptr);
  }
  return out;
}

BurnResult burn_par(const std::vector<Geometry>& geoms, const GridSpec& gs,
                    const BurnOptions& opts, unsigned threads) {
  if (resolve_threads(threads) < 2 || geoms.size() < 2) return burn(geoms, gs, opts);
  validate(gs, opts, geoms.size());
  FullGrid full_grid = make_full_grid(gs);
  BurnResult out;
  bool done = run_parallel(
      geoms.size(), threads, [&](size_t k) { return work(geoms[k], gs); },
      [&](size_t k, BurnResult& o) {
        burn_one(geoms[k], geom_id(opts, k), full_grid, gs, o, opts, nullptr);
      },
      out);
  return done ? out : burn(geoms, gs, opts);
}

BurnResult burn_wkb_par(const std::vector<WKBSpan>& wkb, const GridSpec& gs,
                        const BurnOptions& opts, unsigned threads) {
  if (resolve_threads(threads) < 2 || wkb.size() < 2) return burn_wkb(wkb, gs, opts);
  gs.validate();
  if (opts.id_offset < 0) {
    throw std::invalid_argument("id_offset " + std::to_string(opts.id_offset) + " is negative");
  }
  std::vector<Parsed> parsed(wkb.size());
  for (size_t k = 0; k < wkb.size(); k++) parsed[k] = parse_one(wkb[k]);
  // burn_wkb checks each non-empty blob's id; the last such is enough.
  for (size_t k = wkb.size(); k-- > 0;) {
    if (!parsed[k].skip) { geom_id(opts, k); break; }
  }
  FullGrid full_grid = make_full_grid(gs);
  BurnResult out;
  bool done = run_parallel(
      wkb.size(), threads,
      [&](size_t k) { return parsed[k].skip || parsed[k].failed ? size_t(1) : work(parsed[k].g, gs); },
      [&](size_t k, BurnResult& o) {
        if (parsed[k].skip) return;
        burn_parsed(parsed[k], geom_id(opts, k), full_grid, gs, o, opts, nullptr);
      },
      out);
  if (done) return out;
  BurnResult serial;
  for (size_t k = 0; k < wkb.size(); k++) {
    if (parsed[k].skip) continue;
    burn_parsed(parsed[k], geom_id(opts, k), full_grid, gs, serial, opts, nullptr);
  }
  return serial;
}

BurnResult burn_window(const std::vector<Geometry>& geoms, const GridSpec& gs,
                       uint32_t row0, uint32_t row1, uint32_t col0, uint32_t col1,
                       const BurnOptions& opts) {
  gs.validate();
  check_window(gs, row0, row1, col0, col1);
  validate(gs, opts, geoms.size());
  Window win(gs, row0, row1, col0, col1);
  FullGrid full_grid = make_full_grid(gs);
  BurnResult out;
  for (size_t k = 0; k < geoms.size(); k++) {
    burn_one(geoms[k], geom_id(opts, k), full_grid, gs, out, opts, &win);
  }
  return out;
}

BurnResult burn_wkb_window(const std::vector<WKBSpan>& wkb, const GridSpec& gs,
                           uint32_t row0, uint32_t row1, uint32_t col0, uint32_t col1,
                           const BurnOptions& opts) {
  gs.validate();
  check_window(gs, row0, row1, col0, col1);
  if (opts.id_offset < 0) {
    throw std::invalid_argument("id_offset " + std::to_string(opts.id_offset) + " is negative");
  }
  Window win(gs, row0, row1, col0, col1);
  FullGrid full_grid = make_full_grid(gs);
  BurnResult out;
  for (size_t k = 0; k < wkb.size(); k++) {
    if (!wkb[k].data || wkb[k].size == 0) continue;
    int32_t id = geom_id(opts, k);
    burn_parsed(parse_one(wkb[k]), id, full_grid, gs, out, opts, &win);
  }
  return out;
}

} // namespace controlledburn
