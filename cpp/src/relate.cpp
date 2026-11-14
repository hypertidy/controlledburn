// relate.cpp -- DE-9IM of polygons against grid cells, after the burn.
//
// Port of relate_cells from the controlledburn Rust crate (0.3.1,
// src/relate.rs), step for step, so the matrices are the crate's.
//
// Method, per polygonal geometry:
//
// - Candidates. Each boundary segment is walked row by row; a cell is a
//   candidate when the segment meets its closed box, decided with exact
//   orientation tests. This includes cells the boundary only touches,
//   along a wall or at a corner, which have no record.
// - BI, BB. BI (boundary meets the open cell) is a strict separating axis
//   test. BB is 1 when a segment runs along a wall, 0 when the boundary
//   meets the outline only at points.
// - IB, EB. The contacts split the cell's outline into open arcs, each
//   wholly inside or outside the polygon. Going from one arc to the next
//   past a contact, inside or outside flips once for each boundary branch
//   that leaves the contact into the cell. One known arc per cell then
//   fixes the rest: from a neighbour the boundary does not touch (its
//   record is exact), or from a candidate neighbour already known through
//   the shared wall. A group of candidates with neither gets one exact
//   point-in-polygon test.
// - II, EI. Both 2 when the boundary crosses the open cell; otherwise the
//   open cell is wholly inside or outside, as its arcs are.
// - IE, BE. F only when the polygon lies inside the closed cell.
//
// Copyright (c) 2026 Michael Sumner
// Licensed under Apache License 2.0

#include "no_fma.h"  // first: no fused multiply-add (see header)
#include "controlledburn/relate.hpp"
#include "predicates.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace controlledburn {

namespace {

bool same(const Coord& a, const Coord& b) { return a.x == b.x && a.y == b.y; }

int orient(const Coord& a, const Coord& b, const Coord& c) {
  return detail::orient2d_sign(a.x, a.y, b.x, b.y, c.x, c.y);
}

struct CellBox {
  double x0, y0, x1, y1;
};

// The burn's cell walls.
struct Walls {
  double xmin, xmax, ymin, ymax, dx, dy;
  int64_t ncol, nrow;

  explicit Walls(const GridSpec& g)
      : xmin(g.xmin), xmax(g.xmax), ymin(g.ymin), ymax(g.ymax), dx(g.dx()), dy(g.dy()),
        ncol(g.ncol), nrow(g.nrow) {}

  // As Grid::cell: the last column and row are stretched to the extent.
  CellBox cell(int64_t row, int64_t col) const {
    CellBox b;
    b.x0 = xmin + static_cast<double>(col) * dx;
    b.x1 = col == ncol - 1 ? xmax : xmin + static_cast<double>(col + 1) * dx;
    b.y0 = row == nrow - 1 ? ymin : ymax - static_cast<double>(row + 1) * dy;
    b.y1 = ymax - static_cast<double>(row) * dy;
    return b;
  }

  // Loose index of the column holding x (may be off by one at walls).
  int64_t col_of(double x) const { return static_cast<int64_t>(std::floor((x - xmin) / dx)); }
  int64_t row_of(double y) const { return static_cast<int64_t>(std::floor((ymax - y) / dy)); }
};

// The four walls in counterclockwise order from the bottom-left corner.
// The outline is parameterised by s, its length from that corner.
enum Wall { Bottom = 0, Right = 1, Top = 2, Left = 3 };
constexpr Wall kWalls[4] = {Bottom, Right, Top, Left};

double bw(const CellBox& b) { return b.x1 - b.x0; }
double bh(const CellBox& b) { return b.y1 - b.y0; }

// s at the start of each wall, and the perimeter.
std::array<double, 5> marks(const CellBox& b) {
  double w = bw(b), h = bh(b);
  double c2 = w + h;
  double c3 = c2 + w;
  return {0.0, w, c2, c3, c3 + h};
}

// Corners in counterclockwise order from the bottom left; corner k starts
// wall k.
std::array<Coord, 4> corners(const CellBox& b) {
  return {Coord{b.x0, b.y0}, Coord{b.x1, b.y0}, Coord{b.x1, b.y1}, Coord{b.x0, b.y1}};
}

bool strictly_inside(const CellBox& b, const Coord& p) {
  return p.x > b.x0 && p.x < b.x1 && p.y > b.y0 && p.y < b.y1;
}

// s of a point on `wall`, given its coordinate along the wall.
double s_on(const CellBox& b, Wall wall, double along) {
  auto m = marks(b);
  switch (wall) {
  case Bottom: return along - b.x0;
  case Right: return m[1] + (along - b.y0);
  case Top: return m[2] + (b.x1 - along);
  case Left: return m[3] + (b.y1 - along);
  }
  return 0.0;
}

// The coordinate along `wall` of the point at s.
double along_at(const CellBox& b, Wall wall, double s) {
  auto m = marks(b);
  switch (wall) {
  case Bottom: return b.x0 + s;
  case Right: return b.y0 + (s - m[1]);
  case Top: return b.x1 - (s - m[2]);
  case Left: return b.y1 - (s - m[3]);
  }
  return 0.0;
}

double rem_euclid(double s, double per) {
  double r = std::fmod(s, per);
  return r < 0.0 ? r + std::fabs(per) : r;
}

Coord point_at(const CellBox& b, double s0) {
  auto m = marks(b);
  double s = rem_euclid(s0, m[4]);
  if (s < m[1]) return Coord{along_at(b, Bottom, s), b.y0};
  if (s < m[2]) return Coord{b.x1, along_at(b, Right, s)};
  if (s < m[3]) return Coord{along_at(b, Top, s), b.y1};
  return Coord{b.x0, along_at(b, Left, s)};
}

// Whether a branch from a point on the open `wall` towards o goes into
// the open cell.
bool inward(const CellBox& b, Wall wall, const Coord& o) {
  switch (wall) {
  case Bottom: return o.y > b.y0;
  case Right: return o.x < b.x1;
  case Top: return o.y < b.y1;
  case Left: return o.x > b.x0;
  }
  return false;
}

// Whether a branch from corner k towards o goes into the open cell: from
// a corner, the open cell is the open quadrant on its side.
bool inward_corner(const CellBox& b, size_t k, const Coord& o) {
  bool x_ok = (k == 0 || k == 3) ? o.x > b.x0 : o.x < b.x1;
  bool y_ok = (k == 0 || k == 1) ? o.y > b.y0 : o.y < b.y1;
  return x_ok && y_ok;
}

// Segment ab meets the open box: strict separating axes, with the box's
// corners on the segment's line counting as separated.
bool meets_open(const Coord& a, const Coord& b, const CellBox& bx) {
  if (std::max(a.x, b.x) <= bx.x0 || std::min(a.x, b.x) >= bx.x1 ||
      std::max(a.y, b.y) <= bx.y0 || std::min(a.y, b.y) >= bx.y1) {
    return false;
  }
  auto cs = corners(bx);
  bool all_ge = true, all_le = true;
  for (const auto& c : cs) {
    int v = orient(a, b, c);
    all_ge &= v >= 0;
    all_le &= v <= 0;
  }
  return !(all_ge || all_le);
}

// Segment ab meets the closed box.
bool meets_closed(const Coord& a, const Coord& b, const CellBox& bx) {
  if (std::max(a.x, b.x) < bx.x0 || std::min(a.x, b.x) > bx.x1 ||
      std::max(a.y, b.y) < bx.y0 || std::min(a.y, b.y) > bx.y1) {
    return false;
  }
  auto cs = corners(bx);
  bool all_gt = true, all_lt = true;
  for (const auto& c : cs) {
    int v = orient(a, b, c);
    all_gt &= v > 0;
    all_lt &= v < 0;
  }
  return !(all_gt || all_lt);
}

// A polygonal geometry as one set of boundary segments.
struct Shape {
  std::vector<std::pair<Coord, Coord>> segs;
  double bbox[4];

  // False (no shape) for a non-finite coordinate or no segments.
  bool build(const std::vector<const Polygon*>& polys) {
    const double inf = std::numeric_limits<double>::infinity();
    bbox[0] = inf; bbox[1] = inf; bbox[2] = -inf; bbox[3] = -inf;
    for (const Polygon* p : polys) {
      for (const auto& ring : p->rings) {
        if (ring.size() < 3) continue;
        for (const auto& c : ring) {
          if (!std::isfinite(c.x) || !std::isfinite(c.y)) return false;
          bbox[0] = std::min(bbox[0], c.x);
          bbox[1] = std::min(bbox[1], c.y);
          bbox[2] = std::max(bbox[2], c.x);
          bbox[3] = std::max(bbox[3], c.y);
        }
        bool closed = same(ring.front(), ring.back());
        size_t n = ring.size();
        for (size_t i = 0; i < n - (closed ? 1 : 0); i++) {
          const Coord& a = ring[i];
          const Coord& b = ring[(i + 1) % n];
          if (!same(a, b)) segs.push_back({a, b});
        }
      }
    }
    return !segs.empty();
  }

  // Exact point location by crossing number: true inside. p must not lie
  // on the boundary.
  bool contains(const Coord& p) const {
    bool inside = false;
    for (const auto& s : segs) {
      const Coord& a = s.first;
      const Coord& b = s.second;
      if ((a.y > p.y) != (b.y > p.y)) {
        int o = orient(a, b, p);
        if ((b.y > a.y && o > 0) || (b.y < a.y && o < 0)) inside = !inside;
      }
    }
    return inside;
  }
};

// Runs by (id, row), merged, for "is this cell in a run of this id".
struct RunIndex {
  std::vector<std::array<int32_t, 4>> runs;  // id, row, col_start, col_end

  explicit RunIndex(const BurnResult& r) {
    std::vector<std::array<int32_t, 4>> v;
    v.reserve(r.runs.size());
    for (const auto& x : r.runs) v.push_back({x.id, x.row, x.col_start, x.col_end});
    std::sort(v.begin(), v.end());
    runs.reserve(v.size());
    for (const auto& x : v) {
      if (!runs.empty()) {
        auto& l = runs.back();
        if (l[0] == x[0] && l[1] == x[1] && x[2] <= l[3]) {
          l[3] = std::max(l[3], x[3]);
          continue;
        }
      }
      runs.push_back(x);
    }
  }

  bool covered(int32_t id, int32_t row, int32_t col) const {
    // first run with (id, row, col_start) > (id, row, col)
    auto it = std::partition_point(runs.begin(), runs.end(), [&](const std::array<int32_t, 4>& x) {
      return std::tie(x[0], x[1], x[2]) <= std::tie(id, row, col);
    });
    if (it == runs.begin()) return false;
    const auto& x = *(it - 1);
    return x[0] == id && x[1] == row && col < x[3];
  }
};

// A contact cluster event on the outline: a span lo < hi along a wall, or
// a point lo == hi, with the number of branches leaving it into the open
// cell.
struct Event {
  double lo, hi;
  uint32_t inward;
};

// One candidate cell's outline: its open arcs as (lo, hi) in s, with hi
// past the perimeter for an arc across the bottom-left corner, and each
// arc's inside/outside relative to arc 0.
struct Outline {
  CellBox bx;
  std::vector<std::pair<double, double>> arcs;
  std::vector<bool> parity;
  int base = -1;  // -1 unknown, else 0 or 1

  bool known() const { return base >= 0; }

  bool status(size_t j) const { return known() && ((base != 0) ^ parity[j]); }

  // index of the arc holding s, or -1
  long arc_at(double s) const {
    double per = marks(bx)[4];
    for (size_t j = 0; j < arcs.size(); j++) {
      double lo = arcs[j].first, hi = arcs[j].second;
      if ((lo < s && s < hi) || (lo < s + per && s + per < hi)) return static_cast<long>(j);
    }
    return -1;
  }

  // A point strictly inside an arc on `wall`, as the coordinate along the
  // wall, with the arc's index (-1: none).
  long on_wall(Wall wall, double& along) const {
    auto m = marks(bx);
    size_t k = static_cast<size_t>(wall);
    double ws = m[k], we = m[k + 1];
    for (size_t j = 0; j < arcs.size(); j++) {
      double lo = arcs[j].first, hi = arcs[j].second;
      for (double shift : {0.0, m[4]}) {
        double a = std::max(lo, ws + shift) - shift;
        double b = std::min(hi, we + shift) - shift;
        if (a < b) {
          along = along_at(bx, wall, (a + b) / 2.0);
          return static_cast<long>(j);
        }
      }
    }
    return -1;
  }
};

struct CellFacts {
  bool bi;
  char bb;
  Outline outline;
};

// The open arcs between contact clusters, and each arc's inside/outside
// relative to the first: crossing a cluster flips it once per branch the
// cluster sends into the cell.
void arcs_of(std::vector<Event> ev, double per, std::vector<std::pair<double, double>>& arcs,
             std::vector<bool>& parity) {
  arcs.clear();
  parity.clear();
  if (ev.empty()) {
    arcs.push_back({0.0, per});
    parity.push_back(false);
    return;
  }
  std::stable_sort(ev.begin(), ev.end(), [](const Event& a, const Event& b) { return a.lo < b.lo; });
  std::vector<Event> cl;
  for (const auto& e : ev) {
    if (!cl.empty() && e.lo <= cl.back().hi) {
      cl.back().hi = std::max(cl.back().hi, e.hi);
      cl.back().inward += e.inward;
    } else {
      cl.push_back(e);
    }
  }
  // A cluster touching the end of the outline joins the one at its start.
  if (cl.size() > 1 && cl.back().hi >= per && cl[0].lo <= 0.0) {
    Event last = cl.back();
    cl.pop_back();
    cl[0].lo = last.lo - per;
    cl[0].inward += last.inward;
  }
  if (cl.size() == 1 && cl[0].hi - cl[0].lo >= per) return;
  size_t n = cl.size();
  bool p = false;
  for (size_t i = 0; i < n; i++) {
    double lo = cl[i].hi;
    double hi = i + 1 < n ? cl[i + 1].lo : cl[0].lo + per;
    if (i > 0) p ^= (cl[i].inward % 2 == 1);
    if (lo < hi) {
      arcs.push_back({lo, hi});
      parity.push_back(p);
    }
  }
}

// BI, BB and the outline's arcs for one candidate cell, from the segments
// meeting its closed box.
CellFacts cell_facts(const Shape& shape, const CellBox& bx, const uint32_t* seg_begin,
                     const uint32_t* seg_end) {
  bool bi = false;
  bool line = false, point = false;
  std::vector<Event> ev;
  auto cs = corners(bx);
  auto m = marks(bx);
  for (const uint32_t* sp = seg_begin; sp != seg_end; ++sp) {
    const Coord a = shape.segs[*sp].first;
    const Coord b = shape.segs[*sp].second;
    bi |= meets_open(a, b, bx);
    if (strictly_inside(bx, a) && strictly_inside(bx, b)) continue;
    // Corners: the segment passes through or ends at one.
    for (size_t k = 0; k < 4; k++) {
      const Coord& c = cs[k];
      if (c.x < std::min(a.x, b.x) || c.x > std::max(a.x, b.x) ||
          c.y < std::min(a.y, b.y) || c.y > std::max(a.y, b.y)) {
        continue;
      }
      if (orient(a, b, c) != 0) continue;
      point = true;
      // the boundary leaves the corner towards each endpoint not on it
      uint32_t in = 0;
      for (const Coord* o : {&a, &b}) {
        if (!same(*o, c) && inward_corner(bx, k, *o)) in++;
      }
      ev.push_back({m[k], m[k], in});
    }
    // Along a wall: a span, no branch into the cell.
    bool along_wall = false;
    if (a.y == b.y && (a.y == bx.y0 || a.y == bx.y1)) {
      double lo = std::max(std::min(a.x, b.x), bx.x0), hi = std::min(std::max(a.x, b.x), bx.x1);
      if (lo < hi) {
        Wall w = a.y == bx.y0 ? Bottom : Top;
        double s0 = s_on(bx, w, lo), s1 = s_on(bx, w, hi);
        ev.push_back({std::min(s0, s1), std::max(s0, s1), 0});
        along_wall = true;
      }
    }
    if (a.x == b.x && (a.x == bx.x0 || a.x == bx.x1)) {
      double lo = std::max(std::min(a.y, b.y), bx.y0), hi = std::min(std::max(a.y, b.y), bx.y1);
      if (lo < hi) {
        Wall w = a.x == bx.x1 ? Right : Left;
        double s0 = s_on(bx, w, lo), s1 = s_on(bx, w, hi);
        ev.push_back({std::min(s0, s1), std::max(s0, s1), 0});
        along_wall = true;
      }
    }
    if (along_wall) {
      line = true;
      continue;
    }
    // Open walls: an endpoint on one, or a crossing between its corners.
    for (Wall w : kWalls) {
      bool vertical;
      double at, lo, hi;
      Coord c_lo, c_hi;
      switch (w) {
      case Bottom: vertical = false; at = bx.y0; lo = bx.x0; hi = bx.x1; c_lo = cs[0]; c_hi = cs[1]; break;
      case Right: vertical = true; at = bx.x1; lo = bx.y0; hi = bx.y1; c_lo = cs[1]; c_hi = cs[2]; break;
      case Top: vertical = false; at = bx.y1; lo = bx.x0; hi = bx.x1; c_lo = cs[3]; c_hi = cs[2]; break;
      default: vertical = true; at = bx.x0; lo = bx.y0; hi = bx.y1; c_lo = cs[0]; c_hi = cs[3]; break;
      }
      auto across = [&](const Coord& p) { return vertical ? p.x : p.y; };
      auto along = [&](const Coord& p) { return vertical ? p.y : p.x; };
      const std::pair<Coord, Coord> ends[2] = {{a, b}, {b, a}};
      for (const auto& eo : ends) {
        const Coord& e = eo.first;
        const Coord& o = eo.second;
        if (across(e) == at && along(e) > lo && along(e) < hi) {
          point = true;
          double s = s_on(bx, w, along(e));
          ev.push_back({s, s, inward(bx, w, o) ? 1u : 0u});
        }
      }
      double da = across(a) - at, db = across(b) - at;
      if ((da < 0.0 && db > 0.0) || (da > 0.0 && db < 0.0)) {
        int o1 = orient(a, b, c_lo), o2 = orient(a, b, c_hi);
        if (o1 != 0 && o2 != 0 && o1 != o2) {
          point = true;
          double t = (at - across(a)) / (across(b) - across(a));
          double pos = along(a) + t * (along(b) - along(a));
          pos = std::min(std::max(pos, lo), hi);
          double s = s_on(bx, w, pos);
          ev.push_back({s, s, 1});
        }
      }
    }
  }
  CellFacts f;
  f.bi = bi;
  f.bb = line ? '1' : (point ? '0' : 'F');
  f.outline.bx = bx;
  arcs_of(std::move(ev), m[4], f.outline.arcs, f.outline.parity);
  return f;
}

void relate_one(const Shape& shape, const Walls& walls, const RunIndex& runs, int32_t id,
                std::vector<CellRelation>& out) {
  // candidate (row, col, segment) triples
  struct Hit {
    int32_t row, col;
    uint32_t seg;
    bool operator<(const Hit& o) const { return std::tie(row, col, seg) < std::tie(o.row, o.col, o.seg); }
    bool operator==(const Hit& o) const { return row == o.row && col == o.col && seg == o.seg; }
  };
  std::vector<Hit> hits;
  for (size_t si = 0; si < shape.segs.size(); si++) {
    const Coord a = shape.segs[si].first;
    const Coord b = shape.segs[si].second;
    double miny = std::min(a.y, b.y), maxy = std::max(a.y, b.y);
    int64_t r_lo = std::max<int64_t>(walls.row_of(maxy) - 1, 0);
    int64_t r_hi = std::min<int64_t>(walls.row_of(miny) + 1, walls.nrow - 1);
    for (int64_t row = r_lo; row <= r_hi; row++) {
      CellBox probe = walls.cell(row, 0);
      if (maxy < probe.y0 || miny > probe.y1) continue;
      double xa, xb;
      if (a.y == b.y) {
        xa = std::min(a.x, b.x);
        xb = std::max(a.x, b.x);
      } else {
        auto x_at = [&](double y) { return a.x + (y - a.y) * (b.x - a.x) / (b.y - a.y); };
        double p = x_at(std::max(probe.y0, miny)), q = x_at(std::min(probe.y1, maxy));
        xa = std::min(p, q);
        xb = std::max(p, q);
      }
      int64_t c_lo = std::max<int64_t>(walls.col_of(xa) - 1, 0);
      int64_t c_hi = std::min<int64_t>(walls.col_of(xb) + 1, walls.ncol - 1);
      for (int64_t col = c_lo; col <= c_hi; col++) {
        if (meets_closed(a, b, walls.cell(row, col))) {
          hits.push_back({static_cast<int32_t>(row), static_cast<int32_t>(col),
                          static_cast<uint32_t>(si)});
        }
      }
    }
  }
  if (hits.empty()) return;
  std::sort(hits.begin(), hits.end());
  hits.erase(std::unique(hits.begin(), hits.end()), hits.end());

  std::vector<std::pair<int32_t, int32_t>> cells;
  std::vector<CellFacts> facts;
  std::vector<uint32_t> segs;
  size_t i = 0;
  while (i < hits.size()) {
    int32_t row = hits[i].row, col = hits[i].col;
    size_t j = i;
    segs.clear();
    while (j < hits.size() && hits[j].row == row && hits[j].col == col) {
      segs.push_back(hits[j].seg);
      j++;
    }
    CellBox bx = walls.cell(row, col);
    facts.push_back(cell_facts(shape, bx, segs.data(), segs.data() + segs.size()));
    cells.push_back({row, col});
    i = j;
  }

  // Inside/outside of each cell's arcs: seed from neighbours the boundary
  // does not touch, then across shared walls.
  auto find = [&](int32_t r, int32_t c) -> long {
    auto it = std::lower_bound(cells.begin(), cells.end(), std::make_pair(r, c));
    if (it != cells.end() && it->first == r && it->second == c) return it - cells.begin();
    return -1;
  };
  auto in_grid = [&](int32_t r, int32_t c) {
    return r >= 0 && c >= 0 && static_cast<int64_t>(r) < walls.nrow &&
           static_cast<int64_t>(c) < walls.ncol;
  };
  // neighbour across each wall, and the wall it sees
  auto across = [](Wall w, int32_t& dr, int32_t& dc) -> Wall {
    switch (w) {
    case Bottom: dr = 1; dc = 0; return Top;
    case Right: dr = 0; dc = 1; return Left;
    case Top: dr = -1; dc = 0; return Bottom;
    default: dr = 0; dc = -1; return Right;
    }
  };
  std::vector<size_t> queue;
  for (size_t k = 0; k < cells.size(); k++) {
    int32_t row = cells[k].first, col = cells[k].second;
    Outline& o = facts[k].outline;
    if (o.arcs.empty()) continue;
    for (Wall w : kWalls) {
      int32_t dr, dc;
      across(w, dr, dc);
      int32_t r2 = row + dr, c2 = col + dc;
      if (!in_grid(r2, c2) || find(r2, c2) >= 0) continue;
      // No boundary on this wall: one arc covers its middle.
      auto m = marks(o.bx);
      size_t wi = static_cast<size_t>(w);
      long j = o.arc_at((m[wi] + m[wi + 1]) / 2.0);
      if (j >= 0) {
        o.base = (runs.covered(id, r2, c2) ^ o.parity[static_cast<size_t>(j)]) ? 1 : 0;
        queue.push_back(k);
        break;
      }
    }
  }
  size_t next = 0;
  for (;;) {
    while (next < queue.size()) {
      size_t k = queue[next++];
      int32_t row = cells[k].first, col = cells[k].second;
      for (Wall w : kWalls) {
        int32_t dr, dc;
        Wall w2 = across(w, dr, dc);
        long k2 = find(row + dr, col + dc);
        if (k2 < 0) continue;
        Outline& o2 = facts[static_cast<size_t>(k2)].outline;
        if (o2.known() || o2.arcs.empty()) continue;
        double along = 0.0;
        long j = facts[k].outline.on_wall(w, along);
        if (j < 0) continue;
        bool st = facts[k].outline.status(static_cast<size_t>(j));
        long j2 = o2.arc_at(s_on(o2.bx, w2, along));
        if (j2 >= 0) {
          o2.base = (st ^ o2.parity[static_cast<size_t>(j2)]) ? 1 : 0;
          queue.push_back(static_cast<size_t>(k2));
        }
      }
    }
    // A group no neighbour reaches: one exact point-in-polygon test.
    long k = -1;
    for (size_t q = 0; q < facts.size(); q++) {
      if (!facts[q].outline.known() && !facts[q].outline.arcs.empty()) {
        k = static_cast<long>(q);
        break;
      }
    }
    if (k < 0) break;
    Outline& o = facts[static_cast<size_t>(k)].outline;
    double lo = o.arcs[0].first, hi = o.arcs[0].second;
    o.base = shape.contains(point_at(o.bx, (lo + hi) / 2.0)) ? 1 : 0;
    queue.push_back(static_cast<size_t>(k));
  }

  for (size_t k = 0; k < cells.size(); k++) {
    int32_t row = cells[k].first, col = cells[k].second;
    const CellFacts& f = facts[k];
    const Outline& o = f.outline;
    bool ib = false, eb = false;
    for (size_t j = 0; j < o.arcs.size(); j++) {
      if (o.status(j)) ib = true; else eb = true;
    }
    bool ii, ei;
    if (f.bi) {
      ii = true;
      ei = true;
    } else if (o.arcs.empty()) {
      // The whole outline is boundary (the polygon or a hole is this
      // cell): its coverage is exactly 0 or 1, so the record is exact.
      bool c = runs.covered(id, row, col);
      ii = c;
      ei = !c;
    } else {
      ii = ib;
      ei = eb;
    }
    const CellBox& bx = o.bx;
    bool within = shape.bbox[0] >= bx.x0 && shape.bbox[2] <= bx.x1 &&
                  shape.bbox[1] >= bx.y0 && shape.bbox[3] <= bx.y1;
    auto d = [](bool on, char dim) { return on ? dim : 'F'; };
    IntersectionMatrix mx;
    mx.m = {d(ii, '2'), d(ib, '1'), d(!within, '2'),
            d(f.bi, '1'), f.bb, d(!within, '1'),
            d(ei, '2'), d(eb, '1'), '2'};
    out.push_back({row, col, id, mx});
  }
}

} // namespace

std::vector<CellRelation> relate_cells(const BurnResult& result,
                                       const std::vector<Geometry>& geoms,
                                       const GridSpec& grid,
                                       const BurnOptions& opts) {
  grid.validate();
  if (opts.id_offset < 0) {
    throw std::invalid_argument("id_offset " + std::to_string(opts.id_offset) + " is negative");
  }
  auto id_of = [&](size_t k) -> int32_t {
    if (static_cast<int64_t>(opts.id_offset) + static_cast<int64_t>(k) > INT32_MAX) {
      throw std::invalid_argument("geometry " + std::to_string(k) + " with id_offset " +
                                  std::to_string(opts.id_offset) + ": id beyond i32::MAX");
    }
    return opts.id_offset + static_cast<int32_t>(k);
  };
  if (!geoms.empty()) id_of(geoms.size() - 1);
  RunIndex runs(result);
  Walls walls(grid);
  std::vector<CellRelation> out;
  for (size_t k = 0; k < geoms.size(); k++) {
    int32_t id = id_of(k);
    const Geometry& g = geoms[k];
    if (g.kind != GeomKind::Polygon && g.kind != GeomKind::MultiPolygon) continue;
    std::vector<const Polygon*> polys;
    for (const auto& p : g.polygons) polys.push_back(&p);
    Shape shape;
    if (!shape.build(polys)) continue;
    relate_one(shape, walls, runs, id, out);
  }
  return out;
}

} // namespace controlledburn
