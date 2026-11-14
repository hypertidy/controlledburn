// ops.cpp -- operations on a finished burn: coalesce_runs, merge_cells,
// crop, slice_tiles, slice_tiles_par. Ports of the controlledburn Rust
// crate's (0.3.1), with the same records in the same order.
//
// Copyright (c) 2026 Michael Sumner
// Licensed under Apache License 2.0

#include "no_fma.h"  // first: no fused multiply-add (see header)
#include "controlledburn/ops.hpp"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <exception>
#include <stdexcept>
#include <string>
#include <thread>
#include <tuple>

namespace controlledburn {

namespace {

constexpr float kTol = 1e-6f;

void coalesce_sorted(std::vector<GridRun>& runs) {
  size_t w = 0;
  for (size_t i = 0; i < runs.size(); i++) {
    GridRun cur = runs[i];
    if (w > 0) {
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

} // namespace

void coalesce_runs(std::vector<GridRun>& runs) {
  std::stable_sort(runs.begin(), runs.end(), [](const GridRun& a, const GridRun& b) {
    return std::tie(a.id, a.row, a.col_start) < std::tie(b.id, b.row, b.col_start);
  });
  coalesce_sorted(runs);
}

void merge_cells(BurnResult& result) {
  // stable: a cell's fractions are summed in table order
  auto& edges = result.edges;
  std::stable_sort(edges.begin(), edges.end(), [](const GridEdge& a, const GridEdge& b) {
    return std::tie(a.id, a.row, a.col) < std::tie(b.id, b.row, b.col);
  });
  size_t w = 0;
  for (size_t i = 0; i < edges.size(); i++) {
    GridEdge cur = edges[i];
    if (w > 0) {
      GridEdge& last = edges[w - 1];
      if (last.id == cur.id && last.row == cur.row && last.col == cur.col) {
        last.fraction += cur.fraction;
        continue;
      }
    }
    edges[w++] = cur;
  }
  edges.resize(w);

  w = 0;
  for (size_t i = 0; i < edges.size(); i++) {
    GridEdge e = edges[i];
    if (e.fraction >= 1.0f - kTol) {
      result.runs.push_back({e.row, e.col, e.col + 1, e.id});
    } else {
      edges[w++] = e;
    }
  }
  edges.resize(w);
  coalesce_runs(result.runs);

  // both sorted by (id, row, col): drop the edges a run covers
  const auto& runs = result.runs;
  size_t j = 0;
  w = 0;
  for (size_t i = 0; i < edges.size(); i++) {
    GridEdge e = edges[i];
    while (j < runs.size() &&
           std::make_tuple(runs[j].id, runs[j].row, runs[j].col_end) <=
               std::make_tuple(e.id, e.row, e.col)) {
      j++;
    }
    bool covered = j < runs.size() && runs[j].id == e.id && runs[j].row == e.row &&
                   runs[j].col_start <= e.col;
    if (!covered) edges[w++] = e;
  }
  edges.resize(w);

  auto& lines = result.lines;
  std::stable_sort(lines.begin(), lines.end(), [](const GridLine& a, const GridLine& b) {
    return std::tie(a.id, a.row, a.col) < std::tie(b.id, b.row, b.col);
  });
  w = 0;
  for (size_t i = 0; i < lines.size(); i++) {
    GridLine cur = lines[i];
    if (w > 0) {
      GridLine& last = lines[w - 1];
      if (last.id == cur.id && last.row == cur.row && last.col == cur.col) {
        last.length += cur.length;
        continue;
      }
    }
    lines[w++] = cur;
  }
  lines.resize(w);
}

bool crop(const BurnResult& result, const GridSpec& grid,
          double txmin, double tymin, double txmax, double tymax,
          BurnResult& out, GridSpec& subgrid) {
  const double dx = grid.dx();
  const double dy = grid.dy();
  const int64_t ncol = grid.ncol;
  const int64_t nrow = grid.nrow;

  // 0-based half-open limits, snapped outward. The eps nudge avoids
  // floor/ceil flips when the target aligns exactly with a cell boundary.
  const double eps = 1e-8;
  int64_t col_lo = std::max<int64_t>(static_cast<int64_t>(std::floor((txmin - grid.xmin) / dx + eps)), 0);
  int64_t col_hi = std::min<int64_t>(static_cast<int64_t>(std::ceil((txmax - grid.xmin) / dx - eps)), ncol);
  int64_t row_lo = std::max<int64_t>(static_cast<int64_t>(std::floor((grid.ymax - tymax) / dy + eps)), 0);
  int64_t row_hi = std::min<int64_t>(static_cast<int64_t>(std::ceil((grid.ymax - tymin) / dy - eps)), nrow);
  if (col_lo >= col_hi || row_lo >= row_hi) return false;

  const int32_t c0 = static_cast<int32_t>(col_lo), c1 = static_cast<int32_t>(col_hi);
  const int32_t r0 = static_cast<int32_t>(row_lo), r1 = static_cast<int32_t>(row_hi);

  subgrid = GridSpec{grid.xmin + c0 * dx, grid.ymax - r1 * dy,
                     grid.xmin + c1 * dx, grid.ymax - r0 * dy,
                     c1 - c0, r1 - r0};

  auto in_rows = [&](int32_t r) { return r >= r0 && r < r1; };
  auto in_cols = [&](int32_t c) { return c >= c0 && c < c1; };

  BurnResult o;
  for (const auto& r : result.runs) {
    if (in_rows(r.row) && r.col_end > c0 && r.col_start < c1) {
      o.runs.push_back({r.row - r0, std::max(r.col_start, c0) - c0,
                        std::min(r.col_end, c1) - c0, r.id});
    }
  }
  for (const auto& e : result.edges)
    if (in_rows(e.row) && in_cols(e.col))
      o.edges.push_back({e.row - r0, e.col - c0, e.fraction, e.id});
  for (const auto& l : result.lines)
    if (in_rows(l.row) && in_cols(l.col))
      o.lines.push_back({l.row - r0, l.col - c0, l.length, l.id});
  for (const auto& p : result.points)
    if (in_rows(p.row) && in_cols(p.col))
      o.points.push_back({p.row - r0, p.col - c0, p.id});
  o.notes = result.notes;
  out = std::move(o);
  return true;
}

namespace {

// The tile grid over a grid.
struct Shape {
  const GridSpec& grid;
  uint32_t tile_ncol, tile_nrow;
  size_t ntc, ntr;

  Shape(const GridSpec& g, uint32_t tc, uint32_t tr) : grid(g), tile_ncol(tc), tile_nrow(tr) {
    if (tc == 0 || tr == 0) throw std::invalid_argument("tile shape must be at least 1 x 1");
    ntc = (static_cast<size_t>(g.ncol) + tc - 1) / tc;
    ntr = (static_cast<size_t>(g.nrow) + tr - 1) / tr;
  }

  void check(int32_t row, int32_t col) const {
    if (row < 0 || col < 0 || row >= grid.nrow || col >= grid.ncol) {
      throw std::invalid_argument("record (" + std::to_string(row) + ", " + std::to_string(col) +
                                  ") is outside the " + std::to_string(grid.ncol) + " x " +
                                  std::to_string(grid.nrow) + " grid");
    }
  }

  // The tiles in tile rows [tr0, tr1), from every record (sel null) or the
  // selected ones (indices into each table), all in those tile rows.
  std::vector<Tile> slice(const BurnResult& result, size_t tr0, size_t tr1,
                          const std::vector<uint32_t>* sel) const {
    const int32_t tc = static_cast<int32_t>(tile_ncol), tr = static_cast<int32_t>(tile_nrow);
    auto at = [&](int32_t row, int32_t col) -> size_t {
      check(row, col);
      return (static_cast<size_t>(row / tr) - tr0) * ntc + static_cast<size_t>(col / tc);
    };
    std::vector<BurnResult> out(ntc * (tr1 - tr0));
    auto each = [&](size_t n, const std::vector<uint32_t>* s, auto f) {
      if (!s) {
        for (size_t i = 0; i < n; i++) f(i);
      } else {
        for (uint32_t i : *s) f(i);
      }
    };
    each(result.runs.size(), sel ? &sel[0] : nullptr, [&](size_t i) {
      const GridRun& run = result.runs[i];
      int32_t lr = run.row % tr;
      int32_t c = run.col_start;
      while (c < run.col_end) {
        int32_t c0 = (c / tc) * tc;
        int32_t ce = std::min(c0 + tc, run.col_end);
        out[at(run.row, c)].runs.push_back({lr, c - c0, ce - c0, run.id});
        c = ce;
      }
    });
    each(result.edges.size(), sel ? &sel[1] : nullptr, [&](size_t i) {
      const GridEdge& e = result.edges[i];
      out[at(e.row, e.col)].edges.push_back({e.row % tr, e.col % tc, e.fraction, e.id});
    });
    each(result.lines.size(), sel ? &sel[2] : nullptr, [&](size_t i) {
      const GridLine& l = result.lines[i];
      out[at(l.row, l.col)].lines.push_back({l.row % tr, l.col % tc, l.length, l.id});
    });
    each(result.points.size(), sel ? &sel[3] : nullptr, [&](size_t i) {
      const GridPoint& p = result.points[i];
      out[at(p.row, p.col)].points.push_back({p.row % tr, p.col % tc, p.id});
    });

    const double dx = grid.dx(), dy = grid.dy();
    std::vector<Tile> tiles;
    for (size_t k = 0; k < out.size(); k++) {
      BurnResult& r = out[k];
      if (r.runs.empty() && r.edges.empty() && r.lines.empty() && r.points.empty()) continue;
      uint32_t ti = static_cast<uint32_t>(tr0 + k / ntc), tj = static_cast<uint32_t>(k % ntc);
      uint32_t row0 = ti * tile_nrow, col0 = tj * tile_ncol;
      uint32_t row1 = std::min(row0 + tile_nrow, static_cast<uint32_t>(grid.nrow));
      uint32_t col1 = std::min(col0 + tile_ncol, static_cast<uint32_t>(grid.ncol));
      // Outer edges come from the grid itself so the last tiles end
      // exactly at xmax and ymin.
      double xmax = col1 == static_cast<uint32_t>(grid.ncol) ? grid.xmax : grid.xmin + col1 * dx;
      double ymin = row1 == static_cast<uint32_t>(grid.nrow) ? grid.ymin : grid.ymax - row1 * dy;
      GridSpec tg{grid.xmin + col0 * dx, ymin, xmax, grid.ymax - row0 * dy,
                  static_cast<int32_t>(col1 - col0), static_cast<int32_t>(row1 - row0)};
      tiles.push_back(Tile{ti, tj, row0, col0, tg, std::move(r)});
    }
    return tiles;
  }
};

// slice_tiles_par uses at most one thread per this many records.
constexpr size_t kRecordsPerThread = size_t(1) << 17;
constexpr size_t kBandsPerThread = 8;

} // namespace

std::vector<Tile> slice_tiles(const BurnResult& result, const GridSpec& grid,
                              uint32_t tile_ncol, uint32_t tile_nrow) {
  Shape s(grid, tile_ncol, tile_nrow);
  return s.slice(result, 0, s.ntr, nullptr);
}

std::vector<Tile> slice_tiles_par(const BurnResult& result, const GridSpec& grid,
                                  uint32_t tile_ncol, uint32_t tile_nrow,
                                  unsigned threads_in) {
  Shape s(grid, tile_ncol, tile_nrow);
  size_t threads = threads_in;
  if (threads == 0) {
    unsigned n = std::thread::hardware_concurrency();
    threads = n == 0 ? 1 : n;
  }
  const size_t lens[4] = {result.runs.size(), result.edges.size(), result.lines.size(),
                          result.points.size()};
  size_t total = lens[0] + lens[1] + lens[2] + lens[3];
  threads = std::min(threads, total / kRecordsPerThread);
  bool too_big = false;
  for (size_t n : lens) too_big |= n > UINT32_MAX;
  if (threads < 2 || s.ntr < 2 || too_big) return s.slice(result, 0, s.ntr, nullptr);

  // Band b holds tile rows [band_start(b), band_start(b + 1)). Records are
  // assigned to bands in one serial pass (keeping burn order within each
  // band), then the bands are sliced in parallel. Same tiles, records and
  // order as slice_tiles.
  const size_t ntr = s.ntr;
  const size_t nbands = std::min(std::min(threads * kBandsPerThread, ntr), size_t(65535));
  auto band_start = [&](size_t b) { return (b * ntr + nbands - 1) / nbands; };
  std::vector<size_t> band_of_tile_row(ntr);
  for (size_t t = 0; t < ntr; t++) band_of_tile_row[t] = t * nbands / ntr;
  // band of each record, checked against the grid
  std::vector<std::vector<std::vector<uint32_t>>> sel(nbands,
                                                     std::vector<std::vector<uint32_t>>(4));
  auto place = [&](int32_t row, int32_t col) {
    s.check(row, col);
    return band_of_tile_row[static_cast<size_t>(row) / tile_nrow];
  };
  for (size_t i = 0; i < lens[0]; i++)
    sel[place(result.runs[i].row, result.runs[i].col_start)][0].push_back(static_cast<uint32_t>(i));
  for (size_t i = 0; i < lens[1]; i++)
    sel[place(result.edges[i].row, result.edges[i].col)][1].push_back(static_cast<uint32_t>(i));
  for (size_t i = 0; i < lens[2]; i++)
    sel[place(result.lines[i].row, result.lines[i].col)][2].push_back(static_cast<uint32_t>(i));
  for (size_t i = 0; i < lens[3]; i++)
    sel[place(result.points[i].row, result.points[i].col)][3].push_back(static_cast<uint32_t>(i));

  std::vector<std::vector<Tile>> parts(nbands);
  std::atomic<size_t> next(0);
  std::vector<std::exception_ptr> errors(std::min(threads, nbands));
  std::vector<std::thread> pool;
  for (size_t t = 0; t < errors.size(); t++) {
    pool.emplace_back([&, t]() {
      try {
        for (;;) {
          size_t b = next.fetch_add(1, std::memory_order_relaxed);
          if (b >= nbands) break;
          parts[b] = s.slice(result, band_start(b), band_start(b + 1), sel[b].data());
        }
      } catch (...) {
        errors[t] = std::current_exception();
      }
    });
  }
  for (auto& th : pool) th.join();
  for (auto& e : errors)
    if (e) std::rethrow_exception(e);
  std::vector<Tile> tiles;
  for (auto& p : parts)
    for (auto& t : p) tiles.push_back(std::move(t));
  return tiles;
}

} // namespace controlledburn
