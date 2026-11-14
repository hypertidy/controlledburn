// scanline_shim.cpp -- cpp11 bridge from R to the controlledburn C++ core
//
// The engine lives in the standalone core (canonical sources under cpp/,
// synced here as core_*.cpp by tools/sync-core.sh); this shim is only
// marshalling: raw WKB vectors in, four data.frames out, notes surfaced
// as R warnings, and the R tables back into the core's for the
// operations on a finished burn (crop, tiles, materialise, relate).
//
// The core emits 0-based indices with an exclusive col_end and id = k;
// this shim restores the R package's contract by adding 1 to
// row/col/col_start/id (col_end is left as-is, since a 0-based exclusive
// col_end equals a 1-based inclusive one), and subtracts it again on the
// way back in. R output: 1-based row/col, row 1 at the top, geometry k
// (1-based) as id k + id_offset.
//
// No R API is called while the core runs (it may run on several threads):
// the inputs are copied or viewed before, the outputs built after.
//
// Copyright (c) 2025 Michael Sumner
// Licensed under Apache License 2.0

#include <cpp11.hpp>
#include <cpp11/list.hpp>
#include <cpp11/integers.hpp>
#include <cpp11/doubles.hpp>
#include <cpp11/logicals.hpp>
#include <cpp11/raws.hpp>
#include <cpp11/strings.hpp>

#include <cmath>
#include <limits>
#include <string>
#include <vector>

#include "controlledburn/controlledburn.hpp"

using namespace controlledburn;

namespace {

GridSpec make_grid(double xmin, double ymin, double xmax, double ymax, int ncol, int nrow) {
  // Validate here with the historical messages (the core validates too,
  // with its own wording).
  if (ncol <= 0 || nrow <= 0) {
    cpp11::stop("ncol and nrow must be positive");
  }
  if (!(xmax > xmin) || !(ymax > ymin)) {
    cpp11::stop("Invalid extent: xmax must be > xmin, ymax must be > ymin");
  }
  GridSpec grid{xmin, ymin, xmax, ymax, ncol, nrow};
  try {
    grid.validate();
  } catch (const std::exception& e) {
    cpp11::stop("%s", e.what());
  }
  return grid;
}

BurnOptions make_options(const std::string& mode, bool coalesce_runs, bool merge_parts,
                         int id_offset) {
  BurnOptions opts;
  if (mode == "approx") {
    opts.mode = BurnMode::Approx;
  } else if (mode != "coverage") {
    cpp11::stop("mode must be 'coverage' or 'approx'");
  }
  if (id_offset < 0) cpp11::stop("id_offset must be non-negative");
  opts.coalesce_runs = coalesce_runs;
  opts.merge_parts = merge_parts;
  opts.id_offset = id_offset;
  return opts;
}

// Views of the raw vectors; `keep` holds them alive (WKBSpan does not own).
std::vector<WKBSpan> wkb_spans(const cpp11::list& wkb_list, std::vector<cpp11::raws>& keep) {
  R_xlen_t n = wkb_list.size();
  std::vector<WKBSpan> spans;
  spans.reserve(static_cast<size_t>(n));
  keep.reserve(static_cast<size_t>(n));
  for (R_xlen_t k = 0; k < n; k++) {
    SEXP el = wkb_list[k];
    if (TYPEOF(el) != RAWSXP) {
      spans.push_back({nullptr, 0});
      continue;
    }
    keep.emplace_back(el);
    const cpp11::raws& r = keep.back();
    if (r.size() == 0) {
      spans.push_back({nullptr, 0});
    } else {
      spans.push_back({reinterpret_cast<const uint8_t*>(RAW(r)), static_cast<size_t>(r.size())});
    }
  }
  return spans;
}

void warn_notes(const BurnResult& res) {
  // Surface non-fatal notes (parse failures, GeometryCollection rejection,
  // short rings, per-geometry processing errors) as R warnings, matching
  // the historical warn-and-skip. Note ids are the core's, 0-based.
  for (const auto& note : res.notes) {
    cpp11::warning("geometry %d: %s", note.geom_index + 1, note.message.c_str());
  }
}

cpp11::writable::list data_frame(cpp11::writable::list cols, std::initializer_list<const char*> names,
                                 size_t n) {
  std::vector<std::string> nm(names.begin(), names.end());
  cpp11::writable::strings snames(static_cast<R_xlen_t>(nm.size()));
  for (size_t i = 0; i < nm.size(); i++) snames[static_cast<R_xlen_t>(i)] = nm[i];
  cols.attr("names") = snames;
  cols.attr("class") = "data.frame";
  cols.attr("row.names") = cpp11::writable::integers({NA_INTEGER, -static_cast<int>(n)});
  return cols;
}

// The four tables as R data.frames, 1-based, with row_off / col_off
// (0-based) subtracted from the indices first.
cpp11::writable::list as_r(const BurnResult& res, int row_off = 0, int col_off = 0) {
  size_t n_runs = res.runs.size();
  cpp11::writable::integers runs_row(n_runs), runs_col_start(n_runs), runs_col_end(n_runs),
      runs_id(n_runs);
  for (size_t i = 0; i < n_runs; i++) {
    runs_row[i] = res.runs[i].row - row_off + 1;
    runs_col_start[i] = res.runs[i].col_start - col_off + 1;
    runs_col_end[i] = res.runs[i].col_end - col_off;  // 0-based exclusive == 1-based inclusive
    runs_id[i] = res.runs[i].id + 1;
  }
  cpp11::writable::list runs_cols({runs_row, runs_col_start, runs_col_end, runs_id});
  auto runs_df = data_frame(runs_cols, {"row", "col_start", "col_end", "id"}, n_runs);

  size_t n_edges = res.edges.size();
  cpp11::writable::integers edges_row(n_edges), edges_col(n_edges), edges_id(n_edges);
  cpp11::writable::doubles edges_fraction(n_edges);
  for (size_t i = 0; i < n_edges; i++) {
    edges_row[i] = res.edges[i].row - row_off + 1;
    edges_col[i] = res.edges[i].col - col_off + 1;
    edges_fraction[i] = static_cast<double>(res.edges[i].fraction);
    edges_id[i] = res.edges[i].id + 1;
  }
  cpp11::writable::list edges_cols({edges_row, edges_col, edges_fraction, edges_id});
  auto edges_df = data_frame(edges_cols, {"row", "col", "fraction", "id"}, n_edges);

  size_t n_lines = res.lines.size();
  cpp11::writable::integers lines_row(n_lines), lines_col(n_lines), lines_id(n_lines);
  cpp11::writable::doubles lines_length(n_lines);
  for (size_t i = 0; i < n_lines; i++) {
    lines_row[i] = res.lines[i].row - row_off + 1;
    lines_col[i] = res.lines[i].col - col_off + 1;
    lines_length[i] = static_cast<double>(res.lines[i].length);
    lines_id[i] = res.lines[i].id + 1;
  }
  cpp11::writable::list lines_cols({lines_row, lines_col, lines_length, lines_id});
  auto lines_df = data_frame(lines_cols, {"row", "col", "length", "id"}, n_lines);

  size_t n_points = res.points.size();
  cpp11::writable::integers points_row(n_points), points_col(n_points), points_id(n_points);
  for (size_t i = 0; i < n_points; i++) {
    points_row[i] = res.points[i].row - row_off + 1;
    points_col[i] = res.points[i].col - col_off + 1;
    points_id[i] = res.points[i].id + 1;
  }
  cpp11::writable::list points_cols({points_row, points_col, points_id});
  auto points_df = data_frame(points_cols, {"row", "col", "id"}, n_points);

  cpp11::writable::list result({runs_df, edges_df, lines_df, points_df});
  result.attr("names") = cpp11::writable::strings({"runs", "edges", "lines", "points"});
  return result;
}

// An R table's column, or an empty vector when the table is NULL.
template <typename V>
V column(SEXP df, const char* name) {
  if (Rf_isNull(df)) return V();
  cpp11::list l(df);
  return V(l[name]);
}

// The four R tables (1-based) back into the core's (0-based).
BurnResult from_r(const cpp11::list& x) {
  BurnResult r;
  SEXP runs = x["runs"], edges = x["edges"], lines = x["lines"], points = x["points"];
  {
    auto row = column<cpp11::integers>(runs, "row");
    auto cs = column<cpp11::integers>(runs, "col_start");
    auto ce = column<cpp11::integers>(runs, "col_end");
    auto id = column<cpp11::integers>(runs, "id");
    r.runs.reserve(static_cast<size_t>(row.size()));
    for (R_xlen_t i = 0; i < row.size(); i++)
      r.runs.push_back({row[i] - 1, cs[i] - 1, ce[i], id[i] - 1});
  }
  {
    auto row = column<cpp11::integers>(edges, "row");
    auto col = column<cpp11::integers>(edges, "col");
    auto f = column<cpp11::doubles>(edges, "fraction");
    auto id = column<cpp11::integers>(edges, "id");
    r.edges.reserve(static_cast<size_t>(row.size()));
    for (R_xlen_t i = 0; i < row.size(); i++)
      r.edges.push_back({row[i] - 1, col[i] - 1, static_cast<float>(f[i]), id[i] - 1});
  }
  {
    auto row = column<cpp11::integers>(lines, "row");
    auto col = column<cpp11::integers>(lines, "col");
    auto len = column<cpp11::doubles>(lines, "length");
    auto id = column<cpp11::integers>(lines, "id");
    r.lines.reserve(static_cast<size_t>(row.size()));
    for (R_xlen_t i = 0; i < row.size(); i++)
      r.lines.push_back({row[i] - 1, col[i] - 1, static_cast<float>(len[i]), id[i] - 1});
  }
  {
    auto row = column<cpp11::integers>(points, "row");
    auto col = column<cpp11::integers>(points, "col");
    auto id = column<cpp11::integers>(points, "id");
    r.points.reserve(static_cast<size_t>(row.size()));
    for (R_xlen_t i = 0; i < row.size(); i++) r.points.push_back({row[i] - 1, col[i] - 1, id[i] - 1});
  }
  return r;
}

} // namespace

[[cpp11::register]]
cpp11::writable::list cpp_scanline_burn(
    cpp11::list wkb_list,
    double xmin, double ymin, double xmax, double ymax,
    int ncol, int nrow,
    std::string mode,
    bool coalesce_runs, bool merge_parts, int id_offset, int threads
) {
  GridSpec grid = make_grid(xmin, ymin, xmax, ymax, ncol, nrow);
  BurnOptions opts = make_options(mode, coalesce_runs, merge_parts, id_offset);
  if (threads < 0) cpp11::stop("threads must be non-negative");
  std::vector<cpp11::raws> keep;
  std::vector<WKBSpan> spans = wkb_spans(wkb_list, keep);
  BurnResult res;
  try {
    res = threads == 1 ? burn_wkb(spans, grid, opts)
                       : burn_wkb_par(spans, grid, opts, static_cast<unsigned>(threads));
  } catch (const std::exception& e) {
    cpp11::stop("%s", e.what());
  }
  warn_notes(res);
  return as_r(res);
}

// Burn onto the window rows [row0, row1), cols [col0, col1) (0-based) of
// the grid; the tables are re-based to the window.
[[cpp11::register]]
cpp11::writable::list cpp_scanline_burn_window(
    cpp11::list wkb_list,
    double xmin, double ymin, double xmax, double ymax,
    int ncol, int nrow,
    int row0, int row1, int col0, int col1,
    std::string mode,
    bool coalesce_runs, bool merge_parts, int id_offset
) {
  GridSpec grid = make_grid(xmin, ymin, xmax, ymax, ncol, nrow);
  BurnOptions opts = make_options(mode, coalesce_runs, merge_parts, id_offset);
  if (row0 < 0 || col0 < 0 || row1 < row0 || col1 < col0) cpp11::stop("invalid window");
  std::vector<cpp11::raws> keep;
  std::vector<WKBSpan> spans = wkb_spans(wkb_list, keep);
  BurnResult res;
  try {
    res = burn_wkb_window(spans, grid, static_cast<uint32_t>(row0), static_cast<uint32_t>(row1),
                          static_cast<uint32_t>(col0), static_cast<uint32_t>(col1), opts);
  } catch (const std::exception& e) {
    cpp11::stop("%s", e.what());
  }
  warn_notes(res);
  return as_r(res, row0, col0);
}

// crop: the tables of `x` (a controlledburn object on the grid) cropped
// to the cells overlapping the target (snapped outward), re-based to that
// window, with grid = c(xmin, ymin, xmax, ymax, ncol, nrow) of the window;
// NULL when the target does not overlap the grid.
[[cpp11::register]]
SEXP cpp_crop(cpp11::list x,
              double xmin, double ymin, double xmax, double ymax,
              int ncol, int nrow,
              double txmin, double tymin, double txmax, double tymax) {
  GridSpec grid = make_grid(xmin, ymin, xmax, ymax, ncol, nrow);
  BurnResult r = from_r(x);
  BurnResult o;
  GridSpec sg{};
  if (!crop(r, grid, txmin, tymin, txmax, tymax, o, sg)) return R_NilValue;
  cpp11::writable::list out = as_r(o);
  out.push_back(cpp11::named_arg("grid") = cpp11::writable::doubles(
                    {sg.xmin, sg.ymin, sg.xmax, sg.ymax,
                     static_cast<double>(sg.ncol), static_cast<double>(sg.nrow)}));
  return out;
}

// slice_tiles: a list of tiles, each a list of its four tables plus
// tile = c(tile_row, tile_col) and offset = c(row_off, col_off) (0-based)
// and grid = c(xmin, ymin, xmax, ymax, ncol, nrow).
[[cpp11::register]]
cpp11::writable::list cpp_slice_tiles(
    cpp11::list x,
    double xmin, double ymin, double xmax, double ymax,
    int ncol, int nrow,
    int tile_ncol, int tile_nrow, int threads
) {
  GridSpec grid = make_grid(xmin, ymin, xmax, ymax, ncol, nrow);
  if (tile_ncol <= 0 || tile_nrow <= 0) cpp11::stop("tile dimension must be positive");
  if (threads < 0) cpp11::stop("threads must be non-negative");
  BurnResult r = from_r(x);
  std::vector<Tile> tiles;
  try {
    tiles = threads == 1
                ? slice_tiles(r, grid, static_cast<uint32_t>(tile_ncol), static_cast<uint32_t>(tile_nrow))
                : slice_tiles_par(r, grid, static_cast<uint32_t>(tile_ncol),
                                  static_cast<uint32_t>(tile_nrow), static_cast<unsigned>(threads));
  } catch (const std::exception& e) {
    cpp11::stop("%s", e.what());
  }
  cpp11::writable::list out(static_cast<R_xlen_t>(tiles.size()));
  for (size_t i = 0; i < tiles.size(); i++) {
    const Tile& t = tiles[i];
    cpp11::writable::list tl = as_r(t.result);
    tl.push_back(cpp11::named_arg("tile") =
                     cpp11::writable::integers({static_cast<int>(t.tile_row), static_cast<int>(t.tile_col)}));
    tl.push_back(cpp11::named_arg("offset") =
                     cpp11::writable::integers({static_cast<int>(t.row_off), static_cast<int>(t.col_off)}));
    tl.push_back(cpp11::named_arg("grid") = cpp11::writable::doubles(
                     {t.grid.xmin, t.grid.ymin, t.grid.xmax, t.grid.ymax,
                      static_cast<double>(t.grid.ncol), static_cast<double>(t.grid.nrow)}));
    out[static_cast<R_xlen_t>(i)] = tl;
  }
  return out;
}

// Materialise the four tables into an nrow x ncol column-major R matrix.
// fun "sum": background 0, runs add 1, edges their fraction, lines their
// length, points 1. fun "id": background NA, the id of the last record
// written to a cell (runs, then edges, lines, points).
[[cpp11::register]]
cpp11::writable::doubles_matrix<> cpp_materialise(cpp11::list x, int ncol, int nrow, std::string fun) {
  if (ncol <= 0 || nrow <= 0) cpp11::stop("ncol and nrow must be positive");
  BurnResult r = from_r(x);
  const bool sum = fun == "sum";
  if (!sum && fun != "id") cpp11::stop("fun must be 'sum' or 'id'");
  cpp11::writable::doubles_matrix<> m(nrow, ncol);
  double* px = REAL(m);
  const size_t n = static_cast<size_t>(nrow) * static_cast<size_t>(ncol);
  const double bg = sum ? 0.0 : NA_REAL;
  for (size_t i = 0; i < n; i++) px[i] = bg;
  auto ok = [&](int32_t row, int32_t col) { return row >= 0 && row < nrow && col >= 0 && col < ncol; };
  auto at = [&](int32_t row, int32_t col) -> double& {
    return px[static_cast<size_t>(col) * static_cast<size_t>(nrow) + static_cast<size_t>(row)];
  };
  auto put = [&](int32_t row, int32_t col, double v, int32_t id) {
    if (!ok(row, col)) return;
    if (sum) {
      at(row, col) += v;
    } else {
      at(row, col) = static_cast<double>(id) + 1.0;
    }
  };
  for (const auto& run : r.runs) {
    if (run.row < 0 || run.row >= nrow) continue;
    for (int32_t c = std::max(run.col_start, 0); c < std::min(run.col_end, ncol); c++)
      put(run.row, c, 1.0, run.id);
  }
  for (const auto& e : r.edges) put(e.row, e.col, static_cast<double>(e.fraction), e.id);
  for (const auto& l : r.lines) put(l.row, l.col, static_cast<double>(l.length), l.id);
  for (const auto& p : r.points) put(p.row, p.col, 1.0, p.id);
  return m;
}

// relate_cells: DE-9IM of each polygon against the cells its boundary
// touches, from the burn `x` of the geometries `wkb_list` onto the grid.
[[cpp11::register]]
cpp11::writable::list cpp_relate_cells(
    cpp11::list wkb_list, cpp11::list x,
    double xmin, double ymin, double xmax, double ymax,
    int ncol, int nrow,
    int id_offset
) {
  GridSpec grid = make_grid(xmin, ymin, xmax, ymax, ncol, nrow);
  if (id_offset < 0) cpp11::stop("id_offset must be non-negative");
  BurnOptions opts;
  opts.id_offset = id_offset;
  BurnResult r = from_r(x);
  std::vector<cpp11::raws> keep;
  std::vector<WKBSpan> spans = wkb_spans(wkb_list, keep);
  std::vector<Geometry> geoms(spans.size());
  for (size_t k = 0; k < spans.size(); k++) {
    if (!spans[k].data || spans[k].size == 0) continue;  // empty: no relations
    try {
      geoms[k] = parse_wkb(spans[k].data, spans[k].size);
    } catch (const WKBParseError&) {
      geoms[k] = Geometry();  // the burn noted it; no relations
    }
  }
  std::vector<CellRelation> rel;
  try {
    rel = relate_cells(r, geoms, grid, opts);
  } catch (const std::exception& e) {
    cpp11::stop("%s", e.what());
  }
  size_t n = rel.size();
  cpp11::writable::integers row(n), col(n), id(n);
  cpp11::writable::strings matrix(static_cast<R_xlen_t>(n));
  cpp11::writable::logicals intersects(n), interiors(n), covers(n), touches(n), covered_by(n);
  for (size_t i = 0; i < n; i++) {
    const auto& c = rel[i];
    row[i] = c.row + 1;
    col[i] = c.col + 1;
    id[i] = c.id + 1;
    matrix[static_cast<R_xlen_t>(i)] = c.matrix.str();
    intersects[i] = c.matrix.intersects();
    interiors[i] = c.matrix.interiors_meet();
    covers[i] = c.matrix.covers();
    touches[i] = c.matrix.touches();
    covered_by[i] = c.matrix.covered_by();
  }
  cpp11::writable::list cols({row, col, id, matrix, intersects, interiors, covers, touches, covered_by});
  return data_frame(cols,
                    {"row", "col", "id", "matrix", "intersects", "interiors_meet", "covers",
                     "touches", "covered_by"},
                    n);
}
