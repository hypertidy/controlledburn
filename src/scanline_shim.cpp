// scanline_shim.cpp -- cpp11 bridge from R to the controlledburn C++ core
//
// This file replaces src/scanline_burn.cpp. The engine now lives in the
// standalone core (canonical sources under cpp/, synced here as
// core_scanline.cpp / core_wkb.cpp by tools/sync-core.sh); this shim is
// only marshalling: raw WKB vectors in, four data.frames out, notes
// surfaced as R warnings. WKB parsing is done by the core's own reader,
// not GEOS -- GEOS remains in the package only for the older
// burn_sparse path (controlledburn.cpp).
//
// The output contract is unchanged: runs/edges/lines/points tables with
// 1-based row/col, row 1 at the top, geometry k (1-based) as id k.
//
// Copyright (c) 2025 Michael Sumner
// Licensed under Apache License 2.0

#include <cpp11.hpp>
#include <cpp11/list.hpp>
#include <cpp11/integers.hpp>
#include <cpp11/doubles.hpp>
#include <cpp11/raws.hpp>
#include <cpp11/strings.hpp>

#include <vector>

#include "controlledburn/burn.hpp"

[[cpp11::register]]
cpp11::writable::list cpp_scanline_burn(
    cpp11::list wkb_list,
    double xmin, double ymin, double xmax, double ymax,
    int ncol, int nrow
) {
  // Validate here with the historical messages (the core validates too,
  // with its own wording).
  if (ncol <= 0 || nrow <= 0) {
    cpp11::stop("ncol and nrow must be positive");
  }
  if (xmax <= xmin || ymax <= ymin) {
    cpp11::stop("Invalid extent: xmax must be > xmin, ymax must be > ymin");
  }

  controlledburn::GridSpec grid{xmin, ymin, xmax, ymax, ncol, nrow};

  int n_geoms = wkb_list.size();
  std::vector<controlledburn::WKBSpan> spans;
  spans.reserve(n_geoms);
  // Keep the raws vectors alive for the duration of the burn; WKBSpan is
  // a non-owning view.
  std::vector<cpp11::raws> keep;
  keep.reserve(n_geoms);

  for (int k = 0; k < n_geoms; k++) {
    keep.emplace_back(wkb_list[k]);
    const cpp11::raws& r = keep.back();
    if (r.size() == 0) {
      spans.push_back({nullptr, 0});
    } else {
      spans.push_back({
        reinterpret_cast<const uint8_t*>(RAW(r)),
        static_cast<size_t>(r.size())
      });
    }
  }

  controlledburn::BurnResult res = controlledburn::burn_wkb(spans, grid);

  // Surface non-fatal notes (parse failures, GeometryCollection
  // rejection, per-geometry processing errors) as R warnings, matching
  // the historical behaviour of warn-and-skip.
  for (const auto& note : res.notes) {
    cpp11::warning("geometry %d: %s", note.geom_index, note.message.c_str());
  }

  // Build R data.frames -- four tables, one per geometry kind.
  // Polygon -> runs (interior RLE) + edges (boundary fractions in [0, 1]).
  // Line    -> lines (length in CRS units).
  // Point   -> points (no measure column; implicit 1).

  size_t n_runs = res.runs.size();
  cpp11::writable::integers runs_row(n_runs);
  cpp11::writable::integers runs_col_start(n_runs);
  cpp11::writable::integers runs_col_end(n_runs);
  cpp11::writable::integers runs_id(n_runs);

  for (size_t i = 0; i < n_runs; i++) {
    runs_row[i] = res.runs[i].row;
    runs_col_start[i] = res.runs[i].col_start;
    runs_col_end[i] = res.runs[i].col_end;
    runs_id[i] = res.runs[i].id;
  }

  cpp11::writable::list runs_df(4);
  runs_df[0] = static_cast<SEXP>(runs_row);
  runs_df[1] = static_cast<SEXP>(runs_col_start);
  runs_df[2] = static_cast<SEXP>(runs_col_end);
  runs_df[3] = static_cast<SEXP>(runs_id);
  runs_df.attr("names") = cpp11::writable::strings({"row", "col_start", "col_end", "id"});
  runs_df.attr("class") = "data.frame";
  runs_df.attr("row.names") = cpp11::writable::integers({NA_INTEGER, -static_cast<int>(n_runs)});

  size_t n_edges = res.edges.size();
  cpp11::writable::integers edges_row(n_edges);
  cpp11::writable::integers edges_col(n_edges);
  cpp11::writable::doubles edges_fraction(n_edges);
  cpp11::writable::integers edges_id(n_edges);

  for (size_t i = 0; i < n_edges; i++) {
    edges_row[i] = res.edges[i].row;
    edges_col[i] = res.edges[i].col;
    edges_fraction[i] = static_cast<double>(res.edges[i].fraction);
    edges_id[i] = res.edges[i].id;
  }

  cpp11::writable::list edges_df(4);
  edges_df[0] = static_cast<SEXP>(edges_row);
  edges_df[1] = static_cast<SEXP>(edges_col);
  edges_df[2] = static_cast<SEXP>(edges_fraction);
  edges_df[3] = static_cast<SEXP>(edges_id);
  edges_df.attr("names") = cpp11::writable::strings({"row", "col", "fraction", "id"});
  edges_df.attr("class") = "data.frame";
  edges_df.attr("row.names") = cpp11::writable::integers({NA_INTEGER, -static_cast<int>(n_edges)});

  size_t n_lines = res.lines.size();
  cpp11::writable::integers lines_row(n_lines);
  cpp11::writable::integers lines_col(n_lines);
  cpp11::writable::doubles lines_length(n_lines);
  cpp11::writable::integers lines_id(n_lines);

  for (size_t i = 0; i < n_lines; i++) {
    lines_row[i] = res.lines[i].row;
    lines_col[i] = res.lines[i].col;
    lines_length[i] = static_cast<double>(res.lines[i].length);
    lines_id[i] = res.lines[i].id;
  }

  cpp11::writable::list lines_df(4);
  lines_df[0] = static_cast<SEXP>(lines_row);
  lines_df[1] = static_cast<SEXP>(lines_col);
  lines_df[2] = static_cast<SEXP>(lines_length);
  lines_df[3] = static_cast<SEXP>(lines_id);
  lines_df.attr("names") = cpp11::writable::strings({"row", "col", "length", "id"});
  lines_df.attr("class") = "data.frame";
  lines_df.attr("row.names") = cpp11::writable::integers({NA_INTEGER, -static_cast<int>(n_lines)});

  size_t n_points = res.points.size();
  cpp11::writable::integers points_row(n_points);
  cpp11::writable::integers points_col(n_points);
  cpp11::writable::integers points_id(n_points);

  for (size_t i = 0; i < n_points; i++) {
    points_row[i] = res.points[i].row;
    points_col[i] = res.points[i].col;
    points_id[i] = res.points[i].id;
  }

  cpp11::writable::list points_df(3);
  points_df[0] = static_cast<SEXP>(points_row);
  points_df[1] = static_cast<SEXP>(points_col);
  points_df[2] = static_cast<SEXP>(points_id);
  points_df.attr("names") = cpp11::writable::strings({"row", "col", "id"});
  points_df.attr("class") = "data.frame";
  points_df.attr("row.names") = cpp11::writable::integers({NA_INTEGER, -static_cast<int>(n_points)});

  cpp11::writable::list result(4);
  result[0] = static_cast<SEXP>(runs_df);
  result[1] = static_cast<SEXP>(edges_df);
  result[2] = static_cast<SEXP>(lines_df);
  result[3] = static_cast<SEXP>(points_df);
  result.attr("names") = cpp11::writable::strings({"runs", "edges", "lines", "points"});

  return result;
}
