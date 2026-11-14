// ops.hpp -- operations on a finished burn: run coalescing, merging
// records per cell, cropping to a window, and slicing into tiles.
//
// Ports of the controlledburn Rust crate's coalesce_runs, merge_cells,
// crop, slice_tiles and slice_tiles_par, with the same records in the
// same order.
//
// Copyright (c) 2026 Michael Sumner
// Licensed under Apache License 2.0

#ifndef CONTROLLEDBURN_OPS_HPP
#define CONTROLLEDBURN_OPS_HPP

#include <cstdint>
#include <vector>

#include "burn.hpp"
#include "output.hpp"

namespace controlledburn {

// Merge adjacent or overlapping runs with the same id on the same row.
// The output is sorted by (id, row, col_start).
void coalesce_runs(std::vector<GridRun>& runs);

// One record per (id, row, col), by the rule burn() uses for the parts of
// a multi-part geometry (BurnOptions::merge_parts): edges in one cell with
// one id are summed (float, in table order), a sum reaching 1 - 1e-6
// becomes a one-cell run; runs are sorted and coalesced; an edge in a cell
// a run of the same id covers is dropped; lines in one cell with one id
// are summed. Points and notes are left alone. The tables come out sorted
// by (id, row, col). For records that come together after a burn (a
// front end folding columns back onto a periodic grid, say).
void merge_cells(BurnResult& result);

// Crop a burn to the cells of `grid` overlapping `target` (xmin, ymin,
// xmax, ymax), snapped outward to whole cells and clamped to the grid,
// with indices re-based to the window. Returns false, leaving `out` and
// `subgrid` alone, when the target does not overlap the grid. Notes are
// carried through.
bool crop(const BurnResult& result, const GridSpec& grid,
          double xmin, double ymin, double xmax, double ymax,
          BurnResult& out, GridSpec& subgrid);

// One non-empty tile of slice_tiles.
struct Tile {
    uint32_t tile_row;  // in the tile grid, 0 at the top
    uint32_t tile_col;
    uint32_t row_off;   // first cell of the tile in the whole grid
    uint32_t col_off;
    GridSpec grid;      // the tile's own grid
    BurnResult result;  // indexed on `grid`, in burn order; no notes
};

// Slice a burn into tiles of tile_ncol x tile_nrow cells anchored at the
// grid's top-left cell, in one pass: a record in (row, col) goes to tile
// (row / tile_nrow, col / tile_ncol), runs split at tile column seams.
// Only tiles holding records are returned, in row-major tile order. Each
// equals crop() of the whole burn to its window. Throws
// std::invalid_argument for a zero tile shape or a record outside `grid`.
std::vector<Tile> slice_tiles(const BurnResult& result, const GridSpec& grid,
                              uint32_t tile_ncol, uint32_t tile_nrow);

// slice_tiles on `threads` threads (0: all available), with the same tiles
// and records in the same order. At most one thread per 131,072 records;
// serial with one thread, a single tile row or fewer than 262,144 records.
std::vector<Tile> slice_tiles_par(const BurnResult& result, const GridSpec& grid,
                                  uint32_t tile_ncol, uint32_t tile_nrow,
                                  unsigned threads);

} // namespace controlledburn

#endif // CONTROLLEDBURN_OPS_HPP
