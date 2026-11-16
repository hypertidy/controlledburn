// burn.hpp -- the controlledburn public entry points
//
// Scanline polygon/line/point rasterization with exact coverage
// fractions, O(perimeter) in time and memory, producing sparse output
// (no pixel buffer is ever materialized by the burn itself).
//
// The API follows the controlledburn Rust crate (0.3.2), which is the
// canonical implementation: the same options, defaults and records.
//
// Copyright (c) 2025 Michael Sumner
// Licensed under Apache License 2.0

#ifndef CONTROLLEDBURN_BURN_HPP
#define CONTROLLEDBURN_BURN_HPP

#include <vector>
#include <cmath>
#include <cstdint>
#include <stdexcept>

#include "geometry.hpp"
#include "output.hpp"

namespace controlledburn {

enum class BurnMode {
    Coverage,
    Approx
};

// Options for a burn, as the Rust crate's BurnOptions. The defaults are
// the crate's: coverage mode, runs coalesced, the parts of a multi-part
// geometry merged to one record per cell, ids from 0. A BurnMode
// converts to the default options with that mode.
struct BurnOptions {
    BurnMode mode = BurnMode::Coverage;
    // Merge adjacent runs on a row. Off, the runs table has one run per
    // fully covered boundary cell, as the core emitted until now.
    bool coalesce_runs = true;
    // One record per cell for a multi-part geometry: a cell several parts
    // cross gets one edge with the parts' fractions summed (in part
    // order), or a run once the sum reaches 1; an edge in a cell another
    // part covers fully is dropped; the geometry's runs are sorted and
    // coalesced. Line parts are summed per cell the same way. Overlapping
    // parts are not unioned: the sum is capped at a full cell. Off, each
    // part keeps its own records.
    bool merge_parts = true;
    // Added to each geometry's position to give its id, so burns of
    // consecutive slices of a collection concatenate without renumbering.
    // Must be non-negative, and the last id must fit in int32.
    int32_t id_offset = 0;

    BurnOptions() = default;
    BurnOptions(BurnMode m) : mode(m) {}  // NOLINT: implicit, as Rust's From

    static BurnOptions coverage() { return BurnOptions(); }
    static BurnOptions approx() { return BurnOptions(BurnMode::Approx); }
    // Record-for-record parity with the core before these options (and the Rust
    // crate's parity()): coalescing and part merging off.
    static BurnOptions parity(BurnMode m) {
        BurnOptions o(m);
        o.coalesce_runs = false;
        o.merge_parts = false;
        return o;
    }

    BurnOptions with_mode(BurnMode m) const { BurnOptions o = *this; o.mode = m; return o; }
    BurnOptions with_coalesce_runs(bool on) const { BurnOptions o = *this; o.coalesce_runs = on; return o; }
    BurnOptions with_merge_parts(bool on) const { BurnOptions o = *this; o.merge_parts = on; return o; }
    BurnOptions with_id_offset(int32_t off) const { BurnOptions o = *this; o.id_offset = off; return o; }
};

struct GridSpec {
    double xmin, ymin, xmax, ymax;
    int32_t ncol, nrow;

    double dx() const { return (xmax - xmin) / ncol; }
    double dy() const { return (ymax - ymin) / nrow; }

    void validate() const {
        if (ncol <= 0 || nrow <= 0)
            throw std::invalid_argument("ncol and nrow must be positive");
        // Non-finite bounds (or a span that overflows to inf) give a
        // non-finite cell size, and the walker then indexes far outside
        // the grid. The comparisons are written so NaN fails them.
        if (!std::isfinite(xmin) || !std::isfinite(xmax) ||
            !std::isfinite(ymin) || !std::isfinite(ymax))
            throw std::invalid_argument("invalid extent: bounds must be finite");
        if (!(xmax > xmin) || !(ymax > ymin))
            throw std::invalid_argument(
                "invalid extent: xmax must be > xmin, ymax must be > ymin");
        if (!std::isfinite(dx()) || !std::isfinite(dy()))
            throw std::invalid_argument(
                "invalid extent: extent span overflows double precision");
    }
};

struct WKBSpan {
    const uint8_t* data;
    size_t size;
};

// Burn geometries onto the grid. Geometry k is id = k + opts.id_offset.
// Empty geometries are skipped silently; geometries with non-finite
// coordinates, rings of fewer than 4 coordinates once closed, lines of one
// coordinate and geometries the engine cannot process are skipped with a
// note, the rest of the geometry being burned. Unclosed rings are closed.
// Throws std::invalid_argument for an invalid grid or id_offset.
BurnResult burn(const std::vector<Geometry>& geoms, const GridSpec& grid,
                const BurnOptions& opts = BurnOptions());

// As burn(), from WKB blobs (ISO and EWKB, both byte orders, Z/M
// skipped). Empty or null spans are skipped silently, unparseable blobs
// with a note.
BurnResult burn_wkb(const std::vector<WKBSpan>& wkb, const GridSpec& grid,
                    const BurnOptions& opts = BurnOptions());

// burn() on `threads` threads (0: all available), with the same result
// record for record and in the same order, notes included. Geometries are
// split into contiguous chunks of about equal estimated work (coordinates
// plus grid rows spanned), several per thread; at most one thread per
// 8,192 of that work, and on the calling thread with one thread or one
// geometry. One geometry is never split.
BurnResult burn_par(const std::vector<Geometry>& geoms, const GridSpec& grid,
                    const BurnOptions& opts, unsigned threads);

// burn_wkb() on `threads` threads (0: all available), with the same
// result as burn_wkb's, record for record and in the same order, notes
// included. Each blob is parsed by the thread that burns it, as it burns
// it, so no parsed geometry outlives its burn. Each blob's work, for the
// chunking, is read from its bytes without building the geometry.
// Chunking, the thread cap and the fallback to the calling thread are
// burn_par's; the errors are burn_wkb's.
BurnResult burn_wkb_par(const std::vector<WKBSpan>& wkb, const GridSpec& grid,
                        const BurnOptions& opts, unsigned threads);

namespace detail {
// The threaded burns' work estimate of a geometry, and of a WKB blob read
// without parsing it (1 for an empty or unparseable blob); for tests.
size_t burn_work(const Geometry& g, const GridSpec& grid);
size_t burn_work_wkb(const WKBSpan& s, const GridSpec& grid);
} // namespace detail

// Burn onto a window of the grid: rows [row0, row1), columns [col0,
// col1). The records are those burn() gives in the window's cells, the
// same values in the same order with the full grid's indices: runs are
// cut at the window's columns and nothing outside it is kept. The window
// keeps the full grid's lattice and edge rules, so the windows of a tiling
// put together are the whole burn. Throws std::invalid_argument for an
// empty window or one outside the grid.
BurnResult burn_window(const std::vector<Geometry>& geoms, const GridSpec& grid,
                       uint32_t row0, uint32_t row1, uint32_t col0, uint32_t col1,
                       const BurnOptions& opts = BurnOptions());

BurnResult burn_wkb_window(const std::vector<WKBSpan>& wkb, const GridSpec& grid,
                           uint32_t row0, uint32_t row1, uint32_t col0, uint32_t col1,
                           const BurnOptions& opts = BurnOptions());

struct BBox {
    double xmin, ymin, xmax, ymax;
    bool valid;
};
BBox bbox_wkb(const std::vector<WKBSpan>& wkb);

} // namespace controlledburn

#endif // CONTROLLEDBURN_BURN_HPP
