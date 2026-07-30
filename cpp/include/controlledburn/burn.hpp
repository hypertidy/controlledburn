// burn.hpp -- the controlledburn public entry points
//
// Scanline polygon/line/point rasterization with exact coverage
// fractions, O(perimeter) in time and memory, producing sparse output
// (no pixel buffer is ever materialized by the burn itself).
//
// Copyright (c) 2025 Michael Sumner
// Licensed under Apache License 2.0

#ifndef CONTROLLEDBURN_BURN_HPP
#define CONTROLLEDBURN_BURN_HPP

#include <vector>
#include <cstdint>
#include <stdexcept>

#include "geometry.hpp"
#include "output.hpp"

namespace controlledburn {

// How boundary cells are classified for polygon output.
//
//   Coverage -- exact coverage fractions via analytical traversal.
//               Boundary cells appear in $edges with a fraction in (0, 1).
//               This is the default: exact area conservation.
//
//   Approx   -- cell-center rule (fasterize semantics). A boundary cell
//               is "inside" iff the polygon edge crosses the cell-center
//               scanline to the LEFT of the cell center. Inside cells
//               become runs; outside cells are dropped. No $edges for
//               polygons. Faster (skips traversal math), but approximate.
enum class BurnMode {
    Coverage,
    Approx
};

// Regular grid specification. Cell size is derived: dx = (xmax-xmin)/ncol.
// Row 1 is the TOP row (ymax edge), matching raster convention.
struct GridSpec {
    double xmin, ymin, xmax, ymax;
    int32_t ncol, nrow;

    double dx() const { return (xmax - xmin) / ncol; }
    double dy() const { return (ymax - ymin) / nrow; }

    void validate() const {
        if (ncol <= 0 || nrow <= 0)
            throw std::invalid_argument("ncol and nrow must be positive");
        if (xmax <= xmin || ymax <= ymin)
            throw std::invalid_argument(
                "invalid extent: xmax must be > xmin, ymax must be > ymin");
    }
};

// Burn a set of geometries onto the grid. Geometry k (0-based) is
// assigned id k + 1 in the output tables.
//
// Empty geometries are skipped silently; unsupported inputs are skipped
// with a note in BurnResult::notes. Throws std::invalid_argument for an
// invalid GridSpec only.
BurnResult burn(const std::vector<Geometry>& geoms, const GridSpec& grid,
                BurnMode mode = BurnMode::Coverage);

// Convenience: burn geometries supplied as WKB blobs (one blob per
// geometry). Unparseable blobs are skipped with a note.
struct WKBSpan {
    const uint8_t* data;
    size_t size;
};
BurnResult burn_wkb(const std::vector<WKBSpan>& wkb, const GridSpec& grid,
                    BurnMode mode = BurnMode::Coverage);

} // namespace controlledburn

#endif // CONTROLLEDBURN_BURN_HPP
