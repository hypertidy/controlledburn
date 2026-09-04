// materialize.hpp -- optional consumer: burn a sparse BurnResult into a
// caller-provided pixel buffer with per-pixel reduction functions.
//
// This is the fasterize reconciliation layer. fasterize's engine
// (edgelist + Armadillo raster + pixelfn) materializes as it burns; here
// the burn stays sparse and materialization is a separate, optional pass
// over the runs/edges tables. One engine, two consumption styles.
//
// SEMANTIC NOTE (important for fasterize parity): fasterize burns a cell
// when the CELL CENTER is inside the polygon. controlledburn computes
// EXACT COVERAGE. For interior runs the two agree; for boundary cells
// they differ. The `edge_policy` below controls how boundary cells
// materialize:
//   - Fraction:  write value weighted by coverage fraction (exact mode)
//   - Threshold: include the cell iff fraction >= threshold. With
//                threshold = 0.5 this approximates -- but is not
//                identical to -- the center-in-polygon rule. Exact
//                center-rule parity would require a center-point test
//                per boundary cell and is a candidate for a dedicated
//                sweep mode in the core.
//
// The buffer is row-major, row 0 (top) first, `ncol * nrow` doubles,
// caller-allocated and caller-initialized (typically to NaN background).
//
// Copyright (c) 2025 Michael Sumner
// Licensed under Apache License 2.0

#ifndef CONTROLLEDBURN_MATERIALIZE_HPP
#define CONTROLLEDBURN_MATERIALIZE_HPP

#include <cmath>
#include <cstddef>
#include <cstdint>
#include <vector>
#include <stdexcept>

#include "output.hpp"

namespace controlledburn {

// Reduction applied when multiple geometries touch the same pixel.
// Mirrors fasterize's pixel functions.
enum class PixelFn {
    First,  // keep the first value written
    Last,   // overwrite with each new value
    Sum,    // accumulate
    Min,
    Max,
    Count,  // number of geometries touching the pixel (ignores value)
    Any     // 1 if touched
};

enum class EdgePolicy {
    Fraction,   // value * coverage_fraction, combined per PixelFn
    Threshold   // all-or-nothing at `threshold`
};

struct MaterializeOptions {
    PixelFn fn = PixelFn::Last;
    EdgePolicy edge_policy = EdgePolicy::Threshold;
    double threshold = 0.5;
};

namespace detail {

inline void apply(double* px, double value, PixelFn fn) {
    bool empty = std::isnan(*px);
    switch (fn) {
    case PixelFn::First: if (empty) *px = value; break;
    case PixelFn::Last:  *px = value; break;
    case PixelFn::Sum:   *px = empty ? value : *px + value; break;
    case PixelFn::Min:   *px = empty ? value : (value < *px ? value : *px); break;
    case PixelFn::Max:   *px = empty ? value : (value > *px ? value : *px); break;
    case PixelFn::Count: *px = empty ? 1.0 : *px + 1.0; break;
    case PixelFn::Any:   *px = 1.0; break;
    }
}

} // namespace detail

// Materialize polygon output (runs + edges) into `buffer`.
//
// `values` maps geometry id to burn value: value for id k is
// values[k]. Pass an empty vector to burn the id itself.
// Cells never touched are left untouched (background is whatever the
// caller initialized, conventionally NaN).
inline void materialize(
    const BurnResult& result,
    double* buffer,
    int32_t ncol, int32_t nrow,
    const std::vector<double>& values = {},
    const MaterializeOptions& opts = {}
) {
    if (!buffer || ncol <= 0 || nrow <= 0)
        throw std::invalid_argument("materialize: invalid buffer or dimensions");

    auto value_of = [&](int32_t id) -> double {
        if (values.empty()) return static_cast<double>(id);
        size_t i = static_cast<size_t>(id);
        if (i >= values.size())
            throw std::out_of_range("materialize: geometry id exceeds values size");
        return values[i];
    };

    auto px = [&](int32_t row, int32_t col) -> double* {
        // row/col are 0-based
        return buffer + static_cast<size_t>(row) * ncol + col;
    };

    for (const auto& r : result.runs) {
        if (r.row < 0 || r.row >= nrow) continue;
        int32_t c0 = r.col_start < 0 ? 0 : r.col_start;
        int32_t c1 = r.col_end > ncol ? ncol : r.col_end;  // exclusive
        double v = value_of(r.id);
        for (int32_t c = c0; c < c1; c++) {
            detail::apply(px(r.row, c), v, opts.fn);
        }
    }

    for (const auto& e : result.edges) {
        if (e.row < 0 || e.row >= nrow || e.col < 0 || e.col >= ncol) continue;
        double v = value_of(e.id);
        if (opts.edge_policy == EdgePolicy::Threshold) {
            if (e.fraction >= opts.threshold) {
                detail::apply(px(e.row, e.col), v, opts.fn);
            }
        } else {
            detail::apply(px(e.row, e.col), v * e.fraction, opts.fn);
        }
    }
}

} // namespace controlledburn

#endif // CONTROLLEDBURN_MATERIALIZE_HPP
