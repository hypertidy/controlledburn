// output.hpp -- sparse rasterization output: the four-table contract
//
// Polygon -> runs (interior RLE) + edges (boundary coverage fractions).
// Line    -> lines (length-in-cell, CRS units).
// Point   -> points (no measure column; implicit 1).
//
// All row/col indices are 1-based to match the R package contract.
// Schemas are type-pure: each table's measure column (or absence
// thereof) means exactly one thing.
//
// Copyright (c) 2025 Michael Sumner
// Licensed under Apache License 2.0

#ifndef CONTROLLEDBURN_OUTPUT_HPP
#define CONTROLLEDBURN_OUTPUT_HPP

#include <vector>
#include <string>
#include <cstdint>

namespace controlledburn {

// A single interior run: fully-covered cells [col_start, col_end] on `row`.
struct GridRun {
    int32_t row;
    int32_t col_start;
    int32_t col_end;
    int32_t id;
};

// A single polygon-boundary cell. `fraction` is the dimensionless
// coverage fraction in [0, 1]: area(polygon intersect cell) / area(cell).
struct GridEdge {
    int32_t row;
    int32_t col;
    float fraction;
    int32_t id;
};

// A single line cell. `length` is the absolute length of the line
// within the cell, in CRS units.
struct GridLine {
    int32_t row;
    int32_t col;
    float length;
    int32_t id;
};

// A single point cell.
struct GridPoint {
    int32_t row;
    int32_t col;
    int32_t id;
};

struct BurnResult {
    std::vector<GridRun>   runs;
    std::vector<GridEdge>  edges;
    std::vector<GridLine>  lines;
    std::vector<GridPoint> points;

    // Non-fatal problems encountered per input geometry (parse failures,
    // skipped GeometryCollections, ...). Bindings surface these as
    // warnings; index is the 1-based position in the input.
    struct Note {
        int32_t geom_index;
        std::string message;
    };
    std::vector<Note> notes;
};

} // namespace controlledburn

#endif // CONTROLLEDBURN_OUTPUT_HPP
