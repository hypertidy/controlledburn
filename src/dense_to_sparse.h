// Copyright (c) 2025 Michael Sumner
// Licensed under Apache License 2.0

#ifndef GRIDBURN_DENSE_TO_SPARSE_H
#define GRIDBURN_DENSE_TO_SPARSE_H

#include <vector>
#include <cmath>

// A single interior run: row, col_start, col_end (all 1-based, full raster coords)
struct GridRun {
    int row;
    int col_start;
    int col_end;
    int id;
};

// A single edge cell: row, col, weight (all 1-based, full raster coords)
struct GridEdge {
    int row;
    int col;
    float weight;
    int id;
};

struct SparseResult {
    std::vector<GridRun> runs;
    std::vector<GridEdge> edges;
};

// Convert a dense coverage fraction matrix to sparse two-table format.
//
// mat_data: pointer to row-major float array (nrow x ncol)
// nrow, ncol: dimensions of the coverage matrix (the subgrid)
// row_offset, col_offset: offset of subgrid origin in the full raster (0-based)
// id: polygon identifier
// tol: tolerance for "fully covered" (weight >= 1.0 - tol → interior)
//
// Returns SparseResult with runs (interior RLE) and edges (boundary cells).
inline SparseResult dense_to_sparse(
    const float* mat_data,
    size_t nrow,
    size_t ncol,
    size_t row_offset,
    size_t col_offset,
    int id,
    float tol = 1e-6f
) {
    SparseResult result;

    for (size_t i = 0; i < nrow; i++) {
        int full_row = static_cast<int>(row_offset + i) + 1; // 1-based

        // State for RLE compression of interior cells
        int run_start = -1; // -1 means no active run

        for (size_t j = 0; j < ncol; j++) {
            float w = mat_data[i * ncol + j];

            if (w <= 0.0f) {
                // Outside cell — close any active run
                if (run_start >= 0) {
                    int full_col_end = static_cast<int>(col_offset + j - 1) + 1;
                    result.runs.push_back({full_row, run_start, full_col_end, id});
                    run_start = -1;
                }
                continue;
            }

            int full_col = static_cast<int>(col_offset + j) + 1; // 1-based

            if (w >= (1.0f - tol)) {
                // Interior cell — extend or start a run
                if (run_start < 0) {
                    run_start = full_col;
                }
                // run continues...
            } else {
                // Edge cell (0 < w < 1) — close any active run, emit edge
                if (run_start >= 0) {
                    int full_col_end = static_cast<int>(col_offset + j - 1) + 1;
                    result.runs.push_back({full_row, run_start, full_col_end, id});
                    run_start = -1;
                }
                result.edges.push_back({full_row, full_col, w, id});
            }
        }

        // Close any run that extends to the end of the row
        if (run_start >= 0) {
            int full_col_end = static_cast<int>(col_offset + ncol - 1) + 1;
            result.runs.push_back({full_row, run_start, full_col_end, id});
        }
    }

    return result;
}

#endif // GRIDBURN_DENSE_TO_SPARSE_H
