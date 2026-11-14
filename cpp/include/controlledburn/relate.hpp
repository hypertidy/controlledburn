// relate.hpp -- DE-9IM of polygons against grid cells, after the burn.
//
// Port of the controlledburn Rust crate's relate_cells (0.3.1). For a cell
// the polygon's boundary does not reach, the burn answers already: a run
// cell is COVERED, a cell with no record is DISJOINT. relate_cells gives
// the matrix of every other cell, those whose closed box the boundary
// touches, including cells it only touches along a wall or at a corner
// (which have no record). See the crate's book, "Predicates per cell".
//
// Copyright (c) 2026 Michael Sumner
// Licensed under Apache License 2.0

#ifndef CONTROLLEDBURN_RELATE_HPP
#define CONTROLLEDBURN_RELATE_HPP

#include <array>
#include <cstdint>
#include <string>
#include <vector>

#include "burn.hpp"
#include "geometry.hpp"
#include "output.hpp"

namespace controlledburn {

// A DE-9IM matrix: the polygon's interior, boundary and exterior (rows)
// against the cell's (columns), row-major; each entry 'F', '0', '1', '2'.
struct IntersectionMatrix {
    std::array<char, 9> m;

    static IntersectionMatrix covered() { return from("212FF1FF2"); }
    static IntersectionMatrix disjoint() { return from("FF2FF1212"); }
    static IntersectionMatrix from(const char* s) {
        IntersectionMatrix x{};
        for (int i = 0; i < 9; i++) x.m[i] = s[i];
        return x;
    }

    std::string str() const { return std::string(m.begin(), m.end()); }

    // A DE-9IM pattern of nine characters from T F * 0 1 2.
    bool matches(const std::string& p) const {
        if (p.size() != 9) return false;
        for (int i = 0; i < 9; i++) {
            char c = p[i], v = m[i];
            if (c == '*') continue;
            if (c == 'T' || c == 't') { if (v == 'F') return false; continue; }
            if (c == 'F' || c == 'f') { if (v != 'F') return false; continue; }
            if (c != v) return false;
        }
        return true;
    }
    bool intersects() const { return m[0] != 'F' || m[1] != 'F' || m[3] != 'F' || m[4] != 'F'; }
    bool interiors_meet() const { return m[0] != 'F'; }
    bool covers() const { return m[6] == 'F' && m[7] == 'F'; }
    bool touches() const { return m[0] == 'F' && intersects(); }
    bool covered_by() const { return m[2] == 'F' && m[5] == 'F'; }

    bool operator==(const IntersectionMatrix& o) const { return m == o.m; }
    bool operator!=(const IntersectionMatrix& o) const { return m != o.m; }
};

struct CellRelation {
    int32_t row;
    int32_t col;
    int32_t id;
    IntersectionMatrix matrix;
};

// DE-9IM matrices of polygons against the in-grid cells their boundaries
// touch, from a finished burn of the same geometries onto the same grid
// with the same options (only runs are read, and opts.id_offset gives the
// ids). Sorted by id, row and column. Lines, points, empty geometries and
// geometries with a non-finite coordinate give no relations. Polygons are
// assumed valid (OGC). Throws std::invalid_argument for an invalid grid
// or id_offset.
std::vector<CellRelation> relate_cells(const BurnResult& result,
                                       const std::vector<Geometry>& geoms,
                                       const GridSpec& grid,
                                       const BurnOptions& opts = BurnOptions());

} // namespace controlledburn

#endif // CONTROLLEDBURN_RELATE_HPP
