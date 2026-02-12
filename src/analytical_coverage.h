// analytical_coverage.h — analytical coverage fraction for single-traversal cells
//
// For a single edge traversal through a grid cell, the covered area (to the
// left of the traversal in CCW winding) is a simple polygon: the traversal
// path plus CW cell boundary corners from exit back to entry. This avoids the
// full left_hand_area chain-chasing algorithm in traversal_areas.cpp.
//
// Uses perimeter_distance (same convention as left_hand_area) to determine
// which corners fall in the CW arc from exit to entry.
//
// Copyright (c) 2025 Michael Sumner
// Licensed under Apache License 2.0

#ifndef DENSEBURN_ANALYTICAL_COVERAGE_H
#define DENSEBURN_ANALYTICAL_COVERAGE_H

#include <vector>
#include <cmath>
#include <cstddef>
#include <algorithm>

#include "exactextract/box.h"
#include "exactextract/coordinate.h"
#include "exactextract/side.h"
#include "exactextract/perimeter_distance.h"

namespace controlledburn {

// Compute signed area of a closed polygon ring using the shoelace formula.
// Equivalent to exactextract::area() in measures.cpp.
inline double polygon_signed_area(const std::vector<exactextract::Coordinate>& ring) {
    if (ring.size() < 3) return 0.0;

    double sum = 0.0;
    double x0 = ring[0].x;
    for (size_t i = 1; i < ring.size() - 1; i++) {
        double x = ring[i].x - x0;
        double y1 = ring[i + 1].y;
        double y2 = ring[i - 1].y;
        sum += x * (y2 - y1);
    }
    return sum / 2.0;
}

// Compute exact coverage fraction for a single traversal through a cell.
//
// The traversal enters at coords.front() and exits at coords.back().
// For a CCW-oriented ring, the covered area is to the LEFT of the traversal
// direction. This equals the area of the polygon formed by:
//   1. The traversal path (entry → ... → exit)
//   2. CW cell boundary from exit back to entry (inserting corner points)
//
// The CW direction is used because after exiting a traversal, the left-hand
// region closes by walking BACKWARD along the cell perimeter to the entry.
// This matches the chain-chasing direction in left_hand_area.
//
// Returns coverage fraction in [0, 1].
inline double analytical_covered_fraction(
    const exactextract::Box& box,
    const std::vector<exactextract::Coordinate>& coords,
    exactextract::Side /* entry_side — unused, kept for API compat */,
    exactextract::Side /* exit_side — unused, kept for API compat */
) {
    using exactextract::Coordinate;
    using exactextract::perimeter_distance;

    double cell_area = box.area();
    if (cell_area <= 0.0) return 0.0;
    if (coords.size() < 2) return 0.0;

    double perim = box.perimeter();

    // Perimeter distances of exit and entry points
    double exit_pd  = perimeter_distance(box, coords.back());
    double entry_pd = perimeter_distance(box, coords.front());

    // CW distance from exit to entry (going backward along perimeter)
    // This is the same as perimeter_distance_ccw(exit_pd, entry_pd, perim)
    // which is what left_hand_area uses to find the next chain.
    double arc;
    if (exit_pd > entry_pd + 1e-12) {
        arc = exit_pd - entry_pd;
    } else if (entry_pd > exit_pd + 1e-12) {
        arc = perim - entry_pd + exit_pd;
    } else {
        // Entry ≈ exit (degenerate): traversal starts and ends at same point.
        std::vector<Coordinate> poly(coords.begin(), coords.end());
        if (!(poly.front() == poly.back())) {
            poly.push_back(poly.front());
        }
        return std::abs(polygon_signed_area(poly)) / cell_area;
    }

    // The 4 corners and their perimeter distances (CCW from BL)
    //   BL (xmin,ymin) : 0
    //   TL (xmin,ymax) : h
    //   TR (xmax,ymax) : h + w
    //   BR (xmax,ymin) : 2h + w
    double h = box.height();
    double w = box.width();
    const Coordinate corner_coord[4] = {
        {box.xmin, box.ymin},  // BL
        {box.xmin, box.ymax},  // TL
        {box.xmax, box.ymax},  // TR
        {box.xmax, box.ymin},  // BR
    };
    const double corner_pd[4] = { 0.0, h, h + w, 2.0 * h + w };

    // CW distance from exit_pd to a given perimeter position
    auto cw_from_exit = [&](double pd) -> double {
        double d = exit_pd - pd;
        if (d < 0.0) d += perim;
        return d;
    };

    // Collect corners that lie strictly inside the CW arc from exit to entry,
    // tagged with their CW distance from exit for sorting.
    struct TaggedCorner { Coordinate coord; double dist; };
    TaggedCorner in_arc[4];
    int n_in_arc = 0;

    for (int i = 0; i < 4; i++) {
        double d = cw_from_exit(corner_pd[i]);
        if (d > 1e-12 && d < arc - 1e-12) {
            in_arc[n_in_arc++] = { corner_coord[i], d };
        }
    }

    // Sort by CW distance from exit (nearest first)
    std::sort(in_arc, in_arc + n_in_arc,
        [](const TaggedCorner& a, const TaggedCorner& b) {
            return a.dist < b.dist;
        });

    // Build the left-hand polygon: traversal path + CW corners + close
    std::vector<Coordinate> polygon;
    polygon.reserve(coords.size() + n_in_arc + 1);

    // 1. Traversal coordinates (entry → ... → exit)
    polygon.insert(polygon.end(), coords.begin(), coords.end());

    // 2. Corners in CW order from exit toward entry
    for (int i = 0; i < n_in_arc; i++) {
        polygon.push_back(in_arc[i].coord);
    }

    // 3. Close
    polygon.push_back(polygon.front());

    return std::abs(polygon_signed_area(polygon)) / cell_area;
}

// Compute coverage fraction for a closed ring entirely within one cell.
inline double closed_ring_covered_fraction(
    const exactextract::Box& box,
    const std::vector<exactextract::Coordinate>& ring_coords
) {
    double cell_area = box.area();
    if (cell_area <= 0.0) return 0.0;
    return std::abs(polygon_signed_area(ring_coords)) / cell_area;
}

} // namespace controlledburn

#endif // DENSEBURN_ANALYTICAL_COVERAGE_H
