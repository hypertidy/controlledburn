// geometry.hpp -- minimal planar geometry model for controlledburn
//
// The rasterizer needs only: coordinates, ring structure, and geometry
// kind. No topology, no validity checking, no CRS. Rings are stored as
// flat coordinate vectors; orientation is derived from signed area at
// burn time (rings are normalised to CCW-exterior / CW-hole internally,
// so input orientation does not matter).
//
// Copyright (c) 2025 Michael Sumner
// Licensed under Apache License 2.0

#ifndef CONTROLLEDBURN_GEOMETRY_HPP
#define CONTROLLEDBURN_GEOMETRY_HPP

#include <vector>
#include <cstddef>

namespace controlledburn {

struct Coord {
    double x;
    double y;
};

// A ring or a linestring: an ordered coordinate sequence.
// For polygon rings the sequence should be closed (first == last);
// the burn is tolerant of unclosed input for linestrings.
using CoordSeq = std::vector<Coord>;

enum class GeomKind {
    Point,
    LineString,
    Polygon,
    MultiPoint,
    MultiLineString,
    MultiPolygon,
    // GeometryCollection is intentionally absent: mixed-dimension input
    // produces a sparse table with inconsistent weight semantics
    // (50 m^2 of polygon vs 50 m of line, indistinguishable in output).
    // Callers must split collections into homogeneous groups.
};

// One polygon: exterior ring first, then zero or more hole rings.
struct Polygon {
    std::vector<CoordSeq> rings;
};

// A geometry is a tagged union kept deliberately dumb: exactly one of
// the containers below is populated according to `kind`. Multi types
// use the same container as their singular counterpart with size > 1
// permitted (and size 1 or 0 also fine).
struct Geometry {
    GeomKind kind = GeomKind::Polygon;
    std::vector<CoordSeq> points;      // Point / MultiPoint: each seq has 1 coord
    std::vector<CoordSeq> lines;       // LineString / MultiLineString
    std::vector<Polygon>  polygons;    // Polygon / MultiPolygon

    bool empty() const {
        return points.empty() && lines.empty() && polygons.empty();
    }
};

// Signed area of a closed ring (shoelace). Positive = counter-clockwise
// in a conventional (x east, y north) frame.
double signed_area(const CoordSeq& ring);

inline bool is_ccw(const CoordSeq& ring) { return signed_area(ring) > 0.0; }

} // namespace controlledburn

#endif // CONTROLLEDBURN_GEOMETRY_HPP
