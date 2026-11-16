// wkb.hpp -- minimal WKB reader for the geometry model in geometry.hpp
//
// Supports ISO WKB and EWKB: both byte orders, 2D coordinates with Z
// and/or M ordinates skipped (ISO 1000/2000/3000 offsets and EWKB
// 0x80000000 / 0x40000000 flags), EWKB embedded SRID skipped. Curved
// types and GeometryCollection are rejected -- linearise / split
// upstream, same contract as the R package.
//
// This exists so the core library has zero dependencies. Bindings that
// already have a geometry engine in hand (GEOS in R via libgeos,
// shapely in Python) are free to bypass it and construct Geometry
// directly.
//
// Copyright (c) 2025 Michael Sumner
// Licensed under Apache License 2.0

#ifndef CONTROLLEDBURN_WKB_HPP
#define CONTROLLEDBURN_WKB_HPP

#include <cstdint>
#include <cstddef>
#include <stdexcept>

#include "geometry.hpp"

namespace controlledburn {

struct WKBParseError : std::runtime_error {
    using std::runtime_error::runtime_error;
};

// Parse a single WKB geometry. Throws WKBParseError on malformed input
// or unsupported type (curves, GeometryCollection).
Geometry parse_wkb(const uint8_t* data, size_t size);

namespace detail {

// What the threaded burn's work estimate needs from one blob, read
// without building the geometry: `n`, its number of coordinates, and
// [lo, hi], the y-range of at most 16 coordinates of each part (a
// polygon's exterior ring; a multi-point's points are one part), the
// coordinates the estimate samples from a parsed geometry. False for a
// blob that does not parse or that the burn rejects.
bool wkb_sample_extent(const uint8_t* data, size_t size, size_t& n, double& lo,
                       double& hi) noexcept;

} // namespace detail

} // namespace controlledburn

#endif // CONTROLLEDBURN_WKB_HPP
