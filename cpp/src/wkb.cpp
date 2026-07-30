// wkb.cpp -- minimal WKB reader implementation
//
// Copyright (c) 2025 Michael Sumner
// Licensed under Apache License 2.0

#include "controlledburn/wkb.hpp"

#include <cstring>
#include <string>

namespace controlledburn {

namespace {

class Cursor {
public:
    Cursor(const uint8_t* data, size_t size)
        : data_(data), size_(size), pos_(0), swap_(false) {}

    void set_byte_order(uint8_t order) {
        // 0 = big endian (XDR), 1 = little endian (NDR)
        bool machine_little = is_machine_little_endian();
        swap_ = (order == 1) != machine_little;
    }

    uint8_t read_u8() {
        need(1);
        return data_[pos_++];
    }

    uint32_t read_u32() {
        need(4);
        uint32_t v;
        std::memcpy(&v, data_ + pos_, 4);
        pos_ += 4;
        if (swap_) v = byteswap32(v);
        return v;
    }

    double read_f64() {
        need(8);
        uint64_t v;
        std::memcpy(&v, data_ + pos_, 8);
        pos_ += 8;
        if (swap_) v = byteswap64(v);
        double d;
        std::memcpy(&d, &v, 8);
        return d;
    }

    void skip_f64(size_t n) {
        need(8 * n);
        pos_ += 8 * n;
    }

private:
    void need(size_t n) const {
        if (pos_ + n > size_)
            throw WKBParseError("WKB truncated at byte " + std::to_string(pos_));
    }

    static bool is_machine_little_endian() {
        const uint16_t one = 1;
        uint8_t first;
        std::memcpy(&first, &one, 1);
        return first == 1;
    }

    static uint32_t byteswap32(uint32_t v) {
        return ((v & 0x000000FFu) << 24) |
               ((v & 0x0000FF00u) << 8)  |
               ((v & 0x00FF0000u) >> 8)  |
               ((v & 0xFF000000u) >> 24);
    }

    static uint64_t byteswap64(uint64_t v) {
        return (static_cast<uint64_t>(byteswap32(static_cast<uint32_t>(v))) << 32) |
                static_cast<uint64_t>(byteswap32(static_cast<uint32_t>(v >> 32)));
    }

    const uint8_t* data_;
    size_t size_;
    size_t pos_;
    bool swap_;
};

// EWKB dimension/SRID flags
constexpr uint32_t EWKB_Z    = 0x80000000u;
constexpr uint32_t EWKB_M    = 0x40000000u;
constexpr uint32_t EWKB_SRID = 0x20000000u;

struct TypeInfo {
    uint32_t base;       // 1..7 base geometry type
    size_t extra_dims;   // number of ordinates beyond x, y to skip
    bool has_srid;
};

TypeInfo decode_type(uint32_t raw) {
    TypeInfo t{0, 0, false};
    t.has_srid = (raw & EWKB_SRID) != 0;
    if (raw & EWKB_Z) t.extra_dims++;
    if (raw & EWKB_M) t.extra_dims++;
    uint32_t base = raw & 0x0FFFFFFFu;
    // ISO offsets: 1000 = Z, 2000 = M, 3000 = ZM
    if (base >= 3000 && base < 4000) { t.extra_dims += 2; base -= 3000; }
    else if (base >= 2000 && base < 3000) { t.extra_dims += 1; base -= 2000; }
    else if (base >= 1000 && base < 2000) { t.extra_dims += 1; base -= 1000; }
    t.base = base;
    return t;
}

CoordSeq read_coordseq(Cursor& cur, size_t extra_dims) {
    uint32_t n = cur.read_u32();
    CoordSeq seq;
    seq.reserve(n);
    for (uint32_t i = 0; i < n; i++) {
        double x = cur.read_f64();
        double y = cur.read_f64();
        if (extra_dims) cur.skip_f64(extra_dims);
        seq.push_back({x, y});
    }
    return seq;
}

// Read one point's coordinates as a 1-element sequence.
// WKB POINT EMPTY is conventionally encoded with NaN ordinates; the
// caller filters non-finite coords.
CoordSeq read_point_coords(Cursor& cur, size_t extra_dims) {
    double x = cur.read_f64();
    double y = cur.read_f64();
    if (extra_dims) cur.skip_f64(extra_dims);
    return CoordSeq{{x, y}};
}

Polygon read_polygon_body(Cursor& cur, size_t extra_dims) {
    uint32_t n_rings = cur.read_u32();
    Polygon poly;
    poly.rings.reserve(n_rings);
    for (uint32_t r = 0; r < n_rings; r++) {
        poly.rings.push_back(read_coordseq(cur, extra_dims));
    }
    return poly;
}

// Read the header of a nested geometry (byte order + type + optional
// SRID) and return its TypeInfo, verifying the base type matches
// `expected_base`.
TypeInfo read_nested_header(Cursor& cur, uint32_t expected_base) {
    uint8_t order = cur.read_u8();
    cur.set_byte_order(order);
    TypeInfo t = decode_type(cur.read_u32());
    if (t.has_srid) cur.read_u32();
    if (t.base != expected_base)
        throw WKBParseError("WKB multi-geometry contains mismatched part type " +
                            std::to_string(t.base));
    return t;
}

} // namespace

Geometry parse_wkb(const uint8_t* data, size_t size) {
    if (!data || size < 5)
        throw WKBParseError("WKB too short");

    Cursor cur(data, size);
    uint8_t order = cur.read_u8();
    if (order > 1)
        throw WKBParseError("WKB invalid byte order flag");
    cur.set_byte_order(order);

    TypeInfo t = decode_type(cur.read_u32());
    if (t.has_srid) cur.read_u32(); // skip SRID

    Geometry g;
    switch (t.base) {
    case 1: { // Point
        g.kind = GeomKind::Point;
        g.points.push_back(read_point_coords(cur, t.extra_dims));
        return g;
    }
    case 2: { // LineString
        g.kind = GeomKind::LineString;
        g.lines.push_back(read_coordseq(cur, t.extra_dims));
        return g;
    }
    case 3: { // Polygon
        g.kind = GeomKind::Polygon;
        g.polygons.push_back(read_polygon_body(cur, t.extra_dims));
        return g;
    }
    case 4: { // MultiPoint
        g.kind = GeomKind::MultiPoint;
        uint32_t n = cur.read_u32();
        for (uint32_t i = 0; i < n; i++) {
            TypeInfo pt = read_nested_header(cur, 1);
            g.points.push_back(read_point_coords(cur, pt.extra_dims));
        }
        return g;
    }
    case 5: { // MultiLineString
        g.kind = GeomKind::MultiLineString;
        uint32_t n = cur.read_u32();
        for (uint32_t i = 0; i < n; i++) {
            TypeInfo lt = read_nested_header(cur, 2);
            g.lines.push_back(read_coordseq(cur, lt.extra_dims));
        }
        return g;
    }
    case 6: { // MultiPolygon
        g.kind = GeomKind::MultiPolygon;
        uint32_t n = cur.read_u32();
        for (uint32_t i = 0; i < n; i++) {
            TypeInfo pt = read_nested_header(cur, 3);
            g.polygons.push_back(read_polygon_body(cur, pt.extra_dims));
        }
        return g;
    }
    case 7:
        throw WKBParseError(
            "GeometryCollection is not supported (mixed dimensions break "
            "weight semantics); split into homogeneous groups");
    default:
        throw WKBParseError("WKB unsupported geometry type " +
                            std::to_string(t.base) +
                            " (curved types must be linearised upstream)");
    }
}

} // namespace controlledburn
