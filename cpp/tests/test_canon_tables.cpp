// test_canon_tables.cpp -- record-for-record agreement with the Rust canon.
//
// hypertidy/controlledburn-rs is the canonical implementation of the burn
// algorithm; this C++ core follows it. fixtures/canon-tables.txt is
// written by controlledburn-rs (tests/canon.rs): lattice-biased rings,
// holes, multipolygons, lines and points on non-dyadic grids, with the
// tables the Rust crate emits in parity mode (no run coalescing) for both
// modes. This test burns the same WKB and requires the same records in the
// same order, fractions and lengths compared as bits. See the header of
// tests/canon.rs in controlledburn-rs for the format and for how to
// regenerate the file after a deliberate change.
//
// Copyright (c) 2026 Michael Sumner
// Licensed under Apache License 2.0

#include "controlledburn/controlledburn.hpp"

#include <cstdint>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace controlledburn;

namespace {

struct Mode {
    std::string name;
    std::vector<std::string> records;
};

struct Case {
    std::string name;
    GridSpec grid{};
    std::vector<std::vector<uint8_t>> wkb;
    std::vector<Mode> modes;
};

double from_bits(const std::string& h) {
    uint64_t u = std::stoull(h, nullptr, 16);
    double d;
    std::memcpy(&d, &u, sizeof d);
    return d;
}

std::vector<uint8_t> unhex(const std::string& h) {
    std::vector<uint8_t> out(h.size() / 2);
    for (size_t i = 0; i < out.size(); i++)
        out[i] = static_cast<uint8_t>(std::stoul(h.substr(2 * i, 2), nullptr, 16));
    return out;
}

std::string f32_hex(float f) {
    uint32_t u;
    std::memcpy(&u, &f, sizeof u);
    char buf[16];
    std::snprintf(buf, sizeof buf, "%08x", u);
    return buf;
}

std::vector<std::string> records(const BurnResult& r) {
    std::vector<std::string> out;
    for (const auto& x : r.runs)
        out.push_back("run " + std::to_string(x.row) + " " + std::to_string(x.col_start) + " " +
                      std::to_string(x.col_end) + " " + std::to_string(x.id));
    for (const auto& x : r.edges)
        out.push_back("edge " + std::to_string(x.row) + " " + std::to_string(x.col) + " " +
                      f32_hex(x.fraction) + " " + std::to_string(x.id));
    for (const auto& x : r.lines)
        out.push_back("line " + std::to_string(x.row) + " " + std::to_string(x.col) + " " +
                      f32_hex(x.length) + " " + std::to_string(x.id));
    for (const auto& x : r.points)
        out.push_back("point " + std::to_string(x.row) + " " + std::to_string(x.col) + " " +
                      std::to_string(x.id));
    return out;
}

std::vector<Case> read_cases(const std::string& path) {
    std::ifstream in(path);
    if (!in) {
        std::printf("cannot open %s\n", path.c_str());
        return {};
    }
    std::vector<Case> cases;
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        std::istringstream ss(line);
        std::string tag;
        ss >> tag;
        if (tag == "case") {
            cases.emplace_back();
            ss >> cases.back().name;
        } else if (tag == "grid") {
            std::string a, b, c, d;
            int nc, nr;
            ss >> a >> b >> c >> d >> nc >> nr;
            cases.back().grid = GridSpec{from_bits(a), from_bits(b), from_bits(c), from_bits(d), nc, nr};
        } else if (tag == "wkb") {
            std::string h;
            ss >> h;
            cases.back().wkb.push_back(unhex(h));
        } else if (tag == "mode") {
            cases.back().modes.emplace_back();
            ss >> cases.back().modes.back().name;
        } else if (tag == "end") {
        } else {
            cases.back().modes.back().records.push_back(line);
        }
    }
    return cases;
}

} // namespace

int main(int argc, char** argv) {
    std::string dir = argc > 1 ? argv[1] : "fixtures";
    auto cases = read_cases(dir + "/canon-tables.txt");
    if (cases.size() < 100) {
        std::printf("FAIL: expected the canonical cases, read %zu\n", cases.size());
        return 1;
    }
    int failures = 0;
    size_t nrec = 0;
    for (const auto& c : cases) {
        std::vector<WKBSpan> spans;
        for (const auto& b : c.wkb) spans.push_back(WKBSpan{b.data(), b.size()});
        for (const auto& m : c.modes) {
            BurnMode mode = m.name == "approx" ? BurnMode::Approx : BurnMode::Coverage;
            BurnResult r = burn_wkb(spans, c.grid, mode);
            auto got = records(r);
            nrec += m.records.size();
            if (!r.notes.empty()) {
                std::printf("FAIL %s [%s]: note: %s\n", c.name.c_str(), m.name.c_str(),
                            r.notes[0].message.c_str());
                failures++;
            }
            if (got != m.records) {
                size_t i = 0;
                while (i < got.size() && i < m.records.size() && got[i] == m.records[i]) i++;
                std::printf("FAIL %s [%s]: %zu records, canon %zu; first difference at %zu: '%s' vs canon '%s'\n",
                            c.name.c_str(), m.name.c_str(), got.size(), m.records.size(), i,
                            i < got.size() ? got[i].c_str() : "(none)",
                            i < m.records.size() ? m.records[i].c_str() : "(none)");
                failures++;
            }
        }
    }
    std::printf("%zu cases, %zu canonical records, %d failures\n", cases.size(), nrec, failures);
    return failures == 0 ? 0 : 1;
}
