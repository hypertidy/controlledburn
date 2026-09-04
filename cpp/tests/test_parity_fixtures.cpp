// test_parity_fixtures.cpp -- cross-language parity tests driven by shared
// CSV fixtures in fixtures/. Identical fixtures are read by R and Python
// test suites so any regression is caught across all three surfaces.

#include "controlledburn/controlledburn.hpp"

#include <cmath>
#include <cstdio>
#include <cstring>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>

using namespace controlledburn;

static int failures = 0;

#define CHECK(cond) do { \
    if (!(cond)) { \
        std::printf("FAIL %s:%d: %s\n", __FILE__, __LINE__, #cond); \
        failures++; \
    } \
} while (0)

#define CHECK_NEAR(a, b, tol) do { \
    double aa = (a), bb = (b); \
    if (bb != 0 ? (std::fabs(aa - bb) / std::fabs(bb) > (tol)) \
                : (std::fabs(aa - bb) > (tol))) { \
        std::printf("FAIL %s:%d: %s (%.12g) != %s (%.12g)\n", \
                    __FILE__, __LINE__, #a, aa, #b, bb); \
        failures++; \
    } \
} while (0)

// --- Minimal CSV parser (no quoting edge cases beyond double-quote fields) ---

static std::string unquote(const std::string& s) {
    if (s.size() >= 2 && s.front() == '"' && s.back() == '"')
        return s.substr(1, s.size() - 2);
    return s;
}

static std::vector<std::string> split_csv_line(const std::string& line) {
    std::vector<std::string> fields;
    std::string field;
    bool in_quote = false;
    for (size_t i = 0; i < line.size(); i++) {
        char c = line[i];
        if (c == '"') { in_quote = !in_quote; field += c; }
        else if (c == ',' && !in_quote) { fields.push_back(field); field.clear(); }
        else { field += c; }
    }
    fields.push_back(field);
    return fields;
}

struct FixtureGeom {
    std::string case_name;
    std::string wkb_hex;
    double xmin, ymin, xmax, ymax;
    int ncol, nrow;
};

struct FixtureExpected {
    std::string case_name;
    double covered_area;  // NAN if "NA"
    int edges_empty;      // -1 if "NA"
    double line_length;   // NAN if "NA"
    int n_points;         // -1 if "NA"
    double tol_rel;       // NAN if "NA"
};

static double parse_double_or_nan(const std::string& s) {
    if (s == "NA" || s.empty()) return std::nan("");
    return std::stod(s);
}

static int parse_int_or_neg(const std::string& s) {
    if (s == "NA" || s.empty()) return -1;
    return std::stoi(s);
}

static std::vector<uint8_t> hex_to_bytes(const std::string& hex) {
    std::vector<uint8_t> bytes;
    for (size_t i = 0; i + 1 < hex.size(); i += 2) {
        unsigned val;
        std::sscanf(hex.c_str() + i, "%2x", &val);
        bytes.push_back(static_cast<uint8_t>(val));
    }
    return bytes;
}

static double covered_area(const BurnResult& r, const GridSpec& gs) {
    double cell = gs.dx() * gs.dy();
    double total = 0.0;
    for (const auto& run : r.runs)
        total += cell * (run.col_end - run.col_start);
    for (const auto& e : r.edges)
        total += cell * e.fraction;
    return total;
}

int main(int argc, char** argv) {
    // Determine fixtures directory: either passed as argv[1] or
    // relative to the source tree (when run from build/).
    std::string fixtures_dir = "../../fixtures";
    if (argc > 1) fixtures_dir = argv[1];

    // Load geometries
    std::ifstream gf(fixtures_dir + "/geometries.csv");
    if (!gf.is_open()) {
        std::printf("SKIP: cannot open %s/geometries.csv\n", fixtures_dir.c_str());
        return 0;
    }

    std::string header;
    std::getline(gf, header);
    auto hdr_fields = split_csv_line(header);

    // Find column indices
    int ci_case = -1, ci_wkb = -1, ci_xmin = -1, ci_ymin = -1;
    int ci_xmax = -1, ci_ymax = -1, ci_ncol = -1, ci_nrow = -1;
    for (size_t i = 0; i < hdr_fields.size(); i++) {
        std::string h = unquote(hdr_fields[i]);
        if (h == "case") ci_case = i;
        else if (h == "wkb_hex") ci_wkb = i;
        else if (h == "xmin") ci_xmin = i;
        else if (h == "ymin") ci_ymin = i;
        else if (h == "xmax") ci_xmax = i;
        else if (h == "ymax") ci_ymax = i;
        else if (h == "ncol") ci_ncol = i;
        else if (h == "nrow") ci_nrow = i;
    }

    std::vector<FixtureGeom> geoms;
    std::string line;
    while (std::getline(gf, line)) {
        if (line.empty()) continue;
        auto f = split_csv_line(line);
        FixtureGeom g;
        g.case_name = unquote(f[ci_case]);
        g.wkb_hex = unquote(f[ci_wkb]);
        g.xmin = std::stod(f[ci_xmin]);
        g.ymin = std::stod(f[ci_ymin]);
        g.xmax = std::stod(f[ci_xmax]);
        g.ymax = std::stod(f[ci_ymax]);
        g.ncol = std::stoi(f[ci_ncol]);
        g.nrow = std::stoi(f[ci_nrow]);
        geoms.push_back(g);
    }

    // Load expected
    std::ifstream ef(fixtures_dir + "/expected.csv");
    std::getline(ef, header);
    auto ehdr = split_csv_line(header);
    int ei_case = -1, ei_area = -1, ei_edges_empty = -1;
    int ei_line_len = -1, ei_npoints = -1, ei_tol = -1;
    for (size_t i = 0; i < ehdr.size(); i++) {
        std::string h = unquote(ehdr[i]);
        if (h == "case") ei_case = i;
        else if (h == "covered_area") ei_area = i;
        else if (h == "edges_empty") ei_edges_empty = i;
        else if (h == "line_length") ei_line_len = i;
        else if (h == "n_points") ei_npoints = i;
        else if (h == "tol_rel") ei_tol = i;
    }

    std::vector<FixtureExpected> expected;
    while (std::getline(ef, line)) {
        if (line.empty()) continue;
        auto f = split_csv_line(line);
        FixtureExpected e;
        e.case_name = unquote(f[ei_case]);
        e.covered_area = parse_double_or_nan(f[ei_area]);
        e.edges_empty = parse_int_or_neg(f[ei_edges_empty]);
        e.line_length = parse_double_or_nan(f[ei_line_len]);
        e.n_points = parse_int_or_neg(f[ei_npoints]);
        e.tol_rel = parse_double_or_nan(f[ei_tol]);
        expected.push_back(e);
    }

    // Run tests
    for (const auto& g : geoms) {
        // Find expected
        const FixtureExpected* exp = nullptr;
        for (const auto& e : expected) {
            if (e.case_name == g.case_name) { exp = &e; break; }
        }
        if (!exp) {
            std::printf("WARN: no expected values for case '%s'\n",
                        g.case_name.c_str());
            continue;
        }

        std::printf("  %s ... ", g.case_name.c_str());

        auto wkb_bytes = hex_to_bytes(g.wkb_hex);
        GridSpec gs{g.xmin, g.ymin, g.xmax, g.ymax, g.ncol, g.nrow};
        BurnResult r = burn_wkb({{wkb_bytes.data(), wkb_bytes.size()}}, gs);

        int before = failures;

        if (!std::isnan(exp->covered_area)) {
            CHECK_NEAR(covered_area(r, gs), exp->covered_area, exp->tol_rel);
        }

        if (exp->edges_empty == 1) {
            CHECK(r.edges.empty());
        }

        if (!std::isnan(exp->line_length)) {
            double total = 0;
            for (const auto& l : r.lines) total += l.length;
            CHECK_NEAR(total, exp->line_length, exp->tol_rel);
        }

        if (exp->n_points >= 0) {
            CHECK(static_cast<int>(r.points.size()) == exp->n_points);
        }

        std::printf("%s\n", failures == before ? "ok" : "FAIL");
    }

    if (failures == 0) {
        std::printf("all parity fixture tests passed\n");
        return 0;
    }
    std::printf("%d failure(s)\n", failures);
    return 1;
}
