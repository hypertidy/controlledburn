// burn_triangle.cpp -- minimal usage example
//
// Build (after building the library):
//   c++ -std=c++17 -I../include burn_triangle.cpp -L../build -lcontrolledburn -o burn_triangle
//
// Copyright (c) 2025 Michael Sumner
// Licensed under Apache License 2.0

#include <controlledburn/controlledburn.hpp>
#include <cstdio>

int main() {
    using namespace controlledburn;

    GridSpec grid{0, 0, 10, 10, 10, 10};

    Geometry tri;
    tri.kind = GeomKind::Polygon;
    tri.polygons.push_back({{{{1.2, 1.7}, {8.6, 2.9}, {4.1, 8.3}, {1.2, 1.7}}}});

    BurnResult r = burn({tri}, grid);

    std::printf("runs (interior, fully covered):\n");
    for (const auto& run : r.runs) {
        std::printf("  row %d cols %d..%d id %d\n",
                    run.row, run.col_start, run.col_end, run.id);
    }

    std::printf("edges (boundary, coverage fraction):\n");
    for (const auto& e : r.edges) {
        std::printf("  row %d col %d fraction %.4f id %d\n",
                    e.row, e.col, e.fraction, e.id);
    }

    double cell = grid.dx() * grid.dy();
    double total = 0.0;
    for (const auto& run : r.runs) total += cell * (run.col_end - run.col_start + 1);
    for (const auto& e : r.edges) total += cell * e.fraction;
    double exact = -signed_area(tri.polygons[0].rings[0]);
    std::printf("covered area %.6f, exact polygon area %.6f\n",
                total, exact < 0 ? -exact : exact);
    return 0;
}
