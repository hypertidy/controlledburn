# libcontrolledburn

Pure C++17 core of controlledburn: scanline polygon/line/point
rasterization with exact coverage fractions, O(perimeter) in time and
memory, producing sparse output. No GEOS, no Armadillo, no R, and no
materialized raster.

## Provenance

The walker (`walk_polyline`), ring processing (`walk_ring`), and the
winding sweep are carried over verbatim from the R package's
`src/scanline_burn.cpp`. Only the boundary layer changed:

| R package (src/)                     | cpp/ core                              |
|--------------------------------------|----------------------------------------|
| GEOS WKB reader                      | `wkb.cpp` (~200-line zero-dep reader)  |
| GEOS ring/coord accessors            | `Geometry` structs (`geometry.hpp`)    |
| `geos_is_ccw`                        | shoelace `signed_area` sign            |
| `geos_get_component_boxes`           | coordinate min/max (`ring_envelope`)   |
| `cpp11::stop` / `cpp11::warning`     | exceptions + `BurnResult::notes`       |
| cpp11 data.frame construction        | plain struct vectors (`output.hpp`)    |

The GEOS-free subset of Daniel Baston's exactextract (box, coordinate,
crossing, grid, measures, perimeter_distance, side, traversal,
traversal_areas) is vendored under `src/ee/` together with
`analytical_coverage.h`. The GEOS-coupled parts of exactextract
(geos_utils, cell, floodfill, raster_cell_intersection) are not needed
by this engine and are not vendored.

Note the orientation sign conventions differ deliberately:
`controlledburn::signed_area` (public, `geometry.hpp`) is standard
shoelace with CCW positive; `polygon_signed_area` in
`ee/analytical_coverage.h` uses the opposite convention internal to the
coverage math.

## API

```cpp
#include <controlledburn/controlledburn.hpp>
using namespace controlledburn;

GridSpec grid{xmin, ymin, xmax, ymax, ncol, nrow};

// From geometry structs:
BurnResult r = burn(geometries, grid);

// Or from WKB blobs (ISO and EWKB, Z/M skipped, both byte orders):
BurnResult r = burn_wkb(wkb_spans, grid);

// r.runs   -- interior RLE:      (row, col_start, col_end, id); col_end exclusive
// r.edges  -- boundary cells:    (row, col, fraction, id), fraction in [0, 1]
// r.lines  -- line cells:        (row, col, length, id), CRS units
// r.points -- point cells:       (row, col, id)
// r.notes  -- non-fatal problems (parse failures etc.), per input index
```

All indices are 0-based, row 0 is the top row, and `col_end` is
exclusive: a run covers the half-open column range `[col_start, col_end)`.
Geometry k (0-based input position) gets id = k. A binding that wants a
different convention adds its own offset -- the R package's cpp11 shim
adds 1 to row/col/col_start/id (leaving `col_end`) to restore its
1-based, inclusive contract.

`GeometryCollection` is rejected (mixed dimensions break weight
semantics -- split upstream) and curved types must be linearised by the
caller, same contract as the R package.

## fasterize reconciliation

`materialize.hpp` provides the bridge to fasterize semantics: an
optional pass that burns a sparse `BurnResult` into a caller-provided
pixel buffer with fasterize's pixel functions (first, last, sum, min,
max, count, any).

One engine, two consumption styles: keep the sparse tables
(controlledburn), or materialize them (fasterize).

An honest semantic caveat: fasterize burns a cell when the cell CENTER
is inside the polygon; controlledburn computes exact coverage. Interior
cells agree; boundary cells differ. `MaterializeOptions::edge_policy`
offers `Threshold` (include boundary cell iff fraction >= 0.5, which
approximates but does not equal the center rule) and `Fraction`
(value weighted by coverage -- exact, area-conserving). Exact
center-rule parity would need a center-point winding test per boundary
cell; that is a candidate for a dedicated sweep mode if drop-in
fasterize output equality is ever required.

## Build

```sh
cmake -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build
ctest --test-dir build
```

Requires only a C++17 compiler and CMake >= 3.16.

## Tests

`tests/test_burn.cpp` asserts the core invariant -- total burned
coverage (runs + edge fractions, times cell area) equals exact polygon
area -- across aligned and offset rectangles, triangles on awkward
grids, polygons with holes (both hole orientations), disjoint
multipolygons, beyond-extent polygons (the padding-column winding
case), grid-straddling polygons, line length conservation, point
binning, WKB round trips (including malformed input), materialization
under both edge policies, and degenerate input.

## Bindings (planned)

- R: the existing package becomes a thin cpp11 shim over `burn_wkb`
  (or keeps GEOS-based parsing via libgeos and constructs `Geometry`
  directly -- the core does not care where coordinates come from).
- Python: pybind11/nanobind over the same API, returning numpy arrays
  or Arrow tables of the four-table contract.
