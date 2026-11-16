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

## License

Apache License 2.0 (see [../LICENSE.md](../LICENSE.md)). The files under
`src/ee/` (except `analytical_coverage.h`) are unmodified exactextract
sources, Apache-2.0, copyright ISciences, LLC; the scanline walker and
`materialize.hpp` derive from fasterize (MIT, copyright EcoHealth
Alliance). See [../inst/COPYRIGHTS](../inst/COPYRIGHTS).

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
Geometry k (0-based input position) gets id = k + `id_offset`. A binding that wants a
different convention adds its own offset -- the R package's cpp11 shim
adds 1 to row/col/col_start/id (leaving `col_end`) to restore its
1-based, inclusive contract.

### Options and operations (from the Rust crate)

The API follows the controlledburn Rust crate (0.3.2), the canonical
implementation, with the same defaults and the same records
(`tests/test_canon_ops.cpp` checks every operation against the crate's
own answers in `fixtures/canon-ops.txt`).

```cpp
BurnOptions opts;                    // coverage, coalesce_runs, merge_parts, id_offset 0
opts = BurnOptions::parity(BurnMode::Coverage);  // the core's earlier records
opts = BurnOptions(BurnMode::Approx).with_merge_parts(false).with_id_offset(100);

burn(geoms, grid, opts);                          // or burn_wkb
burn_par(geoms, grid, opts, threads);             // same result, several threads
burn_window(geoms, grid, row0, row1, col0, col1, opts);  // only a window's cells

coalesce_runs(r.runs);               // merge adjacent runs
merge_cells(r);                      // one record per (id, cell)
crop(r, grid, xmin, ymin, xmax, ymax, out, subgrid);
slice_tiles(r, grid, tile_ncol, tile_nrow);       // or slice_tiles_par
relate_cells(r, geoms, grid, opts);  // DE-9IM per (polygon, cell), relate.hpp
```

A `BurnMode` converts to the default options with that mode, so
`burn(geoms, grid, BurnMode::Approx)` coalesces runs and merges parts;
use `BurnOptions::parity(mode)` for the records the core gave before these options.

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

`tests/test_canon_tables.cpp` and `tests/test_canon_ops.cpp` check the
core record for record against the Rust crate (see
[../fixtures/README.md](../fixtures/README.md)).

## Bindings

- R: the package in this repository is a cpp11 shim over this core
  (`src/scanline_shim.cpp`, sources synced by `tools/sync-core.sh`).
- Python: built from the Rust crate itself (`python/`), not this core.

The threaded entry points use `std::thread`; link with your platform's
thread library (CMake's `Threads::Threads`).
