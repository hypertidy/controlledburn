# Vendoring exactextract C++ Core into gridburn

## The Structure

exactextractr (the R package) vendors the upstream [exactextract](https://github.com/isciences/exactextract) C++ library as a git submodule at `src/exactextract/`. The CRAN tarball flattens this to regular files. The upstream library has ~50 files in `src/exactextract/src/`, but the **coverage fraction algorithm** is a self-contained subset.

## What We Need: ~14 Core Files

These files form the computational core — the ring traversal + flood fill algorithm that computes exact cell coverage fractions. They depend only on GEOS C API, nothing else.

### Must vendor (the algorithm)

| File | Role |
|------|------|
| `raster_cell_intersection.h/.cpp` | **THE** entry point. Takes Grid + GEOSGeometry, returns coverage matrix |
| `traversal.h/.cpp` | Ring edge traversal — walks polygon edges across grid cells |
| `traversal_areas.h/.cpp` | Computes partial areas during traversal |
| `floodfill.h/.cpp` | Identifies fully interior cells (coverage = 1.0) |
| `cell.h/.cpp` | Cell geometry representation |
| `coordinate.h/.cpp` | Coordinate types |
| `crossing.h` | Edge crossing types (header only) |
| `side.h/.cpp` | Cell side enumeration (TOP/BOTTOM/LEFT/RIGHT) |
| `box.h/.cpp` | Bounding box |
| `grid.h/.cpp` | Grid definition (extent + dimensions) — this is the "what raster" type |
| `matrix.h` | Dense matrix template (header only) |
| `raster.h` | Raster value template (header only) |
| `geos_utils.h/.cpp` | GEOS helper functions (ring extraction, coordinate access) |
| `utils.h/.cpp` | General utilities |

**Total: ~14 headers, ~11 .cpp files, ~25 files**

### Skip entirely (application layer)

- `exactextract.cpp` — CLI main
- `gdal_*.h/.cpp` — All GDAL I/O (5 files)
- `*_sequential_processor.*` — Processing pipeline (4 files)
- `processor.*` — Abstract processor
- `output_writer.*` — Output abstraction
- `raster_source.h` — Abstract raster source
- `raster_stats.h`, `stats_registry.h` — Statistics framework
- `weighted_quantiles.*`, `variance.h` — Aggregation math
- `measures.*`, `perimeter_distance.*`, `raster_area.h` — Geometry measures
- `operation.h` — Operation abstraction
- `vend/optional.hpp` — std::optional polyfill (we require C++17 anyway)
- All test files

## The Key Interface

The call surface for gridburn is remarkably clean:

```cpp
#include "raster_cell_intersection.h"
#include "grid.h"

// 1. Define a grid from extent + dimensions
exactextract::Grid<exactextract::bounded_extent> grid(
    exactextract::Box(xmin, ymin, xmax, ymax), nx, ny);

// 2. Compute coverage fractions for one GEOSGeometry*
exactextract::RasterCellIntersection rci(grid, geom);

// 3. Access the dense result matrix
// rci covers a sub-grid (clipped to geometry bbox)
// rows/cols map to the sub-grid extent
```

That's it. The class does bbox clipping internally (only processes the subgrid intersecting the geometry), runs the traversal, flood-fills the interior, and gives you back a dense matrix of floats.

## How exactextractr Bridges to R (for reference)

exactextractr's `src/coverage_fraction.cpp`:
- Receives WKB bytes from R (via sf)
- Creates GEOSGeometry from WKB using GEOS C API
- Constructs Grid from R raster extent/resolution
- Calls `RasterCellIntersection(grid, geom)`
- Copies result matrix into an R numeric matrix
- Returns to R

exactextractr's `src/geos_r.h`:
- Manages GEOS context handle lifecycle in R
- Links to system GEOS via `configure` script + `geos-config`

## How gridburn Would Differ

| Aspect | exactextractr | gridburn |
|--------|---------------|----------|
| C++ framework | Rcpp | cpp11 |
| GEOS source | System GEOS via configure | libgeos (LinkingTo) |
| Geometry input | WKB from sf | `geos_geometry` from geos package |
| Output | Dense R matrix (full raster) | Two-table sparse (runs + edges) |
| Scope | Zonal stats engine | Coverage fractions only |

### GEOS Linking: libgeos vs system GEOS

This is the key architectural difference. exactextractr uses system GEOS:
```
SystemRequirements: GEOS (>= 3.5.0)
```
with a `configure` script that finds `geos-config`.

gridburn would use libgeos:
```
LinkingTo: libgeos, cpp11
```
No configure script needed. But we need the GEOS C API function pointers that libgeos exports. The vendored exactextract code calls GEOS functions directly (e.g., `GEOSGetExteriorRing_r()`), which works because libgeos exports the entire C API.

**Potential friction:** The vendored code includes `geos_c.h` directly. With libgeos, we'd include `libgeos.h` instead and call `libgeos_init_api()` at package load. The function signatures are identical — it's the same API — but the include path and init pattern differ. This means light patching of `geos_utils.h` to use `#include "libgeos.h"` instead of `#include <geos_c.h>`.

### Geometry Input

exactextractr gets `GEOSGeometry*` by deserializing WKB bytes that R passes as raw vectors (from sf). We'd get `GEOSGeometry*` from geos package's external pointers. The geos package stores geometries as `GEOSGeometry*` behind R external pointers, accessible via its C API.

This is actually simpler — no WKB serialization/deserialization round-trip.

## Feasibility Assessment

### What makes this tractable

1. **Clean separation.** The coverage algorithm is a self-contained cluster with no GDAL, no R, no stats dependencies. Just geometry + grid → coverage matrix.

2. **Battle-tested.** exactextractr is on CRAN, heavily used, edge cases hammered out over years.

3. **Apache 2.0.** No license friction.

4. **Small surface.** ~25 files, ~14 headers. The code is well-structured C++14.

5. **Precedent.** exactextractr already vendors this code for CRAN. We know it compiles on all CRAN platforms.

### What needs attention

1. **libgeos include path.** Replace `#include <geos_c.h>` with `#include "libgeos.h"` in geos_utils.h. Minor patch.

2. **GEOS context handle.** The vendored code uses `GEOSContextHandle_t` throughout (reentrant API). libgeos provides this. Need to thread the handle through from R → cpp11 → vendored code. exactextractr already does this (see `geos_r.h`).

3. **C++ standard.** The vendored code is C++14. cpp11 requires C++11. No conflict. We could target C++17 for our own code while the vendored code stays C++14.

4. **Makevars.** Need to list all vendored .cpp files in `PKG_CXXFLAGS` include paths and `OBJECTS`. exactextractr's Makevars.in is the template.

5. **Namespace.** The vendored code lives in `namespace exactextract`. No collision risk, but we might want to nest it or leave as-is.

## dense_to_sparse in C++

The second piece we need in C++ (beyond the vendored core) is the dense→sparse conversion. The vendored code gives us a dense matrix of coverage fractions per geometry. We need to:

1. Iterate the dense matrix row by row
2. Separate cells with weight ≈ 1.0 (interior → runs table) from cells with 0 < weight < 1.0 (edges → edges table)
3. RLE-compress interior cells into runs (row, col_start, col_end, id)
4. Collect edge cells as individual entries (row, col, weight, id)

This is straightforward C++ (~100 lines) that we write ourselves. It sits between the vendored core and the R return value.

## Recommended File Layout

```
gridburn/
├── DESCRIPTION          # LinkingTo: libgeos, cpp11
├── R/
│   ├── burn_sparse.R    # Main R entry point
│   └── zzz.R            # .onLoad → init libgeos API
├── src/
│   ├── gridburn.cpp     # cpp11-registered functions
│   ├── dense_to_sparse.h/.cpp  # Our code: matrix → two tables
│   ├── Makevars         # Include paths, object list
│   └── exactextract/    # Vendored subset
│       ├── raster_cell_intersection.h/.cpp
│       ├── traversal.h/.cpp
│       ├── traversal_areas.h/.cpp
│       ├── floodfill.h/.cpp
│       ├── cell.h/.cpp
│       ├── coordinate.h/.cpp
│       ├── crossing.h
│       ├── side.h/.cpp
│       ├── box.h/.cpp
│       ├── grid.h/.cpp
│       ├── matrix.h
│       ├── raster.h
│       ├── geos_utils.h/.cpp  # patched: libgeos.h include
│       └── utils.h/.cpp
└── inst/
    └── COPYRIGHTS       # Apache 2.0 attribution for exactextract
```

## Risk: What If It's Harder Than Expected?

The main risk is **hidden dependencies** in the vendored files — e.g., `raster_cell_intersection.cpp` pulling in something from the stats layer. The mitigation is:

1. Start by copying just the 14 core files
2. Try to compile `raster_cell_intersection.cpp` with libgeos
3. Chase any missing includes — if they pull in GDAL or stats code, we know we have a problem
4. If the dependency graph is wider than expected, we can always fall back to linking against system GEOS and using the GEOS 3.14 C API directly (which is a single function call, no vendoring needed)

## Comparison: Vendor vs Wait for libgeos 3.14

| | Vendor exactextract core | Wait for libgeos 3.14 |
|---|---|---|
| **Works now** | Yes, with GEOS 3.11 | No, need Dewey to bump |
| **CRAN-ready** | Yes, follows exactextractr pattern | Blocked on libgeos release |
| **Maintenance** | Must track upstream fixes | Zero — it's in GEOS |
| **Code size** | ~25 vendored files | 1 function call |
| **GEOS jump** | 3.11→3.14 is 3 major versions, C++17 required | Same issue for libgeos |
| **Future swap** | Drop vendored code when libgeos catches up | N/A |

**Verdict: Vendor now, swap later.** The vendoring surface is clean, the precedent exists (exactextractr already does it), and we can design the R API so that swapping backends is a one-file change in C++.
