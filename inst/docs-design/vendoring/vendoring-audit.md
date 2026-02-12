# exactextract Vendoring Audit — Definitive Results

## Verdict: VERY CLEAN. Go ahead.

The include graph from `raster_cell_intersection.h` is a tight, self-contained
cluster with **zero GDAL**, **zero R**, **zero stats framework** dependencies.

## Files to Vendor: 27 files, 3159 lines

### 12 compilation units (.cpp)

| File | Lines | Role |
|------|-------|------|
| `raster_cell_intersection.cpp` | 433 | **Entry point**: Grid + GEOSGeometry → coverage matrix |
| `cell.cpp` | 177 | Cell traversal and covered fraction computation |
| `traversal.cpp` | 80 | Ring edge traversal state machine |
| `traversal_areas.cpp` | 115 | Area computation during traversal |
| `floodfill.cpp` | 41 | Flood-fill interior cells to coverage = 1.0 |
| `geos_utils.cpp` | 242 | GEOS C API helpers (version-conditional) |
| `measures.cpp` | 53 | Geometric measures (area, length of coord vectors) |
| `perimeter_distance.cpp` | 59 | Perimeter distance computation (used by traversal_areas) |
| `coordinate.cpp` | 22 | Coordinate ostream operator |
| `side.cpp` | 41 | Side enum ostream operator |
| `box.cpp` | 137 | Bounding box operations |
| `grid.cpp` | 100 | Grid cell computation helpers |

### 15 headers

| File | Lines | Kind |
|------|-------|------|
| `raster_cell_intersection.h` | 62 | Main class + free functions |
| `cell.h` | 77 | Cell class declaration |
| `traversal.h` | 65 | Traversal class declaration |
| `traversal_areas.h` | 27 | TraversalAreas declaration |
| `floodfill.h` | 134 | FloodFill class + template impl |
| `geos_utils.h` | 110 | GEOS smart pointers + helpers |
| `grid.h` | 359 | Grid<T> template (header-heavy) |
| `matrix.h` | 148 | Matrix<T> template (header-only) |
| `raster.h` | 333 | Raster<T> template (header-only) |
| `measures.h` | 32 | area(), length() declarations |
| `coordinate.h` | 52 | Coordinate struct {double x, y} |
| `crossing.h` | 43 | Crossing struct (header-only) |
| `side.h` | 33 | Side enum (TOP/BOTTOM/LEFT/RIGHT) |
| `box.h` | 153 | Box class (xmin,ymin,xmax,ymax) |
| `perimeter_distance.h` | 31 | PerimeterDistance declaration |

## Files NOT needed (29 files we skip)

- `exactextract.cpp` — CLI main
- `gdal_dataset_wrapper.h/.cpp` — GDAL dataset I/O
- `gdal_raster_wrapper.h/.cpp` — GDAL raster I/O
- `gdal_writer.h/.cpp` — GDAL output writer
- `processor.h/.cpp` — Abstract processing pipeline
- `raster_sequential_processor.h/.cpp` — Raster processing
- `feature_sequential_processor.h/.cpp` — Feature processing
- `output_writer.h/.cpp` — Output abstraction
- `raster_source.h` — Abstract raster source
- `raster_stats.h` — Statistics framework
- `stats_registry.h` — Stats registry
- `operation.h` — Operation abstraction
- `weighted_quantiles.h/.cpp` — Quantile computation
- `variance.h` — Variance computation
- `raster_area.h` — Raster area computation
- `utils.h/.cpp` — General utilities (only used by CLI + GDAL writer)
- `version.h.in` — Version template

## GEOS Interface: 5 files need `geos_c.h` patching

```
geos_utils.h       — #include <geos_c.h>  →  #include "libgeos.h"
geos_utils.cpp     — (no direct include, gets it via geos_utils.h)
raster_cell_intersection.h  — #include <geos_c.h>
raster_cell_intersection.cpp — #include <geos_c.h>
floodfill.h        — #include <geos_c.h>
floodfill.cpp      — #include <geos_c.h>
```

That's it. 5 lines to change. The `#include <geos_c.h>` becomes `#include "libgeos.h"`.

## GEOS Version Compatibility: Already handled

The code has **built-in version guards** for GEOS 3.7, 3.8, and 3.10:

```cpp
#define HAVE_370 (GEOS_VERSION_MAJOR >= 3 && GEOS_VERSION_MINOR >= 7)
#define HAVE_380 (GEOS_VERSION_MAJOR >= 3 && GEOS_VERSION_MINOR >= 8)
#define HAVE_3100 (GEOS_VERSION_MAJOR >= 3 && GEOS_VERSION_MINOR >= 10)
```

With libgeos 3.11.1, **all three fast paths activate**:
- `GEOSSegmentIntersection_r` (3.7+) — fast segment intersection
- `GEOSGeom_getXMin_r` etc (3.7+) — fast bbox extraction
- `GEOSCoordSeq_isCCW_r` (3.7+) — fast CCW test
- `GEOSGeom_createPointFromXY_r` (3.8+) — fast point creation
- `GEOSCoordSeq_getXY_r` (3.8+) — fast coordinate access
- `GEOSCoordSeq_copyToBuffer_r` (3.10+) — bulk coordinate copy

Each has a fallback for older GEOS versions, so even if libgeos stayed at an older version, it would still work — just slower.

## GEOS Functions Actually Used (42 functions)

All from the reentrant `_r` API. All available in GEOS 3.5+ (with version-conditional
fast paths for newer versions). libgeos exports the complete C API, so no issues.

## Potential CRAN Issues: Minimal

- `<iostream>` included in `coordinate.h`, `side.h`, `traversal_areas.cpp` — for `operator<<` debug overloads. CRAN allows this; it's actual `std::cout`/`std::cerr` writes that are problematic.
- All `std::cout` calls in the code are **already commented out** (debug traces).
- No R API calls, no Rcpp calls in the vendored code.
- Apache 2.0 license — needs `inst/COPYRIGHTS` attribution.

## The Call Site for gridburn

```cpp
#include "libgeos.h"               // from libgeos package
#include "libgeos.c"               // function pointer init (once)
#include "raster_cell_intersection.h"  // from vendored exactextract

// In .onLoad or init:
libgeos_init_api();

// Per geometry:
exactextract::Grid<exactextract::bounded_extent> grid(
    exactextract::Box(xmin, ymin, xmax, ymax),  // raster extent
    dx, dy);                                      // cell size

exactextract::RasterCellIntersection rci(grid, geos_context, geom_ptr);

// rci.results() → Matrix<float> of coverage fractions
// rci.m_geometry_grid → the subgrid (clipped to geometry bbox)
// rci.rows(), rci.cols() → dimensions of the subgrid result

// Then our dense_to_sparse() converts to runs + edges tables
```

Note: Grid constructor takes `(Box, dx, dy)` not `(Box, nx, ny)`. So from R:
- `dx = (xmax - xmin) / ncol`
- `dy = (ymax - ymin) / nrow`

## Accessing the Geometry from geos package

The geos R package stores `GEOSGeometry*` as external pointers. We need to extract
these in C++. The geos package doesn't currently export a C-level API for this, but:

1. We can accept WKB bytes (raw vectors) from R and deserialize in C++ via GEOS — same pattern as exactextractr
2. Or we can look at how geos stores the external pointer and access it directly

Option 1 is safest and doesn't couple us to geos package internals. The WKB
round-trip is fast (no geometry computation, just serialization).

## Summary

This is a textbook vendoring job:
- Clean dependency boundary (zero GDAL/R/stats contamination)
- 5 include path edits (`geos_c.h` → `libgeos.h`)
- All GEOS version compatibility already built in
- 3159 total lines — manageable, well-tested upstream code
- The same 27 files are already vendored in exactextractr on CRAN
