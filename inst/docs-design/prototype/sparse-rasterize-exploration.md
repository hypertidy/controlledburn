# Sparse Rasterization Super-Package: Exploration Design

## Vision

A new R package (working name: **gridburn** or **cellweight** — TBD) that produces a **dense tabular, sparse spatial** representation of geometry-to-grid intersection. For any set of geometries (polygons or lines) and a raster specification (extent + dimension), the output is:

```
row | col_start | col_end | weight | geometry_id
```

Where `weight` is either:

- **1.0** for all cells (cheap/approximate mode — controlledburn equivalent)
- **0 < w ≤ 1** exact coverage fraction (exactextract equivalent)

No raster output, no pixel values, no zonal statistics. Just the sparse intersection table.

## Current Landscape

### controlledburn (hypertidy)

- **Algorithm**: Wylie et al. (1967) scanline rasterization via fasterize internals
- **Output**: `(start_col, end_col, row, geometry_id)` — 0-indexed, returned as list of vectors
- **Input**: sf polygons + extent/dimension specification (no raster object needed)
- **Coverage model**: cell-centre intersection only (binary: in or out)
- **Lines**: partially working, known issues
- **C++ binding**: Rcpp
- **Strengths**: extremely fast, zero-materialisation philosophy, clean separation of raster spec from raster data
- **Weaknesses**: no coverage fractions, lines incomplete, Rcpp not cpp11

### exactextract (isciences, C++ library)

- **Algorithm**: ring traversal tracking entry/exit through grid cells, counter-clockwise bounded area computation within each cell, flood fill for interior cells, point-in-polygon for ambiguous cells
- **Output**: per-cell coverage fractions (0–1) within a bounding box sub-grid (the `Raster` template class)
- **Input**: polygons only (GEOS geometry via WKB, or coordinate sequences)
- **Lines**: **not supported** — the algorithm is fundamentally polygon-oriented (rings, bounded areas, flood fill)
- **Key C++ classes**: `Grid` (defines the raster spec), `Box` (bounding box), `RasterCellIntersection` (the core computation), `CoverageWriter` / various stat operations
- **Dependencies**: GEOS required for geometry handling; GDAL optional (CLI only)
- **License**: Apache-2.0
- **Strengths**: proven, fast, accurate coverage fractions, well-tested
- **Weaknesses**: polygon-only, tightly coupled to its own Grid/Box/Raster abstractions

### exactextractr (isciences, R package)

- **Wraps**: the exactextract C++ library (vendored in `src/exactextract/`)
- **Binding**: Rcpp
- **Input**: sf polygons + raster/terra/stars objects
- **Output**: heavily oriented to zonal statistics (summary values per feature), or `coverage_fraction()` returning full RasterLayer per feature
- **Also has**: `rasterize_polygons()` — assigns cells to polygon with greatest coverage area (added ~0.8.x)
- **Near-miss for our use case**: `exact_extract(rast, polys, fun=NULL)` returns per-cell data frames with `value`, `coverage_fraction`, optionally `cell` and `xy` columns. This is *almost* the sparse table we want, but it requires an actual raster input (not just a spec) and returns one data.frame per feature.
- **Architecture**: deeply integrated with raster I/O patterns — not designed for "just give me the sparse table from a grid spec" use case
- **Lines**: not supported (inherits C++ library limitation)

### fasterize (hypertidy, CRAN)

- **Algorithm**: same Wylie et al. scanline as controlledburn
- **Output**: materialised raster (RasterLayer)
- **Now accepts**: any polygon vector (sfc, wkt, wkb, geos) not just sf — broadened recently
- **Lines**: not supported (open issue #1)

### Relationship to other hypertidy packages

- **vaster**: grid logic without data — `cell_from_row_col()`, `xy_from_cell()` etc. Natural companion for interpreting the sparse output
- **vapour**: GDAL API — could provide geometry input via WKB/WKT without sf dependency
- **whatarelief**: raster data acquisition — upstream of any analysis that would consume the sparse table
- **sds**: spatial data sources — could provide raster spec metadata

## exactextract C++ Core Analysis

The core library in `src/` is well-structured for extraction. Key headers:

- **`grid.h`** — `Grid<bounded_extent>` defines raster spec (extent + resolution). This is the equivalent of controlledburn's `(extent, dimension)` pair. Pure geometry, no data.
- **`raster.h`** — `AbstractRaster<T>` wraps a `Grid` with a `Matrix` for storage. We'd use `Grid` but not `AbstractRaster` (we don't want pixel storage).
- **`matrix.h`** — Simple 2D array. Used by `Raster` for coverage fraction storage.
- **`raster_cell_intersection.h`** — **The key class.** Takes a `Grid` and a geometry (GEOS `GEOSGeometry*`), computes coverage fractions for all intersecting cells. Returns results as a `Raster<float>` over the bounding sub-grid.
- **`box.h`** — Bounding box utilities.
- **`traversal_areas.h`** — The actual ring-traversal algorithm: traces a polygon ring through grid cells, computing the bounded area contributions.
- **`floodfill.h`** — Propagates inside/outside status from ring-touched cells to interior cells.
- **`geos_utils.h`** — GEOS geometry helpers (ring iteration, coordinate access). **This creates a hard GEOS dependency.**

**Build options confirm isolation**: `BUILD_CLI=OFF` removes GDAL dependency. The core library needs only GEOS.

**Vendoring assessment**: The core algorithm (~10-15 source files) could be vendored. However, the GEOS dependency for geometry parsing means we either:
1. Accept GEOS as a system dependency (it's already required for sf, which most users have)
2. Replace GEOS geometry parsing with our own WKB parser feeding coordinate sequences
3. Use the `geos` R package's C API to provide geometry pointers

Option 1 is simplest and most pragmatic for CRAN. Option 3 is cleanest for the hypertidy ecosystem (geos package is lightweight). Option 2 is most independent but most work.

**Alternative to vendoring**: Since the algorithm is well-documented (ring traversal → bounded areas → flood fill), we could reimplement in clean C++ without the GEOS coupling. The algorithm itself isn't patented and is described in the literature. This would give us a clean implementation that takes coordinate arrays directly (consistent with how controlledburn works). The tradeoff is test coverage — exactextract has years of edge-case fixes.

## Critical Finding: GEOSGridIntersectionFractions in GEOS 3.14

**GEOS 3.14** (released August 2025) includes `GEOSGridIntersectionFractions` — Dan Baston's exactextract algorithm ported directly into the GEOS C API ([PR #1295](https://github.com/libgeos/geos/pull/1295)). From the release notes:

> "Grid cell overlay with polygon. Beyond the flood-fill, get actual proportions of each cell covered."

The C API signature (from `geos_c.h`):

```c
int GEOSGridIntersectionFractions_r(
    GEOSContextHandle_t handle,
    const GEOSGeometry* g,
    double xmin, double ymin,
    double xmax, double ymax,
    unsigned nx, unsigned ny,
    float* buf);
```

Takes a geometry + grid extent + grid dimensions, fills a pre-allocated `float[nx*ny]` buffer with coverage fractions. Returns a **dense** buffer — our package converts this to sparse.

However:

- **`libgeos` on CRAN** (March 2025) bundles GEOS 3.11.1 — does **not** have this function yet
- **`geos` R package** (December 2025) links to libgeos — also doesn't have access

**Strategy**: This means we have two paths and should design for both:

1. **Immediate (Phase 1)**: Vendor the exactextract C++ core alongside libgeos for the coverage fraction computation. This works now.

2. **Future (when libgeos updates to GEOS ≥ 3.14)**: Drop the vendored code entirely and call `GEOSGridIntersectionFractions` through the standard libgeos C API. The R-level API stays identical.

The package should be architected so that swapping (1) for (2) is a one-file change in the C++ layer. The R API and output format are identical either way.

**Alternative shortcut**: Could also contribute a PR to `libgeos` to bump to GEOS 3.14, which would unblock path (2) immediately. Dewey Dunnington maintains it.

## Key Architectural Decisions

### 1. New package vs extend exactextractr?

**Recommendation: new package.**

Rationale:
- exactextractr is deeply coupled to zonal statistics workflows and raster I/O
- The sparse table output is a fundamentally different API surface
- Dan Baston (exactextractr author) may not want the scope creep
- A thin, focused package aligns with hypertidy philosophy
- We need lines support which exactextract's algorithm can't provide
- We want cpp11, not Rcpp

However, we should **vendor or link** the exactextract C++ core for the polygon coverage fraction computation. It's Apache-2.0, and the core algorithm classes are well-isolated.

### 2. C++ binding: cpp11 ✓

### 3. Geometry input via {geos} + {libgeos}

The `libgeos` CRAN package bundles a copy of GEOS and exposes the C API for C++ packages. It explicitly supports cpp11:

```cpp
#include <cpp11.hpp>
#include "libgeos.h"
#include "libgeos.c"  // once per package — function pointer init
```

With `LinkingTo: libgeos` in DESCRIPTION — no system GEOS install needed, no configure script.

The `geos` R package stores geometries as external pointers to `GEOSGeometry*`. Our package would:
- `Imports: geos` (for R-level geometry creation/conversion)
- `LinkingTo: libgeos, cpp11` (for C++ GEOS API access)
- Accept `geos_geometry` vectors at R level
- Extract the `GEOSGeometry*` pointers in C++ to pass to the exactextract core

This means users can do `geos::as_geos_geometry(sf_obj)` or `geos::geos_read_wkt(...)` to get input, and our package never touches sf directly.

### 4. Lines: deferred

Lines are a separate problem with a different algorithm (no ring traversal / flood fill). Park it for now. The package API should not preclude adding line support later, but it's not in scope for the initial build.

### 5. Output format

**Approximate mode** (controlledburn equivalent): the established RLE format.

| Column | Type | Description |
|--------|------|-------------|
| `row` | integer | 1-based row index |
| `col_start` | integer | 1-based start column (inclusive) |
| `col_end` | integer | 1-based end column (inclusive) |
| `id` | integer | 1-based geometry index |

This is the format from controlledburn and the earlier RLE work. All cells implicitly have weight = 1.

**Exact mode** (coverage fractions): two record types.

The interior cells (fully covered, weight = 1) still compress beautifully as RLE runs.
The edge cells (partially covered, 0 < weight < 1) are individual cells with a weight.

**Two-table approach** (cleanest):

*Runs table* — interior cells, all weight = 1:

| `row` | `col_start` | `col_end` | `id` |

*Edges table* — boundary cells with fractional coverage:

| `row` | `col` | `weight` | `id` |

Return as a named list: `list(runs = ..., edges = ...)`. The runs table IS the controlledburn output. The edges table adds the exactextract information.

**Single-table shortcut** (if preferred): use a sentinel in `col_end` to distinguish:

| `row` | `col_start` | `col_end` | `id` |

- Run record: `col_end >= col_start` → cells col_start:col_end, weight = 1
- Edge record: `col_end = NA` → single cell at `col_start`, weight stored... somewhere

The "somewhere" is the problem. Options: (a) pack weight into a 5th column that's NA for runs, (b) reinterpret col_end as weight when it's < col_start (hacky but compact), (c) just use two tables.

**Decision: two tables.** It's clean, the runs table works identically to existing controlledburn consumers, and the edges table is a straightforward addition. A utility function can merge them to a fully-expanded per-cell table when needed for downstream work (e.g. the chunk materialisation from the earlier RLE conversation).

## Exploration Tasks

### Phase 1: Extract the exactextract core and get coverage fractions out

**Task 1.1**: Clone exactextract, identify the minimal C++ headers needed for the coverage fraction computation. Key files to isolate: `grid.h`, `raster.h`, `matrix.h`, `raster_cell_intersection.h` (or equivalent entry point), `traversal_areas.h`, `floodfill.h`, `box.h`, `geos_utils.h`. Build standalone with `BUILD_CLI=OFF` (no GDAL).

**Task 1.2**: Write a minimal C++ test that takes a GEOS geometry pointer + Grid spec → vector of `(row, col, coverage_fraction)` tuples. Proves the core can be called without the framework.

**Task 1.3**: Wrap in cpp11 R function. Input: `geos_geometry` vector + extent + dimension. Output: the two-table list (`runs` + `edges`). This requires:
- Accepting geos geometry handles from R (via geos package external pointer)
- Calling the exactextract coverage computation
- Splitting output into runs (weight == 1) and edges (weight < 1)
- Returning as named list of data.frames

This is the MVP.

### Phase 2: Approximate mode (controlledburn equivalent)

**Task 2.1**: Add the "cheap" path. Two options:
- (a) Run exact algorithm, discard fractional weights, RLE-compress → runs table only
- (b) Port the fasterize/controlledburn scanline algorithm into the same codebase for maximum speed

Option (b) is faster but more code. Option (a) gets us there immediately.

**Task 2.2**: Benchmark approximate mode vs controlledburn.

### Phase 3: Polish and hypertidy integration

**Task 3.1**: Ensure output works with vaster for cell ↔ coordinate conversion.

**Task 3.2**: Test exact mode against exactextractr `coverage_fraction()` for consistency.

**Task 3.3**: Utility functions:
- `merge_runs_edges()` → fully-expanded per-cell table with weight column
- `materialise_chunk(runs, edges, row_range, col_range)` → dense matrix block (the chunk-coordinate transform from the earlier RLE conversation)

**Task 3.4**: Package structure, vignette, CRAN checks.

## Risk Assessment

| Risk | Impact | Likelihood | Mitigation |
|------|--------|------------|------------|
| exactextract C++ core too entangled | High | Medium | Phase 1 tests this early. Fallback: reimplement ring-traversal from the published description. |
| cpp11 interop with vendored C++ | Medium | Low | exactextract core is standard C++17, cpp11 is a thin header. |
| geos package C API access from cpp11 | Medium | Low | geos package exposes external pointers; accessing the underlying `GEOSGeometry*` should be straightforward. |
| Competing with exactextractr on CRAN | Low | Low | Different purpose — sparse intermediate representation, not zonal stats. |

## Dependencies (Target)

**Hard requirements:**
- cpp11 (C++ interface)
- geos (geometry input + GEOS C API linkage)

**Suggested:**
- vaster (grid index ↔ coordinate utilities)

**Vendored C++:**
- exactextract core (`Grid`, `RasterCellIntersection`, `traversal_areas`, `floodfill` — Apache-2.0)

**NOT depending on:**
- sf, terra, raster, stars (no spatial framework objects)
- GDAL (no I/O)
- Rcpp
- wk (geos handles geometry; wk conversion is optional convenience)

## Naming Ideas

- **gridburn** — echoes controlledburn, emphasises the grid
- **cellweight** — describes the output directly
- **sparseburn** — the sparse output emphasis
- **burnweight** — controlledburn + weight
- **gridtrace** — traces geometry through a grid

## Next Steps

**Immediate: R prototype using exactextractr as backend**

The `GEOSGridIntersectionFractions_r` signature tells us the GEOS function returns a dense `float[nx*ny]` buffer. This means our package's core value is:
1. Per-geometry bounding box clipping (only compute the subgrid)
2. Dense → sparse conversion (the two-table format)
3. Multi-geometry accumulation

All of this can be prototyped in pure R using exactextractr + terra before writing any C++. See `gridburn-prototype.R` for the implementation.

**Phase 1 revised**:
1. Run the R prototype on real data (NC counties, Antarctic polygons)
2. Validate two-table output format against controlledburn (runs) and exactextractr (coverage fractions)
3. Test the `materialise_chunk()` utility with the chunk-coordinate transform from the earlier RLE work
4. Benchmark: is the R prototype fast enough for your use cases, or is C++ needed?

**Phase 2: C++ package** (when needed):
- If libgeos has updated to GEOS 3.14 by then: call `GEOSGridIntersectionFractions_r` directly via `LinkingTo: libgeos`
- If not: vendor the exactextract C++ core (same algorithm, same output)
- Either way: the `dense_to_sparse` conversion moves to C++ for speed

**Also worth tracking**: Dan Baston's GEOS PR #1232 ("Add grid subdivision / rasterization algorithms") is still open. This may add additional capabilities beyond just coverage fractions.
