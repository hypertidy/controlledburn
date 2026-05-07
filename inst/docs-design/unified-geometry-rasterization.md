# Unified geometry rasterization in controlledburn

Design doc for the next big iteration. Consolidates the framing, decisions, and implementation plan for extending controlledburn from polygon-only rasterization to a unified polygon / line / point rasterizer over a single sparse output schema.

No backward-compatibility constraints — controlledburn has no production users, only sympathetic experimenters. Hard breaking changes are acceptable wherever they improve the design.

---

## Philosophy

Six load-bearing claims. Each justifies the next; pulling any one out weakens the others.

### 1. Measure on a grid is one operation parameterised by dimension

Polygon-area-weighted, line-length-weighted, and point-presence rasterization are not three algorithms — they are one algorithm (`geometry ∩ cell` measure) specialised by the dimension of the input geometry. The 2×2 makes this visible:

|         | binary (touch test) | measure (geometry ∩ cell) |
|---------|---------------------|----------------------------|
| polygon | cell-centre test (legacy fasterize path) | exact area fraction (`walk_ring` + winding sweep) |
| line    | "line intersects cell" | length-in-cell (`walk_polyline` + length sum) |
| point   | cell membership | — (no measure for 0-D) |

The dimension-completeness argument closes the family rigorously: 0, 1, and 2 are the only Hausdorff dimensions a 2-D embedding admits. There is no 3-D primitive in plane geometry, so the family is genuinely complete at three members rather than arbitrarily truncated. MapInfo, Manifold GIS, and ArcInfo coverages historically shipped this insight as a single feature class because, geometrically, it *is* a single feature class.

### 2. The walker is the primitive; post-processing is the variation

`walk_polyline` traversing cells via `Box::crossing()` and recording per-cell coordinate sequences is the shared core. Polygon mode runs winding + coverage on that output. Line mode sums segment lengths on that output. Point mode skips the walker entirely (degenerate case).

The polygon-specific machinery — CCW normalisation, winding sweep, signed area — was never load-bearing for the walk itself. It was post-processing that happened to live inside `walk_ring` for historical reasons. Factoring it out reveals the structure that was already there.

### 3. Emit absolute measures; let downstream compose

The `measure` column carries the primitive quantity in the geometry's coordinate space: square CRS units for polygons, CRS units for lines. Not fractions.

Polygon area fractions appear natural because cell area is a single scalar — but length has no canonical scalar to divide by once cells are non-square (`xres ≠ yres`, common in practice and universal in geographic CRSes). Any choice of denominator (cell width, cell height, diagonal) is arbitrary and silently axis-biased. Absolute measures don't have that problem; downstream normalisation is always one division away.

This matches the silicate-style "primitives in, downstream composes" stance already operating in `vaster` (cell ↔ xy as primitives, no fused transformation), `silicate` (segments and vertices as primitives, no fused topology), and the broader hypertidy ecosystem.

**Coordinate-space, not geodesic.** The measure is `area_of(geometry ∩ cell)` computed in whatever CRS the geometry and grid live in. Geometry in degrees → area in square degrees, length in degrees. For geodesic measures, project first. Documented explicitly; not silently corrected.

### 4. Output schema as the contract

`(row, col, [weight,] id)` is the format. Small enough to be obviously correct, large enough to carry all three primitive types without forking, and aligned with what `vaster` already consumes for cell ↔ xy round-trips and what `polymer2` consumes for sparse overlay.

Points carry no `weight` column. Materialise treats absent-weight as "weight = 1 for every row" — NULL-defaults rather than schema branching.

### 5. API names what's in the column, not how it was computed

User-facing argument: `weight = c("none", "measure")`. Slot-based naming describes the data, not the procedure. Extensible — future modes (`"density"`, `"attribute"`, etc.) slot in without renaming. Matches the reality that downstream callers routinely apply attribute values or other weightings on top of the geometric measure.

`exact = TRUE/FALSE` is retained internally as dispatch / cross-package vocabulary, since "exact coverage" is the established phrase in exactextract / exactextractr discussions. It would not survive the addition of a third weight mode, which is why it's not the public surface.

### 6. Closure is not catholicism

The format closes under the three primitive types and their multipart variants. It explicitly does not close under:

- **GeometryCollection** — mixing dimensions in one input would produce rows with different weight semantics (50 m² of polygon vs 50 m of line, indistinguishable in the output table). Caller responsibility: split into homogeneous groups, run separate burns, label with distinct `id`s.
- **Curved types** (CircularString, CompoundCurve, etc.) — caller responsibility: linearise via `wk::wk_linestring()` or equivalent before passing in.

Knowing exactly what's in and what's out is part of the design, not a limitation of it.

---

## Architecture

### Walker core

`walk_polyline(coords, grid, ...) → vector<CellTraversal>`

Steps a polyline through grid cells using `Box::crossing()` for cell-edge intersections. Records per-cell-visit:

```cpp
struct LightTraversal {
    std::vector<Coordinate> coords;   // entry → intermediates → exit
    Side entry_side = Side::NONE;
    Side exit_side = Side::NONE;
};
```

This is the existing structure inside `walk_ring`, with polygon-specific concerns (CCW normalisation, `is_exterior` flag, coupling to winding sweep) extracted out.

### Geometry-type post-processors

| Type | Walker | Post-processor | Output |
|------|--------|----------------|--------|
| Polygon | `walk_polyline` per ring + CCW normalisation | Winding sweep + coverage fraction × cell area | Runs (interior) + edges (boundary, weighted) |
| LineString | `walk_polyline` | Length sum per cell: `Σ |coords[i+1] − coords[i]|` | One row per non-empty cell, `(row, col, length, id)` |
| Point | None (direct cell-index computation) | None | One row per point, `(row, col, id)` |

The polygon path retains its current two-table output (runs + edges) because interior cells genuinely have constant weight = `cell_area` and RLE compression is a real win. Line and point paths produce a single edges-style table — every touched cell is "boundary" in the sense that there's no interior.

### Geometry-type dispatch

`process_geometry(geom, ...)` switches on `GEOSGeomTypeId_r`:

```
GEOS_POLYGON         → process_polygon
GEOS_MULTIPOLYGON    → recurse on parts
GEOS_LINESTRING      → process_line
GEOS_MULTILINESTRING → recurse on parts
GEOS_POINT           → process_point
GEOS_MULTIPOINT      → recurse on parts
GEOS_GEOMETRYCOLLECTION → reject with clear error
other                → reject with clear error
```

Currently `process_geometry` recurses on `GEOS_GEOMETRYCOLLECTION` (line 406 of `scanline_burn.cpp`). This becomes a rejection.

### R-side API

```r
burn_sparse(x,
            extent = NULL, dimension = NULL, resolution = NULL,
            weight = c("none", "measure"))
```

One function, three geometry types, one sparse output. The `weight` argument is the primary user-facing control. Geometry type is autodetected from WKB.

`burn_scanline()` either becomes the preferred implementation throughout (replacing the `cpp_burn_sparse` exactextract path) or is folded back into `burn_sparse` as a single export. To be settled during implementation — the scanline path is faster and produces the same numbers; the only reason to retain both is if the exactextract path has features the scanline path doesn't.

### Materialise

`materialise_chunk()` accepts any controlledburn output regardless of geometry type. If the input has no `weight` column, treat as `weight = 1` per row (no schema branching, just a NULL default in the chunk-fill loop).

---

## Implementation plan

### Phase 1 — Walker extraction

1. Refactor `walk_ring` → `walk_polyline` core + polygon wrapper. The wrapper does CCW normalisation, calls `walk_polyline`, then runs the existing winding sweep.
2. Confirm polygon output is bit-identical before and after the refactor. This is the load-bearing test for the extraction — if numbers shift, the extraction is wrong.

### Phase 2 — Line case

1. Add `process_line` that calls `walk_polyline` and accumulates segment lengths per cell.
2. Add LINESTRING / MULTILINESTRING dispatch in `process_geometry`.
3. Cross-check against an independent length computation (e.g. clip each segment to each cell using GEOS `intersection`, sum lengths) on a few test geometries.
4. Document that the `weight` column for lines is in CRS units.

### Phase 3 — Point case

1. Add `process_point` (compute cell index, emit row, no walker).
2. Add POINT / MULTIPOINT dispatch.
3. Confirm `materialise_chunk` handles weight-less input.

### Phase 4 — API consolidation

1. Replace `burn_sparse` / `burn_scanline` two-function surface with the unified entry point.
2. Reject GeometryCollection and curved types with clear error messages.
3. Update vignette to walk through all three geometry types with the same data flow.
4. Update README.

### Phase 5 — Cross-package validation

1. `vaster`: confirm `(row, col, weight, id)` and `(row, col, id)` round-trip through `cell_from_row_col` / `xy_from_cell` unchanged.
2. `polymer2`: confirm line-measure output composes through the merge-join overlay path the polygon path uses.
3. `materialise_chunk`: end-to-end test of mixed-id input across all three geometry types in a single sparse table.

### Out of scope for this iteration

- Performance benchmarking (after correctness).
- Migration of existing controlledburn line callers — none exist; the line case is greenfield.
- `density` / `attribute` weight modes — slot is reserved in the API, not implemented now.

---

## Open decisions

These are not blocking; reasonable defaults exist. Worth a quick check during implementation, not a design review.

- [ ] **Whether to retain `cpp_burn_sparse` (exactextract path) alongside `cpp_scanline_burn`.** Default: drop, the scanline path is faster and produces the same numbers. Re-evaluate only if a feature gap appears.
- [ ] **Whether multi-vertex traversals (line vertex falls inside a cell) need any special handling beyond the natural sum.** Default: no. The traversal `[entry, vertex, exit]` produces the correct length under a straight sum.
- [ ] **Self-intersecting linestrings.** Default: caller responsibility — if a line crosses itself within a cell, the recorded length is the geometric truth (sum of segment lengths in cell, including the crossing portion counted once per traversal). No deduplication.
- [ ] **Empty / degenerate geometries.** Default: skip with no warning. Already the polygon path's behaviour for empty inputs.

---

## Cross-references

- Existing design docs: `inst/docs-design/exact-extract-design.md`, `inst/docs-design/sparse-rasterization-exploration.md`, `inst/docs-design/denseburn-refactor/`.
- Source: `src/scanline_burn.cpp` (the walker is `walk_ring` at line 137; dispatch is `process_geometry` at line ~395; `Box::crossing` is in `src/exactextract/box.cpp` line 39).
- Sibling packages: `vaster` (cell ↔ xy primitives, output consumer), `polymer2` (sparse overlay, output consumer), `silicate` (philosophical antecedent — primitives over fused abstractions).
