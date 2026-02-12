# controlledburn: scanline polygon rasterization with exact coverage fractions

## The goal

Compute the exact same sparse two-table output as gridburn (interior runs + weighted boundary edges) without ever allocating a dense matrix. Memory usage is O(perimeter) not O(bounding-box area). The output is a **polygon-grid intersection database** — a pluripotent intermediate from which dense rasters, zonal statistics, and polygon overlays can all be derived on demand.

## Key insight

The exactextract algorithm computes exact coverage fractions per boundary cell by walking polygon edges through the grid, cell by cell. This walk is already O(perimeter). The dense matrix only exists to support the flood fill that classifies interior cells. But interior classification is exactly what a scanline rasterizer does — by counting edge crossings per row, we know which spans are inside without ever touching the interior cells.

So: keep the per-cell coverage fraction computation, replace the flood fill with a scanline winding-count sweep.

## Data structures

```cpp
// A polygon edge: one segment between consecutive ring vertices.
// Preprocessed with grid-row span for the sweep.
struct ScanEdge {
    double x0, y0, x1, y1;   // segment endpoints (y0 <= y1, i.e. y0 is bottom)
    int row_enter;            // first grid row the edge touches (topmost, smallest index)
    int row_exit;             // last grid row the edge touches
    double x_at_row_top;     // x-intersection at top of current row band (maintained during sweep)
    double dx_per_row;        // x increment per row (for non-horizontal edges)
    int ring_id;              // which ring (for winding direction)
    int poly_id;              // which polygon
    bool is_exterior;         // exterior ring or hole
    bool upward;              // edge goes from bottom to top (y increasing)
};

// Active Edge Table entry: an edge currently crossing the active row band.
// Sorted by x_at_row_top for left-to-right processing.
struct ActiveEdge {
    ScanEdge* edge;
    // Crossing information for winding rule:
    // +1 if edge crosses row center going right-to-left (upward edge)
    // -1 if edge crosses row center going left-to-right (downward edge)
    // 0 if edge doesn't cross row center (enters and exits same side)
    int winding_contribution;
    int col_entry;            // leftmost cell this edge touches on this row
    int col_exit;             // rightmost cell this edge touches on this row
};

// Output: same as gridburn
struct GridRun {
    int row, col_start, col_end, id;
};

struct GridEdge {
    int row, col;
    float weight;
    int id;
};
```

## Algorithm overview

```
1. PREPROCESS
   For each polygon ring, decompose into ScanEdge segments.
   Sort all edges into a Global Edge Table (GET) by row_enter.

2. SWEEP rows top to bottom (row 1 to nrow)
   a. Activate: move edges from GET whose row_enter == current_row into AET.
   b. Process: for each active edge, walk it through cells on this row.
      Compute exact coverage fractions for boundary cells.
      Record crossing events for the winding rule.
   c. Classify: sort crossings by column, apply winding rule to determine
      interior spans between boundary cells.
   d. Emit: write runs for interior spans, edges for boundary cells.
   e. Retire: remove edges whose row_exit == current_row from AET.
   f. Advance: update x_at_row_top for remaining edges.
```

## Detailed algorithm

### Step 1: Preprocessing

For each polygon, for each ring, for each consecutive vertex pair (v_i, v_{i+1}):

```cpp
// Skip horizontal edges — they don't contribute crossings and their
// cell coverage is handled by the adjacent non-horizontal edges.
if (y0 == y1) continue;

// Normalise so y0 < y1 (bottom to top in geographic coords).
// In grid coords (row 0 = ymax), this means row_enter < row_exit.
// Track original direction for winding.
bool upward = (original_y1 > original_y0);  // going up in geographic space

ScanEdge e;
e.y0 = min(vy_i, vy_{i+1});  // bottom
e.y1 = max(vy_i, vy_{i+1});  // top
// ... set x0, x1 to match y ordering

// Which grid rows does this edge span?
// Row numbering: row 0 is at ymax, row (nrow-1) is at ymin.
e.row_enter = floor((ymax - e.y1) / dy);  // topmost row
e.row_exit  = floor((ymax - e.y0) / dy);  // bottommost row
// clamp to [0, nrow-1]

e.dx_per_row = (e.x1 - e.x0) / (e.y1 - e.y0) * dy;
e.x_at_row_top = x_at(ymax - e.row_enter * dy);  // interpolate
```

Sort all ScanEdges by `row_enter` into the Global Edge Table.

### Step 2: The row sweep

```cpp
GlobalEdgeTable GET;  // sorted by row_enter
ActiveEdgeTable AET;  // sorted by x_at_row_top

std::vector<GridRun> all_runs;
std::vector<GridEdge> all_edges;

for (int row = 0; row < nrow; row++) {

    double row_ymax = grid_ymax - row * dy;
    double row_ymin = row_ymax - dy;
    double row_ymid = row_ymax - dy / 2;

    // --- 2a. Activate new edges ---
    while (GET has edges with row_enter == row) {
        AET.insert(GET.pop());
    }

    // --- 2b. Process active edges through cells on this row ---
    //
    // For each active edge, we need:
    //   (i)  which cells it passes through on this row
    //   (ii) exact coverage fraction for each such cell
    //   (iii) whether it crosses the row's centre line (for winding)
    //
    // This is the per-row equivalent of exactextract's Cell traversal,
    // but we only instantiate Cell objects for the boundary cells.

    struct CellCrossing {
        int col;
        float coverage_fraction;
        int poly_id;
        int winding_delta;  // +1, -1, or 0
    };

    std::vector<CellCrossing> crossings;

    for (ActiveEdge& ae : AET) {
        ScanEdge* e = ae.edge;

        // Clip the edge segment to this row band [row_ymin, row_ymax]
        Segment seg = clip_to_band(*e, row_ymin, row_ymax);

        // Walk the clipped segment through grid cells left to right.
        // For each cell the segment passes through:
        std::vector<CellHit> hits = walk_segment_through_cells(
            seg, grid_xmin, dx, row_ymax, dy);

        for (auto& hit : hits) {
            // Exact coverage fraction for this cell.
            // This uses exactextract's Cell class:
            //   - create Cell for grid_cell(row, hit.col)
            //   - feed it the entry and exit coordinates
            //   - call covered_fraction()
            //
            // For the controlledburn context, we only need the
            // traversal of THIS edge through THIS cell. Multiple
            // edges through the same cell get combined later.
            float frac = compute_cell_coverage(seg, hit, row, grid);

            // Does this edge cross the row centre line within this cell?
            // An edge that enters from top and exits bottom (or vice versa)
            // crosses the centre. An edge that enters and exits the same
            // horizontal half doesn't.
            int wd = 0;
            if (seg_crosses_y(seg, row_ymid, hit.col_xmin, hit.col_xmax)) {
                wd = e->upward ? +1 : -1;
                // CCW exterior ring: upward edge on right side → entering
                // CW hole ring: opposite
                if (!e->is_exterior) wd = -wd;
            }

            crossings.push_back({hit.col, frac, e->poly_id, wd});
        }
    }

    // --- 2c. Classify interior spans ---
    //
    // Sort crossings by (poly_id, col).
    // Process each polygon independently on this row.

    std::sort(crossings.begin(), crossings.end(),
        [](const CellCrossing& a, const CellCrossing& b) {
            return std::tie(a.poly_id, a.col) < std::tie(b.poly_id, b.col);
        });

    // Group by poly_id
    for (auto [poly_id, poly_crossings] : group_by_poly(crossings)) {

        // Merge crossings in the same cell (multiple edges through one cell).
        // Sum coverage fractions (capped at 1.0).
        // Sum winding deltas.
        auto merged = merge_same_cell(poly_crossings);

        // Sweep left to right with winding count.
        int winding = 0;
        int prev_col = -1;

        for (auto& mc : merged) {

            // Emit interior run between previous boundary cell and this one
            if (winding != 0 && prev_col >= 0 && mc.col > prev_col + 1) {
                all_runs.push_back({
                    row + 1,           // 1-based
                    prev_col + 1 + 1,  // column after previous boundary, 1-based
                    mc.col - 1 + 1,    // column before this boundary, 1-based
                    poly_id
                });
            }

            // Emit boundary cell
            if (mc.coverage_fraction > 0 && mc.coverage_fraction < 1.0f) {
                all_edges.push_back({row + 1, mc.col + 1,
                                     mc.coverage_fraction, poly_id});
            } else if (mc.coverage_fraction >= 1.0f) {
                // Degenerate: boundary cell fully covered.
                // Could be a run of length 1, or include in adjacent run.
                // For simplicity, emit as a single-cell run.
                all_runs.push_back({row + 1, mc.col + 1, mc.col + 1, poly_id});
            }

            // Update winding count AFTER processing this cell
            winding += mc.winding_delta;
            prev_col = mc.col;
        }

        // Note: winding should return to 0 after processing all crossings
        // on a row for a valid polygon. If it doesn't, we have a problem.
    }

    // --- 2e. Retire spent edges ---
    AET.remove_if([row](const ActiveEdge& ae) {
        return ae.edge->row_exit == row;
    });

    // --- 2f. Advance x intercepts ---
    for (ActiveEdge& ae : AET) {
        ae.edge->x_at_row_top += ae.edge->dx_per_row;
    }
}
```

### The boundary cell coverage fraction

This is where exactextract's geometric precision matters. For each cell that a polygon edge passes through, we need the exact area of the polygon within that cell.

The key insight: we don't need the full Cell class machinery for most cases. A polygon edge crossing a grid cell creates a trapezoid (or triangle) whose area can be computed analytically from the entry and exit coordinates.

```
     ┌─────────────────┐
     │           ╱     │
     │         ╱       │
     │       ╱  inside │
     │     ╱           │
     │   ╱             │
     └─────────────────┘

Entry: (x_entry, y_entry) on left edge of cell
Exit:  (x_exit, y_exit) on bottom edge of cell
```

For a single edge crossing a cell, the covered area is:

```
If edge enters from side A and exits side B, the polygon interior
(assuming CCW winding) is to the left of the edge direction.
The covered fraction is the area of the polygon between:
  - the edge segment within the cell
  - the cell boundary segments connecting entry to exit going
    clockwise around the polygon interior
divided by the cell area (dx * dy).
```

exactextract's `Cell` class handles the general case including multiple edges through a single cell. For the scanline algorithm, we could:

**Option A: Reuse Cell class directly.** For each boundary cell, instantiate a Cell, feed it all edge segments that pass through it, call `covered_fraction()`. This is exact and handles all edge cases. Cost: one heap allocation per boundary cell.

**Option B: Analytical computation for single-edge cells.** Most boundary cells are crossed by exactly one edge. The covered fraction is a simple geometric calculation — the area of the polygon slice within the cell. For the rare multi-edge cells (vertices, self-intersections), fall back to the Cell class.

Option B is the performance path for production. Option A is the safe starting point.

### Multi-edge cells

When two edges from the same polygon meet at a vertex that falls inside a grid cell, that cell sees two edge segments. The coverage fraction comes from both:

```
     ┌──────────────────┐
     │    ╲      ╱      │
     │     ╲    ╱       │
     │      ╲  ╱        │
     │       ╲╱ vertex  │
     │      inside      │
     └──────────────────┘
```

The scanline sweep naturally handles this: both edges are active on the same row, both generate a `CellCrossing` for the same cell. The `merge_same_cell` step sums the coverage fractions and winding deltas.

For exact coverage: the two edge segments together with the cell boundary define the covered region. The Cell class computes this correctly by accumulating all traversals. The analytical approach would compute each edge's contribution as a signed area and sum them.

### Holes

A polygon hole (interior ring) contributes edges with opposite winding. In the sweep:

```
Row through a donut:

     exterior edge    interior (run)    hole edge    hole interior    hole edge    interior (run)    exterior edge
col: ──────|───────────────────────────────|──────────────────────────────|──────────────────────────────|──────
winding:   0 → 1                          1 → 0                         0 → 1                          1 → 0

Emits: edge, run, edge, (nothing — winding=0 inside hole), edge, run, edge
```

The winding count naturally drops to 0 inside the hole, so no interior run is emitted there. The hole's boundary cells get their own coverage fractions, and because the hole's winding is opposite to the exterior ring's, the coverage fraction in a cell that's partly inside the hole is correctly reduced.

For a cell that's ON the hole boundary but also fully inside the exterior ring: the exterior ring contributes weight=1.0 for this cell (it's interior to the exterior ring), but the hole edge subtracts the hole's coverage. The net weight is `1.0 - hole_edge_fraction`. This matches exactextract's `add_ring_results` with `factor = -1`.

### Winding vs even-odd

The algorithm uses the winding number rule (nonzero = inside) rather than even-odd. This correctly handles:

- Simple polygons (equivalent to even-odd)
- Polygons with holes (winding drops to 0 inside holes)
- Overlapping multipolygon components (winding accumulates)

For OGC-compliant polygons (exterior CCW, holes CW), the winding count is always 0 or 1. Self-intersecting polygons could produce higher counts — we'd treat any nonzero as "inside."

## Complexity

Let P = total perimeter in grid cells (sum of edge-cell intersections across all polygons).

| Operation | Complexity |
|---|---|
| Preprocessing (sort edges) | O(E log E) where E = number of edge segments |
| Sweep (all rows) | O(nrow × avg_active_edges + P) |
| Per boundary cell | O(1) for single-edge, O(k) for k-edge cells |
| Memory | O(P) for crossings, O(output) for runs+edges |

No dense matrix. No flood fill. No GEOS PreparedGeometry calls.

For Lake Superior at 256k × 128k: the coastline might produce ~500k boundary cells. Memory ≈ a few hundred MB for the crossing structures during sweep, and the output sparse tables. Compare to exactextract's 122 GB dense matrix.

## The traversal strategy question

The algorithm as described is **polygon-sequential within a row-band sweep**. All polygons active on a given row are processed together. This is the natural decomposition.

But the sweep can be chunked:

### Row-band parallelism

Divide the grid into horizontal bands of H rows each. Each band is processed independently — its own GET (edges entering within this band), its own AET (edges that were already active when the band started). Bands are embarrassingly parallel if you precompute the initial AET state for each band.

The initial AET for band B is the set of edges that entered in bands 0..B-1 and haven't yet exited. Computing this requires a single top-to-bottom pass over the edge list, recording which edges are active at each band boundary. O(E) preprocessing, then each band runs independently.

### Tile parallelism

For 2D locality (important if polygons are spatially clustered), divide into tiles. Each tile is a row-band restricted to a column range. Edges are clipped to the tile's column range. Runs that span the tile boundary are emitted as-is (they'll be adjacent in the output, just not merged across tile boundaries).

### Adaptive banding

Most rows in a typical dataset have few active edges (empty space). Rows through a complex coastline have many. Adaptive band heights — wide bands through empty space, narrow through dense edges — could balance work.

## Integration with exactextract's Cell class

The cleanest starting point reuses the Cell class for coverage fractions:

```cpp
float compute_cell_coverage_via_cell(
    const std::vector<Segment>& edge_segments,  // all edges through this cell
    const Box& cell_box,
    bool is_areal)
{
    Cell cell(cell_box);

    for (auto& seg : edge_segments) {
        cell.take(seg.start, &seg.prev);
        if (cell.last_traversal().exited()) {
            cell.force_exit();
        }
    }
    cell.force_exit();

    return static_cast<float>(cell.covered_fraction());
}
```

This produces identical results to exactextract for every boundary cell. The difference is we only instantiate Cells for boundary cells (O(perimeter)) instead of the full bbox.

## Output format

Identical to gridburn:

```
runs:  data.frame(row, col_start, col_end, id)  — interior cells, weight = 1.0
edges: data.frame(row, col, weight, id)          — boundary cells, 0 < weight < 1
```

This is the canonical polygon-grid intersection database. From it:

- **Dense raster**: expand runs, place edges. O(output pixels).
- **Zonal stats**: aggregate over runs (count, area) and edges (weighted). O(sparse records). No pixel data needed.
- **Polygon overlay**: merge-join two sparse databases on (row, col). O(sum of sparse records).
- **Streaming write**: emit rows sequentially to GDAL. O(1) memory per row.
- **Subwindow query**: filter runs/edges by row/col range. O(log n) with sorted index.

## What remains

1. **Prototype the sweep** in R + cpp11, reusing exactextract's Cell class for coverage fractions. Validate against gridburn output.

2. **Analytical single-edge coverage** to avoid Cell allocation for the common case. The geometry is a trapezoid/triangle — closed form.

3. **Benchmark** perimeter-proportional scaling vs gridburn's tiled approach.

4. **Multi-polygon merge**: when two polygons share a boundary, their boundary cells overlap. The output should handle this cleanly (separate edge records with different `id` values, or a merged record with combined weight).

5. **Edge cases**: polygon vertex exactly on a cell boundary, edge exactly along a cell boundary (horizontal or vertical), tiny sliver polygons where the edge enters and exits the same cell side.
