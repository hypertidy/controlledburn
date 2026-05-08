# controlledburn (development version)

Iteration extending controlledburn from polygon-only to a unified
geometry rasterizer over polygon, line, and point input. Output schema
gains type-pure tables; the unification follows the geometric structure
of measure-on-a-grid as one operation parameterised by Hausdorff
dimension. Design philosophy and decisions are recorded in
`inst/docs-design/unified-geometry-rasterization.md`.

## Breaking changes

* Polygon boundary cells: `$edges$weight` renamed to `$edges$fraction`
  to reflect what the column carries (a dimensionless fraction in
  [0, 1]). [`materialise_chunk()`] reads either name for forward
  compatibility with cached results, but new code should use
  `$edges$fraction`.

* `burn_sparse()` now errors on line and point input. Previously it
  fell through to a path that produced silent-empty (points) or
  silently-mistyped (lines) output. The error message redirects to
  `burn_scanline()`.

## What's new

* `burn_scanline()` accepts LINESTRING, MULTILINESTRING, POINT, and
  MULTIPOINT alongside the existing POLYGON / MULTIPOLYGON support.
  Output schema gains `$lines` (`row, col, length, id`; absolute
  length in CRS units) and `$points` (`row, col, id`; no measure
  column — implicit weight = 1).

* `materialise_chunk()` handles all four tables. Line cells contribute
  length, point cells contribute count, polygon cells contribute
  fraction. Mixed-kind input errors rather than silently combining
  units; filter to one kind via `id =` first.

* `print.controlledburn` reports `$lines` and `$points` counts when
  present.

* GeometryCollection input is rejected with a warning. Mixed-dimension
  input would produce a sparse table where rows from different kinds
  carry different units — caller must split into homogeneous groups
  and run separate burns. Curved types (CircularString, CompoundCurve,
  etc.) must be linearised by the caller before passing in.

## Internals

* `walk_polyline()`: cell-stepping walker extracted from `walk_ring()`.
  Geometry-agnostic; takes a `closed` parameter that selects whether
  to push back coords for the start-cell-reentry case (true for rings,
  false for open lines). The polygon path runs the existing
  winding-sweep post-processing on the walker's output; the line path
  sums per-cell segment lengths from the same `CellMap`.

* `process_line()`: walks lines on the full grid (no `shrink_to_fit`),
  because horizontal/vertical lines have a degenerate bbox that
  `Box::empty()` reports as true. Walker cost is O(cells touched), so
  the bbox-shrinking optimisation is dead weight for lines anyway.

* `process_point()`: trivial. Computes cell index via
  `Grid::get_row` / `get_column`, drops out-of-extent silently.

* Edge zoo (`tests/testthat/test-edge-zoo.R`) pins canonical
  rasterizer edge cases by category (Horizontal, Sub-pixel,
  Alignment, Precision, Topology, Collinear, CRS-boundary). Each
  test pins a *convention*, not just a numerical expectation, so
  the test surface doubles as API contract documentation.

# controlledburn 0.1.0

Complete rewrite of controlledburn using the exactextract algorithm (Daniel
Baston, vendored from exactextractr) for exact polygon-grid coverage fractions.

## What's new

* `burn_scanline()`: O(perimeter) scanline sweep with winding-number interior
  classification and exact boundary coverage fractions. No dense matrix
  allocation — output is sparse runs + edges tables.

* `burn_sparse()`: Reference implementation using the exactextract dense
  algorithm, compressed to the same sparse output format.

* `materialise_chunk()`: Opt-in expansion to dense matrix or vector, with
  per-polygon-id filtering.

* Default grid parameters: extent derived from geometry bbox via `wk::wk_bbox()`,
  dimension auto-fitted to 256 cells on the long axis preserving aspect ratio,
  or specified as `resolution`.

* Geometry input via `wk::wkb()`, `geos_geometry`, `sfc`, `blob`, or raw WKB
  list (compatible with vapour/gdalraster output).

* Moved from Rcpp to cpp11, using libgeos for GEOS access.

## Internals

* Scanline algorithm: lightweight walk using `Box::crossing()` directly (no
 `Cell` class allocation), winding-count interior classification, analytical
  single-traversal coverage via `perimeter_distance()`.

* Validated against burn_sparse across 52 test cases: simple shapes, NC
  counties, shared boundaries, edge cases (grid-aligned edges, slivers,
  extent clipping, degenerate shapes).

* O(perimeter) scaling confirmed by benchmark: 17× faster than dense at
  3200×3200 resolution for complex shapes. Memory: sparse output ~50 MB vs
  ~2 GB dense for real-world 32K×16K grids.

## History

* Derived from fasterize (Noam Ross, EcoHealth Alliance) — original scanline
  edge logic for polygon rasterization.

* Previous controlledburn (0.0.2) returned binary in/out run-length indexes
  without coverage fractions. That code is archived at
  `archive0-2026-02-13` branch.

# controlledburn 0.0.2

* Basic function working, return a list of triplets of zero-index
  start,end,line index.

* Raster package objects and use removed.

* Converted from fasterize edge scanline logic, to not materialize any
  raster data. Armadillo usage removed.
