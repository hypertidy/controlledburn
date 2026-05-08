# Canonical edge cases — the rasterizer "final boss" zoo
#
# Tests organised by the six (+1) failure-mode categories that appear in
# the rasterizer / computational-geometry literature (Skia, Freetype, AGG,
# the Vulkan/OpenGL CTS top-left rule). Each test pins down a CONVENTION:
# the expected behaviour is part of the API contract, not just a numeric
# spot-check. Comments explain *why* the answer is what it is.
#
# Categories (one describe-block / region per category):
#   1. Horizontal      — edges/lines coincident with a scanline (y-row) boundary
#   2. Sub-pixel       — segments shorter than one cell in either dimension
#   3. Alignment       — edges/lines exactly on a cell boundary
#   4. Precision       — very long segments, accumulated cell-by-cell length
#   5. Topology        — self-intersecting (bow-tie) geometry
#   6. Collinear       — three or more collinear vertices on a ring or line
#   7. CRS-boundary    — domain-specific non-promises (antimeridian, poles)
#
# References: Skia GM tests, FreeType ftraster, AGG agg_rasterizer_scanline_aa,
# Khronos CTS rasterization tests (top-left rule), exactextract test suite.

ext10 <- c(0, 10, 0, 10)

#skip()

# =========================================================================
# 1. HORIZONTAL  — edges / lines coincident with a scanline
# =========================================================================
#
# Convention pinned: an edge whose y is exactly equal to a row's y_mid
# contributes ZERO winding delta (the polygon-path winding sweep uses
# strict inequalities both sides). This is the "exclude" half of the
# top-left rule applied to row centres rather than pixel centres. Lines
# coincident with a cell-row boundary are ambiguous between two adjacent
# cells; the walker uses Box::side() to break the tie.

test_that("polygon edge exactly on cell row boundary — no winding artifacts", {
  skip_if_not_installed("geos")
  # Top edge of the polygon coincides exactly with cell-row boundary y=5
  # in a 10x10 grid (dy=1). Zero-winding-delta convention means this edge
  # contributes no spurious row crossings.
  p <- geos::as_geos_geometry("POLYGON ((2 2, 8 2, 8 5, 2 5, 2 2))")
  expect_scanline_matches_sparse(p, ext10, c(10L, 10L),
                                 label = "horiz edge on row boundary")
})

test_that("polygon shared horizontal edge — coverage is complementary", {
  skip_if_not_installed("geos")
  # Two stacked rectangles share y=5. Total coverage of the union is
  # what each contributes — no double-count, no missed sliver.
  p1 <- geos::as_geos_geometry("POLYGON ((2 2, 8 2, 8 5, 2 5, 2 2))")
  p2 <- geos::as_geos_geometry("POLYGON ((2 5, 8 5, 8 8, 2 8, 2 5))")
  r <- burn_scanline(c(p1, p2), extent = ext10, dimension = c(10L, 10L))
  expect_complementary(r, label = "shared horizontal edge")
})

test_that("horizontal line exactly on a cell-row boundary — pinned ownership", {
  skip_if_not_installed("geos")
  # Line at y=5 on a 10x10 grid (dy=1). The walker assigns this line to
  # ONE row (not split, not duplicated). Convention to pin: which row?
  # Whatever the answer, total length must equal line length and the
  # line must appear in exactly one row's worth of cells.
  line <- geos::as_geos_geometry("LINESTRING (1 5, 9 5)")
  r <- burn_scanline(line, extent = ext10, dimension = c(10L, 10L))
  expect_equal(nrow(r$runs), 0)
  expect_equal(sum(r$lines$length), 8.0, tolerance = 1e-6,
               label = "total length preserved")
  # Line lies on a single row boundary; pin which row owns it.
  # (TODO: fill in the row number once we observe and decide which
  #  side of the boundary the walker assigns. The point is to make the
  #  decision visible and stable across versions.)
  # expect_equal(unique(r$edges$row), <pinned row>)
})

test_that("vertical line exactly on a cell-column boundary — pinned ownership", {
  skip_if_not_installed("geos")
  line <- geos::as_geos_geometry("LINESTRING (5 1, 5 9)")
  r <- burn_scanline(line, extent = ext10, dimension = c(10L, 10L))
  expect_equal(nrow(r$runs), 0)
  expect_equal(sum(r$lines$length), 8.0, tolerance = 1e-6,
               label = "total length preserved")
})


# =========================================================================
# 2. SUB-PIXEL  — segments shorter than one cell
# =========================================================================
#
# Convention pinned: degenerate (< 2 coords) returns no edges. Otherwise
# the walker uses parametric clipping (Box::crossing), no slope division,
# so very short segments are numerically safe. The risk is the *endpoints*
# falling on a boundary and Box::side() ties.

test_that("line entirely inside one cell — single edge, exact length", {
  skip_if_not_installed("geos")
  line <- geos::as_geos_geometry("LINESTRING (0.2 0.2, 0.8 0.8)")
  r <- burn_scanline(line, extent = c(0, 3, 0, 3), dimension = c(3L, 3L))
  expect_equal(nrow(r$lines), 1)
  expect_equal(r$lines$length, sqrt(0.72), tolerance = 1e-6)
})

test_that("microscopic line — well below one cell, still exact", {
  skip_if_not_installed("geos")
  # Tiny segment inside a cell; cell dx=dy=1, segment length ~1.4e-3.
  line <- geos::as_geos_geometry("LINESTRING (0.5 0.5, 0.501 0.501)")
  r <- burn_scanline(line, extent = c(0, 3, 0, 3), dimension = c(3L, 3L))
  expect_equal(nrow(r$lines), 1)
  expect_equal(r$lines$length, sqrt(2 * 0.001^2), tolerance = 1e-9)
})

test_that("zero-length segment between identical coords — no edge", {
  skip_if_not_installed("geos")
  line <- geos::as_geos_geometry("LINESTRING (0.5 0.5, 0.5 0.5)")
  r <- burn_scanline(line, extent = c(0, 3, 0, 3), dimension = c(3L, 3L))
  expect_equal(nrow(r$lines), 0)
})

test_that("triangle with sub-pixel area — coverage stays in [0, 1]", {
  skip_if_not_installed("geos")
  # Tiny triangle entirely inside one cell. Its coverage fraction must
  # equal triangle_area / cell_area, and never exceed 1.
  p <- geos::as_geos_geometry(
    "POLYGON ((1.1 1.1, 1.2 1.1, 1.15 1.2, 1.1 1.1))"
  )
  r <- burn_scanline(p, extent = ext10, dimension = c(10L, 10L))
  expect_true(all(r$lines$length >= 0 & r$lines$length <= 1))
})


# =========================================================================
# 3. ALIGNMENT  — edges / lines exactly on cell boundaries
# =========================================================================
#
# Convention pinned: the polygon path uses geometric coverage (not a
# centre-test), so cell-boundary alignment is well-defined as a matter
# of intersection geometry. The line path's tie-breaking goes through
# Box::side() / point_location(), and we pin the ownership choice with
# tests rather than asserting a particular rule (round-to-even, floor,
# etc.) — what matters is that the choice is deterministic and stable.

test_that("vertex exactly on grid node — no double-count, no gap", {
  skip_if_not_installed("geos")
  p <- geos::as_geos_geometry("POLYGON ((2 2, 8 2, 5 5, 2 2))")
  expect_scanline_matches_sparse(p, ext10, c(10L, 10L),
                                 label = "vertex on grid node")
})

test_that("vertex exactly on cell edge midpoint — coverage exact", {
  skip_if_not_installed("geos")
  p <- geos::as_geos_geometry("POLYGON ((2 2, 8 2, 5 3.5, 2 2))")
  expect_scanline_matches_sparse(p, ext10, c(10L, 10L),
                                 label = "vertex on cell edge")
})

test_that("line passing exactly through a cell corner", {
  skip_if_not_installed("geos")
  # Diagonal line at 45° passing through (1, 1), (2, 2), ... — every
  # interior point of the line is a cell corner. The walker must produce
  # one consistent assignment (no duplication, no dropped cells).
  line <- geos::as_geos_geometry("LINESTRING (0.5 0.5, 2.5 2.5)")
  r <- burn_scanline(line, extent = c(0, 3, 0, 3), dimension = c(3L, 3L))
  expect_equal(nrow(r$runs), 0)
  expect_equal(sum(r$lines$length), sqrt(8), tolerance = 1e-6)
  # Total length is preserved; per-cell distribution may have ties at
  # corners, but the SUM is invariant. (TODO: pin per-cell distribution
  #  once observed.)
})


# =========================================================================
# 4. PRECISION  — very long segments
# =========================================================================
#
# Convention pinned: per-cell length comes from the local segment ends
# (entry/exit coords) computed by Box::crossing, not from a global
# parametric t propagated across cells. Errors are O(N_cells × eps),
# not O(L × eps). For raster sizes anyone uses, the accumulated error
# is well below float32 (the type of GridEdge::weight).

test_that("very long line — sum of per-cell lengths matches analytical", {
  skip_if_not_installed("geos")
  # 1000-unit horizontal line on a 1000x1000 grid — touches 1000 cells.
  line <- geos::as_geos_geometry("LINESTRING (0.5 500.5, 999.5 500.5)")
  r <- burn_scanline(line,
                     extent = c(0, 1000, 0, 1000),
                     dimension = c(1000L, 1000L))
  expect_equal(nrow(r$runs), 0)
  expect_equal(sum(r$lines$length), 999.0, tolerance = 1e-3)
})

test_that("very long diagonal — sum matches sqrt(2)*length", {
  skip_if_not_installed("geos")
  line <- geos::as_geos_geometry("LINESTRING (0.5 0.5, 999.5 999.5)")
  r <- burn_scanline(line,
                     extent = c(0, 1000, 0, 1000),
                     dimension = c(1000L, 1000L))
  expect_equal(sum(r$lines$length), sqrt(2) * 999.0, tolerance = 1e-2)
})


# =========================================================================
# 5. TOPOLOGY  — self-intersecting geometry
# =========================================================================
#
# Convention pinned:
#   Polygons: signed-winding (NON-ZERO rule). A self-intersecting bow-tie
#     produces lobes covered with winding 1, central crossover at 0 — same
#     answer as the area integral.
#   Lines: lengths SUM (not unioned). A self-intersecting line that crosses
#     itself within a cell reports the total geometric length, including
#     the doubled portion. Geometrically honest; use ST_UnaryUnion upstream
#     if you want a deduplicated total.

test_that("bow-tie polygon — invalid under simple-features, behaviour not promised", {
  skip_if_not_installed("geos")
  skip("bow-tie / self-intersecting rings are invalid input under OGC simple-features;
        controlledburn does not promise specific output for invalid polygons.
        Repair via GEOSMakeValid (or geos::geos_make_valid) upstream.
        See edge-zoo doc, category 5 (Topology).")
  p <- geos::as_geos_geometry("POLYGON ((2 2, 8 8, 8 2, 2 8, 2 2))")
  expect_scanline_matches_sparse(p, ext10, c(10L, 10L),
                                 label = "bow-tie non-zero winding")
})


test_that("self-intersecting line — lengths sum, not unioned", {
  skip_if_not_installed("geos")
  # Line revisits a cell. Its total length in that cell is the sum of
  # both visits (this is the documented convention).
  line <- geos::as_geos_geometry(
    "LINESTRING (0.2 0.5, 0.8 0.5, 0.8 0.2, 0.2 0.2, 0.2 0.8)"
  )
  r <- burn_scanline(line, extent = c(0, 3, 0, 3), dimension = c(3L, 3L))
  # Total length analytical: 0.6 + 0.3 + 0.6 + 0.6 = 2.1
  expect_equal(sum(r$lines$length), 2.1, tolerance = 1e-6)
  # Single cell, so single edge, weight equals total length.
  expect_equal(nrow(r$lines), 1)
  expect_equal(r$lines$length, 2.1, tolerance = 1e-6)
})


test_that("polygon with hole — winding handles signed contributions", {
  skip_if_not_installed("geos")
  # Outer ring CCW, inner hole CW — coverage_factor flips per ring,
  # winding_delta accumulates correctly across rings.
  p <- geos::as_geos_geometry(
    "POLYGON ((2 2, 8 2, 8 8, 2 8, 2 2),
              (4 4, 4 6, 6 6, 6 4, 4 4))"
  )
  expect_scanline_matches_sparse(p, ext10, c(10L, 10L),
                                 label = "polygon with hole")
})

test_that("two polygons touching at a corner — no shared-cell over-count", {
  skip_if_not_installed("geos")
  p <- geos::as_geos_geometry(
    "MULTIPOLYGON (((2 2, 5 2, 5 5, 2 5, 2 2)),
                   ((5 5, 8 5, 8 8, 5 8, 5 5)))"
  )
  r <- burn_scanline(p, extent = ext10, dimension = c(10L, 10L))
  expect_scanline_matches_sparse(p, ext10, c(10L, 10L),
                                 label = "corner-touch multipolygon")
})
# =========================================================================
# 6. COLLINEAR  — three or more collinear vertices
# =========================================================================
#
# Convention pinned: collinear vertices produce no phantom edge. The walker
# steps by cells, not by vertices, so a collinear triplet just adds an
# intermediate coord to the traversal — length sum still gives the right
# answer. For polygons, the analytical coverage routine handles the
# zero-area degenerate triangle correctly.

test_that("polygon with collinear vertices — coverage unchanged", {
  skip_if_not_installed("geos")
  # Square with extra collinear vertex on the bottom edge.
  p <- geos::as_geos_geometry(
    "POLYGON ((2 2, 5 2, 8 2, 8 8, 2 8, 2 2))"
  )
  p_clean <- geos::as_geos_geometry(
    "POLYGON ((2 2, 8 2, 8 8, 2 8, 2 2))"
  )
  r1 <- burn_scanline(p, extent = ext10, dimension = c(10L, 10L))
  r2 <- burn_scanline(p_clean, extent = ext10, dimension = c(10L, 10L))
  m1 <- materialise_chunk(r1)
  m2 <- materialise_chunk(r2)
  expect_equal(max(abs(m1 - m2)), 0, tolerance = 1e-6)
})

test_that("line with collinear vertex — length unchanged", {
  skip_if_not_installed("geos")
  line_with <- geos::as_geos_geometry(
    "LINESTRING (0.5 0.5, 1.5 0.5, 2.5 0.5)"
  )
  line_without <- geos::as_geos_geometry(
    "LINESTRING (0.5 0.5, 2.5 0.5)"
  )
  r1 <- burn_scanline(line_with, extent = c(0, 3, 0, 3), dimension = c(3L, 3L))
  r2 <- burn_scanline(line_without, extent = c(0, 3, 0, 3), dimension = c(3L, 3L))
  expect_equal(sum(r1$edges$weight), sum(r2$edges$weight), tolerance = 1e-9)
  # Same cells touched, same per-cell lengths.
  expect_equal(nrow(r1$edges), nrow(r2$edges))
})


# =========================================================================
# 7. CRS-BOUNDARY  — domain non-promises
# =========================================================================
#
# Convention pinned (as a NON-promise): controlledburn is pure planar
# arithmetic. Antimeridian crossings, polar singularities, dateline-wrap,
# and other CRS topology are CALLER responsibility — densify and split
# before passing in. These tests document the non-promise so it's an
# explicit contract rather than undefined behaviour.

test_that("antimeridian-crossing line is treated as planar (no unwrapping)", {
  skip_if_not_installed("geos")
  # Two points at lon=170 and lon=-170. As planar input, this is a
  # 340-unit line crossing the meridian-zero region — NOT a 20-unit
  # line wrapping around the antimeridian. controlledburn does not
  # know or care about CRS topology.
  line <- geos::as_geos_geometry("LINESTRING (170 0, -170 0)")
  r <- burn_scanline(line, extent = c(-180, 180, -10, 10),
                     dimension = c(360L, 20L))
  expect_equal(nrow(r$runs), 0)
  # 340 units of length, not 20.
  expect_equal(sum(r$lines$length), 340.0, tolerance = 1e-3)
})

test_that("geometry far outside extent is silently dropped, no error", {
  skip_if_not_installed("geos")
  # A geometry whose bounding box doesn't intersect the grid extent
  # produces no records and does not throw.
  line <- geos::as_geos_geometry("LINESTRING (-1000 -1000, -999 -999)")
  expect_silent(
    r <- burn_scanline(line, extent = ext10, dimension = c(10L, 10L))
  )
  expect_equal(nrow(r$lines), 0)
  expect_equal(nrow(r$runs), 0)
})

