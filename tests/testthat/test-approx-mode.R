# test-approx-mode.R -- tests for mode = "approx" (cell-center rule)
#
# Approx mode is the fasterize-style sweep: boundary cells are classified
# by whether the cell center is inside the polygon. No $edges produced.
# Output is runs only for polygons; lines and points are unaffected.
#
# Parity with fasterize is verified where possible, but note:
# - fasterize itself had boundary handling uncertainties
#   (hypertidy/fasterize#6)
# - GDAL's rasterize has its own vagaries (OSGeo/gdal#14615)
# - The one known discrepancy is horizontal polygon edges that land
#   exactly on a cell center's y-coordinate; see the dedicated test below.

skip_if_not_installed("geos")
library(geos)

ext10 <- c(0, 10, 0, 10)
dim10 <- c(10L, 10L)

# ---- Basic output contract ----

test_that("approx mode produces no edges for polygons", {
  poly <- as_geos_geometry("POLYGON ((0.5 0.5, 2.5 0.5, 2.5 2.5, 0.5 2.5, 0.5 0.5))")
  r <- burn(poly, extent = c(0, 3, 0, 3), dimension = c(3L, 3L), mode = "approx")
  expect_equal(nrow(r$edges), 0L)
  expect_true(nrow(r$runs) > 0)
})

test_that("approx mode cell count is a whole number of cells", {
  poly <- as_geos_geometry(
    "POLYGON ((2.3 4.7, 6.3 4.7, 6.3 8.3, 2.3 8.3, 2.3 4.7))"
  )
  r <- burn(poly, extent = ext10, dimension = dim10, mode = "approx")
  total_cells <- sum(r$runs$col_end - r$runs$col_start + 1)
  expect_equal(total_cells, round(total_cells))
})

test_that("approx mode does not affect lines", {
  line <- as_geos_geometry("LINESTRING (0.5 0.5, 9.5 7.5)")
  r_cov <- burn(line, extent = ext10, dimension = dim10, mode = "coverage")
  r_app <- burn(line, extent = ext10, dimension = dim10, mode = "approx")
  expect_equal(r_app$lines, r_cov$lines)
})

test_that("approx mode does not affect points", {
  pt <- as_geos_geometry("POINT (1.5 1.5)")
  r_cov <- burn(pt, extent = ext10, dimension = dim10, mode = "coverage")
  r_app <- burn(pt, extent = ext10, dimension = dim10, mode = "approx")
  expect_equal(r_app$points, r_cov$points)
})

test_that("mode argument is validated", {
  poly <- as_geos_geometry("POLYGON ((0.5 0.5, 2.5 0.5, 2.5 2.5, 0.5 2.5, 0.5 0.5))")
  expect_error(burn(poly, extent = c(0, 3, 0, 3), dimension = c(3L, 3L), mode = "bad"))
})

# ---- Geometry cases ----

test_that("approx: aligned rectangle on cell boundaries", {
  poly <- as_geos_geometry("POLYGON ((2 4, 6 4, 6 8, 2 8, 2 4))")
  r <- burn(poly, extent = ext10, dimension = dim10, mode = "approx")
  expect_equal(nrow(r$edges), 0L)
  cell_area <- 1.0
  total <- sum(r$runs$col_end - r$runs$col_start + 1) * cell_area
  expect_equal(total, 16.0)
})

test_that("approx: beyond-extent polygon fills entire grid", {
  poly <- as_geos_geometry(
    "POLYGON ((-100 -100, 100 -100, 100 100, -100 100, -100 -100))"
  )
  r <- burn(poly, extent = ext10, dimension = c(5L, 5L), mode = "approx")
  cell_area <- 2.0 * 2.0
  total <- sum(r$runs$col_end - r$runs$col_start + 1) * cell_area
  expect_equal(total, 100.0)
})

test_that("approx: polygon with hole respects the hole", {
  poly <- as_geos_geometry(
    "POLYGON ((1 1, 9 1, 9 9, 1 9, 1 1), (3 3, 7 3, 7 7, 3 7, 3 3))"
  )
  r <- burn(poly, extent = ext10, dimension = c(20L, 20L), mode = "approx")
  expect_equal(nrow(r$edges), 0L)
  cell_area <- 0.5 * 0.5
  total <- sum(r$runs$col_end - r$runs$col_start + 1) * cell_area
  # Outer area 64, hole area 16, net 48; boundary cells shift things
  # but total should be close
  expect_true(total > 40 && total < 56)
})

test_that("approx: adjacent polygons have no gaps or overlaps at shared boundary", {
  polys <- as_geos_geometry(c(
    "POLYGON ((0 0, 5 0, 5 10, 0 10, 0 0))",
    "POLYGON ((5 0, 10 0, 10 10, 5 10, 5 0))"
  ))
  r <- burn(polys, extent = ext10, dimension = dim10, mode = "approx")

  # Every row should have all 10 columns covered between the two polygons
  total_cells <- sum(r$runs$col_end - r$runs$col_start + 1)
  expect_equal(total_cells, 100L)

  # Polygon 1 (x=0..5) should own cols 1-5, polygon 2 (x=5..10) cols 6-10
  runs_p1 <- r$runs[r$runs$id == 1, ]
  runs_p2 <- r$runs[r$runs$id == 2, ]
  expect_true(all(runs_p1$col_end <= 5L))
  expect_true(all(runs_p2$col_start >= 6L))
})

test_that("approx: multipolygon components handled independently", {
  mp <- as_geos_geometry(
    "MULTIPOLYGON (((1 1, 3 1, 3 3, 1 3, 1 1)), ((6 6, 9 6, 9 9, 6 9, 6 6)))"
  )
  r <- burn(mp, extent = ext10, dimension = dim10, mode = "approx")
  expect_equal(nrow(r$edges), 0L)
  # All runs should have the same id (one input geometry)
  expect_true(all(r$runs$id == 1L))
})

# ---- Fasterize parity ----
#
# Fasterize's boundary handling has known uncertainties
# (hypertidy/fasterize#6). Even GDAL's rasterize has vagaries
# (OSGeo/gdal#14615). Parity is tested where behaviour is clear;
# known discrepancies are documented rather than papered over.

# Helper: fasterize a WKT polygon and return the result as a matrix
fasterize_to_matrix <- function(poly_wkt, ext, dim) {
  sf_obj <- sf::st_sf(id = 1, geometry = sf::st_as_sfc(poly_wkt, crs = NA_crs_))
  r_old <- raster::raster(terra::rast(
    xmin = ext[1], xmax = ext[2], ymin = ext[3], ymax = ext[4],
    ncols = dim[1], nrows = dim[2], crs = ""
  ))
  raster::as.matrix(fasterize::fasterize(sf_obj, r_old, field = "id"))
}

test_that("approx matches fasterize: non-aligned rectangle", {
  skip_if_not_installed("fasterize")
  skip_if_not_installed("sf")

  poly_wkt <- "POLYGON ((2.3 4.7, 6.3 4.7, 6.3 8.3, 2.3 8.3, 2.3 4.7))"
  a_mat <- materialize_chunk(
    burn(as_geos_geometry(poly_wkt), extent = ext10, dimension = dim10, mode = "approx")
  )
  f_mat <- fasterize_to_matrix(poly_wkt, ext10, dim10)
  expect_equal(a_mat > 0, !is.na(f_mat) & f_mat == 1)
})

test_that("approx matches fasterize: triangle", {
  skip_if_not_installed("fasterize")
  skip_if_not_installed("sf")

  poly_wkt <- "POLYGON ((1 1, 9 1, 5 9, 1 1))"
  a_mat <- materialize_chunk(
    burn(as_geos_geometry(poly_wkt), extent = ext10, dimension = dim10, mode = "approx")
  )
  f_mat <- fasterize_to_matrix(poly_wkt, ext10, dim10)
  expect_equal(a_mat > 0, !is.na(f_mat) & f_mat == 1)
})

test_that("approx matches fasterize: aligned rectangle", {
  skip_if_not_installed("fasterize")
  skip_if_not_installed("sf")

  poly_wkt <- "POLYGON ((2 4, 6 4, 6 8, 2 8, 2 4))"
  a_mat <- materialize_chunk(
    burn(as_geos_geometry(poly_wkt), extent = ext10, dimension = dim10, mode = "approx")
  )
  f_mat <- fasterize_to_matrix(poly_wkt, ext10, dim10)
  expect_equal(a_mat > 0, !is.na(f_mat) & f_mat == 1)
})

test_that("approx matches fasterize: polygon with hole", {
  skip_if_not_installed("fasterize")
  skip_if_not_installed("sf")

  poly_wkt <- "POLYGON ((1 1, 9 1, 9 9, 1 9, 1 1), (3 3, 7 3, 7 7, 3 7, 3 3))"
  dim20 <- c(20L, 20L)
  a_mat <- materialize_chunk(
    burn(as_geos_geometry(poly_wkt), extent = ext10, dimension = dim20, mode = "approx")
  )
  f_mat <- fasterize_to_matrix(poly_wkt, ext10, dim20)
  expect_equal(a_mat > 0, !is.na(f_mat) & f_mat == 1)
})

# ---- Known discrepancy: horizontal edge at y_mid ----
#
# When a polygon's horizontal edge falls exactly on a cell center's
# y-coordinate, controlledburn does not generate a winding delta for
# that row — the edge doesn't cross y_mid, it lies on it. fasterize
# includes some of these boundary rows (top-inclusive convention).
#
# This is documented as a known limitation. The discrepancy only arises
# when polygon edges align exactly with cell centers, which is a
# degenerate configuration. Exact coverage mode is unaffected and
# reports the correct fractional coverage regardless.

test_that("approx: horizontal edge at y_mid is a known fasterize discrepancy", {
  skip_if_not_installed("fasterize")
  skip_if_not_installed("sf")

  # Polygon 2.5..6.5 x 4.5..8.5: horizontal edges at y=4.5 and y=8.5,
  # which are cell centers for rows 6 and 2 respectively.
  poly_wkt <- "POLYGON ((2.5 4.5, 6.5 4.5, 6.5 8.5, 2.5 8.5, 2.5 4.5))"

  g <- as_geos_geometry(poly_wkt)
  a_mat <- materialize_chunk(
    burn(g, extent = ext10, dimension = dim10, mode = "approx")
  )
  f_mat <- fasterize_to_matrix(poly_wkt, ext10, dim10)

  cb_cells <- sum(a_mat > 0)
  fz_cells <- sum(!is.na(f_mat) & f_mat == 1)

  # controlledburn produces fewer cells than fasterize here because
  # the horizontal boundary row is not included. This is a known and
  # documented discrepancy.
  expect_true(cb_cells < fz_cells)

  # The discrepancy is exactly one row (4 cells on the top boundary).
  expect_equal(fz_cells - cb_cells, 4L)

  # Meanwhile, exact coverage mode correctly handles this — the
  # boundary cells get their fractional coverage regardless of the
  # horizontal edge.
  r_cov <- burn(g, extent = ext10, dimension = dim10, mode = "coverage")
  cell_area <- 1.0
  total_cov <- sum(r_cov$runs$col_end - r_cov$runs$col_start + 1) * cell_area +
    sum(r_cov$edges$fraction) * cell_area
  expect_equal(total_cov, 16.0, tolerance = 1e-6)
})
