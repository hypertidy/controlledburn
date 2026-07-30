# Edge case geometry tests — verify area conservation and no crashes.

skip_if_not_installed("geos")

ext10 <- c(0, 10, 0, 10)

# Helper: total coverage from a burn result
total_coverage <- function(r) {
  cell_area <- ((r$extent[2] - r$extent[1]) / r$dimension[1]) *
               ((r$extent[4] - r$extent[3]) / r$dimension[2])
  run_cells <- if (nrow(r$runs) > 0)
    sum(r$runs$col_end - r$runs$col_start + 1) else 0
  edge_frac <- if (nrow(r$edges) > 0) sum(r$edges$fraction) else 0
  (run_cells + edge_frac) * cell_area
}

# Helper: expected polygon area via shoelace formula
poly_area <- function(wkt) {
  g <- geos::as_geos_geometry(wkt)
  geos::geos_area(g)
}

# Helper: burn and check area conservation.
# For polygons extending beyond the grid, expected area is clipped.
expect_area_conserved <- function(wkt, ext, dim, tol = 1e-4, label = NULL) {
  g <- geos::as_geos_geometry(wkt)
  r <- burn(g, extent = ext, dimension = dim)
  total <- total_coverage(r)
  # Clip polygon to grid extent for the expected area
  grid_box <- geos::as_geos_geometry(sprintf(
    "POLYGON ((%s %s, %s %s, %s %s, %s %s, %s %s))",
    ext[1], ext[3], ext[2], ext[3], ext[2], ext[4], ext[1], ext[4], ext[1], ext[3]
  ))
  clipped <- geos::geos_intersection(g, grid_box)
  expected <- geos::geos_area(clipped)
  expect_equal(total, expected, tolerance = tol, label = label)
}

test_that("vertex on grid node", {
  expect_area_conserved("POLYGON ((2 2, 8 2, 5 5, 2 2))",
                        ext10, c(10L, 10L), label = "grid node 10")
  expect_area_conserved("POLYGON ((2 2, 8 2, 5 5, 2 2))",
                        ext10, c(20L, 20L), label = "grid node 20")
})

test_that("vertex on cell edge", {
  expect_area_conserved("POLYGON ((2 2, 8 2, 5 3.5, 2 2))",
                        ext10, c(10L, 10L), label = "cell edge")
})

test_that("horizontal edges", {
  expect_area_conserved("POLYGON ((2 3, 8 3, 8 7, 2 7, 2 3))",
                        ext10, c(10L, 10L), label = "horiz 10")
  expect_area_conserved("POLYGON ((2 3, 8 3, 8 7, 2 7, 2 3))",
                        ext10, c(30L, 30L), label = "horiz 30")
  expect_area_conserved("POLYGON ((2 5, 8 5, 8 8, 2 8, 2 5))",
                        ext10, c(10L, 10L), label = "horiz on grid")
})

test_that("vertical edges", {
  expect_area_conserved("POLYGON ((3 2, 7 2, 7 8, 3 8, 3 2))",
                        ext10, c(10L, 10L), label = "vert 10")
  expect_area_conserved("POLYGON ((5 2, 8 2, 8 8, 5 8, 5 2))",
                        ext10, c(10L, 10L), label = "vert on grid")
})

test_that("thin slivers", {
  expect_area_conserved("POLYGON ((1 4.9, 9 4.9, 9 5.1, 1 5.1, 1 4.9))",
                        ext10, c(10L, 10L), label = "thin horiz")
  expect_area_conserved("POLYGON ((4.9 1, 5.1 1, 5.1 9, 4.9 9, 4.9 1))",
                        ext10, c(10L, 10L), label = "thin vert")
  expect_area_conserved("POLYGON ((1 1, 9 9, 8.9 9.1, 0.9 1.1, 1 1))",
                        ext10, c(20L, 20L), label = "thin diag")
  expect_area_conserved("POLYGON ((2 4.95, 8 4.95, 8 5.05, 2 5.05, 2 4.95))",
                        ext10, c(10L, 10L), label = "sub-cell")
})

test_that("polygon within one cell", {
  expect_area_conserved("POLYGON ((4.2 4.2, 4.8 4.2, 4.8 4.8, 4.2 4.8, 4.2 4.2))",
                        ext10, c(10L, 10L), label = "rect in cell")
  expect_area_conserved("POLYGON ((4.3 4.3, 4.7 4.3, 4.5 4.7, 4.3 4.3))",
                        ext10, c(10L, 10L), label = "tri in cell")
})

test_that("edge along grid line", {
  expect_area_conserved("POLYGON ((0 0, 10 0, 10 5, 0 5, 0 0))",
                        ext10, c(10L, 10L), label = "edge on grid")
})

test_that("near-degenerate shapes", {
  expect_area_conserved("POLYGON ((1 5, 9 4.99, 9 5.01, 1 5))",
                        ext10, c(20L, 20L), label = "acute")
  expect_area_conserved("POLYGON ((1 5, 9 5.001, 9 4.999, 1 5))",
                        ext10, c(20L, 20L), label = "needle")
})

test_that("extreme resolution ratios", {
  wkt <- "POLYGON ((1 1, 9 1, 9 9, 1 9, 1 1))"
  expect_area_conserved(wkt, ext10, c(1L, 1L), label = "1x1")
  expect_area_conserved(wkt, ext10, c(2L, 2L), label = "2x2")
  expect_area_conserved(wkt, ext10, c(500L, 500L), label = "500x500")
})

test_that("polygon at grid extent boundary", {
  expect_area_conserved("POLYGON ((0 0, 10 0, 10 10, 0 10, 0 0))",
                        ext10, c(10L, 10L), label = "at extent")
  expect_area_conserved("POLYGON ((-1 -1, 11 -1, 11 11, -1 11, -1 -1))",
                        ext10, c(10L, 10L), label = "beyond extent")
  expect_area_conserved("POLYGON ((5 5, 15 5, 15 15, 5 15, 5 5))",
                        ext10, c(10L, 10L), label = "partial outside")
})

test_that("collinear vertices", {
  expect_area_conserved(
    "POLYGON ((2 2, 5 2, 8 2, 8 5, 8 8, 5 8, 2 8, 2 5, 2 2))",
    ext10, c(10L, 10L), label = "collinear")
})
