ext10 <- c(0, 10, 0, 10)

test_that("vertex on grid node", {
  skip_if_not_installed("geos")
  p <- geos::as_geos_geometry("POLYGON ((2 2, 8 2, 5 5, 2 2))")
  expect_scanline_matches_sparse(p, ext10, c(10L,10L), label = "grid node 10")
  expect_scanline_matches_sparse(p, ext10, c(20L,20L), label = "grid node 20")
})

test_that("vertex on cell edge", {
  skip_if_not_installed("geos")
  p <- geos::as_geos_geometry("POLYGON ((2 2, 8 2, 5 3.5, 2 2))")
  expect_scanline_matches_sparse(p, ext10, c(10L,10L), label = "cell edge")
})

test_that("horizontal edges", {
  skip_if_not_installed("geos")
  p <- geos::as_geos_geometry("POLYGON ((2 3, 8 3, 8 7, 2 7, 2 3))")
  expect_scanline_matches_sparse(p, ext10, c(10L,10L), label = "horiz 10")
  expect_scanline_matches_sparse(p, ext10, c(30L,30L), label = "horiz 30")

  # Horizontal edge ON grid line
  p2 <- geos::as_geos_geometry("POLYGON ((2 5, 8 5, 8 8, 2 8, 2 5))")
  expect_scanline_matches_sparse(p2, ext10, c(10L,10L), label = "horiz on grid")
})

test_that("vertical edges", {
  skip_if_not_installed("geos")
  p <- geos::as_geos_geometry("POLYGON ((3 2, 7 2, 7 8, 3 8, 3 2))")
  expect_scanline_matches_sparse(p, ext10, c(10L,10L), label = "vert 10")

  p2 <- geos::as_geos_geometry("POLYGON ((5 2, 8 2, 8 8, 5 8, 5 2))")
  expect_scanline_matches_sparse(p2, ext10, c(10L,10L), label = "vert on grid")
})

test_that("thin slivers", {
  skip_if_not_installed("geos")
  # Horizontal
  p1 <- geos::as_geos_geometry("POLYGON ((1 4.9, 9 4.9, 9 5.1, 1 5.1, 1 4.9))")
  expect_scanline_matches_sparse(p1, ext10, c(10L,10L), label = "thin horiz")

  # Vertical
  p2 <- geos::as_geos_geometry("POLYGON ((4.9 1, 5.1 1, 5.1 9, 4.9 9, 4.9 1))")
  expect_scanline_matches_sparse(p2, ext10, c(10L,10L), label = "thin vert")

  # Diagonal
  p3 <- geos::as_geos_geometry("POLYGON ((1 1, 9 9, 8.9 9.1, 0.9 1.1, 1 1))")
  expect_scanline_matches_sparse(p3, ext10, c(20L,20L), label = "thin diag")

  # Sub-cell
  p4 <- geos::as_geos_geometry("POLYGON ((2 4.95, 8 4.95, 8 5.05, 2 5.05, 2 4.95))")
  expect_scanline_matches_sparse(p4, ext10, c(10L,10L), label = "sub-cell")
})

test_that("polygon within one cell", {
  skip_if_not_installed("geos")
  p1 <- geos::as_geos_geometry("POLYGON ((4.2 4.2, 4.8 4.2, 4.8 4.8, 4.2 4.8, 4.2 4.2))")
  expect_scanline_matches_sparse(p1, ext10, c(10L,10L), label = "rect in cell")

  p2 <- geos::as_geos_geometry("POLYGON ((4.3 4.3, 4.7 4.3, 4.5 4.7, 4.3 4.3))")
  expect_scanline_matches_sparse(p2, ext10, c(10L,10L), label = "tri in cell")
})

test_that("edge along grid line", {
  skip_if_not_installed("geos")
  p <- geos::as_geos_geometry("POLYGON ((0 0, 10 0, 10 5, 0 5, 0 0))")
  expect_scanline_matches_sparse(p, ext10, c(10L,10L), label = "edge on grid")
})

test_that("near-degenerate shapes", {
  skip_if_not_installed("geos")
  # Acute triangle
  p1 <- geos::as_geos_geometry("POLYGON ((1 5, 9 4.99, 9 5.01, 1 5))")
  expect_scanline_matches_sparse(p1, ext10, c(20L,20L), label = "acute")

  # Needle
  p2 <- geos::as_geos_geometry("POLYGON ((1 5, 9 5.001, 9 4.999, 1 5))")
  expect_scanline_matches_sparse(p2, ext10, c(20L,20L), label = "needle")
})

test_that("extreme resolution ratios", {
  skip_if_not_installed("geos")
  p <- geos::as_geos_geometry("POLYGON ((1 1, 9 1, 9 9, 1 9, 1 1))")
  expect_scanline_matches_sparse(p, ext10, c(1L,1L), label = "1x1")
  expect_scanline_matches_sparse(p, ext10, c(2L,2L), label = "2x2")
  expect_scanline_matches_sparse(p, ext10, c(500L,500L), label = "500x500")
})

test_that("polygon at grid extent boundary", {
  skip_if_not_installed("geos")
  # Exact match
  p1 <- geos::as_geos_geometry("POLYGON ((0 0, 10 0, 10 10, 0 10, 0 0))")
  expect_scanline_matches_sparse(p1, ext10, c(10L,10L), label = "at extent")

  # Beyond extent
  p2 <- geos::as_geos_geometry("POLYGON ((-1 -1, 11 -1, 11 11, -1 11, -1 -1))")
  expect_scanline_matches_sparse(p2, ext10, c(10L,10L), label = "beyond extent")

  # Partial outside
  p3 <- geos::as_geos_geometry("POLYGON ((5 5, 15 5, 15 15, 5 15, 5 5))")
  expect_scanline_matches_sparse(p3, ext10, c(10L,10L), label = "partial outside")
})

test_that("collinear vertices", {
  skip_if_not_installed("geos")
  p <- geos::as_geos_geometry(
    "POLYGON ((2 2, 5 2, 8 2, 8 5, 8 8, 5 8, 2 8, 2 5, 2 2))")
  expect_scanline_matches_sparse(p, ext10, c(10L,10L), label = "collinear")
})
