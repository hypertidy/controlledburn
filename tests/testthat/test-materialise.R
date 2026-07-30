test_that("materialise_chunk produces correct dimensions", {
  skip_if_not_installed("geos")
  poly <- geos::as_geos_geometry(
    "POLYGON ((0.5 0.5, 2.5 0.5, 2.5 2.5, 0.5 2.5, 0.5 0.5))")
  r <- burn(poly, extent = c(0, 3, 0, 3), dimension = c(3, 3))
  mat <- materialise_chunk(r)
  expect_equal(dim(mat), c(3, 3))
})

test_that("materialise_chunk sums runs and edges", {
  skip_if_not_installed("geos")
  poly <- geos::as_geos_geometry(
    "POLYGON ((0.5 0.5, 2.5 0.5, 2.5 2.5, 0.5 2.5, 0.5 0.5))")
  r <- burn(poly, extent = c(0, 3, 0, 3), dimension = c(3, 3))
  mat <- materialise_chunk(r)
  # Total coverage should equal polygon area (2x2 = 4)
  cell_area <- 1.0
  expect_equal(sum(mat) * cell_area, 4.0, tolerance = 1e-6)
})

test_that("materialise_chunk handles line output", {
  skip_if_not_installed("geos")
  line <- geos::as_geos_geometry("LINESTRING (0.5 0.5, 2.5 2.5)")
  r <- burn(line, extent = c(0, 3, 0, 3), dimension = c(3, 3))
  mat <- materialise_chunk(r)
  expect_equal(dim(mat), c(3, 3))
  expect_true(sum(mat) > 0)
})

test_that("materialise_chunk handles point output", {
  skip_if_not_installed("geos")
  pt <- geos::as_geos_geometry("POINT (1.5 1.5)")
  r <- burn(pt, extent = c(0, 3, 0, 3), dimension = c(3, 3))
  mat <- materialise_chunk(r)
  expect_equal(sum(mat), 1.0)
})

test_that("materialize_chunk is an alias for materialise_chunk", {
  skip_if_not_installed("geos")
  poly <- geos::as_geos_geometry(
    "POLYGON ((0.5 0.5, 2.5 0.5, 2.5 2.5, 0.5 2.5, 0.5 0.5))")
  r <- burn(poly, extent = c(0, 3, 0, 3), dimension = c(3, 3))
  expect_equal(materialise_chunk(r), materialize_chunk(r))
})

test_that("print.controlledburn works", {
  skip_if_not_installed("geos")
  poly <- geos::as_geos_geometry(
    "POLYGON ((0.5 0.5, 2.5 0.5, 2.5 2.5, 0.5 2.5, 0.5 0.5))")
  r <- burn(poly, extent = c(0, 3, 0, 3), dimension = c(3, 3))
  expect_output(print(r), "controlledburn")
})
