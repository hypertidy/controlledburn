test_that("default extent from geometry bbox", {
  skip_if_not_installed("geos")
  poly <- geos::as_geos_geometry(
    "POLYGON ((0.5 0.5, 2.5 0.5, 2.5 2.5, 0.5 2.5, 0.5 0.5))")
  r <- burn(poly)
  expect_equal(r$extent, c(0.5, 2.5, 0.5, 2.5))
  expect_equal(r$dimension, c(256L, 256L))
})

test_that("default dimension preserves aspect ratio", {
  skip_if_not_installed("geos")
  poly <- geos::as_geos_geometry(
    "POLYGON ((0 0, 10 0, 10 5, 0 5, 0 0))")
  r <- burn(poly)
  # wider than tall: ncol=256, nrow=128
  expect_equal(r$dimension[1], 256L)
  expect_equal(r$dimension[2], 128L)
})

test_that("resolution parameter computes dimension", {
  skip_if_not_installed("geos")
  poly <- geos::as_geos_geometry(
    "POLYGON ((0.5 0.5, 2.5 0.5, 2.5 2.5, 0.5 2.5, 0.5 0.5))")
  r <- burn(poly, resolution = 0.1)
  expect_equal(r$dimension, c(20L, 20L))
})

test_that("resolution with explicit extent", {
  skip_if_not_installed("geos")
  poly <- geos::as_geos_geometry(
    "POLYGON ((0.5 0.5, 2.5 0.5, 2.5 2.5, 0.5 2.5, 0.5 0.5))")
  r <- burn(poly, extent = c(0, 10, 0, 5), resolution = 0.5)
  expect_equal(r$dimension, c(20L, 10L))
})

test_that("dimension and resolution are mutually exclusive", {
  skip_if_not_installed("geos")
  poly <- geos::as_geos_geometry(
    "POLYGON ((0.5 0.5, 2.5 0.5, 2.5 2.5, 0.5 2.5, 0.5 0.5))")
  expect_error(
    burn(poly, dimension = c(10, 10), resolution = 0.1),
    "not both")
})

test_that("burn_sparse is deprecated", {
  skip_if_not_installed("geos")
  poly <- geos::as_geos_geometry(
    "POLYGON ((0.5 0.5, 2.5 0.5, 2.5 2.5, 0.5 2.5, 0.5 0.5))")
  expect_error(burn_sparse(poly), "no longer available")
})

test_that("non-finite extent is rejected", {
  g <- wk::as_wkb("POLYGON ((1 1, 3 1, 3 3, 1 3, 1 1))")
  ## Inf used to reach the core with an infinite cell size and crash
  expect_error(burn(g, extent = c(0, Inf, 0, 4), dimension = c(4L, 4L)))
  expect_error(burn(g, extent = c(-Inf, 4, 0, 4), dimension = c(4L, 4L)))
  expect_error(burn(g, extent = c(0, 4, NA, 4), dimension = c(4L, 4L)))
  expect_error(burn(g, extent = c(-1e308, 1e308, 0, 4), dimension = c(4L, 4L)),
               "overflow")
})
