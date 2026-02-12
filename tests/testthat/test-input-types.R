test_that("geos geometry input works", {
  skip_if_not_installed("geos")
  poly <- geos::as_geos_geometry(
    "POLYGON ((0.5 0.5, 2.5 0.5, 2.5 2.5, 0.5 2.5, 0.5 0.5))")
  r <- burn_scanline(poly, extent = c(0, 3, 0, 3), dimension = c(3, 3))
  expect_s3_class(r, "controlledburn")
  expect_true(nrow(r$edges) > 0 || nrow(r$runs) > 0)
})

test_that("sf sfc input works", {
  skip_if_not_installed("sf")
  poly <- sf::st_sfc(
    sf::st_polygon(list(matrix(c(0.5,0.5, 2.5,0.5, 2.5,2.5, 0.5,2.5, 0.5,0.5),
                                ncol = 2, byrow = TRUE))))
  r <- burn_scanline(poly, extent = c(0, 3, 0, 3), dimension = c(3, 3))
  expect_s3_class(r, "controlledburn")
})

test_that("wk_wkb input works", {
  skip_if_not_installed("geos")
  poly <- geos::as_geos_geometry(
    "POLYGON ((0.5 0.5, 2.5 0.5, 2.5 2.5, 0.5 2.5, 0.5 0.5))")
  wkb <- wk::wkb(geos::geos_write_wkb(poly))
  r <- burn_scanline(wkb, extent = c(0, 3, 0, 3), dimension = c(3, 3))
  expect_s3_class(r, "controlledburn")
})

test_that("raw WKB list input works", {
  skip_if_not_installed("geos")
  poly <- geos::as_geos_geometry(
    "POLYGON ((0.5 0.5, 2.5 0.5, 2.5 2.5, 0.5 2.5, 0.5 0.5))")
  wkb <- unclass(geos::geos_write_wkb(poly))
  r <- burn_scanline(wkb, extent = c(0, 3, 0, 3), dimension = c(3, 3))
  expect_s3_class(r, "controlledburn")
})

test_that("scanline and sparse agree across input types", {
  skip_if_not_installed("geos")
  skip_if_not_installed("sf")

  wkt <- "POLYGON ((0.5 0.5, 2.5 0.5, 2.5 2.5, 0.5 2.5, 0.5 0.5))"
  ext <- c(0, 3, 0, 3)
  dim <- c(10L, 10L)

  g <- geos::as_geos_geometry(wkt)
  s <- sf::st_as_sfc(wkt)
  w <- geos::geos_write_wkb(g)

  r_geos <- materialise_chunk(burn_scanline(g, ext, dim))
  r_sf   <- materialise_chunk(burn_scanline(s, ext, dim))
  r_wkb  <- materialise_chunk(burn_scanline(w, ext, dim))

  expect_equal(r_geos, r_sf)
  expect_equal(r_geos, r_wkb)
})
