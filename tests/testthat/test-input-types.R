test_that("geos geometry input works", {
  skip_if_not_installed("geos")
  poly <- geos::as_geos_geometry(
    "POLYGON ((0.5 0.5, 2.5 0.5, 2.5 2.5, 0.5 2.5, 0.5 0.5))")
  r <- burn(poly, extent = c(0, 3, 0, 3), dimension = c(3, 3))
  expect_s3_class(r, "controlledburn")
  expect_true(nrow(r$edges) > 0 || nrow(r$runs) > 0)
})

# test_that("sf sfc input works", {
# poly <- structure(list(structure(list(structure(c(0.5, 2.5, 2.5, 0.5,
#                                                   0.5, 0.5, 0.5, 2.5, 2.5, 0.5), dim = c(5L, 2L))), class = c("XY",
#                                                                                                               "POLYGON", "sfg"))), class = c("sfc_POLYGON", "sfc"), precision = 0, bbox = structure(c(xmin = 0.5,
#                                                                                                                                                                                                       ymin = 0.5, xmax = 2.5, ymax = 2.5), class = "bbox"), crs = structure(list(
#                                                                                                                                                                                                         input = NA_character_, wkt = NA_character_), class = "crs"), n_empty = 0L)
#   r <- burn(poly, extent = c(0, 3, 0, 3), dimension = c(3, 3))
#   expect_s3_class(r, "controlledburn")
# })

test_that("wk_wkb input works", {
  skip_if_not_installed("geos")
  poly <- geos::as_geos_geometry(
    "POLYGON ((0.5 0.5, 2.5 0.5, 2.5 2.5, 0.5 2.5, 0.5 0.5))")
  wkb <- wk::wkb(geos::geos_write_wkb(poly))
  r <- burn(wkb, extent = c(0, 3, 0, 3), dimension = c(3, 3))
  expect_s3_class(r, "controlledburn")
})

test_that("raw WKB list input works", {
  skip_if_not_installed("geos")
  poly <- geos::as_geos_geometry(
    "POLYGON ((0.5 0.5, 2.5 0.5, 2.5 2.5, 0.5 2.5, 0.5 0.5))")
  wkb <- unclass(geos::geos_write_wkb(poly))
  r <- burn(wkb, extent = c(0, 3, 0, 3), dimension = c(3, 3))
  expect_s3_class(r, "controlledburn")
})

test_that("burn agrees across input types (geos, wk_wkt, wk_wkb)", {
  skip_if_not_installed("geos")

  wkt <- "POLYGON ((0.5 0.5, 2.5 0.5, 2.5 2.5, 0.5 2.5, 0.5 0.5))"
  ext <- c(0, 3, 0, 3)
  dim <- c(10L, 10L)

  g <- geos::as_geos_geometry(wkt)
  w <- geos::geos_write_wkb(g)
  # wk_wkt -> wk_wkb via wk_handle round-trip
  w2 <- wk::wk_handle(wk::wkt(wkt), wk::wkb_writer())

  r_geos <- materialise_chunk(burn(g, ext, dim))
  r_wkb  <- materialise_chunk(burn(w, ext, dim))
  r_wk   <- materialise_chunk(burn(w2, ext, dim))

  expect_equal(r_geos, r_wkb)
  expect_equal(r_geos, r_wk)
})
