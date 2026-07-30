test_that("scanline matches sparse: simple rectangle", {
  skip_if_not_installed("geos")
  poly <- geos::as_geos_geometry(
    "POLYGON ((0.5 0.5, 2.5 0.5, 2.5 2.5, 0.5 2.5, 0.5 0.5))")
  expect_scanline_matches_sparse(poly, c(0,3,0,3), c(3L,3L), label = "rect 3x3")
  expect_scanline_matches_sparse(poly, c(0,3,0,3), c(30L,30L), label = "rect 30x30")
})

test_that("scanline matches sparse: triangle", {
  skip_if_not_installed("geos")
  poly <- geos::as_geos_geometry(
    "POLYGON ((1 1, 9 1, 5 9, 1 1))")
  expect_scanline_matches_sparse(poly, c(0,10,0,10), c(20L,20L), label = "triangle")
})

test_that("scanline matches sparse: polygon with hole", {
  skip_if_not_installed("geos")
  poly <- geos::as_geos_geometry(
    "POLYGON ((1 1, 9 1, 9 9, 1 9, 1 1), (3 3, 7 3, 7 7, 3 7, 3 3))")
  expect_scanline_matches_sparse(poly, c(0,10,0,10), c(20L,20L), label = "donut")
})

test_that("scanline matches sparse: L-shaped polygon", {
  skip_if_not_installed("geos")
  poly <- geos::as_geos_geometry(
    "POLYGON ((0 0, 6 0, 6 4, 4 4, 4 10, 0 10, 0 0))")
  expect_scanline_matches_sparse(poly, c(-1,11,-1,11), c(20L,20L), label = "L-shape")
})

test_that("scanline matches sparse: multiple polygons", {
  skip_if_not_installed("geos")
  polys <- geos::as_geos_geometry(c(
    "POLYGON ((1 1, 4 1, 4 4, 1 4, 1 1))",
    "POLYGON ((6 6, 9 6, 9 9, 6 9, 6 6))"))
  expect_scanline_matches_sparse(polys, c(0,10,0,10), c(20L,20L), label = "multi")
})

test_that("scanline matches sparse: diamond (all diagonal)", {
  skip_if_not_installed("geos")
  poly <- geos::as_geos_geometry(
    "POLYGON ((5 0, 10 5, 5 10, 0 5, 5 0))")
  expect_scanline_matches_sparse(poly, c(-1,11,-1,11), c(40L,40L), label = "diamond")
})

test_that("scanline matches sparse: NC counties", {
  skip_if_not_installed("vapour")
  nc_path <- system.file("shape/nc.shp", package = "sf")
  skip_if_not(file.exists(nc_path), "NC shapefile not available")
  nc_g <- geos::as_geos_geometry(wk::wkb(vapour::vapour_read_geometry(nc_path)))
  nc5 <- nc_g[1:5]
  bb <- geos::geos_extent(geos::geos_make_collection(nc5))
  ext <- c(bb$xmin, bb$xmax, bb$ymin, bb$ymax)
  expect_scanline_matches_sparse(nc5, ext, c(500L, 200L),
                                 tol = 1e-4, label = "nc5")
})
