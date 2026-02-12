test_that("adjacent rectangles: vertical shared edge", {
  skip_if_not_installed("geos")
  polys <- geos::as_geos_geometry(c(
    "POLYGON ((0 0, 5 0, 5 10, 0 10, 0 0))",
    "POLYGON ((5 0, 10 0, 10 10, 5 10, 5 0))"))
  r <- burn_scanline(polys, extent = c(0,10,0,10), dimension = c(20L,20L))
  expect_complementary(r, label = "vert on-grid")

  # Boundary between grid lines
  r2 <- burn_scanline(polys, extent = c(0,10,0,10), dimension = c(12L,12L))
  expect_complementary(r2, label = "vert off-grid")
})

test_that("adjacent rectangles: horizontal shared edge", {
  skip_if_not_installed("geos")
  polys <- geos::as_geos_geometry(c(
    "POLYGON ((0 0, 10 0, 10 5, 0 5, 0 0))",
    "POLYGON ((0 5, 10 5, 10 10, 0 10, 0 5))"))
  r <- burn_scanline(polys, extent = c(0,10,0,10), dimension = c(15L,15L))
  expect_complementary(r, label = "horiz")
})

test_that("four quadrants", {
  skip_if_not_installed("geos")
  polys <- geos::as_geos_geometry(c(
    "POLYGON ((0 0, 5 0, 5 5, 0 5, 0 0))",
    "POLYGON ((5 0, 10 0, 10 5, 5 5, 5 0))",
    "POLYGON ((0 5, 5 5, 5 10, 0 10, 0 5))",
    "POLYGON ((5 5, 10 5, 10 10, 5 10, 5 5))"))
  r <- burn_scanline(polys, extent = c(0,10,0,10), dimension = c(16L,16L))
  expect_complementary(r, label = "quadrants")
})

test_that("diagonal shared edge", {
  skip_if_not_installed("geos")
  polys <- geos::as_geos_geometry(c(
    "POLYGON ((0 0, 10 0, 10 10, 0 0))",
    "POLYGON ((0 0, 10 10, 0 10, 0 0))"))
  r <- burn_scanline(polys, extent = c(0,10,0,10), dimension = c(20L,20L))
  expect_complementary(r, label = "diagonal 20x20")

  r2 <- burn_scanline(polys, extent = c(0,10,0,10), dimension = c(100L,100L))
  expect_complementary(r2, label = "diagonal 100x100")
})

test_that("irregular tiling: L-shape + complement", {
  skip_if_not_installed("geos")
  polys <- geos::as_geos_geometry(c(
    "POLYGON ((0 0, 6 0, 6 4, 4 4, 4 10, 0 10, 0 0))",
    "POLYGON ((6 0, 10 0, 10 10, 4 10, 4 4, 6 4, 6 0))"))
  r <- burn_scanline(polys, extent = c(0,10,0,10), dimension = c(20L,20L))
  expect_complementary(r, label = "irregular")
})

test_that("donut filled with inner polygon", {
  skip_if_not_installed("geos")
  polys <- geos::as_geos_geometry(c(
    "POLYGON ((1 1, 9 1, 9 9, 1 9, 1 1), (3 3, 7 3, 7 7, 3 7, 3 3))",
    "POLYGON ((3 3, 7 3, 7 7, 3 7, 3 3))"))
  r <- burn_scanline(polys, extent = c(0,10,0,10), dimension = c(20L,20L))
  expect_complementary(r, label = "donut filled")
})

test_that("NC counties: no cell exceeds 1.0", {
  skip_if_not_installed("sf")
  nc <- sf::st_read(system.file("shape/nc.shp", package = "sf"), quiet = TRUE)
  nc_geom <- sf::st_geometry(nc[1:10, ])
  ext <- as.numeric(sf::st_bbox(nc_geom))[c(1, 3, 2, 4)]
  r <- burn_scanline(nc_geom, extent = ext, dimension = c(500L, 200L))

  ids <- unique(c(r$runs$id, r$edges$id))
  total <- matrix(0, nrow = 200, ncol = 500)
  for (id in ids) {
    total <- total + materialise_chunk(r, id = id)
  }
  expect_lte(max(total), 1.0 + 1e-5)
})
