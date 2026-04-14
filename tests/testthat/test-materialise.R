test_that("materialise_chunk produces correct dimensions", {
  skip_if_not_installed("geos")
  poly <- geos::as_geos_geometry(
    "POLYGON ((0.5 0.5, 2.5 0.5, 2.5 2.5, 0.5 2.5, 0.5 0.5))")
  r <- burn_scanline(poly, extent = c(0, 3, 0, 3), dimension = c(3, 3))
  mat <- materialise_chunk(r)
  expect_equal(dim(mat), c(3, 3))  # nrow x ncol
})

test_that("materialise_chunk by id", {
  skip_if_not_installed("geos")
  polys <- geos::as_geos_geometry(c(
    "POLYGON ((1 1, 4 1, 4 4, 1 4, 1 1))",
    "POLYGON ((6 6, 9 6, 9 9, 6 9, 6 6))"))
  r <- burn_scanline(polys, extent = c(0, 10, 0, 10), dimension = c(10L, 10L))

  mat1 <- materialise_chunk(r, id = 1)
  mat2 <- materialise_chunk(r, id = 2)
  mat_all <- materialise_chunk(r)

  # Each polygon contributes independently
  expect_true(all(mat1[1:4, 7:10] == 0))  # polygon 1 not in top-right
  expect_true(all(mat2[7:10, 1:4] == 0))  # polygon 2 not in bottom-left
  # Combined equals sum
  expect_equal(mat_all, mat1 + mat2)
})

test_that("materialise_chunk vector output", {
  skip_if_not_installed("geos")
  poly <- geos::as_geos_geometry(
    "POLYGON ((0.5 0.5, 2.5 0.5, 2.5 2.5, 0.5 2.5, 0.5 0.5))")
  r <- burn_scanline(poly, extent = c(0, 3, 0, 3), dimension = c(3, 3))
  v <- materialise_chunk(r, type = "vector")
  mat <- materialise_chunk(r, type = "matrix")
  expect_equal(length(v), 9L)
  expect_equal(v, as.vector(t(mat)))  # row-major
})

test_that("print.controlledburn works", {
  skip_if_not_installed("geos")
  poly <- geos::as_geos_geometry(
    "POLYGON ((0.5 0.5, 2.5 0.5, 2.5 2.5, 0.5 2.5, 0.5 0.5))")
  r <- burn_scanline(poly, extent = c(0, 3, 0, 3), dimension = c(3, 3))
  expect_output(print(r), "controlledburn")
  expect_output(print(r), "3 x 3")
})

test_that("materialise_chunk with target extent (subwindow)", {
  skip_if_not_installed("geos")
  poly <- geos::as_geos_geometry(
    "POLYGON ((1 1, 9 1, 9 9, 1 9, 1 1))")
  r <- burn_scanline(poly, extent = c(0, 10, 0, 10), dimension = c(10L, 10L))
  mat_full <- materialise_chunk(r)

  # Request a subwindow that snaps to cell boundaries
  mat_sub <- materialise_chunk(r, target = c(3, 7, 3, 7))
  expect_equal(dim(mat_sub), c(4, 4))
  expect_equal(attr(mat_sub, "extent"), c(3, 7, 3, 7))

  # Content should match the corresponding full-grid submatrix
  # rows 4:7 (y 3-7 in a 10x10 grid with ymax=10), cols 4:7
  expect_equal(mat_sub, mat_full[4:7, 4:7])
})

test_that("materialise_chunk target snaps outward", {
  skip_if_not_installed("geos")
  poly <- geos::as_geos_geometry(
    "POLYGON ((1 1, 9 1, 9 9, 1 9, 1 1))")
  r <- burn_scanline(poly, extent = c(0, 10, 0, 10), dimension = c(10L, 10L))

  # Request extent that doesn't align — should snap out
  mat_sub <- materialise_chunk(r, target = c(3.3, 6.7, 3.3, 6.7))
  # snap out: xmin 3.3 -> 3, xmax 6.7 -> 7, same for y
  expect_equal(dim(mat_sub), c(4, 4))
  expect_equal(attr(mat_sub, "extent"), c(3, 7, 3, 7))
})

test_that("materialise_chunk target clamps to source", {
  skip_if_not_installed("geos")
  poly <- geos::as_geos_geometry(
    "POLYGON ((1 1, 9 1, 9 9, 1 9, 1 1))")
  r <- burn_scanline(poly, extent = c(0, 10, 0, 10), dimension = c(10L, 10L))

  # Request extends beyond source
  mat_sub <- materialise_chunk(r, target = c(-5, 15, 2, 8))
  # clamp to source: xmin 0, xmax 10, snap y: 2, 8
  expect_equal(dim(mat_sub), c(6, 10))
  expect_equal(attr(mat_sub, "extent"), c(0, 10, 2, 8))
})

test_that("materialise_chunk max_cells safety", {
  skip_if_not_installed("geos")
  poly <- geos::as_geos_geometry(
    "POLYGON ((1 1, 9 1, 9 9, 1 9, 1 1))")
  r <- burn_scanline(poly, extent = c(0, 10, 0, 10), dimension = c(10L, 10L))

  expect_error(materialise_chunk(r, max_cells = 10), "max_cells")
})

test_that("materialise_chunk target outside source warns", {
  skip_if_not_installed("geos")
  poly <- geos::as_geos_geometry(
    "POLYGON ((1 1, 9 1, 9 9, 1 9, 1 1))")
  r <- burn_scanline(poly, extent = c(0, 10, 0, 10), dimension = c(10L, 10L))

  expect_warning(materialise_chunk(r, target = c(20, 30, 20, 30)), "does not intersect")
})

test_that("materialise_chunk target with id filter", {
  skip_if_not_installed("geos")
  polys <- geos::as_geos_geometry(c(
    "POLYGON ((1 1, 4 1, 4 4, 1 4, 1 1))",
    "POLYGON ((6 6, 9 6, 9 9, 6 9, 6 6))"))
  r <- burn_scanline(polys, extent = c(0, 10, 0, 10), dimension = c(10L, 10L))

  # Window covering only polygon 2, filtered to id 2
  mat <- materialise_chunk(r, target = c(5, 10, 5, 10), id = 2)
  expect_equal(dim(mat), c(5, 5))
  expect_true(sum(mat) > 0)

  # Same window, filtered to id 1 — polygon 1 is outside this window
  mat1 <- materialise_chunk(r, target = c(5, 10, 5, 10), id = 1)
  expect_equal(sum(mat1), 0)
})
