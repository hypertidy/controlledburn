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
