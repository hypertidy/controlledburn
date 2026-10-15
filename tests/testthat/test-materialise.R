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

## ---- target: chunked materialisation ----

.wk_poly <- function(wkt) wk::as_wkb(wk::wkt(wkt))

test_that("materialise_chunk(target) equals the sub-matrix of the full materialisation", {
  poly <- .wk_poly(c("POLYGON ((1 1, 9 1, 9 9, 1 9, 1 1))",
                     "POLYGON ((0.35 0.35, 4.75 0.35, 2.55 6.25, 0.35 0.35))"))
  r <- burn(poly, extent = c(0, 10, 0, 10), dimension = c(100L, 100L))
  for (fun in c("sum", "id")) {
    full <- materialise_chunk(r, fun = fun)
    target <- c(2.04, 5.13, 3.3, 7.01)
    chunk <- materialise_chunk(r, fun = fun, target = target)
    sub <- crop_burn(r, target)
    expect_equal(dim(chunk), rev(sub$dimension))
    ## snapped outward: covers at least the target
    expect_true(sub$extent[1] <= target[1] && sub$extent[2] >= target[2])
    expect_true(sub$extent[3] <= target[3] && sub$extent[4] >= target[4])
    ## index of the window within the full grid (0-based offsets)
    col0 <- round((sub$extent[1] - r$extent[1]) / 0.1)
    row0 <- round((r$extent[4] - sub$extent[4]) / 0.1)
    expect_equal(chunk,
                 full[row0 + seq_len(sub$dimension[2]), col0 + seq_len(sub$dimension[1])])
  }
})

test_that("materialise_chunk(target) on an exact cell boundary does not drift", {
  poly <- .wk_poly("POLYGON ((1 1, 9 1, 9 9, 1 9, 1 1))")
  r <- burn(poly, extent = c(0, 10, 0, 10), dimension = c(100L, 100L))
  full <- materialise_chunk(r)
  chunk <- materialise_chunk(r, target = c(2, 5, 3, 7))
  expect_equal(dim(chunk), c(40L, 30L))
  expect_equal(chunk, full[31:70, 21:50])
})

test_that("materialise_chunk(target) covering the whole grid equals the full materialisation", {
  poly <- .wk_poly("POLYGON ((1 1, 9 1, 9 9, 1 9, 1 1))")
  r <- burn(poly, extent = c(0, 10, 0, 10), dimension = c(100L, 100L))
  expect_equal(materialise_chunk(r, target = r$extent), materialise_chunk(r))
  expect_equal(materialise_chunk(r, target = c(-5, 15, -5, 15)), materialise_chunk(r))
})

test_that("materialise_chunk(target) outside the grid errors", {
  poly <- .wk_poly("POLYGON ((1 1, 9 1, 9 9, 1 9, 1 1))")
  r <- burn(poly, extent = c(0, 10, 0, 10), dimension = c(100L, 100L))
  expect_error(suppressWarnings(materialise_chunk(r, target = c(20, 30, 20, 30))),
               "does not overlap")
})

test_that("materialise_chunk(target) with lines and points", {
  geom <- .wk_poly(c("LINESTRING (0.5 0.5, 9.5 9.5)", "POINT (3.33 6.66)"))
  r <- burn(geom, extent = c(0, 10, 0, 10), dimension = c(50L, 50L))
  full <- materialise_chunk(r)
  chunk <- materialise_chunk(r, target = c(3, 4, 6, 7))
  expect_equal(chunk, full[16:20, 16:20])
})

test_that("materialise_chunk sum accumulates overlapping interiors", {
  ## Interior runs used to be assigned 1 rather than added, so where two
  ## polygons' interiors overlapped the sum undercounted, and disagreed with
  ## extract_burn(fun = "sum"). No geos needed.
  g <- wk::as_wkb(wk::wkt(c(
    "POLYGON ((0 0, 6 0, 6 6, 0 6, 0 0))",
    "POLYGON ((4 4, 10 4, 10 10, 4 10, 4 4))",
    "POLYGON ((3.5 3.5, 7.5 3.5, 7.5 7.5, 3.5 7.5, 3.5 3.5))"
  )))
  r <- burn(g, extent = c(0, 10, 0, 10), dimension = c(10L, 10L))
  m <- materialise_chunk(r)
  expect_equal(sum(m), 36 + 36 + 16, tolerance = 1e-6)
  ## cell x 4-5, y 5-6 is inside all three
  expect_equal(m[5, 5], 3)
  ## extract_burn reproduces materialise_chunk at every cell centre
  xy <- expand.grid(x = seq(0.5, 9.5), y = seq(9.5, 0.5))
  rc <- cbind(rep(1:10, each = 10), rep(1:10, times = 10))
  expect_equal(extract_burn(r, xy, "sum"), m[rc], tolerance = 1e-6)
})

