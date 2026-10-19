# Vertices within a few ulp of a cell edge are snapped onto it before the
# coverage walk (see snap_to_lattice in src/core_scanline.cpp). No geos needed.

test_that("vertex a few ulp off a cell edge does not derail the coverage walk", {
  ## 0.3 parsed from input vs 3 * 0.1 computed from the grid differ by ~5e-17.
  ## Before snapping, the walk left the grid and ran until bad_alloc
  ## (exactextractr crashes the session on the same input).
  tri <- function(a, b, c) abs((b[1] - a[1]) * (c[2] - a[2]) - (c[1] - a[1]) * (b[2] - a[2])) / 2
  cases <- list(
    c("POLYGON ((0.3 0.3, 4.7 0.3, 2.5 6.2, 0.3 0.3))", tri(c(0.3, 0.3), c(4.7, 0.3), c(2.5, 6.2))),
    c("POLYGON ((4.7 0.3, 2.5 6.2, 0.3 0.3, 4.7 0.3))", tri(c(0.3, 0.3), c(4.7, 0.3), c(2.5, 6.2))),
    c("POLYGON ((0.3 0.3, 4.7 0.3, 2.5 6.25, 0.3 0.3))", tri(c(0.3, 0.3), c(4.7, 0.3), c(2.5, 6.25)))
  )
  for (cs in cases) {
    g <- wk::as_wkb(wk::wkt(cs[[1]]))
    r <- burn(g, extent = c(0, 10, 0, 10), dimension = c(100L, 100L))
    expect_gt(nrow(r$runs), 0)
    expect_equal(sum(materialise_chunk(r)) * 0.01, as.numeric(cs[[2]]), tolerance = 1e-6)
  }
})

test_that("a row-centre crossing on a column wall keeps its winding", {
  ## Falkland Islands from Natural Earth 110m, shifted so every vertex sits
  ## on a 0.05 multiple. At 20 columns (res 0.5) the first edge crosses the
  ## row centre y = 1.25 at x = 1.0, exactly on a column wall; with a strict
  ## crossing test neither adjacent cell registered it and the row lost two
  ## interior cells. Exact at every resolution now.
  wkt <- "POLYGON ((0.8 1.15, 2.0 1.75, 2.85 1.5, 3.45 1.9, 4.25 1.45, 3.95 1.1, 2.6 0.8, 2.15 1.15, 1.3 0.7, 0.8 1.15))"
  g <- wk::as_wkb(wk::wkt(wkt))
  for (n in c(20L, 100L, 200L, 1000L)) {
    r <- burn(g, extent = c(0, 10, 0, 10), dimension = c(n, n))
    cell <- (10 / n)^2
    a <- (sum(r$runs$col_end - r$runs$col_start + 1) + sum(r$edges$fraction)) * cell
    expect_equal(a, 2.12875, tolerance = 1e-6, label = sprintf("ncol %d", n))
  }
})

test_that("approx: a vertex on a row centre line does not fill the row to the grid edge", {
  ## Found by the controlledburn-rs property suite. The middle vertex sits on
  ## the centre line of row 24 as a user would compute it (ymax - 24.5 / 3),
  ## one ulp below the same line recomputed from the extent-derived dy. The
  ## row-range test admitted the downward edge starting there; the sign test
  ## (y0 >= y_mid) called it upward; the two intercepts then had the same
  ## sign and the row filled from the polygon to the right edge of the grid.
  wkt <- "POLYGON ((-179.89654756364214 -88.83333333333334, -180.94076559960791 -92.3680897750401, -179.6 -88.518152273978, -179.89654756364214 -88.83333333333334))"
  g <- wk::as_wkb(wk::wkt(wkt))
  r <- burn(g, extent = c(-180, -175, -90, -80.66666666666667), dimension = c(15L, 28L), mode = "approx")
  expect_true(all(r$runs$col_end <= 2L))
  expect_equal(nrow(r$runs), 2L)
})

test_that("a cell with edges only along its walls stays a boundary cell", {
  ## From a raster-derived waterbody (DEA Waterbodies uid q9g5b3eyr_v1): a
  ## lattice-aligned notch. Cell [4, 5] x [2, 3] has polygon edges along
  ## both walls, one crossing the row centre up and the other down, so it
  ## has zero coverage and zero net winding. It used to get no boundary
  ## record and was filled as part of the interior run: 3.5 cells, not 2.5.
  wkt <- "POLYGON ((4 2, 2 7, 6 2, 5 2, 5 3, 3 5.5, 4 2.5, 4 2))"
  r <- burn(wk::as_wkb(wk::wkt(wkt)), extent = c(0, 10, 0, 10), dimension = c(10L, 10L))
  a <- sum(r$runs$col_end - r$runs$col_start + 1) + sum(r$edges$fraction)
  expect_equal(a, 2.5, tolerance = 1e-6)
  ## the notch cell is row 8, col 5 (1-based): outside
  expect_false(any(r$runs$row == 8 & r$runs$col_start <= 5 & r$runs$col_end >= 5))
  expect_false(any(r$edges$row == 8 & r$edges$col == 5))
})
