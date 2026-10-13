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
