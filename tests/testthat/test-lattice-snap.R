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

test_that("a zero-width spike into a covered cell keeps it covered", {
  ## Found by the controlledburn-rs property suite: the ring runs in through
  ## the top-left corner of cell [4.8, 5.1] x [1.2, 1.5] to its centre and
  ## straight back out. That traversal leaves where it entered and encloses
  ## nothing; the cell is inside the polygon and must stay covered.
  wkt <- paste0(
    "POLYGON ((2.2955346720684804 -0.28831957426678434, 4.2871984006926445 -0.6968703145658899, ",
    "6.183801073287196 -1.4711665488619694, 5.24135385306532 1.4474746273177665, ",
    "5.275730794077829 1.5426958381891471, 5.824891418164306 2.019485421020327, ",
    "5.151353853065319 1.696766993828962, 3.5999999999999996 2.6999999999999997, ",
    "4.95 1.3499999999999996, 4.8 1.4999999999999996, 2.2955346720684804 -0.28831957426678434))"
  )
  r <- burn(wk::as_wkb(wk::wkt(wkt)), extent = c(0, 6, 0, 5.3999999999999995), dimension = c(20L, 18L))
  a <- (sum(r$runs$col_end - r$runs$col_start + 1) + sum(r$edges$fraction)) * 0.09
  expect_equal(a, 3.026275470389206, tolerance = 1e-6)
  ## the spike cell is row 14, col 17 (1-based): fully covered
  expect_true(any(r$runs$row == 14 & r$runs$col_start <= 17 & r$runs$col_end >= 17))
})

test_that("an edge clipping a cell corner does not add the whole cell", {
  ## The corner clip is recorded as a traversal with entry and exit an ulp
  ## or two apart; compared exactly it counted as movement, and the cell
  ## gained a whole cell of coverage (controlledburn#25,
  ## controlledburn-rs#6). Exact values are the polygon's intersection with
  ## the cell (GEOS).
  cell_value <- function(r, row, col) {
    sum(r$runs$row == row & r$runs$col_start <= col & r$runs$col_end >= col) +
      sum(r$edges$fraction[r$edges$row == row & r$edges$col == col])
  }
  ## #25: cell (13, 19) came out 1 instead of 0.0943
  w25 <- paste0(
    "POLYGON ((13.49307281327689 0.280951025273552, 13.774064548591133 0.1953615710514231, ",
    "14.150447247739818 0.19536157105142304, 14.454947247739819 0.2945346261241548, ",
    "14.571255898165475 0.455, 14.454947247739819 0.6154653738758452, ",
    "17.207761779641103 4.932666216410308, 12.6 1.4000000000000004, ",
    "13.649999999999999 0.34999999999999964, 6.0614012054372814 0.45500000000000046, ",
    "13.49307281327689 0.280951025273552))"
  )
  r <- burn(wk::as_wkb(wk::wkt(w25)), extent = c(0, 20.299999999999997, 0, 9.1), dimension = c(29L, 13L))
  expect_equal(cell_value(r, 13, 19), 0.09430835509369219, tolerance = 1e-5)
  expect_equal((sum(r$runs$col_end - r$runs$col_start + 1) + sum(r$edges$fraction)) * 0.49,
               6.157834037584925, tolerance = 1e-6)

  ## controlledburn-rs#6: cell (5, 5) came out 0.739 instead of 0.455
  wo <- paste0(
    "POLYGON ((2.9404541633227357 6.235726447802992, 7.140490483291759 -0.6880522184656641, ",
    "12.79922599805937 6.649999999999999, 13.222799970968294 7.999704574445243, ",
    "12.79922599805937 8.351528741555857, 8.3 10.399999999999999, 12.05 6.649999999999999, ",
    "2.9404541633227357 6.235726447802992))"
  )
  r <- burn(wk::as_wkb(wk::wkt(wo)), extent = c(-1.7, 35.8, 2.9, 17.9), dimension = c(15L, 6L))
  expect_equal(cell_value(r, 5, 5), 0.4545231387355179, tolerance = 1e-5)

  ## controlledburn-rs#6: cell (4, 1) came out 0.839 instead of 0.339, and
  ## depended on the ring's start vertex
  x <- c(0.22857142857142856, 0.3, 0.5385327833222915, 0.3826046772369731, 0.44285714285714284,
         0.6141351860119378, 0.3, 0.37142857142857144)
  y <- c(0.657142857142857, 0.7285714285714284, 0.5669218770762636, 0.8968780211084503,
         0.9150607669103901, 1.2110801232422874, 0.9407750526246758, 0.7999999999999999)
  for (k in seq_along(x)) {
    i <- c(k:length(x), seq_len(k - 1), k)
    w <- sprintf("POLYGON ((%s))", paste(sprintf("%.17g %.17g", x[i], y[i]), collapse = ", "))
    r <- burn(wk::as_wkb(wk::wkt(w)), extent = c(0.3, 1.0142857142857142, 0.3, 1.157142857142857),
              dimension = c(5L, 6L))
    expect_equal(cell_value(r, 4, 1), 0.338841372753266, tolerance = 1e-5)
  }
})
