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

## ---- corner passes and near-coincident edges (#25, controlledburn-rs#6) ----

## value of one 1-based cell from the tables
cell_value <- function(r, row, col) {
  sum(r$runs$row == row & r$runs$col_start <= col & r$runs$col_end >= col) +
    sum(r$edges$fraction[r$edges$row == row & r$edges$col == col])
}

## WKT of every start vertex, both orientations, of a closed ring
all_starts <- function(x, y) {
  n <- length(x) - 1L
  out <- character()
  for (rev in c(FALSE, TRUE)) for (k in seq_len(n)) {
    i <- ((seq_len(n) + k - 2L) %% n) + 1L
    if (rev) i <- rev(i)
    i <- c(i, i[1L])
    out <- c(out, sprintf("POLYGON ((%s))", paste(sprintf("%.17g %.17g", x[i], y[i]), collapse = ", ")))
  }
  out
}

test_that("an edge exactly through a cell corner is not overcounted", {
  ## (12.6, 1.4) -> (13.65, 0.35) passes through the top-right corner of
  ## 0-based cell (12, 18). It burned 1.0; exact 0.0943.
  x <- c(13.49307281327689, 13.774064548591133, 14.150447247739818, 14.454947247739819,
         14.571255898165475, 14.454947247739819, 17.207761779641103, 12.6,
         13.649999999999999, 6.0614012054372814, 13.49307281327689)
  y <- c(0.280951025273552, 0.1953615710514231, 0.19536157105142304, 0.2945346261241548,
         0.455, 0.6154653738758452, 4.932666216410308, 1.4000000000000004,
         0.34999999999999964, 0.45500000000000046, 0.280951025273552)
  for (w in all_starts(x, y)) {
    r <- burn(wk::as_wkb(wk::wkt(w)), extent = c(0, 20.299999999999997, 0, 9.1), dimension = c(29L, 13L))
    expect_equal(cell_value(r, 13, 19), 0.09430835509369358, tolerance = 1e-5)
    expect_equal((sum(r$runs$col_end - r$runs$col_start + 1) + sum(r$edges$fraction)) * 0.49,
                 6.157834037584929, tolerance = 1e-6)
  }
})

test_that("a 4e-13 sliver with near-coincident edges covers nothing", {
  x <- c(1008.7666666666667, 1003.7666666666667, 1002.1, 1008.7666666666667)
  y <- c(-2000.7, -1999.7, -1999.3666666666666, -2000.7)
  for (k in 1:4) for (w in all_starts(x, y)) {
    r <- burn(wk::as_wkb(wk::wkt(w)), extent = c(1000.1, 1009.7666666666667, -2000.7, -1998.0333333333333),
              dimension = c(29L * k, 8L * k))
    expect_lt(sum(r$runs$col_end - r$runs$col_start + 1) + sum(r$edges$fraction), 1e-5)
  }
})

test_that("a spike to a corner, and a wall traced both ways, stay exact", {
  ## zero-width spike to the corner of 0-based cell (25, 4): fully covered
  x <- c(0.09000464930083663, 0.225, 0.2, 0.275, 0.25, 0.2279343842256219, 0.22617971754901423,
         0.22173674555812387, 0.2166843842256219, 0.09000464930083663)
  y <- c(0.2951258728043881, 0.2749999999999999, 0.3500000000000001, 0.125, 0.44999999999999996,
         0.4548286815756889, 0.4856824868406428, 0.5020994537202749, 0.49639790095734193, 0.2951258728043881)
  for (w in all_starts(x, y)) {
    r <- burn(wk::as_wkb(wk::wkt(w)), extent = c(0, 0.25, 0, 1.6), dimension = c(5L, 32L))
    expect_equal(cell_value(r, 26, 5), 1)
  }
  ## a triangle below the grid traces its bottom edge out and back: nothing
  x <- c(2.837265806357478, 3.014285714285714, 3.0402714541578892, 2.9374143113007465, 2.837265806357478)
  y <- c(0.30000000000000004, 0.30000000000000004, 0.22857142857142854, 0.3, 0.30000000000000004)
  for (w in all_starts(x, y)) {
    r <- burn(wk::as_wkb(wk::wkt(w)), extent = c(0.3, 3.728571428571428, 0.3, 1.0142857142857142),
              dimension = c(24L, 5L))
    expect_equal(nrow(r$runs) + nrow(r$edges), 0L)
  }
})

test_that("an approx centre tie is decided exactly", {
  ## two edges along one line both pass 1e-13 right of the centre of 0-based
  ## cell (13, 2); a rounded intercept landed on it and the cell burned.
  ## Hex literals: R parses 17-digit decimal literals with its own strtod,
  ## which is an ulp off on platforms without an extended long double
  ## (macOS arm64), and an ulp here moves the edges across the centre.
  ## Decimal: x 1000.4 1001.3000000000001 1001.293897355882 1001.3511730051581
  ## 1002.0812267040778 1002.398382455067 1001.293897355882 1001.45 1000.4;
  ## y -2000.7 -1999.8000000000002 -1999.95 -1999.8000000000002
  ## -1999.6101141954598 -1997.7693056968183 -1999.4751141954598 -1999.65
  ## -2000.7; extent 1000.1 1002.8000000000001 -2000.7 -1996.2.
  x <- c(0x1.f433333333333p+9, 0x1.f4a6666666667p+9, 0x1.f4a59e6db5f27p+9, 0x1.f4acf33cae323p+9,
         0x1.f50a65a2fac9cp+9, 0x1.f532fe323fe81p+9, 0x1.f4a59e6db5f27p+9, 0x1.f4b999999999ap+9,
         0x1.f433333333333p+9)
  y <- c(-0x1.f42cccccccccdp+10, -0x1.f3f3333333334p+10, -0x1.f3fcccccccccdp+10, -0x1.f3f3333333334p+10,
         -0x1.f3e70c1c6914dp+10, -0x1.f3713c4df61d8p+10, -0x1.f3de68455ed76p+10, -0x1.f3e999999999ap+10,
         -0x1.f42cccccccccdp+10)
  ext <- c(0x1.f40cccccccccdp+9, 0x1.f566666666667p+9, -0x1.f42cccccccccdp+10, -0x1.f30cccccccccdp+10)
  for (w in all_starts(x, y)) {
    r <- burn(wk::as_wkb(wk::wkt(w)), extent = ext, dimension = c(9L, 15L), mode = "approx")
    expect_equal(cell_value(r, 14, 3), 0)
  }
})
