## A cell crossed only by a hole's boundary, inside the exterior, keeps the
## exterior's cover (1 minus the hole's fraction). It used to get the hole's
## negative fraction alone and be dropped, so holes removed too much area.
## No geos needed.

covered <- function(r) {
  cell <- prod((r$extent[c(2, 4)] - r$extent[c(1, 3)]) / r$dimension)
  (sum(r$runs$col_end - r$runs$col_start + 1) + sum(r$edges$fraction)) * cell
}

shoelace <- function(x, y) {
  n <- length(x)
  abs(sum(x[-n] * y[-1] - x[-1] * y[-n])) / 2
}

ring_wkt <- function(x, y) paste0("(", paste(x, y, collapse = ", "), ")")

test_that("a hole inside the exterior removes exactly its area", {
  ex <- c(0.5, 9.5, 9.5, 0.5, 0.5); ey <- c(0.5, 0.5, 9.5, 9.5, 0.5)
  holes <- list(
    list(x = c(4.2, 5.8, 5.8, 4.2, 4.2), y = c(4.2, 4.2, 5.8, 5.8, 4.2)), # across four cells
    list(x = c(4.2, 4.8, 4.8, 4.2, 4.2), y = c(4.2, 4.2, 4.8, 4.8, 4.2)), # within one cell
    list(x = c(4.2, 5.8, 5.0, 4.2),      y = c(4.6, 4.6, 5.4, 4.6)),      # off the row centres
    list(x = c(4.2, 4.2, 5.8, 5.8, 4.2), y = c(4.2, 5.8, 5.8, 4.2, 4.2))  # the other winding
  )
  for (h in holes) {
    wkt <- sprintf("POLYGON (%s, %s)", ring_wkt(ex, ey), ring_wkt(h$x, h$y))
    r <- burn(wk::as_wkb(wk::wkt(wkt)), extent = c(0, 10, 0, 10), dimension = c(10L, 10L))
    expect_equal(covered(r), shoelace(ex, ey) - shoelace(h$x, h$y), tolerance = 1e-6)
    ## cell (row 6, col 5) holds part of every hole: partly covered
    e <- r$edges[r$edges$row == 6 & r$edges$col == 5, ]
    expect_equal(nrow(e), 1L)
    expect_true(e$fraction > 0 && e$fraction < 1)
  }
})

test_that("two holes, and many holes, cover their area", {
  ex <- c(0, 10, 10, 0, 0); ey <- c(0, 0, 10, 10, 0)
  h1x <- c(1.3, 3.7, 3.7, 1.3, 1.3); h1y <- c(1.3, 1.3, 3.7, 3.7, 1.3)
  h2x <- c(6.1, 8.9, 7.5, 6.1);      h2y <- c(6.6, 6.6, 8.8, 6.6)
  wkt <- sprintf("POLYGON (%s, %s, %s)", ring_wkt(ex, ey), ring_wkt(h1x, h1y), ring_wkt(h2x, h2y))
  r <- burn(wk::as_wkb(wk::wkt(wkt)), extent = c(0, 10, 0, 10), dimension = c(10L, 10L))
  expect_equal(covered(r), 100 - shoelace(h1x, h1y) - shoelace(h2x, h2y), tolerance = 1e-6)

  ## 20 x 20 diamond holes, none aligned with the grid
  rings <- ring_wkt(c(0, 100, 100, 0, 0), c(0, 0, 100, 100, 0))
  want <- 10000
  for (i in 0:19) for (j in 0:19) {
    cx <- 2.5 + 5 * i + 0.013; cy <- 2.5 + 5 * j + 0.029; rr <- 1.7
    rings <- c(rings, ring_wkt(c(cx - rr, cx, cx + rr, cx, cx - rr), c(cy, cy - rr, cy, cy + rr, cy)))
    want <- want - 2 * rr^2
  }
  g <- wk::as_wkb(wk::wkt(sprintf("POLYGON (%s)", paste(rings, collapse = ", "))))
  r <- burn(g, extent = c(0, 100, 0, 100), dimension = c(333L, 333L))
  expect_equal(covered(r), want, tolerance = 1e-6)
  ## materialise_chunk agrees
  expect_equal(sum(materialise_chunk(r)) * (100 / 333)^2, want, tolerance = 1e-6)
})
