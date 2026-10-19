two_polys <- function(dimension = c(40L, 40L)) {
  skip_if_not_installed("geos")
  poly <- geos::as_geos_geometry(c(
    "POLYGON ((1 1, 9 1, 9 9, 1 9, 1 1), (3 3, 3 6, 6 6, 6 3, 3 3))",
    "POLYGON ((2.1 7.2, 4.3 8.1, 3.3 9.6, 2.1 7.2))"
  ))
  burn(poly, extent = c(0, 10, 0, 10), dimension = dimension)
}

test_that("summary marginals and areas match the dense materialisation", {
  r <- two_polys()
  s <- summary(r)
  m <- materialise_chunk(r)
  expect_equal(s$row_marginal, rowSums(m))
  expect_equal(s$col_marginal, colSums(m))
  expect_equal(s$n_geom, 2L)
  expect_equal(sum(s$area$interior + s$area$boundary), sum(m))
  ## cell area is 0.25 * 0.25: polygon 1 is 64 - 9, polygon 2 its triangle
  cell <- 0.25^2
  expect_equal((s$area$interior + s$area$boundary)[1] * cell, 55, tolerance = 1e-9)
  expect_equal(s$records[["runs"]], nrow(r$runs))
  expect_output(print(s), "controlledburn summary")
  expect_output(print(s, ascii = TRUE), "overview")
})

test_that("content extent is the bounding box of touched cells", {
  r <- two_polys()
  ce <- summary(r)$content
  expect_equal(ce$extent, c(1, 9, 1, 9.75))
  expect_equal(ce$cols, c(5L, 36L))
})

test_that("overview with block 1 is the exact coverage", {
  r <- two_polys()
  m <- materialise_chunk(r)
  m[m == 0] <- NA
  ov <- overview_burn(r, block = 1L)
  expect_equal(unname(ov[, ]), m)
  expect_equal(attr(ov, "extent"), r$extent)
})

test_that("overview tile means match block averages of the dense grid", {
  r <- two_polys(c(37L, 23L))  # blocks do not divide the grid
  m <- materialise_chunk(r)
  b <- 4L
  ov <- overview_burn(r, block = b)
  ri <- (seq_len(nrow(m)) - 1L) %/% b
  ci <- (seq_len(ncol(m)) - 1L) %/% b
  dense <- t(rowsum(t(rowsum(m, ri)), ci)) / outer(tabulate(ri + 1L), tabulate(ci + 1L))
  dense[dense == 0] <- NA
  expect_equal(unname(ov[, ]), unname(dense))
  ext <- attr(ov, "extent")
  expect_equal(ext[1], 0)
  expect_equal(ext[4], 10)
  expect_gte(ext[2], 10)
  expect_lte(ext[3], 0)
})

test_that("summary works for lines, points and an empty burn", {
  skip_if_not_installed("geos")
  g <- geos::as_geos_geometry(c("LINESTRING (0.5 0.5, 9.5 9.5)", "POINT (5 2)"))
  r <- burn(g, extent = c(0, 10, 0, 10), dimension = c(10L, 10L))
  s <- summary(r)
  expect_equal(sum(s$area$lines), sqrt(2) * 9, tolerance = 1e-6)
  expect_equal(sum(s$area$points), 1)
  e <- burn(geos::as_geos_geometry("POLYGON ((20 20, 21 20, 21 21, 20 20))"),
            extent = c(0, 10, 0, 10), dimension = c(10L, 10L))
  se <- summary(e)
  expect_null(se$content)
  expect_true(all(is.na(se$overview)))
  expect_output(print(se), "empty")
})

test_that("plot draws exact and overview without error", {
  r <- two_polys()
  f <- tempfile(fileext = ".png")
  grDevices::png(f, width = 600, height = 500)
  on.exit({grDevices::dev.off(); unlink(f)})
  expect_invisible(plot(r))
  expect_s3_class(attr(plot(r, by_id = TRUE, content = TRUE, main = "id"), "summary"),
                  "summary.controlledburn")
  expect_silent(plot(r, what = "overview", tiles = 8, marginals = FALSE, main = "ov"))
  expect_silent(plot(crop_burn(r, c(2, 5, 3, 7))))
})
