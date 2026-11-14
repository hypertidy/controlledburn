# Features aligned with the controlledburn Rust crate (0.3.1): burn options
# (merge_parts, coalesce_runs, id_offset, threads), burning onto a window
# (target), slice_burn() and relate_cells(). The core's records are checked
# against the crate itself in cpp/tests/test_canon_ops.cpp; these tests
# check the R interface, and relate_cells() against GEOS.

skip_if_not_installed("geos")

wkt <- function(...) geos::as_geos_geometry(c(...))

cells_of <- function(r) {
  m <- materialise_chunk(r)
  m
}

# A ring of k vertices around (cx, cy), radius wobbled.
ring_wkt <- function(cx, cy, r, k, phase) {
  a <- seq(0, 2 * pi, length.out = k + 1)[-(k + 1)]
  rr <- r * (0.6 + 0.4 * abs(sin(phase + 1.7 * seq_len(k))))
  x <- cx + rr * cos(a)
  y <- cy + rr * sin(a)
  sprintf("POLYGON ((%s))", paste(sprintf("%.6f %.6f", c(x, x[1]), c(y, y[1])), collapse = ", "))
}

many_polygons <- function(n) {
  geos::as_geos_geometry(vapply(seq_len(n), function(i) {
    ring_wkt((i * 37.31) %% 370, (i * 17.77) %% 230, 3, 40, i)
  }, ""))
}

test_that("merge_parts: two halves of a cell either side of a border make a full cell", {
  mp <- wkt("MULTIPOLYGON (((0 0, 1.5 0, 1.5 3, 0 3, 0 0)), ((1.5 0, 3 0, 3 3, 1.5 3, 1.5 0)))")
  r <- burn(mp, extent = c(0, 3, 0, 3), dimension = c(3L, 3L))
  expect_equal(nrow(r$edges), 0L)
  expect_equal(sum(r$runs$col_end - r$runs$col_start + 1L), 9)
  expect_equal(nrow(r$runs), 3L)  # one coalesced run per row

  old <- burn(mp, extent = c(0, 3, 0, 3), dimension = c(3L, 3L),
              merge_parts = FALSE, coalesce_runs = FALSE)
  expect_equal(nrow(old$edges), 6L)  # 0.5 from each part in the middle column
  expect_equal(materialise_chunk(old), materialise_chunk(r))
})

test_that("coalesce_runs changes the runs, not the cells", {
  # sides on column walls: their boundary cells are full, one-cell runs
  p <- wkt("POLYGON ((1 0.5, 9 0.5, 9 9.5, 1 9.5, 1 0.5))")
  a <- burn(p, extent = c(0, 10, 0, 10), dimension = c(10L, 10L))
  b <- burn(p, extent = c(0, 10, 0, 10), dimension = c(10L, 10L), coalesce_runs = FALSE)
  expect_gt(nrow(b$runs), nrow(a$runs))
  expect_equal(materialise_chunk(a), materialise_chunk(b))
  expect_identical(a$edges, b$edges)
})

test_that("id_offset shifts ids so burns of slices bind together", {
  g <- many_polygons(20)
  ext <- c(0, 370, 0, 230)
  whole <- burn(g, extent = ext, dimension = c(370L, 230L))
  a <- burn(g[1:12], extent = ext, dimension = c(370L, 230L))
  b <- burn(g[13:20], extent = ext, dimension = c(370L, 230L), id_offset = 12L)
  expect_identical(rbind(a$runs, b$runs), whole$runs)
  expect_identical(rbind(a$edges, b$edges), whole$edges)
  expect_equal(range(b$runs$id), c(13L, 20L))
  expect_identical(b$options$id_offset, 12L)
  expect_error(burn(g, id_offset = -1), "non-negative")
})

test_that("threads give the same result as one thread", {
  g <- many_polygons(3000)
  ext <- c(0, 370, 0, 230)
  for (mode in c("coverage", "approx")) {
    one <- burn(g, extent = ext, dimension = c(1530L, 910L), mode = mode)
    for (t in c(2L, 0L)) {
      expect_identical(burn(g, extent = ext, dimension = c(1530L, 910L), mode = mode,
                            threads = t), one)
    }
  }
  expect_error(burn(g, threads = -1), "non-negative")
})

test_that("target burns the window: crop_burn() of the whole burn, record for record", {
  g <- c(
    wkt("MULTIPOLYGON (((1 1, 5 1.2, 4.7 6, 1 5.5, 1 1), (2 2, 3 2, 3 3, 2 3, 2 2)), ((5 5, 9 5, 9 9, 5 9, 5 5)))"),
    wkt("LINESTRING (0.2 9.7, 3 3, 9.9 0.1)", "MULTIPOINT ((1.5 1.5), (6 6), (5 5))")
  )
  ext <- c(0, 10, 0, 10)
  targets <- list(c(2.2, 6.1, 3.3, 7.9), c(0, 10, 4.5, 5.5), c(5, 10, 0, 5), c(9.5, 10, 0, 0.5))
  for (mode in c("coverage", "approx")) {
    whole <- burn(g, extent = ext, dimension = c(37L, 23L), mode = mode)
    for (tg in targets) {
      w <- burn(g, extent = ext, dimension = c(37L, 23L), mode = mode, target = tg)
      cr <- crop_burn(whole, tg)
      expect_identical(w, cr)
    }
  }
  expect_error(burn(g, extent = ext, dimension = c(37L, 23L), target = c(20, 30, 20, 30)),
               "does not overlap")
})

test_that("slice_burn tiles are crops of the whole burn, and reassemble it", {
  g <- c(many_polygons(40), wkt("LINESTRING (0 0, 370 230)"))
  ext <- c(0, 370, 0, 230)
  r <- burn(g, extent = ext, dimension = c(370L, 230L))
  tiles <- slice_burn(r, c(64L, 50L))
  expect_true(length(tiles) > 1)
  full <- materialise_chunk(r)
  back <- matrix(0, nrow(full), ncol(full))
  for (t in tiles) {
    cr <- crop_burn(r, t$extent)
    expect_identical(t$runs, cr$runs)
    expect_identical(t$edges, cr$edges)
    expect_identical(t$lines, cr$lines)
    expect_identical(t$dimension, cr$dimension)
    m <- materialise_chunk(t)
    rows <- t$offset["row"] + seq_len(t$dimension[2])
    cols <- t$offset["col"] + seq_len(t$dimension[1])
    back[rows, cols] <- m
    expect_equal(t$tile, c(row = t$offset[["row"]] %/% 50L + 1L,
                           col = t$offset[["col"]] %/% 64L + 1L))
  }
  expect_equal(back, full)
  expect_identical(slice_burn(r, c(64L, 50L), threads = 3L), tiles)
  expect_identical(length(slice_burn(r, 1000L)), 1L)
  expect_error(slice_burn(r, 0L), "positive")
})

# The DE-9IM matrix of a polygon against each cell, from GEOS.
geos_relate_cells <- function(poly, extent, dimension, rows, cols) {
  dx <- (extent[2] - extent[1]) / dimension[1]
  dy <- (extent[4] - extent[3]) / dimension[2]
  x0 <- extent[1] + (cols - 1) * dx
  x1 <- ifelse(cols == dimension[1], extent[2], extent[1] + cols * dx)
  y1 <- extent[4] - (rows - 1) * dy
  y0 <- ifelse(rows == dimension[2], extent[3], extent[4] - rows * dy)
  boxes <- geos::geos_make_polygon(
    c(rbind(x0, x1, x1, x0)), c(rbind(y0, y0, y1, y1)),
    feature_id = rep(seq_along(rows), each = 4)
  )
  geos::geos_relate(poly, boxes)
}

test_that("relate_cells: a square on cell walls", {
  sq <- wkt("POLYGON ((1 1, 3 1, 3 3, 1 3, 1 1))")
  r <- burn(sq, extent = c(0, 4, 0, 4), dimension = c(4L, 4L))
  rel <- relate_cells(sq, r)
  expect_equal(sum(rel$covers), 4L)
  expect_equal(sum(rel$touches), 12L)
  expect_true(all(rel$intersects))
  expect_named(rel, c("row", "col", "id", "matrix", "intersects", "interiors_meet",
                      "covers", "touches", "covered_by"))
})

test_that("relate_cells agrees with GEOS on every cell, and the burn gives the rest", {
  polys <- c(
    wkt("POLYGON ((0.5 0.5, 6.2 0.7, 7 6, 3 7.5, 0.4 5, 0.5 0.5), (2 2, 4 2, 4 4, 2 4, 2 2))"),
    wkt("POLYGON ((1 1, 5 1, 5 3, 3 3, 3 5, 1 5, 1 1))"),          # on walls and corners
    wkt("POLYGON ((4 0, 8 4, 4 8, 0 4, 4 0))"),                    # through corners
    wkt("POLYGON ((2.2 2.2, 2.8 2.3, 2.5 2.9, 2.2 2.2))"),         # inside one cell
    wkt("MULTIPOLYGON (((0 0, 2 0, 2 2, 0 2, 0 0)), ((2 2, 4 2, 4 4, 2 4, 2 2)))"),
    wkt("POLYGON ((-1 -1, 9 -1, 9 9, -1 9, -1 -1), (3 3, 4 3, 4 4, 3 4, 3 3))")  # hole is a cell
  )
  ext <- c(0, 8, 0, 8)
  dim <- c(8L, 8L)
  for (k in seq_along(polys)) {
    p <- polys[k]
    for (mode in c("coverage", "approx")) {
      r <- burn(p, extent = ext, dimension = dim, mode = mode)
      rel <- relate_cells(p, r)
      grid <- expand.grid(row = seq_len(dim[2]), col = seq_len(dim[1]))
      want <- geos_relate_cells(p, ext, dim, grid$row, grid$col)
      runs <- r$runs
      in_run <- vapply(seq_len(nrow(grid)), function(i) {
        any(runs$row == grid$row[i] & runs$col_start <= grid$col[i] & runs$col_end >= grid$col[i])
      }, NA)
      got <- ifelse(in_run, "212FF1FF2", "FF2FF1212")
      key <- paste(grid$row, grid$col)
      got[match(paste(rel$row, rel$col), key)] <- rel$matrix
      expect_identical(got, want, info = paste("polygon", k, mode))
    }
  }
})

test_that("relate_cells uses the burn's id_offset and skips lines and points", {
  g <- c(wkt("POLYGON ((1 1, 3 1, 3 3, 1 3, 1 1))"), wkt("LINESTRING (0 0, 4 4)", "POINT (2 2)"))
  r <- burn(g, extent = c(0, 4, 0, 4), dimension = c(4L, 4L), id_offset = 10L)
  rel <- relate_cells(g, r)
  expect_true(all(rel$id == 11L))
  expect_equal(nrow(rel), 16L)
})
