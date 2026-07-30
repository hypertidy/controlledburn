# Tests for point-geometry input to burn.
#
# Points are the 0-D member of the unified geometry rasterization family.
# A point is either in a cell or it isn't — there is no fractional
# intersection — so points carry no `weight` column. The output schema
# is a separate $points data.frame with columns (row, col, id).
# Materialise treats absent-weight as implicit weight = 1.

test_that("single point lands in expected cell", {
  skip_if_not_installed("geos")
  pt <- geos::as_geos_geometry("POINT (1.5 1.5)")
  r <- burn(pt, extent = c(0, 3, 0, 3), dimension = c(3L, 3L))

  expect_s3_class(r, "controlledburn")
  expect_equal(nrow(r$runs), 0)
  expect_equal(nrow(r$edges), 0)
  expect_equal(nrow(r$points), 1)

  # Point at (1.5, 1.5) on a 3x3 grid with extent (0,0,3,3) is in the
  # middle cell. With top-down 1-based row indexing: y=1.5 in a 3-row
  # grid (dy=1) lands in row 2 (the middle); x=1.5 in a 3-col grid is
  # col 2 (the middle).
  expect_equal(r$points$row, 2L)
  expect_equal(r$points$col, 2L)
})

test_that("multiple points produce one record each", {
  skip_if_not_installed("geos")
  pts <- geos::as_geos_geometry(c(
    "POINT (0.5 0.5)",   # bottom-left cell (row 3, col 1)
    "POINT (1.5 1.5)",   # middle cell      (row 2, col 2)
    "POINT (2.5 2.5)"    # top-right cell   (row 1, col 3)
  ))
  r <- burn(pts, extent = c(0, 3, 0, 3), dimension = c(3L, 3L))

  expect_equal(nrow(r$points), 3)
  # Three distinct geometry ids
  expect_equal(sort(r$points$id), 1:3)
  # Diagonal placement — sorted by row gives 1, 2, 3 with cols 3, 2, 1
  expect_setequal(paste(r$points$row, r$points$col),
                  c("3 1", "2 2", "1 3"))
})

test_that("MULTIPOINT shares one id across components", {
  skip_if_not_installed("geos")
  mp <- geos::as_geos_geometry("MULTIPOINT ((0.5 0.5), (2.5 2.5))")
  r <- burn(mp, extent = c(0, 3, 0, 3), dimension = c(3L, 3L))

  expect_equal(nrow(r$points), 2)
  # Both records share the same MULTIPOINT's id (= 1)
  expect_equal(unique(r$points$id), 1L)
})

test_that("points outside the grid extent are dropped silently", {
  skip_if_not_installed("geos")
  pts <- geos::as_geos_geometry(c(
    "POINT (1.5 1.5)",       # in
    "POINT (-1 -1)",         # out (lower-left)
    "POINT (10 10)",         # out (upper-right)
    "POINT (0.5 5)"          # out (above)
  ))
  expect_silent(
    r <- burn(pts, extent = c(0, 3, 0, 3), dimension = c(3L, 3L))
  )
  expect_equal(nrow(r$points), 1L)
  expect_equal(r$points$row, 2L)
  expect_equal(r$points$col, 2L)
  # Note: out-of-extent geometries still consume an id slot (geom 1, 3, 4
  # produced no records). Currently the surviving record is geom 1's id.
  expect_equal(r$points$id, 1L)
})

test_that("point exactly on cell-corner / boundary lands deterministically", {
  skip_if_not_installed("geos")
  # Point at the shared corner of four cells. The Grid::get_row/get_column
  # tie-breaking is deterministic; this test pins which cell wins so the
  # convention is stable across versions.
  #
  # Derivation: get_row(y=1) on a bounded grid with ymax=3, dy=1 returns
  # floor((3 - 1) / 1) = 2 (0-based) → row 3 (1-based). get_column(x=1)
  # returns floor((1 - 0) / 1) = 1 (0-based) → col 2 (1-based).
  # Pinning: a point on a horizontal cell boundary belongs to the cell
  # BELOW it (closer to ymin); a point on a vertical cell boundary
  # belongs to the cell to the RIGHT (closer to xmax).
  pt <- geos::as_geos_geometry("POINT (1 1)")
  r <- burn(pt, extent = c(0, 3, 0, 3), dimension = c(3L, 3L))
  expect_equal(nrow(r$points), 1)
  # If this assertion ever needs to change, the change is an API contract
  # change, not a routine fix.
  expect_equal(r$points$row, 3L)
  expect_equal(r$points$col, 2L)
})

test_that("point at exactly xmin / xmax / ymin / ymax — extent-edge inclusion", {
  skip_if_not_installed("geos")
  # All four boundary-edge points are INCLUDED in the grid (not dropped).
  # The cell each one lands in is determined by Grid::get_row / get_column's
  # special-case handling for boundary coords.
  pts <- geos::as_geos_geometry(c(
    "POINT (0 0)",   # xmin, ymin — bottom-left corner
    "POINT (3 3)",   # xmax, ymax — top-right corner
    "POINT (0 3)",   # xmin, ymax — top-left corner
    "POINT (3 0)"    # xmax, ymin — bottom-right corner
  ))
  r <- burn(pts, extent = c(0, 3, 0, 3), dimension = c(3L, 3L))
  expect_equal(nrow(r$points), 4)
})

test_that("multiple points in same cell accumulate in materialise (count)", {
  skip_if_not_installed("geos")
  pts <- geos::as_geos_geometry(c(
    "POINT (1.2 1.2)",
    "POINT (1.5 1.5)",
    "POINT (1.8 1.8)"
  ))
  r <- burn(pts, extent = c(0, 3, 0, 3), dimension = c(3L, 3L))
  expect_equal(nrow(r$points), 3)

  m <- materialise_chunk(r)
  # All three points land in the middle cell — cell shows count = 3.
  expect_equal(m[2, 2], 3)
  # All other cells are empty.
  expect_equal(sum(m), 3)
})

test_that("materialise_chunk fills point cell with weight 1", {
  skip_if_not_installed("geos")
  pt <- geos::as_geos_geometry("POINT (1.5 1.5)")
  r <- burn(pt, extent = c(0, 3, 0, 3), dimension = c(3L, 3L))
  m <- materialise_chunk(r)
  expect_equal(dim(m), c(3, 3))
  expect_equal(m[2, 2], 1)
  expect_equal(sum(m), 1)
})

test_that("materialise_chunk id filter applies to points", {
  skip_if_not_installed("geos")
  pts <- geos::as_geos_geometry(c(
    "POINT (0.5 0.5)",   # id 1, cell (3, 1)
    "POINT (2.5 2.5)"    # id 2, cell (1, 3)
  ))
  r <- burn(pts, extent = c(0, 3, 0, 3), dimension = c(3L, 3L))

  m1 <- materialise_chunk(r, id = 1L)
  expect_equal(m1[3, 1], 1)
  expect_equal(sum(m1), 1)

  m2 <- materialise_chunk(r, id = 2L)
  expect_equal(m2[1, 3], 1)
  expect_equal(sum(m2), 1)
})

test_that("MULTIPOINT counts each component", {
  skip_if_not_installed("geos")
  mp <- geos::as_geos_geometry(
    "MULTIPOINT ((1.5 1.5), (1.5 1.5), (2.5 2.5))"
  )
  r <- burn(mp, extent = c(0, 3, 0, 3), dimension = c(3L, 3L))
  m <- materialise_chunk(r)
  # Two points at (1.5, 1.5) → cell (2, 2) gets 2; one at (2.5, 2.5) → (1, 3) gets 1.
  expect_equal(m[2, 2], 2)
  expect_equal(m[1, 3], 1)
  expect_equal(sum(m), 3)
})

test_that("point output has the expected schema (no weight column)", {
  skip_if_not_installed("geos")
  pt <- geos::as_geos_geometry("POINT (1.5 1.5)")
  r <- burn(pt, extent = c(0, 3, 0, 3), dimension = c(3L, 3L))
  expect_named(r$points, c("row", "col", "id"))
  expect_false("weight" %in% names(r$points))
})

test_that("empty input — no points table issues", {
  skip_if_not_installed("geos")
  pt <- geos::as_geos_geometry("POINT EMPTY")
  r <- burn(pt, extent = c(0, 3, 0, 3), dimension = c(3L, 3L))
  expect_equal(nrow(r$points), 0)
  expect_equal(nrow(r$edges), 0)
  expect_equal(nrow(r$runs), 0)
})

test_that("mixed polygon + line + point requires separate burns", {
  # Same documented workflow as the line tests: split mixed-dimension
  # input into homogeneous groups, run separate burns. Each output has
  # the same shape but different weight semantics.
  skip_if_not_installed("geos")
  poly <- geos::as_geos_geometry(
    "POLYGON ((0.5 0.5, 1.5 0.5, 1.5 1.5, 0.5 1.5, 0.5 0.5))"
  )
  line <- geos::as_geos_geometry("LINESTRING (1.5 1.5, 2.5 2.5)")
  pt   <- geos::as_geos_geometry("POINT (2.5 0.5)")

  ext <- c(0, 3, 0, 3); dim <- c(3L, 3L)
  rp <- burn(poly, extent = ext, dimension = dim)
  rl <- burn(line, extent = ext, dimension = dim)
  rt <- burn(pt,   extent = ext, dimension = dim)

  expect_true(nrow(rp$edges) + nrow(rp$runs) > 0)
  expect_true(nrow(rl$lines) > 0)
  expect_equal(nrow(rt$points), 1)
})

test_that("burn_sparse() errors on point input rather than silently returning empty", {
  # Previous behaviour: raster_cell_intersection threw 'Unsupported geometry
  # type.' which cpp_burn_sparse's try/catch demoted to a per-geometry
  # warning, leaving the user with a controlledburn object holding zero
  # records and no clear signal that their input was wrong. The new
  # behaviour errors up front with a message that names burn()
  # as the working alternative.
  skip_if_not_installed("geos")
  pt <- geos::as_geos_geometry("POINT (1.5 1.5)")
  expect_error(
    burn_sparse(pt, extent = c(0, 3, 0, 3), dimension = c(3L, 3L)),
    "burn"
  )

  mp <- geos::as_geos_geometry("MULTIPOINT ((0.5 0.5), (2.5 2.5))")
  expect_error(
    burn_sparse(mp, extent = c(0, 3, 0, 3), dimension = c(3L, 3L)),
    "burn"
  )

  # Mixed input: the error fires on the first point encountered.
  poly <- geos::as_geos_geometry(
    "POLYGON ((0.5 0.5, 1.5 0.5, 1.5 1.5, 0.5 1.5, 0.5 0.5))"
  )
  expect_error(
    burn_sparse(c(poly, pt), extent = c(0, 3, 0, 3), dimension = c(3L, 3L)),
    "burn"
  )
})
