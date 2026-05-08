# Tests for line-geometry input to burn_scanline.
#
# Lines produce only an `edges` table (no interior `runs`), with `weight`
# carrying the absolute length of the line within each cell, in CRS units
# (not a fraction). This is the line-row of the unified geometry rasterization
# scheme; see inst/docs-design/unified-geometry-rasterization.md.

test_that("simple horizontal line crossing three cells emits length per cell", {
  skip_if_not_installed("geos")
  # Line from x=0.5 to x=2.5 at y=0.5 crosses three 1x1 cells in a 3x3 grid.
  # Cell (3, 1): x in [0, 1], length = 0.5 (from 0.5 to 1.0)
  # Cell (3, 2): x in [1, 2], length = 1.0 (full cell width)
  # Cell (3, 3): x in [2, 3], length = 0.5 (from 2.0 to 2.5)
  line <- geos::as_geos_geometry("LINESTRING (0.5 0.5, 2.5 0.5)")
  r <- burn_scanline(line, extent = c(0, 3, 0, 3), dimension = c(3, 3))
  expect_s3_class(r, "controlledburn")

  # Lines have no interior runs
  expect_equal(nrow(r$runs), 0)

  # Three edge cells with the expected lengths
  expect_equal(nrow(r$edges), 3)
  expect_equal(sort(r$edges$weight), c(0.5, 0.5, 1.0))

  # Sum of per-cell lengths equals total line length
  expect_equal(sum(r$edges$weight), 2.0)
})

test_that("vertical line crossing three cells emits length per cell", {
  skip_if_not_installed("geos")
  line <- geos::as_geos_geometry("LINESTRING (0.5 0.5, 0.5 2.5)")
  r <- burn_scanline(line, extent = c(0, 3, 0, 3), dimension = c(3, 3))
  expect_equal(nrow(r$runs), 0)
  expect_equal(nrow(r$edges), 3)
  expect_equal(sum(r$edges$weight), 2.0)
})

test_that("diagonal line: per-cell lengths sum to total length", {
  skip_if_not_installed("geos")
  # Diagonal from (0.5, 0.5) to (2.5, 2.5). Total length = sqrt(8) ≈ 2.828.
  line <- geos::as_geos_geometry("LINESTRING (0.5 0.5, 2.5 2.5)")
  r <- burn_scanline(line, extent = c(0, 3, 0, 3), dimension = c(3, 3))
  expect_equal(nrow(r$runs), 0)
  expect_true(nrow(r$edges) > 0)
  expect_equal(sum(r$edges$weight), sqrt(8), tolerance = 1e-6)
})

test_that("line with vertex inside a cell sums sub-segment lengths", {
  skip_if_not_installed("geos")
  # L-shape: enter cell (3, 2) at (1.0, 0.5), bend at vertex (1.5, 0.5),
  # then up to (1.5, 1.0). The cell's traversal records [entry, vertex, exit]
  # and the length sum is 0.5 + 0.5 = 1.0.
  line <- geos::as_geos_geometry("LINESTRING (0.5 0.5, 1.5 0.5, 1.5 2.5)")
  r <- burn_scanline(line, extent = c(0, 3, 0, 3), dimension = c(3, 3))
  expect_equal(nrow(r$runs), 0)
  expect_true(nrow(r$edges) > 0)
  # Total length: 1.0 (horizontal) + 2.0 (vertical) = 3.0
  expect_equal(sum(r$edges$weight), 3.0, tolerance = 1e-6)
})

test_that("line entirely within a single cell", {
  skip_if_not_installed("geos")
  # Both endpoints inside cell (3, 1) of a 3x3 grid on extent (0,3)x(0,3).
  line <- geos::as_geos_geometry("LINESTRING (0.2 0.2, 0.8 0.8)")
  r <- burn_scanline(line, extent = c(0, 3, 0, 3), dimension = c(3, 3))
  expect_equal(nrow(r$runs), 0)
  expect_equal(nrow(r$edges), 1)
  expect_equal(r$edges$weight, sqrt(0.72), tolerance = 1e-6) # sqrt(0.6^2 + 0.6^2)
})

test_that("MULTILINESTRING accumulates lengths across components", {
  skip_if_not_installed("geos")
  # Two parallel horizontal lines. Each contributes length 2.0, total 4.0.
  ml <- geos::as_geos_geometry(
    "MULTILINESTRING ((0.5 0.5, 2.5 0.5), (0.5 2.5, 2.5 2.5))"
  )
  r <- burn_scanline(ml, extent = c(0, 3, 0, 3), dimension = c(3, 3))
  expect_equal(nrow(r$runs), 0)
  expect_equal(sum(r$edges$weight), 4.0, tolerance = 1e-6)
})

test_that("degenerate line (single point repeated) emits no edges", {
  skip_if_not_installed("geos")
  line <- geos::as_geos_geometry("LINESTRING (0.5 0.5, 0.5 0.5)")
  r <- burn_scanline(line, extent = c(0, 3, 0, 3), dimension = c(3, 3))
  expect_equal(nrow(r$runs), 0)
  expect_equal(nrow(r$edges), 0)
})

test_that("line outside grid extent emits no edges", {
  skip_if_not_installed("geos")
  line <- geos::as_geos_geometry("LINESTRING (-2 -2, -1 -1)")
  r <- burn_scanline(line, extent = c(0, 3, 0, 3), dimension = c(3, 3))
  expect_equal(nrow(r$runs), 0)
  expect_equal(nrow(r$edges), 0)
})

test_that("GeometryCollection input is rejected with a warning", {
  skip_if_not_installed("geos")
  gc <- geos::as_geos_geometry(
    "GEOMETRYCOLLECTION (POINT (0.5 0.5), LINESTRING (0.5 0.5, 1.5 1.5))"
  )
  expect_warning(
    r <- burn_scanline(gc, extent = c(0, 3, 0, 3), dimension = c(3, 3)),
    "GeometryCollection"
  )
  expect_equal(nrow(r$runs), 0)
  expect_equal(nrow(r$edges), 0)
})

test_that("mixed POLYGON + LINESTRING input via separate burns", {
  # Demonstrates the documented workflow: split mixed-dimension inputs
  # into homogeneous groups and run separate burns. The output tables
  # both have the same shape (row, col, weight, id) but the `weight`
  # column carries area for polygons and length for lines.
  skip_if_not_installed("geos")
  poly <- geos::as_geos_geometry(
    "POLYGON ((0.5 0.5, 1.5 0.5, 1.5 1.5, 0.5 1.5, 0.5 0.5))"
  )
  line <- geos::as_geos_geometry("LINESTRING (1.5 1.5, 2.5 2.5)")

  rp <- burn_scanline(poly, extent = c(0, 3, 0, 3), dimension = c(3, 3))
  rl <- burn_scanline(line, extent = c(0, 3, 0, 3), dimension = c(3, 3))

  expect_true(nrow(rp$edges) > 0 || nrow(rp$runs) > 0)
  expect_true(nrow(rl$edges) > 0)
  expect_equal(nrow(rl$runs), 0)
})
