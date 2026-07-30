# test-parity-fixtures.R -- cross-language parity tests driven by shared
# fixtures in fixtures/. Identical fixtures are read by C++ and Python
# test suites so any regression is caught across all three surfaces.

fixtures_dir <- file.path(testthat::test_path(), "..", "..", "fixtures")
skip_if_not(dir.exists(fixtures_dir), "fixtures/ directory not found")

geoms  <- read.csv(file.path(fixtures_dir, "geometries.csv"),
                    stringsAsFactors = FALSE)
expect <- read.csv(file.path(fixtures_dir, "expected.csv"),
                   stringsAsFactors = FALSE)
fixtures <- merge(geoms, expect, by = "case")

for (i in seq_len(nrow(fixtures))) {
  f <- fixtures[i, ]

  test_that(paste0("parity fixture: ", f$case), {
    wkb_raw <- wk::wk_handle(wk::wkt(f$wkt), wk::wkb_writer())
    gs <- controlledburn:::grid_params_from_ext(
      c(f$xmin, f$xmax, f$ymin, f$ymax), c(f$ncol, f$nrow)
    )
    result <- controlledburn::burn_scanline(wkb_raw, gs$extent, gs$dimension)

    if (!is.na(f$covered_area)) {
      cell_area <- ((f$xmax - f$xmin) / f$ncol) * ((f$ymax - f$ymin) / f$nrow)
      total <- 0
      if (nrow(result$runs) > 0) {
        total <- total + sum(result$runs$col_end - result$runs$col_start + 1) * cell_area
      }
      if (nrow(result$edges) > 0) {
        total <- total + sum(result$edges$fraction) * cell_area
      }
      expect_equal(total, f$covered_area, tolerance = f$tol_rel)
    }

    if (!is.na(f$edges_empty) && f$edges_empty == 1) {
      expect_equal(nrow(result$edges), 0L)
    }

    if (!is.na(f$line_length)) {
      expect_equal(sum(result$lines$length), f$line_length,
                   tolerance = f$tol_rel)
    }

    if (!is.na(f$n_points)) {
      expect_equal(nrow(result$points), as.integer(f$n_points))
    }
  })
}
