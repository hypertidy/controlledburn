# Helper: check coverage complementarity across polygon ids.
# At every cell touched by any polygon, the sum of coverage fractions
# should equal 1.0 (no gaps, no overlaps).
expect_complementary <- function(result, tol = 1e-5, label = NULL) {
  nc <- result$dimension[1]
  nr <- result$dimension[2]

  total <- matrix(0, nrow = nr, ncol = nc)

  # Runs contribute 1.0 per cell
  if (nrow(result$runs) > 0) {
    for (i in seq_len(nrow(result$runs))) {
      r <- result$runs$row[i]
      cs <- result$runs$col_start[i]
      ce <- result$runs$col_end[i]
      total[r, cs:ce] <- total[r, cs:ce] + 1.0
    }
  }

  # Edges contribute their fraction
  if (nrow(result$edges) > 0) {
    for (i in seq_len(nrow(result$edges))) {
      r <- result$edges$row[i]
      c <- result$edges$col[i]
      total[r, c] <- total[r, c] + result$edges$fraction[i]
    }
  }

  any_coverage <- total > tol
  bad <- any_coverage & abs(total - 1.0) > tol
  expect_equal(sum(bad), 0L, label = label)
}
