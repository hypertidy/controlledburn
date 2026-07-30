#' Crop a controlledburn result to a target extent
#'
#' Filters and clips the sparse tables (runs, edges, lines, points)
#' to a sub-window defined by a target extent. Returns a new
#' `controlledburn` object with offset row/col indices and snapped
#' extent. No dense allocation — just data frame filtering.
#'
#' @param x A `controlledburn` object returned by [burn()].
#' @param target numeric vector `c(xmin, xmax, ymin, ymax)` defining
#'   the target window. Snapped outward to cell boundaries.
#'
#' @return A `controlledburn` object covering only the target window,
#'   with row/col indices re-based to 1.
#'
#' @export
#' @examples
#' if (requireNamespace("geos", quietly = TRUE)) {
#'   library(geos)
#'   poly <- as_geos_geometry("POLYGON ((1 1, 9 1, 9 9, 1 9, 1 1))")
#'   r <- burn(poly, extent = c(0, 10, 0, 10), dimension = c(100L, 100L))
#'   sub <- crop_burn(r, target = c(2, 5, 3, 7))
#'   sub
#' }
crop_burn <- function(x, target) {
  stopifnot(inherits(x, "controlledburn"))
  stopifnot(length(target) == 4)

  ext <- x$extent
  dim <- x$dimension
  dx <- (ext[2] - ext[1]) / dim[1]
  dy <- (ext[4] - ext[3]) / dim[2]

  # Target extent to row/col bounds (snap outward, 1-based)
  col_lo <- max(1L, as.integer(floor((target[1] - ext[1]) / dx)) + 1L)
  col_hi <- min(dim[1], as.integer(ceiling((target[2] - ext[1]) / dx)))
  row_hi <- min(dim[2], as.integer(ceiling((ext[4] - target[3]) / dy)))
  row_lo <- max(1L, as.integer(floor((ext[4] - target[4]) / dy)) + 1L)

  if (col_lo > col_hi || row_lo > row_hi) {
    warning("target extent does not overlap the grid")
    new_ext <- target
    new_dim <- c(0L, 0L)
  } else {
    new_ext <- c(
      ext[1] + (col_lo - 1L) * dx,
      ext[1] + col_hi * dx,
      ext[4] - row_hi * dy,
      ext[4] - (row_lo - 1L) * dy
    )
    new_dim <- c(col_hi - col_lo + 1L, row_hi - row_lo + 1L)
  }

  crop_table <- function(df, row_col, col_cols = NULL) {
    if (is.null(df) || nrow(df) == 0) return(df)
    keep <- df$row >= row_lo & df$row <= row_hi
    if (!is.null(col_cols)) {
      if (length(col_cols) == 1) {
        keep <- keep & df[[col_cols]] >= col_lo & df[[col_cols]] <= col_hi
      } else {
        # Runs: col_end >= col_lo and col_start <= col_hi
        keep <- keep & df[[col_cols[2]]] >= col_lo & df[[col_cols[1]]] <= col_hi
      }
    }
    df <- df[keep, , drop = FALSE]
    df$row <- df$row - row_lo + 1L
    if (!is.null(col_cols)) {
      if (length(col_cols) == 1) {
        df[[col_cols]] <- df[[col_cols]] - col_lo + 1L
      } else {
        df[[col_cols[1]]] <- pmax(df[[col_cols[1]]], col_lo) - col_lo + 1L
        df[[col_cols[2]]] <- pmin(df[[col_cols[2]]], col_hi) - col_lo + 1L
      }
    }
    rownames(df) <- NULL
    df
  }

  out <- list(
    runs   = crop_table(x$runs, "row", c("col_start", "col_end")),
    edges  = crop_table(x$edges, "row", "col"),
    lines  = crop_table(x$lines, "row", "col"),
    points = crop_table(x$points, "row", "col"),
    extent = new_ext,
    dimension = new_dim
  )
  class(out) <- "controlledburn"
  out
}
