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
  if (!is.numeric(target) || length(target) != 4L || anyNA(target)) {
    stop("'target' must be numeric c(xmin, xmax, ymin, ymax)")
  }
  ext <- x$extent
  dim <- x$dimension
  res <- cpp_crop(x, ext[1], ext[3], ext[2], ext[4], dim[1], dim[2],
                  target[1], target[3], target[2], target[4])
  if (is.null(res)) {
    warning("target extent does not overlap the grid")
    res <- .empty_tables()
    new_ext <- target
    new_dim <- c(0L, 0L)
  } else {
    g <- res$grid
    new_ext <- g[c(1, 3, 2, 4)]
    new_dim <- as.integer(g[5:6])
  }
  out <- list(
    runs   = res$runs,
    edges  = res$edges,
    lines  = res$lines,
    points = res$points,
    extent = new_ext,
    dimension = new_dim
  )
  if (!is.null(x$options)) out$options <- x$options
  class(out) <- "controlledburn"
  out
}

# The four tables with no records.
.empty_tables <- function() {
  list(
    runs = data.frame(row = integer(), col_start = integer(),
                      col_end = integer(), id = integer()),
    edges = data.frame(row = integer(), col = integer(),
                       fraction = numeric(), id = integer()),
    lines = data.frame(row = integer(), col = integer(),
                       length = numeric(), id = integer()),
    points = data.frame(row = integer(), col = integer(), id = integer())
  )
}
