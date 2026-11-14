#' Slice a burn into tiles
#'
#' Cuts a [burn()] result into a regular grid of tiles in one pass over its
#' records, anchored at the top-left cell of the grid. A record in cell
#' `(row, col)` goes to the tile that holds it, with indices local to that
#' tile, and a run is split where it crosses a tile seam. Nothing is
#' recomputed, so each tile is [crop_burn()] of the whole result to the
#' tile's extent, record for record, and the materialised tiles put back
#' together are the materialised whole grid.
#'
#' Calling [crop_burn()] once per tile gives the same tables but scans
#' every record for every tile; `slice_burn()` costs about one pass over
#' the records whatever the number of tiles.
#'
#' @param x A `controlledburn` object returned by [burn()].
#' @param dimension integer, the tile size in cells, `c(ncol, nrow)`, or
#'   one number for square tiles. Tiles on the right and bottom edges are
#'   smaller when this does not divide the grid.
#' @param threads integer, the number of threads (default `1`; `0` for all
#'   available). The tiles are the same. Results of fewer than 262,144
#'   records are sliced on one thread whatever `threads` is.
#'
#' @return A list of `controlledburn` objects, one per tile holding at
#'   least one record, in row-major tile order (tiles without records are
#'   left out). Each has the usual tables, its own `extent` and
#'   `dimension`, and also
#'   \describe{
#'     \item{`tile`}{integer `c(row = , col = )`, the tile's 1-based
#'       position in the tile grid, row 1 at the top.}
#'     \item{`offset`}{integer `c(row = , col = )`, the cells before the
#'       tile in the whole grid: cell `(r, c)` of the tile is cell
#'       `(r + offset["row"], c + offset["col"])` of `x`.}
#'   }
#' @export
#' @examples
#' if (requireNamespace("geos", quietly = TRUE)) {
#'   poly <- geos::as_geos_geometry("POLYGON ((1 1, 9 1, 9 9, 1 9, 1 1))")
#'   r <- burn(poly, extent = c(0, 10, 0, 10), dimension = c(100L, 100L))
#'   tiles <- slice_burn(r, c(32L, 32L))
#'   length(tiles)
#'   tiles[[1]]$tile
#'   tiles[[1]]$extent
#'   ## each tile is the crop of the whole burn to the tile
#'   identical(materialise_chunk(tiles[[1]]),
#'             materialise_chunk(r, target = tiles[[1]]$extent))
#' }
slice_burn <- function(x, dimension = c(256L, 256L), threads = 1L) {
  stopifnot(inherits(x, "controlledburn"))
  if (!is.numeric(dimension) || !length(dimension) %in% 1:2 || anyNA(dimension) ||
      any(dimension < 1) || any(dimension != floor(dimension))) {
    stop("'dimension' must be one or two positive whole numbers, c(ncol, nrow)")
  }
  dimension <- rep_len(as.integer(dimension), 2L)
  threads <- .count(threads, "threads")
  ext <- x$extent
  dim <- x$dimension
  tiles <- cpp_slice_tiles(x, ext[1], ext[3], ext[2], ext[4], dim[1], dim[2],
                           dimension[1], dimension[2], threads)
  lapply(tiles, function(t) {
    g <- t$grid
    out <- list(
      runs = t$runs, edges = t$edges, lines = t$lines, points = t$points,
      extent = g[c(1, 3, 2, 4)],
      dimension = as.integer(g[5:6]),
      tile = c(row = t$tile[1] + 1L, col = t$tile[2] + 1L),
      offset = c(row = t$offset[1], col = t$offset[2])
    )
    if (!is.null(x$options)) out$options <- x$options
    class(out) <- "controlledburn"
    out
  })
}
