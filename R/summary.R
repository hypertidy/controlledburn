#' Summarise a burn without materialising it
#'
#' Compute a compact description of a [burn()] result in time proportional
#' to the number of records in its sparse tables, never to the number of
#' grid cells, so it is cheap even for a national grid of billions of cells.
#'
#' The pieces are:
#'
#' \describe{
#'   \item{`overview`}{a coarse matrix of mean coverage per tile, see
#'     [overview_burn()]. When the grid is small enough that one tile is one
#'     cell, this is the exact coverage.}
#'   \item{`row_marginal`, `col_marginal`}{exact covered-cell totals per grid
#'     row and per grid column, at full resolution.}
#'   \item{`content`}{the extent and 1-based row/col range of touched cells,
#'     or `NULL` for an empty burn.}
#'   \item{`area`}{per-geometry totals in cells: `interior` (run cells),
#'     `boundary` (summed edge fractions), `lines` (summed length, in
#'     extent units), `points` (count).}
#'   \item{`fractions`}{quantiles of the edge coverage fractions.}
#'   \item{`records`}{row counts of the sparse tables.}
#' }
#'
#' Coverage here means: run cells 1 each, edge cells their fraction, line
#' and point cells 1 each (presence).
#'
#' The number of geometries is taken as the largest id present, so
#' geometries after the last one that touched the grid are not counted.
#'
#' @param object A `controlledburn` object returned by [burn()].
#' @param tiles Target number of tiles along the longer grid axis for the
#'   overview, passed to [overview_burn()].
#' @param ... Ignored.
#'
#' @return A list of class `summary.controlledburn`.
#' @export
#' @examples
#' if (requireNamespace("geos", quietly = TRUE)) {
#'   poly <- geos::as_geos_geometry(c("POLYGON ((1 1, 9 1, 9 9, 1 9, 1 1))",
#'                                    "POLYGON ((2 2, 4 3, 3 5, 2 2))"))
#'   r <- burn(poly, extent = c(0, 10, 0, 10), dimension = c(100L, 100L))
#'   s <- summary(r)
#'   s
#'   s$area
#' }
summary.controlledburn <- function(object, tiles = 256L, ...) {
  x <- object
  dm <- as.integer(x$dimension)
  ncol <- dm[1L]
  nrow <- dm[2L]
  runs <- x$runs
  edges <- x$edges
  lines <- .burn_table(x, "lines")
  points <- .burn_table(x, "points")
  len <- runs$col_end - runs$col_start + 1L

  ## row marginal: rowsum by row, O(records)
  row_marginal <- numeric(nrow)
  row_marginal <- .accumulate_into(row_marginal, runs$row, len)
  row_marginal <- .accumulate_into(row_marginal, edges$row, edges$fraction)
  row_marginal <- .accumulate_into(row_marginal, lines$row, 1)
  row_marginal <- .accumulate_into(row_marginal, points$row, 1)

  ## column marginal: difference array for runs, O(records + ncol)
  d <- numeric(ncol + 1L)
  d <- .accumulate_into(d, runs$col_start, 1)
  d <- .accumulate_into(d, runs$col_end + 1L, -1)
  col_marginal <- cumsum(d)[seq_len(ncol)]
  col_marginal <- .accumulate_into(col_marginal, edges$col, edges$fraction)
  col_marginal <- .accumulate_into(col_marginal, lines$col, 1)
  col_marginal <- .accumulate_into(col_marginal, points$col, 1)

  n <- max(0L, runs$id, edges$id, lines$id, points$id)
  area <- data.frame(
    id = seq_len(n),
    interior = .accumulate_into(numeric(n), runs$id, len),
    boundary = .accumulate_into(numeric(n), edges$id, edges$fraction),
    lines = .accumulate_into(numeric(n), lines$id, lines$length),
    points = .accumulate_into(numeric(n), points$id, 1)
  )

  fractions <- if (nrow(edges) > 0) {
    stats::quantile(edges$fraction, c(0, 0.25, 0.5, 0.75, 1), names = FALSE)
  } else {
    rep(NA_real_, 5L)
  }

  structure(
    list(
      extent = x$extent,
      dimension = dm,
      records = c(runs = nrow(runs), edges = nrow(edges),
                  lines = nrow(lines), points = nrow(points)),
      n_geom = n,
      content = .content_extent(x),
      area = area,
      fractions = fractions,
      row_marginal = row_marginal,
      col_marginal = col_marginal,
      overview = overview_burn(x, tiles = tiles)
    ),
    class = "summary.controlledburn"
  )
}

#' @param x A `summary.controlledburn` object.
#' @param ascii Also print the overview as text.
#' @rdname summary.controlledburn
#' @export
print.summary.controlledburn <- function(x, ascii = FALSE, ...) {
  ext <- x$extent
  dm <- x$dimension
  cat("<controlledburn summary>\n")
  cat(sprintf("  grid:     %d col x %d row\n", dm[1L], dm[2L]))
  cat(sprintf("  extent:   x[%g, %g] y[%g, %g]\n", ext[1L], ext[2L], ext[3L], ext[4L]))
  if (is.null(x$content)) {
    cat("  content:  (empty)\n")
  } else {
    ce <- x$content$extent
    cat(sprintf("  content:  cols %d..%d, rows %d..%d  x[%g, %g] y[%g, %g]\n",
                x$content$cols[1L], x$content$cols[2L],
                x$content$rows[1L], x$content$rows[2L],
                ce[1L], ce[2L], ce[3L], ce[4L]))
  }
  cat("  records: ", paste(sprintf("%s=%d", names(x$records), x$records),
                           collapse = " "), "\n")
  touched <- with(x$area, interior > 0 | boundary > 0 | lines > 0 | points > 0)
  cat(sprintf("  geoms:    %d (largest id), %d burned, %d empty\n",
              x$n_geom, sum(touched), sum(!touched)))
  cat(sprintf("  area:     interior %s cells, boundary %s cells",
              .fmt_num(sum(x$area$interior)), .fmt_num(sum(x$area$boundary))))
  if (!anyNA(x$fractions)) {
    cat(sprintf(" (fraction min %.3g, median %.3g, max %.3g)",
                x$fractions[1L], x$fractions[3L], x$fractions[5L]))
  }
  cat("\n")
  if (sum(x$area$lines) > 0) {
    cat(sprintf("  lines:    %s length units\n", .fmt_num(sum(x$area$lines))))
  }
  if (sum(x$area$points) > 0) {
    cat(sprintf("  points:   %s in grid\n", .fmt_num(sum(x$area$points))))
  }
  ov <- dim(x$overview)
  cat(sprintf("  overview: %d x %d tiles of %d cells\n", ov[2L], ov[1L],
              attr(x$overview, "block")))
  if (ascii) cat(.ascii_overview(x$overview, ...), sep = "\n")
  invisible(x)
}

#' Coarse overview of a burn
#'
#' Aggregate a [burn()] result onto a grid of square tiles of `block` cells
#' and return the mean coverage per tile (in `[0, 1]` for non-overlapping
#' polygons). Tiles cut off by the grid boundary are normalised by their
#' actual cell count. Cost is proportional to the number of records plus
#' the number of tiles: a run that spans several tiles contributes partial
#' lengths to its two end tiles, and the fully covered tiles between them
#' are filled with a difference array, so nothing is expanded to cells
#' whatever the block size.
#'
#' The result is a matrix in raster order (row 1 is the top), ready for
#' `ximage::ximage(ov, extent = attr(ov, "extent"))`.
#'
#' @param x A `controlledburn` object returned by [burn()].
#' @param tiles Target number of tiles along the longer grid axis; ignored
#'   if `block` is given.
#' @param block Tile size in cells. `block = 1` gives the exact coverage,
#'   the same as [materialise_chunk()] for polygons, with `NA` for 0.
#'
#' @return A numeric matrix (tile rows by tile columns) with attributes
#'   `extent` (`c(xmin, xmax, ymin, ymax)` of the tile grid, which extends
#'   past the burn grid on the right and bottom when `block` does not
#'   divide its dimension) and `block`. Untouched tiles are `NA`.
#' @export
#' @examples
#' if (requireNamespace("geos", quietly = TRUE)) {
#'   poly <- geos::as_geos_geometry("POLYGON ((1 1, 9 1, 9 9, 1 9, 1 1))")
#'   r <- burn(poly, extent = c(0, 10, 0, 10), dimension = c(1000L, 1000L))
#'   ov <- overview_burn(r, tiles = 20)
#'   dim(ov)
#'   attr(ov, "extent")
#' }
overview_burn <- function(x, tiles = 256L, block = NULL) {
  stopifnot(inherits(x, "controlledburn"))
  ext <- x$extent
  dm <- as.integer(x$dimension)
  ncol <- dm[1L]
  nrow <- dm[2L]
  b <- if (is.null(block)) {
    max(1L, as.integer(ceiling(max(dm) / tiles)))
  } else {
    as.integer(block)
  }
  ntx <- as.integer(ceiling(ncol / b))
  nty <- as.integer(ceiling(nrow / b))
  edges <- x$edges
  lines <- .burn_table(x, "lines")
  points <- .burn_table(x, "points")

  ## tiles are numbered row-major from the top left
  tile_index <- function(row, col) {
    ((row - 1L) %/% b) * ntx + ((col - 1L) %/% b) + 1L
  }

  ## a run touches tile columns tc0..tc1 of tile row tr; its two end tiles
  ## get partial lengths and every tile strictly between them is full
  runs <- x$runs
  tr <- (runs$row - 1L) %/% b
  tc0 <- (runs$col_start - 1L) %/% b
  tc1 <- (runs$col_end - 1L) %/% b
  base <- tr * ntx + 1L
  single <- tc0 == tc1
  multi <- !single

  acc <- numeric(ntx * nty)
  acc <- .accumulate_into(acc, (base + tc0)[single],
                          (runs$col_end - runs$col_start + 1L)[single])
  acc <- .accumulate_into(acc, (base + tc0)[multi],
                          ((tc0 + 1L) * b - runs$col_start + 1L)[multi])
  acc <- .accumulate_into(acc, (base + tc1)[multi],
                          (runs$col_end - tc1 * b)[multi])
  span <- multi & (tc1 - tc0 > 1L)
  if (any(span)) {
    d <- numeric(ntx * nty + 1L)
    d <- .accumulate_into(d, (base + tc0 + 1L)[span], b)
    d <- .accumulate_into(d, (base + tc1)[span], -b)
    ## cumsum within each tile row: a tile row is a column of matrix(d, ntx)
    full <- apply(matrix(d[seq_len(ntx * nty)], ntx, nty), 2L, cumsum)
    acc <- acc + as.vector(full)
  }
  acc <- .accumulate_into(acc, tile_index(edges$row, edges$col), edges$fraction)
  acc <- .accumulate_into(acc, tile_index(lines$row, lines$col), 1)
  acc <- .accumulate_into(acc, tile_index(points$row, points$col), 1)

  ## actual cell count per tile (tiles on the right and bottom are smaller)
  tw <- pmin(b, ncol - (seq_len(ntx) - 1L) * b)
  th <- pmin(b, nrow - (seq_len(nty) - 1L) * b)
  cells <- as.vector(t(outer(th, tw)))

  out <- acc / cells
  out[acc == 0] <- NA_real_
  out <- matrix(out, nrow = nty, ncol = ntx, byrow = TRUE)

  resx <- vaster::x_res(dm, ext)
  resy <- vaster::y_res(dm, ext)
  attr(out, "extent") <- c(ext[1L], ext[1L] + ntx * b * resx,
                           ext[4L] - nty * b * resy, ext[4L])
  attr(out, "block") <- b
  out
}

## lines and points tables may be absent from older or cropped objects
.burn_table <- function(x, name) {
  t <- x[[name]]
  if (is.null(t)) {
    t <- data.frame(row = integer(), col = integer(), id = integer())
    if (name == "lines") t$length <- numeric()
  }
  t
}

## extent and row/col range of touched cells, or NULL for an empty burn
.content_extent <- function(x) {
  ext <- x$extent
  dm <- x$dimension
  lines <- .burn_table(x, "lines")
  points <- .burn_table(x, "points")
  rows <- c(x$runs$row, x$edges$row, lines$row, points$row)
  if (!length(rows)) return(NULL)
  cols <- c(x$runs$col_start, x$runs$col_end, x$edges$col, lines$col, points$col)
  rr <- range(rows)
  cc <- range(cols)
  resx <- vaster::x_res(dm, ext)
  resy <- vaster::y_res(dm, ext)
  list(
    rows = rr, cols = cc,
    extent = c(ext[1L] + (cc[1L] - 1L) * resx, ext[1L] + cc[2L] * resx,
               ext[4L] - rr[2L] * resy, ext[4L] - (rr[1L] - 1L) * resy)
  )
}

## out[idx] += w, summing duplicates in idx; w recycled to length(idx)
.accumulate_into <- function(out, idx, w) {
  if (!length(idx)) return(out)
  s <- rowsum(rep_len(as.numeric(w), length(idx)), idx, reorder = FALSE)
  at <- as.integer(rownames(s))
  out[at] <- out[at] + s[, 1L]
  out
}

.fmt_num <- function(v) {
  format(signif(v, 4), big.mark = ",", scientific = FALSE, trim = TRUE)
}

## text rendering of an overview at most `width` characters wide
.ascii_overview <- function(m, width = 72L, ramp = " .:-=+*#%@") {
  step <- max(1L, ceiling(ncol(m) / width))
  ## characters are about twice as tall as wide: pool 2 * step rows per
  ## step columns
  ci <- (seq_len(ncol(m)) - 1L) %/% step
  ri <- (seq_len(nrow(m)) - 1L) %/% (2L * step)
  m <- t(rowsum(t(rowsum(m, ri, na.rm = TRUE)), ci, na.rm = TRUE)) /
    outer(tabulate(ri + 1L), tabulate(ci + 1L))
  chars <- strsplit(ramp, "")[[1]]
  lev <- pmin(length(chars), 1L + floor(pmin(m, 1) * (length(chars) - 1L) + 0.5))
  lev[is.na(lev)] <- 1L
  apply(matrix(chars[lev], nrow(m)), 1L, paste, collapse = "")
}
