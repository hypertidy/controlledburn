#' Plot a burn
#'
#' Draw a [burn()] result in map coordinates (`asp = 1`) with the grid
#' extent as the frame. Two renderings are available:
#'
#' \describe{
#'   \item{`"exact"`}{every record is drawn: runs as rectangles, edge cells
#'     shaded by coverage fraction, line cells shaded by length (relative
#'     to the cell diagonal), points as cells. Cost is one rectangle per
#'     record, so this is the picture of the sparse representation itself.}
#'   \item{`"overview"`}{the coarse mean-coverage matrix from
#'     [overview_burn()] drawn with [ximage::ximage()]. Cost does not
#'     depend on grid size.}
#' }
#'
#' `"auto"` (default) picks `"exact"` when the total record count is at
#' most `threshold`, else `"overview"`. Either way the optional marginal
#' strips are exact: covered cells per grid row (right) and per grid
#' column (top), from [summary.controlledburn()].
#'
#' @param x A `controlledburn` object returned by [burn()].
#' @param what `"auto"`, `"exact"`, or `"overview"`.
#' @param marginals Draw the row and column marginal strips. The layout is
#'   sized from the extent's aspect ratio so the strips span exactly the
#'   map; leftover device space is left blank.
#' @param col Colour ramp (low to high coverage) for `"overview"` and for
#'   edge shading in `"exact"`.
#' @param by_id In `"exact"` mode, colour cells by geometry `id` instead of
#'   by coverage (edge cells take their fraction as opacity).
#' @param tiles Passed to [overview_burn()] for the `"overview"` rendering.
#' @param threshold Record count above which `"auto"` switches to
#'   `"overview"`.
#' @param content Draw the content extent (bounding box of touched cells)
#'   as a dashed rectangle.
#' @param main Optional title, drawn in the outer margin above the strips.
#' @param ... Passed to [graphics::plot()] for the main panel.
#'
#' @return `x`, invisibly, with the summary computed for the marginals
#'   attached as attribute `summary`.
#' @export
#' @examples
#' if (requireNamespace("geos", quietly = TRUE)) {
#'   poly <- geos::as_geos_geometry(c("POLYGON ((1 1, 9 1, 9 9, 1 9, 1 1))",
#'                                    "POLYGON ((2 2, 4 3, 3 5, 2 2))"))
#'   r <- burn(poly, extent = c(0, 10, 0, 10), dimension = c(40L, 40L))
#'   plot(r)
#'   plot(r, by_id = TRUE, marginals = FALSE)
#'   plot(r, what = "overview", tiles = 10)
#' }
plot.controlledburn <- function(x, what = c("auto", "exact", "overview"),
                                marginals = TRUE,
                                col = grDevices::hcl.colors(64, "YlOrRd", rev = TRUE),
                                by_id = FALSE, tiles = 256L, threshold = 5e4,
                                content = FALSE, main = NULL, ...) {
  what <- match.arg(what)
  ext <- x$extent
  dm <- x$dimension
  nrec <- sum(vapply(c("runs", "edges", "lines", "points"),
                     function(k) nrow(.burn_table(x, k)), integer(1)))
  if (what == "auto") what <- if (nrec <= threshold) "exact" else "overview"

  s <- summary(x, tiles = tiles)

  if (marginals) {
    op <- graphics::par(no.readonly = TRUE)
    on.exit(graphics::par(op), add = TRUE)
    .marginal_layout(ext, strip = 1 / 6, main = main)
    graphics::par(mar = c(3.5, 3.5, 0.5, 0.5), mgp = c(2, 0.6, 0))
  } else if (!is.null(main)) {
    op <- graphics::par(mar = c(3.5, 3.5, 2.5, 0.5))
    on.exit(graphics::par(op), add = TRUE)
  }

  graphics::plot(NA, xlim = ext[1:2], ylim = ext[3:4], asp = 1,
                 xlab = "", ylab = "", xaxs = "i", yaxs = "i", bty = "n", ...)
  ## asp = 1 widens one axis to fill the panel: the strips must share these
  ## limits, not the extent, to stay aligned with the map
  usr <- graphics::par("usr")

  if (what == "overview") {
    ov <- s$overview
    ximage::ximage(ov, extent = attr(ov, "extent"), add = TRUE, col = col,
                   zlim = c(0, max(1, ov, na.rm = TRUE)))
  } else {
    .draw_exact(x, col = col, by_id = by_id, n_geom = s$n_geom)
  }
  graphics::rect(ext[1L], ext[3L], ext[2L], ext[4L], border = "grey30")
  if (content && !is.null(s$content)) {
    ce <- s$content$extent
    graphics::rect(ce[1L], ce[3L], ce[2L], ce[4L], border = "grey30", lty = 2)
  }

  if (marginals) {
    xc <- vaster::x_from_col(dm, ext, seq_len(dm[1L]))
    yc <- vaster::y_from_row(dm, ext, seq_len(dm[2L]))
    cmax <- max(1, s$col_marginal)
    rmax <- max(1, s$row_marginal)

    ## top strip: column marginal against x
    graphics::par(mar = c(0.5, 3.5, 0.5, 0.5))
    graphics::plot(NA, xlim = usr[1:2], ylim = c(0, cmax),
                   xaxs = "i", yaxs = "i", axes = FALSE, xlab = "", ylab = "")
    graphics::polygon(c(ext[1L], xc, ext[2L]), c(0, s$col_marginal, 0),
                      col = "grey70", border = NA)
    graphics::axis(2, at = c(0, cmax), las = 1, cex.axis = 0.7,
                   labels = c("", .fmt_num(max(s$col_marginal))))

    ## right strip: row marginal against y
    graphics::par(mar = c(3.5, 0.5, 0.5, 0.5))
    graphics::plot(NA, ylim = usr[3:4], xlim = c(0, rmax),
                   xaxs = "i", yaxs = "i", axes = FALSE, xlab = "", ylab = "")
    graphics::polygon(c(0, s$row_marginal, 0), c(ext[4L], yc, ext[3L]),
                      col = "grey70", border = NA)
    graphics::axis(1, at = c(0, rmax), cex.axis = 0.7,
                   labels = c("", .fmt_num(max(s$row_marginal))))
  }

  if (!is.null(main)) graphics::title(main = main, outer = marginals)

  attr(x, "summary") <- s
  invisible(x)
}

## Layout sized to the extent's aspect ratio: main panel bottom left, column
## strip above it, row strip to its right, and spacer cells for whatever the
## device has left over so the strips are never longer than the map.
.marginal_layout <- function(ext, strip = 1 / 6, main = NULL) {
  ## the title lives in the outer margin; layout() divides what is left
  oma_top <- if (is.null(main)) 0 else 2
  graphics::par(oma = c(0, 0, oma_top, 0))
  din <- grDevices::dev.size("in")
  ## shave a little so rounding never asks for more than the device has
  W <- din[1L] * 0.99
  H <- (din[2L] - oma_top * graphics::par("csi")) * 0.99
  aspect <- (ext[4L] - ext[3L]) / (ext[2L] - ext[1L])
  sw <- strip * min(W, H)
  ## the main panel's margins (mar = c(3.5, 3.5, 0.5, 0.5) lines) are part
  ## of its layout cell, so size the plot region inside them to the
  ## extent's aspect, as large as fits
  mw <- 4 * graphics::par("csi")
  wm <- min(W - sw - mw, (H - sw - mw) / aspect)
  hm <- wm * aspect
  widths <- c(wm + mw, sw, max(0.01, W - wm - mw - sw))
  heights <- c(sw, hm + mw, max(0.01, H - hm - mw - sw))
  ## panels: 1 main, 2 column strip, 3 row strip; 0 empty
  mat <- rbind(c(2L, 0L, 0L),
               c(1L, 3L, 0L),
               c(0L, 0L, 0L))
  ## relative sizes, not lcm(): a recorded plot replayed on a device of
  ## another size (knitr, pkgdown) then scales instead of failing
  graphics::layout(mat, widths = widths, heights = heights)
  invisible(NULL)
}

## one rectangle per record
.draw_exact <- function(x, col, by_id, n_geom) {
  ext <- x$extent
  dm <- x$dimension
  resx <- vaster::x_res(dm, ext)
  resy <- vaster::y_res(dm, ext)
  ## cbind(xmin, xmax, ymin, ymax) for cells in row, from col0 to col1
  cells <- function(row, col0, col1 = col0) {
    cbind(ext[1L] + (col0 - 1L) * resx, ext[1L] + col1 * resx,
          ext[4L] - row * resy, ext[4L] - (row - 1L) * resy)
  }
  ramp <- function(v) col[pmax(1L, pmin(length(col), ceiling(v * length(col))))]
  idcol <- function(id) grDevices::hcl.colors(max(1L, n_geom), "Dark 3")[id]
  draw <- function(xy, fill) ximage::xrect(xy, add = TRUE, col = fill, border = NA)

  r <- x$runs
  if (nrow(r) > 0) {
    draw(cells(r$row, r$col_start, r$col_end),
         if (by_id) idcol(r$id) else col[length(col)])
  }
  e <- x$edges
  if (nrow(e) > 0) {
    draw(cells(e$row, e$col),
         if (by_id) .with_alpha(idcol(e$id), e$fraction) else ramp(e$fraction))
  }
  l <- .burn_table(x, "lines")
  if (nrow(l) > 0) {
    draw(cells(l$row, l$col),
         if (by_id) idcol(l$id) else ramp(l$length / sqrt(resx^2 + resy^2)))
  }
  p <- .burn_table(x, "points")
  if (nrow(p) > 0) {
    draw(cells(p$row, p$col), if (by_id) idcol(p$id) else col[length(col)])
  }
  invisible(NULL)
}

## per-cell alpha on a vector of colours
.with_alpha <- function(cols, alpha) {
  m <- grDevices::col2rgb(cols)
  grDevices::rgb(m[1L, ], m[2L, ], m[3L, ],
                 alpha = pmin(1, pmax(0, alpha)) * 255, maxColorValue = 255)
}
