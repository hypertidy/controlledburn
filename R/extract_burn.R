#' Extract values from a controlledburn result at point locations
#'
#' Queries the sparse tables directly; no dense matrix is created.
#' For points at cell centres this reproduces [materialize_chunk()]
#' exactly:
#' `extract_burn(x, xy, fun) == materialize_chunk(x, fun)[cbind(row, col)]`
#'
#' Interior runs are matched by row equality plus column interval
#' containment; edges, lines, and points by exact (row, col) match.
#' Cost is O(points + hits), independent of grid dimension.
#'
#' @param x A `controlledburn` object returned by [burn()].
#' @param xy A two-column matrix or data.frame of x, y coordinates in
#'   the CRS of the burn.
#' @param fun character, `"id"` (last geometry written wins, `NA`
#'   background) or `"sum"` (accumulated coverage/length/count, `0`
#'   background). Same semantics as [materialize_chunk()].
#' @return A numeric vector of length `nrow(xy)`. Points outside the
#'   grid extent are `NA` for both funs.
#' @export
extract_burn <- function(x, xy, fun = c("id", "sum")) {
  stopifnot(inherits(x, "controlledburn"))
  fun <- match.arg(fun)
  xy <- as.matrix(xy)[, 1:2, drop = FALSE]
  storage.mode(xy) <- "double"

  if (fun == "sum") .check_single_kind(x)

  rc <- .cell_from_xy(x, xy)
  n <- nrow(rc)
  inside <- !is.na(rc[, 1L])

  out <- rep(if (fun == "id") NA_real_ else 0, n)
  out[!inside] <- NA_real_
  if (!any(inside)) return(out)

  q <- data.frame(qi = which(inside),
                  row = rc[inside, 1L],
                  col = rc[inside, 2L])

  ## interior runs: row equality + column interval containment
  run_hits <- NULL
  if (!is.null(x$runs) && nrow(x$runs) > 0) {
    r <- x$runs
    r$ti <- seq_len(nrow(r))
    m <- merge(q, r, by = "row")
    run_hits <- m[m$col >= m$col_start & m$col <= m$col_end, ,
                  drop = FALSE]
  }

  ## single-cell tables: exact (row, col) join
  cell_hits <- function(df) {
    if (is.null(df) || nrow(df) == 0) return(NULL)
    df$ti <- seq_len(nrow(df))
    merge(q, df, by = c("row", "col"))
  }
  edge_hits  <- cell_hits(x$edges)
  line_hits  <- cell_hits(x$lines)
  point_hits <- cell_hits(x$points)

  if (fun == "id") {
    ## replicate materialize_chunk write order: runs, edges, lines,
    ## points, each ascending table index; last write wins
    apply_id <- function(out, hits) {
      if (is.null(hits) || nrow(hits) == 0) return(out)
      hits <- hits[order(hits$ti), , drop = FALSE]
      out[hits$qi] <- hits$id
      out
    }
    out <- apply_id(out, run_hits)
    out <- apply_id(out, edge_hits)
    out <- apply_id(out, line_hits)
    out <- apply_id(out, point_hits)
  } else {
    add <- function(out, hits, val) {
      if (is.null(hits) || nrow(hits) == 0) return(out)
      s <- rowsum(val, hits$qi)
      idx <- as.integer(rownames(s))
      out[idx] <- out[idx] + s[, 1L]
      out
    }
    out <- add(out, run_hits, rep(1, if (is.null(run_hits)) 0 else nrow(run_hits)))
    out <- add(out, edge_hits, edge_hits$fraction)
    out <- add(out, line_hits, line_hits$length)
    out <- add(out, point_hits, rep(1, if (is.null(point_hits)) 0 else nrow(point_hits)))
  }
  out
}

## map coordinates to 1-based (row, col); row 1 is the top (ymax) row;
## points exactly on the xmax/ymin edge belong to the last col/row
.cell_from_xy <- function(x, xy) {
  ext <- x$extent
  dim <- x$dimension
  dx <- (ext[2] - ext[1]) / dim[1]
  dy <- (ext[4] - ext[3]) / dim[2]
  col <- floor((xy[, 1] - ext[1]) / dx) + 1
  row <- floor((ext[4] - xy[, 2]) / dy) + 1
  col[!is.na(col) & xy[, 1] == ext[2]] <- dim[1]
  row[!is.na(row) & xy[, 2] == ext[3]] <- dim[2]
  bad <- is.na(col) | is.na(row) |
    col < 1 | col > dim[1] | row < 1 | row > dim[2]
  row[bad] <- NA
  col[bad] <- NA
  cbind(as.integer(row), as.integer(col))
}

## fun = "sum" mixes units across kinds (fraction vs CRS length vs
## count), mirroring the materialize_chunk() refusal
.check_single_kind <- function(x) {
  kinds <- c(
    polygon = (!is.null(x$runs)  && nrow(x$runs)  > 0) ||
      (!is.null(x$edges) && nrow(x$edges) > 0),
    line    =  !is.null(x$lines)  && nrow(x$lines)  > 0,
    point   =  !is.null(x$points) && nrow(x$points) > 0
  )
  if (sum(kinds) > 1) {
    stop("fun = \"sum\" on mixed geometry kinds would mix units (",
         paste(names(kinds)[kinds], collapse = ", "),
         "); filter to one kind first", call. = FALSE)
  }
  invisible(TRUE)
}
