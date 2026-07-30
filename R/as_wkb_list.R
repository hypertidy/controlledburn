# Internal generic: convert geometry input to a wk_wkb vector.
# Methods handle sf (sfc), geos, wk_wkb, blob, and raw-list input.

as_wkb_list <- function(x) {
  UseMethod("as_wkb_list")
}

#' @export
as_wkb_list.list <- function(x) {
  stopifnot(all(vapply(x, is.raw, logical(1))))
  wk::wkb(x)
}

#' @export
as_wkb_list.blob <- function(x) {
  wk::wkb(x)
}

#' @export
as_wkb_list.wk_wkb <- function(x) {
  x
}

#' @export
as_wkb_list.sfc <- function(x) {
  wk::as_wkb(x)
}

#' @export
as_wkb_list.geos_geometry <- function(x) {
  wk::as_wkb(x)
}
