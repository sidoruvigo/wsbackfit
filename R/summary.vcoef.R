#' @name summary.vcoef
#' @aliases summary.vcoef
#' @title Summary for a \code{vcoef} fitted object.
#' @description Takes a fitted \code{vcoef} object produced by \code{vcoef()} and produces various useful summaries from it.
#' @param object an object of class \code{vcoef} as produced by \code{vcoef()}.
#' @param \dots other arguments (not implemented).
#' @export
summary.vcoef <- function(object, ...) {
  class(object) <- "summary.vcoef"
  object
}
