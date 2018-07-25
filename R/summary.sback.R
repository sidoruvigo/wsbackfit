#' @name summary.sback
#' @aliases summary.sback
#' @title Summary for a \code{sback} fitted object.
#' @description Takes a fitted \code{sback} object produced by \code{sback()} and produces various useful summaries from it.
#' @param object an object of class \code{sback} as produced by \code{sback()}.
#' @param \dots other arguments (not implemented).
#' @examples
#' ##---- Should be DIRECTLY executable !! ----
#' ##-- ==>  Define data, use random,
#' ##--	or do  help(data=index)  for the standard data sets.
#'
#' ## The function is currently defined as
#' function (object, ...)
#' {
#'   class(object) <- "summary.sback"
#'   object
#' }
#' @export
summary.sback <- function(object, ...) {
  class(object) <- "summary.sback"
  object
}
