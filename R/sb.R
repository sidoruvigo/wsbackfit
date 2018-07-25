#' @name sb
#' @aliases sb
#' @title Specify a nonparametric and/or a varying coefficient term in a wsbackfit formula
#' @description Function used to indicate a nonparametric term in a \code{sback} formula and a varying coefficient term in a \code{vcoef} formula.
#' @param x1 the univariate predictor.
#' @param by numeric predictor of the same dimension as x1. If present, the coefficients of this predictor depend, nonparametrically, on x1, i.e., a varying coefficient term. Only applied for the \code{vcoef()} function.
#' @param h bandwidth for this term. If h = -1, the bandwidth is automatically selected using cross-validation. For the \code{sback()} fuction, h = 0 indicates a parametric/linear fit. Defaults to -1.
#' @export
sb <- function(x1 = NULL, by = NULL, h = -1) {
  args <- match.call()
  if (!is.null(args$x1) & is.null(args$by)) {
    cov <- c("ONE", deparse(args$x1, backtick = TRUE, width.cutoff = 500))
  } else if (!is.null(args$x1) & !is.null(args$by)) {
    cov <- c(deparse(args$by, backtick = TRUE, width.cutoff = 500),
            deparse(args$x1, backtick = TRUE, width.cutoff = 500))
  } else {
    stop("Invalid expression")
  }
  res <- list(cov = cov, h = h)
  res
}


