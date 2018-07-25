#' @name print.sback
#' @aliases print.sback
#' @title Print a \code{sback} object.
#' @description The default print method for a \code{sback} object.
#' @param x an object of class \code{sback} as produced by \code{sback()}.
#' @param \dots other arguments, ignored.
#' @export
print.sback <- function(x, ...) {
  foo0 <- cbind(colnames(x$effects), round(x$h, 4))
  colnames(foo0) <- c("Effect", "h")
  rownames(foo0) <- rep("", dim(foo0)[1])
  cat("Generalized Smooth Backfitting / wsbackfit:\n\n")
  cat("Call: ")
  print(x$call)
  cat("\nSample size:", length(x$fitted.values), "\n\nBandwidths used in model:\n")
  print(foo0, quote = FALSE)
  cat("\nLinear components:\n")
  print(x$coeff, quote = FALSE)
  cat("\n")
}
