#' @name print.vcoef
#' @aliases print.vcoef
#' @title Print a \code{vcoef} object.
#' @description The default print method for a \code{vcoef} object.
#' @param x an object of class \code{vcoef} as produced by \code{vcoef()}.
#' @param \dots other arguments, ignored.
#' @export
print.vcoef <- function(x, ...) {
  foo0 <- cbind(colnames(x$effects), round(x$h, 4))
  colnames(foo0) <- c("Effect", "h")
  rownames(foo0) <- rep("", dim(foo0)[1])
  cat("Generalized Varying Coefficent Model / wsbackfit:\n\n")
  cat("Call: ")
  print(x$call)
  cat("\nSample size:", length(x$fitted.values), "\n\nBandwidths used in model:\n")
  print(foo0, quote = FALSE)
  cat("\nLinear components:\n")
  print(x$coeff, quote = FALSE)
  cat("\n")
}
