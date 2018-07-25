#' @name plot.vcoef
#' @aliases plot.vcoef
#' @title Default vcoef plotting
#' @description Takes a fitted \code{vcoef} object produced by \code{vcoef()} and plot the estimated partial functions on the scale of the linear predictor.
#' @param x an object of class \code{vcoef} as produced by \code{vcoef()}.
#' @param ask a logical value. If \code{TRUE}, the default, the user is asked for confirmation, before a new figure is drawn.
#' @param select Allows the plot for a single model term to be selected for printing. e.g. if you just want the plot for the second smooth term set select=2.
#' @param \dots other graphics parameters to pass on to plotting commands.
#' @export
plot.vcoef <- function(x, ask = TRUE, select = NULL, ...) {
  dots <- list(...)
	grph.opt <- names(dots)
	p.functions <- colnames(x$effects)[x$h != 0]
	if (is.null(select)) {
	  ind <- 1:length(p.functions)
	} else if (!is.numeric(select) | !all(select %in% (1:length(p.functions)))) {
	  stop("The model terms selected for printing do not exists")
	} else {
	  ind <- select
	}
	j <- 0
	for(i in ind) {
		j <- j + 1
		if(j > 1 & length(ind) > 1) {
			if(ask) readline("Press return for next page....")
		}

		# Partial effect
		int.f <- interpret.sbformula(stats::as.formula(paste("~", p.functions[i])), method = "vcoef")
		var <- int.f$II[2,]
		ord <- order(x$pdata[, var])
		x.data <- x$pdata[ord, var]
		y.data <- x$effects[ord, i]
		stub <- paste(ifelse("xlab" %in% grph.opt, "", paste(",xlab = \"", var, "\"",sep = "")),
                ifelse("ylab" %in% grph.opt, "", paste(",ylab = \"", p.functions[i], "\"", sep = "")),
                ifelse("main" %in% grph.opt, "", paste(",main = \"Partial effect for ", var, "\"", sep = "")),
                ifelse("type" %in% grph.opt, "", ",type = \"l\""), ",...)", sep = "")
		plot <- paste("plot(x.data, y.data", stub, sep = "")
		eval(parse(text = plot))
		graphics::rug(x$data[, var])

		# Interaction surface
		if(int.f$II[1,] != "ONE") {
			x.data <- seq(min(x$pdata[, var]), max(x$pdata[, var]), l = 50)
			y.data <- seq(min(x$pdata[, int.f$II[1,]]), max(x$pdata[, int.f$II[1, ]]), l = 50)
			z.data <- outer(stats::approxfun(x$pdata[ord,var], x$effects[ord, i])(x.data), y.data, '*')
			if(ask) readline("Press return for next page....")
			graphics::persp(x.data, y.data, z.data, xlab = var, ylab = int.f$II[1, ], zlab = p.functions[i],
				main = paste("Estimated surface for ", var, " and ", int.f$II[1, ], sep = ""),
				theta = ifelse(is.null(dots$theta), 45, dots$theta),
				phi = ifelse(is.null(dots$phi), 45, dots$phi),
				ticktype = "detailed",
				shade = ifelse(is.null(dots$shade), 0.5, dots$shade),
				cex.axis = dots$cex.axis, cex.lab = dots$cex.lab, cex.sub = dots$cex.sub, cex = dots$cex)
		}
	}
}
