#' @name plot.sback
#' @aliases plot.sback
#' @title Default sback plotting
#' @description Takes a fitted \code{sback} object produced by \code{sback()} and plot the estimated partial functions on the scale of the linear predictor
#' @param x an object of class \code{sback} as produced by \code{sback()}.
#' @param ask a logical value. If \code{TRUE}, the default, the user is asked for confirmation, before a new figure is drawn.
#' @param select Allows the plot for a single model term to be selected for printing. e.g. if you just want the plot for the second smooth term set select = 2.
#' @param \dots other graphics parameters to pass on to plotting commands.
#' @export
plot.sback <- function(x, ask = TRUE, select = NULL, ...) {
	p.functions <- colnames(x$effects)
	grph.opt <- names(list(...))
	if(is.null(select)) {
		ind <- 1:length(p.functions)
	} else if(!is.numeric(select) | !all(select %in% (1:length(p.functions)))) {
		stop("The model terms selected for printing do not exists")
	} else {
		ind <- select
	}
	j = 0
	for(i in ind) {
		j = j + 1
		if(j > 1 & length(ind) > 1) {
			if(ask) readline("Press return for next page....")
		}
		var <- interpret.sbformula(stats::as.formula(paste("~", p.functions[i])), method = "sback")$II[2, ]
		ord <- order(x$pdata[,var])
		x.data <- x$pdata[ord,var]
		y.data.nl <- x$effects[ord,i]
		y.data.l <- x$coeff[var] * x.data
		range <- max(c(y.data.nl, y.data.l)) - min(c(y.data.nl, y.data.l))
		min.ylim <- min(c(y.data.nl, y.data.l)) - 0.1*range
		max.ylim <- max(c(y.data.nl, y.data.l)) + 0.1*range
		stub <- paste(ifelse("xlab" %in% grph.opt, "", paste(",xlab = \"", var, "\"",sep = "")),
                ifelse("ylab" %in% grph.opt, "", paste(",ylab = \"", p.functions[i], "\"", sep = "")),
                ifelse("main" %in% grph.opt, "", paste(",main = \"Partial effect for ", var, "\"", sep = "")),
                ifelse("type" %in% grph.opt, "", ",type = \"l\""),
                ifelse("ylim" %in% grph.opt, "", paste(",ylim = c(", min.ylim,",", max.ylim,")", sep = "")), ",...)", sep = "")
		plot <- paste("plot(x.data, y.data.nl", stub, sep = "")
		eval(parse(text = plot))

		stub <- paste(ifelse("lty" %in% grph.opt, "", ",lty = 2"), ",...)", sep = "")
		lines <- paste("lines(x.data, y.data.l", stub, sep = "")
		eval(parse(text = lines))

		graphics::rug(x$data[,var])
		stub <- paste(ifelse("lty" %in% grph.opt, "", "lty = 1:2"), ", legend = c(\"Non linear\",\"Linear\"), bty = \"n\"", ",...)", sep = "")
		legend <- paste("legend(\"topleft\",", stub, sep = "")
		eval(parse(text = legend))
	}
}
