#' @name sback
#' @aliases sbackMain
#' @title Generalized additive and partially linear models.
#' @description Main function for fitting generalized additive and partially linear models by using smooth backfitting.
#' @param formula TODO
#' @param data TODO
#' @param weights TODO
#' @param kbin TODO
#' @param family TODO
#' @param newdata TODO
#' @export
sback <- function(formula, data, weights = NULL, kbin = 15, family = 2, newdata = NULL) {
		if(missing(formula)) {
		stop("Argument \"formula\" is missing, with no default")
	}
	if(missing(formula)) {
		stop("Argument \"data\" is missing, with no default")
	}
	if(!(family %in% 1:4)) {
		stop("Family not supported")
	}
	fsb <- interpret.sbformula(formula, method = "sback")
	if(is.null(fsb$response)) {
		stop("Response variable should be specified in argument \"formula\"")
	}
	varnames <- fsb$II[2,]
	if(any(is.na(match(c(fsb$response, varnames), names(data))))) {
		stop("Not all needed variables are supplied in data")
	}
	if(!is.null(newdata)) {
		if(any(is.na(match(varnames, names(data))))) {
			stop("Not all needed variables are supplied in newdata")
		}
	} else {
		newdata = data
	}
	data <- na.omit(data[, c(fsb$response, varnames)])
	newdata <- na.omit(newdata[, varnames, drop = FALSE])

	n <- nrow(data)
	n0 <- nrow(newdata)
	if(is.null(weights)) {
		weights <- rep(1, n)
	} else {
		if(sum(weights) <=0 || any(weights < 0) || length(weights) != n)
			stop("The specified weights are not correct")
	}
 	sback.fit  <- .Fortran("sbackMain",
               x      = matrix(as.double(as.matrix(data[, varnames])), ncol = length(varnames)),
               y      = as.double(data[, fsb$response]),
               w      = as.double(weights),
               n      = as.integer(n),
               npar   = as.integer(length(varnames)),
               kbin   = as.integer(kbin),
               h      = as.double(fsb$h),
               m      = matrix(as.double(rep(0.0, n * length(varnames))), nrow = n, ncol = length(varnames)),
               muhat  = as.double(rep(0.0, n)),
               family = as.double(family),
	  		       x0	    = matrix(as.double(as.matrix(newdata[, varnames])), ncol = length(varnames)),
	  		       m0     = matrix(as.double(rep(0.0, n0 * length(varnames))), nrow = n0, ncol = length(varnames)),
	  		       muhat0 = as.double(rep(0.0, n0)),
	   		       n0	    = as.integer(n0),
	   		       B      = as.double(rep(0.0, as.integer(length(varnames) + 1))), PACKAGE = "wsbackfit")

	effects <- sback.fit$m0
	colnames(effects) <- fsb$partial
	names(sback.fit$B) <- c("Intercept", varnames)

	res <- list(call = match.call(), data = data, pdata = newdata, effects = effects,
	            fitted.values = sback.fit$muhat0, h = sback.fit$h, fit = sback.fit, coeff = sback.fit$B)
  	class(res) <- "sback"
  	res
}
