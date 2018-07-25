#' @name vcoef
#' @aliases vcoefMain
#' @title Generalized varying coefficient models.
#' @description Main function for fitting generalized varying coefficient models.
#' @param formula TODO
#' @param data TODO
#' @param weights TODO
#' @param kbin TODO
#' @param family TODO
#' @param newdata TODO
#' @export
vcoef <- function(formula, data, weights = NULL, kbin = 15, family = 2, newdata = NULL) {
  if (missing(formula)) {
    stop("Argument \"formula\" is missing, with no default")
  }
  if (missing(formula)) {
    stop("Argument \"data\" is missing, with no default")
  }
  if (!(family %in% 1:4)) {
    stop("Family not supported")
	}
  data[, "ONE"] <- 1
  fsb <- interpret.sbformula(formula, method = "vcoef")
  if (is.null(fsb$response)) {
    stop("Response variable should be specified in argument \"formula\"")
  }
  z.varnames <- fsb$II[1, ]
  x.varnames <- fsb$II[2, ]
  if (anyNA(match(c(fsb$response, c(x.varnames, z.varnames)), names(data)))) {
    stop("Not all needed variables are supplied in data")
  }
  if (!is.null(newdata)) {
    newdata[, "ONE"] <- 1
    if (anyNA(match(c(x.varnames, z.varnames), names(data)))) {
      stop("Not all needed variables are supplied in newdata")
    }
  } else {
    newdata <- data
  }
	data <- stats::na.omit(data[, c(fsb$response, unique(c(x.varnames, z.varnames)))])
	newdata <- stats::na.omit(newdata[, unique(c(x.varnames, z.varnames))])

	n <- nrow(data)
	n0 <- nrow(newdata)
	# ifelse(is.null(weights), weights <- rep(1, n), if(sum(weights) <=0 || any(weights < 0) || length(weights) != n)	stop("The specified weights are not correct"))
	if (is.null(weights)) {
	  weights <- rep(1, n)
	} else {
	  if (sum(weights) <= 0 || any(weights < 0) || length(weights) != n)
	    stop("The specified weights are not correct")
	}
	ind.z.l <- z.varnames != "ONE"
	ind.l <- fsb$h == 0
	if(any(ind.z.l) | any (ind.l)) {
	  z.varnames.l <- unique(c(z.varnames[ind.z.l], x.varnames[ind.l]))
		zl <- matrix(as.double(as.matrix(data[, z.varnames.l])), ncol = length(z.varnames.l))
		z0l <- matrix(as.double(as.matrix(newdata[, z.varnames.l])), ncol = length(z.varnames.l))
	} else {
	  zl <- as.double(0)
	  z0l <- as.double(0)
	  z.varnames.l <- NULL
	}

	vcoef.fit  <- .Fortran("vcoefMain",
		x 	   = matrix(as.double(as.matrix(data[, x.varnames])), ncol = length(x.varnames)),
		z 	   = matrix(as.double(as.matrix(data[, z.varnames])), ncol = length(z.varnames)),
		y 	   = as.double(data[, fsb$response]),
		w 	   = as.double(weights),
		n 	   = as.integer(n),
		npar   = as.integer(length(x.varnames)),
		zl 	   = zl,
		nparl  = as.integer(length(z.varnames.l)),
		kbin   = as.integer(kbin),
		h 	   = as.double(fsb$h),
		m      = matrix(as.double(rep(0.0, n*length(z.varnames))), nrow = n, ncol = length(z.varnames)),
		mx	   = matrix(as.double(rep(0.0, n*length(x.varnames))), nrow = n, ncol = length(x.varnames)),
		muhat  = as.double(rep(0.0, n)),
		family = as.double(family),
		x0 	   = matrix(as.double(as.matrix(newdata[, x.varnames])), ncol = length(x.varnames)),
		z0 	   = matrix(as.double(as.matrix(newdata[, z.varnames])), ncol = length(z.varnames)),
		z0l    = z0l,
		m0 	   = matrix(as.double(rep(0.0, n0*length(z.varnames))), nrow = n0, ncol = length(z.varnames)),
		mx0    = matrix(as.double(rep(0.0, n0*length(x.varnames))), nrow = n0, ncol = length(x.varnames)),
		muhat0 = as.double(rep(0.0, n0)),
		n0     = as.integer(n0),
		B      = as.double(rep(0.0, as.integer(length(z.varnames.l) + 1))), PACKAGE = "wsbackfit")

	effects <- vcoef.fit$mx0
	colnames(effects) <- fsb$partial
	names(vcoef.fit$B) <- c("Intercept", z.varnames.l)

	res <- list(call = match.call(), data = data, pdata = newdata, effects = effects,
	            fitted.values = vcoef.fit$muhat0, h = vcoef.fit$h, fit = vcoef.fit,
	            coeff = vcoef.fit$B)
  	class(res) <- "vcoef"
  	res
}


