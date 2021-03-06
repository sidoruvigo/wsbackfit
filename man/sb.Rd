% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sb.R
\name{sb}
\alias{sb}
\title{Specify a nonparametric and/or a varying coefficient term in a wsbackfit formula}
\usage{
sb(x1 = NULL, by = NULL, h = -1)
}
\arguments{
\item{x1}{the univariate predictor.}

\item{by}{numeric predictor of the same dimension as x1. If present, the coefficients of this predictor depend, nonparametrically, on x1, i.e., a varying coefficient term. Only applied for the \code{vcoef()} function.}

\item{h}{bandwidth for this term. If h = -1, the bandwidth is automatically selected using cross-validation. For the \code{sback()} fuction, h = 0 indicates a parametric/linear fit. Defaults to -1.}
}
\description{
Function used to indicate a nonparametric term in a \code{sback} formula and a varying coefficient term in a \code{vcoef} formula.
}
