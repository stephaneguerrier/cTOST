% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internal_average_univariate.R
\name{size_TOST}
\alias{size_TOST}
\title{Size of Univariate TOST}
\usage{
size_TOST(alpha, sigma, nu, delta, ...)
}
\arguments{
\item{alpha}{A \code{numeric} value specifying the significance level.}

\item{sigma}{A \code{numeric} value representing the estimated standard error of \code{theta}.}

\item{nu}{A \code{numeric} value specifying the degrees of freedom.}

\item{delta}{A \code{numeric} value defining the (bio)equivalence margin. The procedure assumes symmetry, i.e., the (bio)equivalence region is \eqn{(-\delta, \delta)}.}

\item{...}{Additional parameters.}
}
\value{
A \code{numeric} value corresponding to the probability (size) of the TOST procedure.
}
\description{
Computes the size (type I error rate) of the univariate Two One-Sided Tests (TOST) procedure.
}
\keyword{internal}
