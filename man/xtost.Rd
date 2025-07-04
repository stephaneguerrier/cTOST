% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internal.R
\name{xtost}
\alias{xtost}
\title{TO BE DOCUMENTED}
\usage{
xtost(
  theta_hat,
  sig_hat,
  nu,
  alpha,
  delta,
  correction = "no",
  B = 10^4,
  seed = 85
)
}
\arguments{
\item{theta_hat}{The estimated mean.}

\item{sig_hat}{The estimated standard deviation.}

\item{nu}{The degrees of freedom parameter.}

\item{alpha}{The significance level for the test.}

\item{delta}{The equivalence bound used for the TOST decision.}

\item{correction}{TO BE DOCUMENTED}

\item{B}{TO BE DOCUMENTED (numb of bootstrap for correction based on bootstrap)}

\item{seed}{TO BE DOCUMENTED (also for bootstrap)}
}
\value{
TO BE DOCUMENTED
}
\description{
TO BE DOCUMENTED
}
\keyword{internal}
