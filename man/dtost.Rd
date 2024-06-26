% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/equiv.R
\name{dtost}
\alias{dtost}
\title{The delta-TOST Corrective Procedure for (Bio)Equivalence Testing}
\usage{
dtost(theta, sigma, nu, alpha, delta)
}
\arguments{
\item{theta}{A \code{numeric} value corresponding to the difference of means.}

\item{sigma}{A \code{numeric} value corresponding to the standard error.}

\item{nu}{A \code{numeric} value corresponding to the number of degrees of freedom.}

\item{alpha}{A \code{numeric} value specifying the significance level.}

\item{delta}{A \code{numeric} value corresponding to (bio)equivalence limit. We assume symmetry, i.e, the (bio)equivalence interval corresponds to (-delta,delta)}
}
\value{
A \code{tost} object with the structure:
\itemize{
 \item decision:    A boolean variable indicating whether (bio)equivalence is accepted or not.
 \item ci:          Confidence interval at the 1 - 2*alpha level.
 \item theta:       The difference of means used for the test.
 \item sigma:       The standard error used for the test.
 \item nu:          The number of degrees of freedom used for the test.
 \item alpha:       The significance level used for the test.
 \item delta:       The (bio)equivalence limits used for the test.
 \item method:      The method used for the test (here the "delta-TOST").
}
}
\description{
This functions is used to compute the delta-TOST, a corrective procedure of the (bio)equivalence bounds applied to the Two One-Sided Test (TOST) for (bio)equivalence testing in the univariate framework.
}
\examples{
data(skin)

theta_hat = diff(apply(skin,2,mean))
nu = nrow(skin) - 1
sig_hat = sd(apply(skin,1,diff))/sqrt(nu)
res_dtost = cTOST:::dtost(theta = theta_hat, sigma = sig_hat, nu = nu,
              alpha = 0.05, delta = log(1.25))
compare_to_tost(res_dtost)

}
\author{
Younes Boulaguiem, Stéphane Guerrier, Dominique-Laurent Couturier
}
\keyword{internal}
