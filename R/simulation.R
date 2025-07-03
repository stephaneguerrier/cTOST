#' Simulate data from the canonical form of the univariate average equivalence problem
#'
#' @param mu Numeric value indicating the true mean of the data.
#' @param sigma Numeric value indicating the variance of the data.
#' @param nu Numeric value indicating the degrees of freedom parameter of the considered setting.
#' @param seed Numeric value (optional) used to set the seed for reproducibility of random number generation. Default value is 18.
#'
#' @return A list containing the estimated mean and estimated standard deviation of the simulated data.
#'
#' @examples
#' simulate_data(mu = 10, sigma = 2, nu = 10)
#'
#' @importFrom stats rnorm rchisq
simulate_data = function(mu, sigma, nu, seed = 18){
  set.seed(seed)
  theta_hat = rnorm(1, mean = mu, sd = sigma)
  sig_hat = sqrt(rchisq(1, df = nu)/nu*sigma)
  list(theta_hat = theta_hat, sig_hat = sig_hat)
}
