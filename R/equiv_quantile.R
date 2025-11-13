#' Core quantile equivalence testing procedure
#'
#' @description
#' This is the core the qTOST procedure to obtain CI, etc to test both a single
#' or multiple quantiles. It takes the test statistic (`theta`), its standard
#' error (`sigma`), and other parameters to perform the two one-sided tests for
#' equivalence. This is a low-level function; users should typically use the
#' main `qtost` wrapper function.
#'
#' @param theta The calculated test statistic, $\theta$.
#' @param sigma The standard error of the test statistic, $\sigma_\theta$.
#' @param pi_x A numeric scalar or vector representing the quantile(s) of interest
#'   in the reference population ($X$).
#' @param delta_l A numeric scalar or vector for the lower equivalence margin(s).
#' @param delta_u A numeric scalar or vector for the upper equivalence margin(s).
#' @param alpha The nominal significance level for the test (e.g., 0.05).
#' @param corrected_alpha An optional corrected significance level for `alpha-qTOST`.
#'
#' @return
#' An object of class `qtost` or `mqtost` containing the results of the test.
#' This list includes:
#' \itemize{
#'   \item `decision`: A boolean value (`TRUE` for equivalence, `FALSE` otherwise).
#'   \item `ci`: The confidence interval for the estimated quantile in population Y, $\hat{\pi}_y$.
#'   \item `pi_y_hat`: The point estimate of the quantile in population Y, calculated as $\Phi(\theta)$.
#'   \item `theta`: The value of the test statistic $\theta$.
#'   \item `sigma`: The standard error of $\theta$.
#'   \item `ci_theta`: The confidence interval for $\theta$.
#'   \item `alpha`: The significance level used for the test.
#'   \item `pi_x`: The reference quantile(s).
#'   \item `delta`: The original equivalence margin.
#'   \item `eq_region`: The equivalence region defined by `delta_l` and `delta_u`.
#'   \item `method`: The name of the method used ("qTOST" or "alpha-qTOST").
#'   \item `setting`: The setting used ("single" or "multiple" quantiles).
#' }
#'
#' @keywords internal
qtost_core <- function(theta, sigma, pi_x, delta_l, delta_u, alpha = 0.05, corrected_alpha = NULL) {
  if (!is.null(corrected_alpha)) {
    alpha0 = alpha
    alpha = corrected_alpha
  }
  qnt <- qnorm(1 - alpha)
  up_theta <- theta + qnt * sigma
  dw_theta <- theta - qnt * sigma
  upper <- pnorm(up_theta)
  lower <- pnorm(dw_theta)
  decision <- upper < delta_u & lower > delta_l
  out <- list(
    decision = decision,
    ci = cbind(lower, upper),
    pi_y_hat = pnorm(theta),
    theta = theta,
    sigma = sigma,
    ci_theta = cbind(lower = dw_theta, upper = up_theta),
    alpha = alpha,
    corrected_alpha = corrected_alpha,
    pi_x = pi_x,
    delta = delta_u[1]-pi_x[1],
    eq_region = cbind(delta_l, pi_x, delta_u),
    method = "qTOST",
    setting = if (length(pi_x) > 1) "multiple" else "single"
  )
  if (is.null(corrected_alpha)) {
    out = out[out != "corrected_alpha"]
  } else {
    out$alpha = alpha0
    out$method = "alpha-qTOST"
  }
  class(out) <- if (length(pi_x) > 1) "mqtost" else "qtost"
  return(out)
}

#' @title Generate distribution of estimated parameter of interest
#'
#' @description This function is used to generate estimated parameter of interest and its standard deviation
#' using Monte Carlo simulation.
#'
#' @param theta   A \code{numeric} value corresponding to the parameter of interest.
#' @param gamma   A \code{numeric} value corresponding to the ratio of sample variance between the target and reference group.
#' @param n_x     An \code{integer} value representing the sample size of the reference group.
#' @param n_y     An \code{integer} value representing the sample size of the target group.
#' @param pi_x    A \code{numeric} value defining the quantile of interest.
#' @param B       An optional \code{integer} value specifying the number of Monte Carlo replications (default: \code{B = 10^5}).
#' @param seed    An optional \code{integer} value specifying the random seed for reproducibility (default: \code{seed = 12345}).
#'
#' @keywords internal
#'
#' @return A list with the structure:
#' \itemize{
#'  \item theta_hat:         B \code{numeric} values corresponding to the \eqn{\hat{\theta}} used in the test.
#'  \item sigma_theta_hat:   B \code{numeric} values corresponding to the \eqn{\hat{\sigma}} of \eqn{\hat{\theta}} used in the test.
#' }
#'
get_qTOST_rvs = function(theta, gamma, n_x, n_y, pi_x, B = 10^4, seed = 12345) {
  nu_x = n_x - 1
  nu_y = n_y - 1
  set.seed(seed)
  W1 = sqrt(rchisq(n = B, df = nu_x))
  W2 = sqrt(rchisq(n = B, df = nu_y))
  Z = rnorm(B)
  sgamma = sqrt(gamma)
  snu_y = sqrt(nu_y)
  snu_x = sqrt(nu_x)
  q_x = qnorm(pi_x)
  b1 = snu_y / sgamma * (theta * sgamma - q_x)
  b2 = snu_y / sgamma * sqrt(1 / n_x + gamma / n_y)
  b3 = snu_y * q_x / (snu_x * sgamma)
  theta_hat = b1 / W2 + b2 * (Z / W2) + b3 * W1 / W2
  gamma_hat = W2^2 / W1^2 * gamma * nu_x / nu_y
  sigma_theta_hat = sqrt((1 / n_y) * (1 + theta_hat^2 / 2) +
                           (1 / n_x) / gamma_hat * (1 + qnorm(pi_x)^2 / 2))
  out = list(theta_hat = theta_hat, sigma_theta_hat = sigma_theta_hat)
}

#' Power calculation for quantile TOST using Monte Carlo simulation
#'
#' @description
#' This function estimates the statistical power of the quantile-based
#' Two One-Sided Tests (qTOST) procedure using Monte Carlo simulations.
#'
#' @param theta     A \code{numeric} value corresponding to the parameter of interest.
#' @param gamma     A \code{numeric} value corresponding to the ratio of sample variance between the target and reference group.
#' @param n_x       An \code{integer} value representing the sample size of the reference group.
#' @param n_y       An \code{integer} value representing the sample size of the target group.
#' @param pi_x      A \code{numeric} value defining the quantile of interest.
#' @param delta_l   A \code{numeric} value defining the equivalence lower limit.
#' @param delta_u   A \code{numeric} value defining the equivalence upper limit.
#' @param alpha     A \code{numeric} value specifying the significance level.
#' @param B         An \code{integer} value specifying the number of Monte Carlo replications.
#' @param seed      An \code{integer} value specifying the random seed for reproducibility.
#' @param ...       Additional parameters
#'
#' @details
#' The estimated power is the proportion of Monte Carlo replicates in which equivalence is concluded.
#'
#' @return The function returns a \code{numeric} value that corresponds to a probability.
#'
#' @keywords internal
#'
power_qTOST_MC = function(theta, gamma = gamma, n_x, n_y, pi_x, delta_l, delta_u, alpha, B, seed, ...) {
  sol = get_qTOST_rvs(theta = theta, gamma = gamma, n_x = n_x, n_y = n_y, pi_x = pi_x, B = B, seed = seed)
  theta_hat = sol$theta_hat
  sigma_theta_hat = sol$sigma_theta_hat
  up_theta = theta_hat + qnorm(1 - alpha) * sigma_theta_hat
  dw_theta = theta_hat - qnorm(1 - alpha) * sigma_theta_hat
  tmp = up_theta < qnorm(delta_u) & dw_theta > qnorm(delta_l)
  res = sum(tmp) / B
  res
}

#' @title Objective function to optimize for quantile TOST
#'
#' @param test      A \code{numeric} value corresponding to the significance level to optimize.
#' @param theta     A \code{numeric} value corresponding to the parameter of interest.
#' @param gamma     A \code{numeric} value corresponding to the ratio of sample variance between the target and reference group.
#' @param n_x       An \code{integer} value representing the sample size of the reference group.
#' @param n_y       An \code{integer} value representing the sample size of the target group.
#' @param pi_x      A \code{numeric} value defining the quantile of interest.
#' @param delta_l   A \code{numeric} value defining the equivalence lower limit.
#' @param delta_u   A \code{numeric} value defining the equivalence upper limit.
#' @param alpha     A \code{numeric} value specifying the significance level.
#' @param B         An \code{integer} value specifying the number of Monte Carlo replications.
#' @param seed      An \code{integer} value specifying the random seed for reproducibility.
#'
#' @keywords internal
#'
#' @return The function returns a \code{numeric} value for the objective function.
#'
obj_fun_qTOST_MC = function(test, theta, gamma, n_x, n_y, pi_x, delta_l, delta_u, alpha, B, seed) {
  size = power_qTOST_MC(theta, gamma, n_x, n_y, pi_x, delta_l, delta_u, alpha=test, B, seed)
  10^5 * (size - alpha)
}


#' @title Get alpha star for quantile TOST
#'
#' @param gamma     A \code{numeric} value corresponding to the ratio of sample variance between the target and reference group.
#' @param n_x       An \code{integer} value representing the sample size of the reference group.
#' @param n_y       An \code{integer} value representing the sample size of the target group.
#' @param pi_x      A \code{numeric} value defining the quantile of interest.
#' @param delta_l   A \code{numeric} value defining the equivalence lower limit.
#' @param delta_u   A \code{numeric} value defining the equivalence upper limit.
#' @param alpha     A \code{numeric} value specifying the significance level.
#' @param B         An \code{integer} value specifying the number of Monte Carlo replications.
#' @param seed      An \code{integer} value specifying the random seed for reproducibility.
#' @param tol       An optional \code{numeric} value defining the tolerance for root-finding (default: \code{tol = .Machine$double.eps^0.5}).
#' @param ...       Additional parameters
#'
#' @keywords internal
#'
#' @return A list with at least two components:
#' \itemize{
#'  \item \code{root}:       A \code{numeric} value corresponding to the location of the root.
#'  \item \code{f.root}:     A \code{numeric} value corresponding to the function evaluated at \code{root}.
#' }
#'
get_alpha_qTOST_MC = function(gamma, n_x, n_y, pi_x, delta_l, delta_u, alpha, B, seed, tol = .Machine$double.eps, ...) {
  # check sup
  thetas = c(qnorm(delta_l), qnorm(delta_u))
  pw_l = power_qTOST_MC(thetas[1], gamma, n_x, n_y, pi_x, delta_l, delta_u, alpha, B, seed)
  pw_u = power_qTOST_MC(thetas[2], gamma, n_x, n_y, pi_x, delta_l, delta_u, alpha, B, seed)
  lambda = thetas[which.max(c(pw_l, pw_u))]
  # compute adjustment
  obj_func_a1 = obj_fun_qTOST_MC(alpha, lambda, gamma, n_x, n_y, pi_x, delta_l, delta_u, alpha, B, seed)
  if (obj_func_a1 > tol) return(list(root = alpha, f.root = obj_func_a1)) # size close to alpha
  obj_func_a2 = obj_fun_qTOST_MC(1/2, lambda, gamma, n_x, n_y, pi_x, delta_l, delta_u, alpha, B, seed)
  if (abs(obj_func_a2) < tol) return(list(root = 1/2, f.root = obj_func_a2)) # size close to 0.5
  out = uniroot(obj_fun_qTOST_MC, interval = c(alpha, 1/2),
                lambda, gamma, n_x, n_y, pi_x, delta_l, delta_u,
                alpha, B, seed)
  out
}
