#' One Sample DP-TOST for Proportions
#'
#' @description
#' Internal function implementing the DP-TOST procedure from Section 3.1 of
#' Pareek et al. (2025) for one-sample proportion equivalence testing under
#' differential privacy.
#'
#' @param p_hat Privatized sample proportion (with noise already added)
#' @param n Sample size (public information)
#' @param lower Lower equivalence bound
#' @param upper Upper equivalence bound
#' @param epsilon Privacy budget. Default is 1
#' @param alpha Significance level. Default is 0.05
#' @param B Number of Monte Carlo replications. Default is 10000
#' @param max_resample Maximum resampling attempts for invalid solutions. Default is 100
#' @param seed Random seed for reproducibility. Default is 1337
#' @param ... Additional arguments (currently unused)
#'
#' @return A list containing:
#'   \item{conf.int}{(1-2α) confidence interval}
#'   \item{lower}{Lower equivalence bound}
#'   \item{upper}{Upper equivalence bound}
#'   \item{epsilon}{Privacy budget}
#'   \item{alpha}{Significance level}
#'   \item{decision}{Logical, TRUE if equivalence established}
#'   \item{p_hat}{Privatized proportion}
#'   \item{n}{Sample size}
#'   \item{seed}{Random seed used}
#'
#' @keywords internal
#' @noRd
prop_test_dp_one_sample <- function(p_hat, n, lower, upper,
                                     epsilon = 1, alpha = 0.05, B = 10^3,
                                     max_resample = 100, seed = 1337, ...) {

  # Check seed validity
  if (is.na(seed) || !is.numeric(seed) || length(seed) != 1) {
    stop("seed must be a single numeric value")
  }

  # Monte Carlo simulation
  nu_samples <- rep(NA, B)

  for (b in 1:B) {
    # Reconstruct p1
    p_check <- solve_proportion_matching(
      p_hat = p_hat,
      n = n,
      epsilon = epsilon,
      max_resample = max_resample,
      seed = seed + b
    )

    # Difference from reference
    nu_samples[b] <- p_check
  }

  # Confidence interval (percentile method)
  ci_lower <- quantile(nu_samples, probs = alpha, na.rm = TRUE)
  ci_upper <- quantile(nu_samples, probs = 1 - alpha, na.rm = TRUE)

  # Equivalence decision: CI must be entirely within [lower, upper]
  decision <- (ci_lower > lower) && (ci_upper < upper)

  # Return result
  result <- list(
    conf.int = c(ci_lower, ci_upper),
    lower = lower,
    upper = upper,
    epsilon = epsilon,
    alpha = alpha,
    decision = decision,
    p_hat = p_hat,
    n = n,
    seed = seed
  )

  return(result)
}


#' Solve Proportion Matching Problem
#'
#' @description
#' Internal function that solves the quadratic equation (5) from Pareek et al.
#' (2025) to reconstruct proportion values from privatized statistics via
#' moment matching.
#'
#' @param p_hat Privatized proportion estimate
#' @param n Sample size
#' @param epsilon Privacy budget
#' @param max_resample Maximum number of resampling attempts
#' @param seed Random seed. Default is 1337
#'
#' @return Reconstructed proportion value in [0,1], or NA if no valid solution found
#'
#' @keywords internal
#' @noRd
solve_proportion_matching <- function(p_hat, n, epsilon, max_resample, seed = 1337) {

  # Check seed validity
  if (is.na(seed) || !is.numeric(seed) || length(seed) != 1) {
    stop("seed must be a single numeric value")
  }

  # Set seed
  set.seed(seed)

  valid_solution <- FALSE
  attempts <- 0

  while (!valid_solution && attempts <= max_resample) {
    attempts <- attempts + 1

    # Generate noise and normal random variable
    Z_star <- rnorm(1)
    U_star <- rlaplace_custom(1, scale = 1 / (n * epsilon))

    # Compute quadratic coefficients
    delta_val <- Z_star / sqrt(n)
    gamma_val <- delta_val^2

    Lambda <- -4 * p_hat^2 + 4 * p_hat + gamma_val +
              8 * p_hat * U_star - 4 * U_star^2 - 4 * U_star

    # Check if solutions exist
    if (Lambda < 0) next

    # Two candidate roots (equation 5)
    numerator1 <- -delta_val * sqrt(Lambda) + (2 * p_hat + gamma_val - 2 * U_star)
    numerator2 <- delta_val * sqrt(Lambda) + (2 * p_hat + gamma_val - 2 * U_star)
    denominator <- 2 * (gamma_val + 1)

    p_check_1 <- numerator1 / denominator
    p_check_2 <- numerator2 / denominator

    # Check validity: must be in [0, 1]
    candidates <- c()
    losses <- c()

    if (!is.na(p_check_1) && p_check_1 >= 0 && p_check_1 <= 1) {
      candidates <- c(candidates, p_check_1)
      loss1 <- abs(p_hat - p_check_1 - sqrt(p_check_1 * (1 - p_check_1) / n) * Z_star - U_star)
      losses <- c(losses, loss1)
    }

    if (!is.na(p_check_2) && p_check_2 >= 0 && p_check_2 <= 1) {
      candidates <- c(candidates, p_check_2)
      loss2 <- abs(p_hat - p_check_2 - sqrt(p_check_2 * (1 - p_check_2) / n) * Z_star - U_star)
      losses <- c(losses, loss2)
    }

    # Select solution with smallest loss
    if (length(candidates) > 0) {
      p_check <- candidates[which.min(losses)]
      valid_solution <- TRUE
      return(p_check)
    }
  }

  # No valid solution found
  return(NA)
}


#' Laplace Random Number Generator
#'
#' @description
#' Internal function generating random samples from the Laplace distribution
#' using the inverse CDF method with uniform random variables.
#'
#' @param n Number of samples to generate
#' @param location Location parameter (μ). Default is 0
#' @param scale Scale parameter (b > 0). Default is 1
#'
#' @return Numeric vector of length n with samples from Laplace(location, scale)
#'
#' @details
#' Uses the transformation: X = μ - b × sign(U - 0.5) × log(1 - 2|U - 0.5|)
#' where U ~ Uniform(0,1)
#'
#' @keywords internal
#' @noRd
rlaplace_custom <- function(n, location = 0, scale = 1) {

  u <- runif(n)

  # Inverse CDF: X = location - scale * sign(U - 0.5) * log(1 - 2|U - 0.5|)
  centered <- u - 0.5
  samples <- location - scale * sign(centered) * log(1 - 2 * abs(centered))

  return(samples)
}
