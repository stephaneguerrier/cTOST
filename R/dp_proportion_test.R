#' Differentially Private TOST for Proportion Equivalence Testing
#'
#' @description
#' Performs equivalence testing for proportions under differential privacy using
#' the DP-TOST procedure. Supports both one-sample and two-sample tests.
#'
#' @details
#' This function implements the DP-TOST (Differentially Private Two One-Sided Tests)
#' procedure for testing equivalence of proportions while maintaining differential
#' privacy guarantees. The test uses Monte Carlo simulation to reconstruct the
#' sampling distribution from privatized statistics and constructs confidence
#' intervals using the percentile method.
#'
#' For the one-sample case, the test evaluates whether a proportion is equivalent
#' to a reference range [lower, upper]. For the two-sample case, it tests whether
#' the difference between two proportions (p1 - p2) falls within the equivalence
#' bounds [lower, upper].
#'
#' The privacy mechanism adds Laplace noise to the observed proportions. The privacy
#' budget epsilon controls the amount of noise: larger epsilon values provide less
#' privacy but more accurate inference.
#'
#' @param p_hat Privatized sample proportion for the first (or only) group. This
#'   should be the observed proportion with differential privacy noise already added.
#' @param n Sample size for the first (or only) group.
#' @param p_hat2 Privatized sample proportion for the second group (optional).
#'   If NULL, a one-sample test is performed. Default is NULL.
#' @param n2 Sample size for the second group (optional). Required if p_hat2 is
#'   provided. Default is NULL.
#' @param lower Lower equivalence bound. For one-sample tests, this is the lower
#'   bound for the proportion. For two-sample tests, this is the lower bound for
#'   the difference (p1 - p2). Typical value: -delta.
#' @param upper Upper equivalence bound. For one-sample tests, this is the upper
#'   bound for the proportion. For two-sample tests, this is the upper bound for
#'   the difference (p1 - p2). Typical value: delta.
#' @param epsilon Privacy budget. Can be a single value (used for all samples) or
#'   a vector of length 2 (epsilon[1] for first group, epsilon[2] for second group).
#'   Larger values provide less privacy but more statistical power. Default is 1.
#' @param alpha Significance level for the test. The confidence interval will have
#'   level (1 - 2*alpha). Default is 0.05.
#' @param B Number of Monte Carlo replications for reconstructing the sampling
#'   distribution. Larger values provide more accurate results but take longer.
#'   Default is 10000.
#' @param max_resample Maximum number of resampling attempts when solving the
#'   moment matching equation. Default is 100.
#' @param seed Random seed for reproducibility. Default is 1337.
#' @param ... Additional arguments (currently unused).
#'
#' @return A list of class "dp_tost_prop" containing:
#'   \item{decision}{Logical. TRUE if equivalence is established, FALSE otherwise.}
#'   \item{conf.int}{Confidence interval for the parameter of interest (proportion
#'     for one-sample, difference for two-sample).}
#'   \item{lower}{Lower equivalence bound.}
#'   \item{upper}{Upper equivalence bound.}
#'   \item{alpha}{Significance level used.}
#'   \item{epsilon}{Privacy budget(s) used.}
#'   \item{estimate}{Point estimate(s) of the privatized proportion(s).}
#'   \item{sample_size}{Sample size(s).}
#'   \item{test_type}{Character string indicating "one-sample" or "two-sample".}
#'   \item{seed}{Random seed used.}
#'   \item{B}{Number of Monte Carlo replications.}
#'
#' @examples
#' \dontrun{
#' # One-sample test: Is proportion equivalent to range [0.3, 0.5]?
#' set.seed(123)
#' n <- 500
#' p_true <- 0.42
#' epsilon <- 1
#'
#' # Generate Bernoulli data
#' x <- rbinom(n, 1, p_true)
#' p_obs <- mean(x)
#'
#' # Add Laplace noise using inverse CDF method
#' u <- runif(1)
#' scale <- 1 / (n * epsilon)
#' laplace_noise <- -scale * sign(u - 0.5) * log(1 - 2 * abs(u - 0.5))
#' p_hat <- p_obs + laplace_noise
#'
#' result <- prop_test_equiv_dp(
#'   p_hat = p_hat,
#'   n = n,
#'   lower = 0.3,
#'   upper = 0.5,
#'   epsilon = 1,
#'   B = 1000
#' )
#' result
#'
#' # Two-sample test: Is difference (p1-p2) equivalent to [-0.1, 0.1]?
#' set.seed(456)
#' n1 <- 300
#' n2 <- 300
#' p1_true <- 0.52
#' p2_true <- 0.48
#' epsilon1 <- 1
#' epsilon2 <- 1
#'
#' # Generate Bernoulli data for both groups
#' x1 <- rbinom(n1, 1, p1_true)
#' x2 <- rbinom(n2, 1, p2_true)
#' p1_obs <- mean(x1)
#' p2_obs <- mean(x2)
#'
#' # Add Laplace noise to each proportion
#' u1 <- runif(1)
#' u2 <- runif(1)
#' scale1 <- 1 / (n1 * epsilon1)
#' scale2 <- 1 / (n2 * epsilon2)
#' p1_hat <- p1_obs - scale1 * sign(u1 - 0.5) * log(1 - 2 * abs(u1 - 0.5))
#' p2_hat <- p2_obs - scale2 * sign(u2 - 0.5) * log(1 - 2 * abs(u2 - 0.5))
#'
#' result <- prop_test_equiv_dp(
#'   p_hat = p1_hat,
#'   n = n1,
#'   p_hat2 = p2_hat,
#'   n2 = n2,
#'   lower = -0.1,
#'   upper = 0.1,
#'   epsilon = c(1, 1),
#'   B = 1000
#' )
#' result
#' }
#'
#' @references
#' Pareek, K., Shang, H. L., & Guerrier, S. (2025). Differentially Private
#' Quantile-based TOST for Equivalence Testing. Manuscript in preparation.
#'
#' @export
prop_test_equiv_dp <- function(p_hat, n, p_hat2 = NULL, n2 = NULL,
                               lower, upper, epsilon = 1, alpha = 0.05,
                               B = 10^4, max_resample = 100, seed = 1337, ...) {

  # Determine test type
  two_sample <- !is.null(p_hat2)

  # Validate inputs
  if (two_sample && is.null(n2)) {
    stop("n2 must be provided when p_hat2 is specified (two-sample test)")
  }

  if (!two_sample && !is.null(n2)) {
    warning("n2 is ignored for one-sample test")
  }

  # Handle epsilon: convert to vector of length 2 if needed
  if (length(epsilon) == 1) {
    epsilon <- rep(epsilon, 2)
  } else if (length(epsilon) != 2) {
    stop("epsilon must be a single value or a vector of length 2")
  }

  # Perform the appropriate test
  if (two_sample) {
    # Two-sample test
    result <- prop_test_dp_two_sample(
      p1_hat = p_hat,
      p2_hat = p_hat2,
      n1 = n,
      n2 = n2,
      lower = lower,
      upper = upper,
      epsilon1 = epsilon[1],
      epsilon2 = epsilon[2],
      alpha = alpha,
      B = B,
      max_resample = max_resample,
      seed = seed,
      ...
    )

    # Format output
    out <- list(
      decision = result$decision,
      conf.int = result$conf.int,
      lower = result$lower,
      upper = result$upper,
      alpha = result$alpha,
      epsilon = c(result$epsilon1, result$epsilon2),
      estimate = c(p1 = result$p1_hat, p2 = result$p2_hat),
      sample_size = c(n1 = result$n1, n2 = result$n2),
      test_type = "two-sample",
      seed = result$seed,
      B = B
    )

  } else {
    # One-sample test
    result <- prop_test_dp_one_sample(
      p_hat = p_hat,
      n = n,
      lower = lower,
      upper = upper,
      epsilon = epsilon[1],
      alpha = alpha,
      B = B,
      max_resample = max_resample,
      seed = seed,
      ...
    )

    # Format output
    out <- list(
      decision = result$decision,
      conf.int = result$conf.int,
      lower = result$lower,
      upper = result$upper,
      alpha = result$alpha,
      epsilon = result$epsilon,
      estimate = result$p_hat,
      sample_size = result$n,
      test_type = "one-sample",
      seed = result$seed,
      B = B
    )
  }

  class(out) <- "dp_tost_prop"
  return(out)
}


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


#' Two Sample DP-TOST for Proportions
#'
#' @description
#' Internal function implementing the DP-TOST procedure for two-sample proportion
#' difference equivalence testing under differential privacy.
#'
#' @param p1_hat Privatized sample proportion for group 1 (with noise already added)
#' @param p2_hat Privatized sample proportion for group 2 (with noise already added)
#' @param n1 Sample size for group 1 (public information)
#' @param n2 Sample size for group 2 (public information)
#' @param lower Lower equivalence bound for difference (p1 - p2)
#' @param upper Upper equivalence bound for difference (p1 - p2)
#' @param epsilon1 Privacy budget for group 1. Default is 1
#' @param epsilon2 Privacy budget for group 2. Default is 1
#' @param alpha Significance level. Default is 0.05
#' @param B Number of Monte Carlo replications. Default is 10000
#' @param max_resample Maximum resampling attempts for invalid solutions. Default is 100
#' @param seed Random seed for reproducibility. Default is 1337
#' @param ... Additional arguments (currently unused)
#'
#' @return A list containing:
#'   \item{conf.int}{(1-2α) confidence interval for difference}
#'   \item{lower}{Lower equivalence bound}
#'   \item{upper}{Upper equivalence bound}
#'   \item{epsilon1}{Privacy budget for group 1}
#'   \item{epsilon2}{Privacy budget for group 2}
#'   \item{alpha}{Significance level}
#'   \item{decision}{Logical, TRUE if equivalence established}
#'   \item{p1_hat}{Privatized proportion for group 1}
#'   \item{p2_hat}{Privatized proportion for group 2}
#'   \item{n1}{Sample size for group 1}
#'   \item{n2}{Sample size for group 2}
#'   \item{seed}{Random seed used}
#'
#' @keywords internal
#' @noRd
prop_test_dp_two_sample <- function(p1_hat, p2_hat, n1, n2, lower, upper,
                                     epsilon1 = 1, epsilon2 = 1, alpha = 0.05,
                                     B = 10^3, max_resample = 100, seed = 1337, ...) {

  # Check seed validity
  if (is.na(seed) || !is.numeric(seed) || length(seed) != 1) {
    stop("seed must be a single numeric value")
  }

  # Monte Carlo simulation for both proportions
  nu_samples_1 <- rep(NA, B)
  nu_samples_2 <- rep(NA, B)

  for (b in 1:B) {
    # Reconstruct p1
    nu_samples_1[b] <- solve_proportion_matching(
      p_hat = p1_hat,
      n = n1,
      epsilon = epsilon1,
      max_resample = max_resample,
      seed = seed + b
    )

    # Reconstruct p2
    nu_samples_2[b] <- solve_proportion_matching(
      p_hat = p2_hat,
      n = n2,
      epsilon = epsilon2,
      max_resample = max_resample,
      seed = seed + B + b
    )
  }

  # Compute difference distribution
  diff_samples <- nu_samples_1 - nu_samples_2

  # Confidence interval (percentile method)
  ci_lower <- quantile(diff_samples, probs = alpha, na.rm = TRUE)
  ci_upper <- quantile(diff_samples, probs = 1 - alpha, na.rm = TRUE)

  # Equivalence decision: CI must be entirely within [lower, upper]
  decision <- (ci_lower > lower) && (ci_upper < upper)

  # Return result
  result <- list(
    conf.int = c(ci_lower, ci_upper),
    lower = lower,
    upper = upper,
    epsilon1 = epsilon1,
    epsilon2 = epsilon2,
    alpha = alpha,
    decision = decision,
    p1_hat = p1_hat,
    p2_hat = p2_hat,
    n1 = n1,
    n2 = n2,
    seed = seed
  )

  return(result)
}


#' Print Results of DP-TOST for Proportions
#'
#' @param x A \code{dp_tost_prop} object returned by \code{prop_test_equiv_dp}
#' @param ticks Number of ticks to print the confidence interval in the console. Default is 30.
#' @param rn Number of digits to consider when printing the results. Default is 5.
#' @param ... Further arguments to be passed to or from methods.
#' @return Prints object invisibly.
#' @importFrom cli cli_text col_green col_red symbol
#'
#' @export
print.dp_tost_prop <- function(x, ticks = 30, rn = 5, ...) {

  # Decision message
  if (x$decision) {
    cli_text(col_green("{symbol$tick} Accept equivalence"))
  } else {
    cli_text(col_red("{symbol$cross} Can't accept equivalence"))
  }

  # Compute visual representation
  lower_be <- x$conf.int[1] > x$lower
  upper_be <- x$conf.int[2] < x$upper

  rg <- range(c(x$conf.int, x$lower, x$upper))
  rg_delta <- rg[2] - rg[1]
  std_be_interval <- round(ticks * (c(x$lower, x$upper) - rg[1]) / rg_delta) + 1
  std_zero <- round(-ticks * rg[1] / rg_delta) + 1
  std_fit_interval <- round(ticks * (x$conf.int - rg[1]) / rg_delta) + 1
  std_fit_interval_center <- round(ticks * (sum(x$conf.int) / 2 - rg[1]) / rg_delta) + 1

  # Print equivalence region
  cat("Equiv. Region:  ")
  for (i in 1:(ticks + 1)) {
    if (i >= std_be_interval[1] && i <= std_be_interval[2]) {
      if (i == std_be_interval[1]) {
        cat("|-")
      } else if (i == std_be_interval[2]) {
        cat("-|")
      } else if (i == std_zero) {
        cat("-0-")
      } else {
        cat("-")
      }
    } else {
      cat(" ")
    }
  }
  cat("\n")

  # Print estimated interval
  cat("Estim. Inter.:  ")
  for (i in 1:(ticks + 1)) {
    if (i >= std_fit_interval[1] && i <= std_fit_interval[2]) {
      if (i == std_fit_interval[1]) {
        if (i > std_be_interval[1] && i < std_be_interval[2]) {
          cat(col_green("(-"))
        } else if (lower_be) {
          cat(col_green("(-"))
        } else {
          cat(col_red("(-"))
        }
      } else if (i == std_fit_interval[2]) {
        if (i > std_be_interval[1] && i < std_be_interval[2]) {
          cat(col_green("-)"))
        } else if (upper_be) {
          cat(col_green("-)"))
        } else {
          cat(col_red("-)"))
        }
      } else if (i == std_fit_interval_center) {
        if (i >= std_be_interval[1] && i <= std_be_interval[2]) {
          cat(col_green("-x-"))
        } else {
          cat(col_red("-x-"))
        }
      } else {
        if (i >= std_be_interval[1] && i <= std_be_interval[2]) {
          cat(col_green("-"))
        } else {
          cat(col_red("-"))
        }
      }
    } else {
      cat(" ")
    }
  }
  cat("\n")

  # Print confidence interval
  cat("CI = (")
  cat(format(round(x$conf.int[1], rn), nsmall = rn))
  cat(" ; ")
  cat(format(round(x$conf.int[2], rn), nsmall = rn))
  cat(")\n\n")

  # Print method and parameters
  if (x$test_type == "one-sample") {
    cat("Method: DP-TOST (proportion)\n")
  } else {
    cat("Method: DP-TOST (two-sample proportion)\n")
  }

  cat("alpha = ")
  cat(x$alpha)
  cat("; Equiv. bounds = [")
  cat(format(round(x$lower, rn)))
  cat(", ")
  cat(format(round(x$upper, rn)))
  cat("]\n")

  # Privacy parameters
  if (x$test_type == "one-sample") {
    cat("epsilon = ")
    cat(x$epsilon)
  } else {
    cat("epsilon = (")
    cat(x$epsilon[1])
    cat(", ")
    cat(x$epsilon[2])
    cat(")")
  }
  cat("; B = ")
  cat(x$B)
  cat("\n")

  # Sample size and estimates
  if (x$test_type == "one-sample") {
    cat("n = ")
    cat(x$sample_size)
    cat("; p_hat = ")
    cat(format(round(x$estimate, rn), nsmall = rn))
  } else {
    cat("n = (")
    cat(x$sample_size[1])
    cat(", ")
    cat(x$sample_size[2])
    cat("); p_hat = (")
    cat(format(round(x$estimate[1], rn), nsmall = rn))
    cat(", ")
    cat(format(round(x$estimate[2], rn), nsmall = rn))
    cat(")")
  }
  cat("\n")

  invisible(x)
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
