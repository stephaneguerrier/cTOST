# DP-TOST for Bounded Means (One Sample)
# Based on Section 3.2 of DP_equivalence_test.pdf

#' Differentially Private Mean Equivalence Test (Internal)
#'
#' @description
#' Internal unified interface for DP-TOST that automatically selects the best implementation.
#' Uses ultra-fast C++ by default (15-250x speedup), with option to use pure R.
#'
#' @param a Lower bound for data truncation
#' @param b Upper bound for data truncation
#' @param n Sample size
#' @param epsilon Privacy budget
#' @param mean_private_obs Observed privatized mean
#' @param sd_private_obs Observed privatized standard deviation
#' @param lower Lower equivalence bound for TOST
#' @param upper Upper equivalence bound for TOST
#' @param B Number of bootstrap iterations (default: 1000)
#' @param alpha Significance level (default: 0.05)
#' @param seed Random seed for reproducibility (default: NULL)
#' @param method Implementation to use: "cpp" (default, fast) or "r" (exact R)
#'
#' @return List with components:
#'   \item{decision}{1 if equivalence established, 0 otherwise}
#'   \item{ci_lower}{Lower bound of confidence interval}
#'   \item{ci_upper}{Upper bound of confidence interval}
#'   \item{mean_private_obs}{Observed privatized mean (input)}
#'   \item{sd_private_obs}{Observed privatized SD (input)}
#'   \item{mu_estimates}{Bootstrap estimates of mu}
#'   \item{sigma_estimates}{Bootstrap estimates of sigma}
#'   \item{method}{Method used ("cpp" or "r")}
#'
#' @details
#' This function performs equivalence testing under differential privacy.
#'
#' **Method Selection:**
#' - `"cpp"` (default): Ultra-fast C++ implementation using Nelder-Mead
#'   - 15-250x faster depending on problem size
#'   - Uses native C++ RNG (different from R but statistically equivalent)
#'   - Recommended for production and large-scale simulations
#'
#' - `"r"`: Pure R implementation using L-BFGS-B
#'   - Uses R's optim() with gradient-based optimization
#'   - Exact reproducibility with R workflows
#'   - Useful for validation and debugging
#'
#' **Performance Guide:**
#' - Small problems (n~100, B~1000): C++ is 60-250x faster
#' - Large problems (n~1000, B~10000): C++ is 15-20x faster
#'
#' @keywords internal
dp_mean_test <- function(a, b, n, epsilon, mean_private_obs, sd_private_obs,
                         lower, upper, B = 1000, alpha = 0.05, seed = NULL,
                         method = c("cpp", "r")) {

  # Match method argument
  method <- match.arg(method)

  # Dispatch to appropriate implementation
  if (method == "cpp") {
    # Load C++ implementation if not already loaded
    if (!exists("tost_dp_one_sample_ultra_fast")) {
      stop("C++ implementation not loaded. Please run: sourceCpp('src/dp_mean_ultra_fast.cpp')")
    }

    result <- tost_dp_one_sample_ultra_fast(
      a = a, b = b, n = n, epsilon = epsilon,
      mean_private_obs = mean_private_obs,
      sd_private_obs = sd_private_obs,
      lower = lower, upper = upper,
      B = B, alpha = alpha, seed = seed
    )

  } else {
    # Use pure R implementation
    result <- tost_dp_one_sample(
      a = a, b = b, n = n, epsilon = epsilon,
      mean_private_obs = mean_private_obs,
      sd_private_obs = sd_private_obs,
      lower = lower, upper = upper,
      B = B, alpha = alpha, seed = seed,
      use_cpp = FALSE
    )
  }

  # Add method info to result
  result$method <- method

  return(result)
}

# Compute privatized statistics T_X for bounded means
# Inputs:
#   a, b: bounds for data truncation
#   mu, sigma: true population parameters
#   n: sample size
#   scale_mean, scale_sd: scale parameters for Laplace noise
#   z: vector of n standard normal N(0,1) random variables
#   u1, u2: two Uniform(0,1) random variables for Laplace noise
# Returns: c(mean_private, sd_private)
compute_tx_private <- function(a, b, mu, sigma, n, scale_mean, scale_sd, z, u1, u2) {

  # Step 1: Generate data from N(mu, sigma^2) using standard normals
  x <- mu + sigma * z

  # Step 2: Truncate data to [a, b]
  x_truncated <- pmin(pmax(x, a), b)

  # Step 3: Compute sample statistics (non-private)
  mean_obs <- mean(x_truncated)
  sd_obs <- sd(x_truncated)

  # Step 4: Add Laplace noise for DP
  # Laplace(0, lambda) can be generated using: -lambda * sign(u - 0.5) * log(1 - 2 * |u - 0.5|)
  noise_mean <- -scale_mean * sign(u1 - 0.5) * log(1 - 2 * abs(u1 - 0.5))
  noise_sd <- -scale_sd * sign(u2 - 0.5) * log(1 - 2 * abs(u2 - 0.5))

  # Privatized statistics
  mean_private <- mean_obs + noise_mean
  sd_private <- sd_obs + noise_sd

  # Return only privatized statistics
  c(mean_private, sd_private)
}

# Objective function for moment-matching optimization (Section 3.2, equation 6)
# Inputs:
#   theta: c(mu, sigma) - parameters to optimize
#   a, b, n, scale_mean, scale_sd, z, u1, u2: inputs for compute_tx_private
#   mean_private_obs, sd_private_obs: observed privatized statistics
#   use_cpp: if TRUE, use C++ version of compute_tx_private
# Returns: squared L2 norm between observed and simulated privatized statistics
obj_fun_one_sample <- function(theta, a, b, n, scale_mean, scale_sd, z, u1, u2,
                                mean_private_obs, sd_private_obs, use_cpp = FALSE) {

  # Compute privatized statistics for given theta = c(mu, sigma)
  if (use_cpp) {
    stat_private <- compute_tx_private_cpp(a, b, theta[1], theta[2], n, scale_mean, scale_sd, z, u1, u2)
  } else {
    stat_private <- compute_tx_private(a, b, theta[1], theta[2], n, scale_mean, scale_sd, z, u1, u2)
  }

  # Observed privatized statistics
  stat_private_obs <- c(mean_private_obs, sd_private_obs)

  # Compute squared L2 norm
  sum((stat_private_obs - stat_private)^2)
}


# DP-TOST for one-sample bounded means (R implementation for reference)
# Inputs:
#   a, b: bounds for data truncation
#   n: sample size
#   epsilon: privacy budget
#   mean_private_obs, sd_private_obs: observed privatized statistics
#   lower, upper: equivalence bounds for TOST
#   B: number of bootstrap iterations
#   alpha: significance level (default 0.05)
#   seed: random seed for reproducibility
#   use_cpp: if TRUE, use C++ compute_tx_private for speedup (default FALSE)
# Returns: list with decision, ci_lower, ci_upper, observed stats, and estimates
# NOTE: Use tost_dp_one_sample_cpp() for full C++ version (faster but different optimizer)
tost_dp_one_sample <- function(a, b, n, epsilon, mean_private_obs, sd_private_obs,
                               lower, upper, B = 1000, alpha = 0.05, seed = NULL, use_cpp = FALSE) {

  # Set seed
  if (!is.null(seed)) {
    set.seed(seed)
  }

  # Compute scale parameters for Laplace noise
  delta_f_mean <- (b - a) / n
  delta_f_sd <- (b - a) / sqrt(n - 1)
  scale_mean <- delta_f_mean / (epsilon / 2)
  scale_sd <- delta_f_sd / (epsilon / 2)

  # Generate random variables: B times n matrix of N(0,1) and B times 2 matrix of Unif(0,1)
  z_matrix <- matrix(rnorm(B * n), nrow = B, ncol = n)
  u_matrix <- matrix(runif(B * 2), nrow = B, ncol = 2)

  # Store optimization results
  mu_estimates <- numeric(B)
  sigma_estimates <- numeric(B)

  # Initial values for optimization: use observed privatized statistics
  theta_init <- c(mean_private_obs, sd_private_obs)

  # For each bootstrap iteration
  for (i in 1:B) {
    # Extract i-th row of random variables
    z_i <- z_matrix[i, ]
    u1_i <- u_matrix[i, 1]
    u2_i <- u_matrix[i, 2]

    # Optimize objective function using L-BFGS-B with box constraints
    opt_result <- optim(
      par = theta_init,
      fn = obj_fun_one_sample,
      method = "L-BFGS-B",
      lower = c(a, 0.000001),           # mu >= a, sigma > 0
      upper = c(b, (b - a) / 2 + 5 * scale_sd),  # mu <= b, sigma <= (b-a)/2 + 5*scale_sd
      a = a, b = b, n = n,
      scale_mean = scale_mean,
      scale_sd = scale_sd,
      z = z_i,
      u1 = u1_i,
      u2 = u2_i,
      mean_private_obs = mean_private_obs,
      sd_private_obs = sd_private_obs,
      use_cpp = use_cpp
    )

    mu_estimates[i] <- opt_result$par[1]
    sigma_estimates[i] <- opt_result$par[2]
  }

  # Compute quantiles for TOST
  ci_lower <- quantile(mu_estimates, alpha, names = FALSE)
  ci_upper <- quantile(mu_estimates, 1 - alpha, names = FALSE)

  # TOST decision: accept equivalence if ci_lower > lower AND ci_upper < upper
  decision <- as.numeric(ci_lower > lower & ci_upper < upper)

  # Return results
  return(list(
    decision = decision,
    ci_lower = ci_lower,
    ci_upper = ci_upper,
    mean_private_obs = mean_private_obs,
    sd_private_obs = sd_private_obs,
    mu_estimates = mu_estimates,
    sigma_estimates = sigma_estimates
  ))
}


#' Differentially Private TOST for Mean Equivalence Testing
#'
#' @description
#' Performs equivalence testing for bounded means under differential privacy using
#' the DP-TOST procedure. Supports both one-sample and two-sample tests.
#'
#' @details
#' This function implements the DP-TOST (Differentially Private Two One-Sided Tests)
#' procedure for testing equivalence of means while maintaining differential privacy
#' guarantees. The test uses Monte Carlo simulation with moment-matching optimization
#' to reconstruct the sampling distribution from privatized statistics.
#'
#' For the one-sample case, the test evaluates whether a mean is equivalent to a
#' reference range [lower, upper]. For the two-sample case, it tests whether the
#' difference between two means (μ1 - μ2) falls within the equivalence bounds
#' [lower, upper].
#'
#' The privacy mechanism adds Laplace noise to both the observed mean and standard
#' deviation. The privacy budget epsilon controls the amount of noise: larger epsilon
#' values provide less privacy but more accurate inference.
#'
#' @param mean_private_obs Privatized sample mean for the first (or only) group.
#'   This should be the observed mean with differential privacy noise already added.
#' @param sd_private_obs Privatized sample standard deviation for the first (or only)
#'   group with differential privacy noise already added.
#' @param a Lower bound for data truncation in the first (or only) group.
#' @param b Upper bound for data truncation in the first (or only) group.
#' @param n Sample size for the first (or only) group.
#' @param mean_private_obs2 Privatized sample mean for the second group (optional).
#'   If NULL, a one-sample test is performed. Default is NULL.
#' @param sd_private_obs2 Privatized sample standard deviation for the second group
#'   (optional). Required if mean_private_obs2 is provided. Default is NULL.
#' @param a2 Lower bound for data truncation in the second group (optional).
#'   If NULL but mean_private_obs2 is provided, defaults to a. Default is NULL.
#' @param b2 Upper bound for data truncation in the second group (optional).
#'   If NULL but mean_private_obs2 is provided, defaults to b. Default is NULL.
#' @param n2 Sample size for the second group (optional). Required if mean_private_obs2
#'   is provided. Default is NULL.
#' @param lower Lower equivalence bound. For one-sample tests, this is the lower bound
#'   for the mean. For two-sample tests, this is the lower bound for the difference
#'   (μ1 - μ2). Typical value: -delta.
#' @param upper Upper equivalence bound. For one-sample tests, this is the upper bound
#'   for the mean. For two-sample tests, this is the upper bound for the difference
#'   (μ1 - μ2). Typical value: delta.
#' @param epsilon Privacy budget. Can be a single value (used for all samples) or a
#'   vector of length 2 (epsilon[1] for first group, epsilon[2] for second group).
#'   Larger values provide less privacy but more statistical power.
#' @param alpha Significance level for the test. The confidence interval will have
#'   level (1 - 2*alpha). Default is 0.05.
#' @param B Number of Monte Carlo replications for reconstructing the sampling
#'   distribution. Larger values provide more accurate results but take longer.
#'   Default is 10000.
#' @param seed Random seed for reproducibility. Default is 1337.
#' @param method Implementation to use: "cpp" (default, ultra-fast) or "r" (pure R).
#'   The C++ implementation is 15-250x faster depending on problem size.
#'
#' @return A list of class "tost_dp" containing:
#'   \item{decision}{Logical. TRUE if equivalence is established, FALSE otherwise.}
#'   \item{conf.int}{Confidence interval for the parameter of interest (mean for
#'     one-sample, difference for two-sample).}
#'   \item{lower}{Lower equivalence bound.}
#'   \item{upper}{Upper equivalence bound.}
#'   \item{alpha}{Significance level used.}
#'   \item{epsilon}{Privacy budget(s) used.}
#'   \item{estimate}{Point estimate(s) of the privatized mean(s).}
#'   \item{sample_size}{Sample size(s).}
#'   \item{test_type}{Character string indicating "one-sample" or "two-sample".}
#'   \item{method}{Method used ("cpp" or "r").}
#'   \item{seed}{Random seed used.}
#'   \item{B}{Number of Monte Carlo replications.}
#'
#' @examples
#' \dontrun{
#' # One-sample test: Is mean equivalent to range [1.5, 3.5]?
#' set.seed(123)
#' n <- 100
#' a <- 0
#' b <- 5
#' mu_true <- 2.5
#' sigma_true <- 1.0
#' epsilon <- 1
#'
#' # Generate truncated normal data
#' z <- rnorm(n)
#' x <- pmin(pmax(mu_true + sigma_true * z, a), b)
#' mean_obs <- mean(x)
#' sd_obs <- sd(x)
#'
#' # Add Laplace noise
#' delta_f_mean <- (b - a) / n
#' delta_f_sd <- (b - a) / sqrt(n - 1)
#' scale_mean <- delta_f_mean / (epsilon / 2)
#' scale_sd <- delta_f_sd / (epsilon / 2)
#'
#' u <- runif(2)
#' mean_private <- mean_obs - scale_mean * sign(u[1] - 0.5) * log(1 - 2 * abs(u[1] - 0.5))
#' sd_private <- sd_obs - scale_sd * sign(u[2] - 0.5) * log(1 - 2 * abs(u[2] - 0.5))
#'
#' result <- tost_equiv_dp(
#'   mean_private_obs = mean_private,
#'   sd_private_obs = sd_private,
#'   a = a, b = b, n = n,
#'   lower = 1.5,
#'   upper = 3.5,
#'   epsilon = 1,
#'   B = 1000
#' )
#' result
#'
#' # Two-sample test: Is difference (μ1 - μ2) equivalent to [-1, 1]?
#' set.seed(456)
#' n1 <- 100
#' n2 <- 100
#' a1 <- 0; b1 <- 5
#' a2 <- 0; b2 <- 5
#' mu1_true <- 2.5
#' mu2_true <- 2.3
#' sigma1_true <- 1.0
#' sigma2_true <- 1.2
#' epsilon <- 1
#'
#' # Generate data for both groups
#' z1 <- rnorm(n1)
#' z2 <- rnorm(n2)
#' x1 <- pmin(pmax(mu1_true + sigma1_true * z1, a1), b1)
#' x2 <- pmin(pmax(mu2_true + sigma2_true * z2, a2), b2)
#'
#' # Add Laplace noise to both groups
#' scale_mean1 <- ((b1 - a1) / n1) / (epsilon / 2)
#' scale_sd1 <- ((b1 - a1) / sqrt(n1 - 1)) / (epsilon / 2)
#' scale_mean2 <- ((b2 - a2) / n2) / (epsilon / 2)
#' scale_sd2 <- ((b2 - a2) / sqrt(n2 - 1)) / (epsilon / 2)
#'
#' u1 <- runif(2)
#' u2 <- runif(2)
#' mean_private1 <- mean(x1) - scale_mean1 * sign(u1[1] - 0.5) * log(1 - 2 * abs(u1[1] - 0.5))
#' sd_private1 <- sd(x1) - scale_sd1 * sign(u1[2] - 0.5) * log(1 - 2 * abs(u1[2] - 0.5))
#' mean_private2 <- mean(x2) - scale_mean2 * sign(u2[1] - 0.5) * log(1 - 2 * abs(u2[1] - 0.5))
#' sd_private2 <- sd(x2) - scale_sd2 * sign(u2[2] - 0.5) * log(1 - 2 * abs(u2[2] - 0.5))
#'
#' result <- tost_equiv_dp(
#'   mean_private_obs = mean_private1,
#'   sd_private_obs = sd_private1,
#'   a = a1, b = b1, n = n1,
#'   mean_private_obs2 = mean_private2,
#'   sd_private_obs2 = sd_private2,
#'   a2 = a2, b2 = b2, n2 = n2,
#'   lower = -1,
#'   upper = 1,
#'   epsilon = 1,
#'   B = 1000
#' )
#' result
#' }
#'
#'
#' @export
tost_equiv_dp <- function(mean_private_obs, sd_private_obs, a, b, n,
                          mean_private_obs2 = NULL, sd_private_obs2 = NULL,
                          a2 = NULL, b2 = NULL, n2 = NULL,
                          lower, upper, epsilon, alpha = 0.05,
                          B = 10^4, seed = 1337, method = "cpp", ...) {

  # Validate method argument
  if (!method %in% c("cpp", "r")) {
    stop("method must be either 'cpp' or 'r'")
  }

  # Determine test type
  two_sample <- !is.null(mean_private_obs2)

  # Validate inputs
  if (two_sample && (is.null(sd_private_obs2) || is.null(n2))) {
    stop("sd_private_obs2 and n2 must be provided when mean_private_obs2 is specified (two-sample test)")
  }

  if (!two_sample && (!is.null(sd_private_obs2) || !is.null(n2))) {
    warning("sd_private_obs2 and n2 are ignored for one-sample test")
  }

  # Set default bounds for second sample if not provided
  if (two_sample) {
    if (is.null(a2)) a2 <- a
    if (is.null(b2)) b2 <- b
  }

  # Handle epsilon: convert to vector of length 2 if needed
  if (length(epsilon) == 1) {
    epsilon_vec <- rep(epsilon, 2)
  } else if (length(epsilon) == 2) {
    epsilon_vec <- epsilon
  } else {
    stop("epsilon must be a single value or a vector of length 2")
  }

  # Perform the appropriate test
  if (two_sample) {
    # Two-sample test: Run DP test on both samples, then take difference

    # Run first sample
    result1 <- dp_mean_test(
      a = a, b = b, n = n,
      epsilon = epsilon_vec[1],
      mean_private_obs = mean_private_obs,
      sd_private_obs = sd_private_obs,
      lower = -Inf, upper = Inf,  # Dummy bounds, we only need estimates
      B = B,
      alpha = alpha,
      seed = seed,
      method = method
    )

    # Run second sample (use different seed to ensure independence)
    result2 <- dp_mean_test(
      a = a2, b = b2, n = n2,
      epsilon = epsilon_vec[2],
      mean_private_obs = mean_private_obs2,
      sd_private_obs = sd_private_obs2,
      lower = -Inf, upper = Inf,  # Dummy bounds, we only need estimates
      B = B,
      alpha = alpha,
      seed = seed + 999,  # Different seed for independence
      method = method
    )

    # Compute difference distribution: mu1 - mu2
    diff_estimates <- result1$mu_estimates - result2$mu_estimates

    # Compute confidence interval from difference distribution
    ci_lower <- quantile(diff_estimates, alpha, names = FALSE)
    ci_upper <- quantile(diff_estimates, 1 - alpha, names = FALSE)

    # TOST decision
    decision <- (ci_lower > lower) && (ci_upper < upper)

    # Format output
    out <- list(
      decision = decision,
      conf.int = c(ci_lower, ci_upper),
      lower = lower,
      upper = upper,
      alpha = alpha,
      epsilon = epsilon_vec,
      estimate = c(mu1 = mean_private_obs, mu2 = mean_private_obs2),
      sd_estimate = c(sigma1 = sd_private_obs, sigma2 = sd_private_obs2),
      sample_size = c(n1 = n, n2 = n2),
      bounds = list(group1 = c(a, b), group2 = c(a2, b2)),
      test_type = "two-sample",
      method = method,
      seed = seed,
      B = B
    )

  } else {
    # One-sample test
    result <- dp_mean_test(
      a = a, b = b, n = n,
      epsilon = epsilon_vec[1],
      mean_private_obs = mean_private_obs,
      sd_private_obs = sd_private_obs,
      lower = lower,
      upper = upper,
      B = B,
      alpha = alpha,
      seed = seed,
      method = method
    )

    # Extract confidence interval from result
    ci_lower <- result$ci_lower
    ci_upper <- result$ci_upper
    decision <- result$decision == 1

    # Format output
    out <- list(
      decision = decision,
      conf.int = c(ci_lower, ci_upper),
      lower = lower,
      upper = upper,
      alpha = alpha,
      epsilon = epsilon_vec[1],
      estimate = mean_private_obs,
      sd_estimate = sd_private_obs,
      sample_size = n,
      bounds = c(a, b),
      test_type = "one-sample",
      method = method,
      seed = seed,
      B = B
    )
  }

  class(out) <- "tost_dp"
  return(out)
}

#' Print method for DP-TOST mean test
#'
#' @param x An object of class \code{tost_dp}
#' @param ticks Number of characters for visual representation
#' @param rn Number of decimal places for rounding
#' @param ... Additional arguments (not used)
#'
#' @return Invisibly returns the input object
#' @export
#' @importFrom cli cli_text col_green col_red symbol
print.tost_dp <- function(x, ticks = 30, rn = 5, ...) {

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
    cat("Method: DP-TOST (mean)\n")
  } else {
    cat("Method: DP-TOST (two-sample mean)\n")
  }

  cat("alpha = ")
  cat(x$alpha)
  cat("; Equiv. bounds = [")
  cat(format(round(x$lower, rn)))
  cat(", ")
  cat(format(round(x$upper, rn)))
  cat("]\n")

  # Print DP parameters
  cat("epsilon = ")
  if (length(x$epsilon) == 1) {
    cat(x$epsilon)
  } else {
    cat("(")
    cat(x$epsilon[1])
    cat(", ")
    cat(x$epsilon[2])
    cat(")")
  }
  cat("; B = ")
  cat(x$B)
  cat("\n")

  # Print truncation bounds
  if (x$test_type == "one-sample") {
    cat("Truncation: (a, b) = (")
    cat(x$bounds[1])
    cat(", ")
    cat(x$bounds[2])
    cat(")\n")
  } else {
    cat("Truncation: (a1, a2) = (")
    cat(x$bounds$group1[1])
    cat(", ")
    cat(x$bounds$group2[1])
    cat("); (b1, b2) = (")
    cat(x$bounds$group1[2])
    cat(", ")
    cat(x$bounds$group2[2])
    cat(")\n")
  }

  invisible(x)
}
