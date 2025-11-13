#' Quantile equivalence testing procedures
#'
#' @description
#' Performs a Two One-Sided Test (TOST) to assess the equivalence of a quantile
#' from a test population ($Y$) with the corresponding quantile from a reference
#' population ($X$), assuming the data are normally distributed. The test evaluates
#' if the true quantile $\pi_y$ is within a pre-specified equivalence margin $\delta$
#' around the reference quantile $\pi_x$.
#'
#' The null hypotheses for the two one-sided tests are:
#' $H_{01}: \pi_y \ge \pi_x - \delta$ and $H_{02}: \pi_y \le \pi_x + \delta$.
#' Equivalence is concluded if both null hypotheses are rejected.
#'
#' @param x A \code{numeric} vector of data for the reference group ($X$), or a  \code{list} containing `mean`, `sd`, and `n`.
#' @param y A  \code{numeric} vector of data for the test group ($Y$), or a  \code{list} containing `mean`, `sd`, and `n`.
#' @param pi_x A  \code{numeric} scalar or vector specifying the quantile(s) of interest in the reference group $X$ (e.g., 0.9 for the 90th percentile).
#' @param delta A \code{numeric} value or vector defining the (bio)equivalence margin(s). The procedure assumes symmetry, i.e., the (bio)equivalence region is \eqn{(\pi_x-\delta, \pi_x+\delta)}. In the multivariate case, it is assumed to be the same across all dimensions.
#' @param alpha A \code{numeric} value specifying the significance level, which must be between 0 and 0.5 (default: \code{alpha = 0.05}).
#' @param method A \code{character} string specifying the finite sample adjustment method. Available methods are: \code{"unadjusted"} (standard unadjusted qTOST), \code{"alpha"} (alpha-qTOST). Default: \code{method = "alpha"}.
#' @param B A \code{numeric} value specifying the number of Monte Carlo replications, required for the `"alpha"` method (default: \code{B = 10^4}).
#' @param seed A \code{numeric} value specifying a seed for reproducibility (default: \code{seed = 12345}).
#' @param tol A \code{numeric} value specifying a tolerance level (default: \code{tol = .Machine$double.eps}).
#' @param max_iter A \code{numeric} value specifying a maximum number of iteration to compute the supremum at which the size is assessed (default: \code{max_iter = 10}).
#' @param ... Additional arguments (not currently used).
#'
#' @return
#' An object of class `qtost` (for a single quantile) or `mqtost` (for multiple quantiles) with the following components:
#' \itemize{
#'   \item `decision`: The (component-wise) equivalence decision (`TRUE` or `FALSE`).
#'   \item `method`: The method used (`"qTOST"` or `"alpha-qTOST"`).
#'   \item `ci`: The $(1 - 2\alpha)$ confidence interval for the estimated quantile $\hat{\pi}_y$.
#'   \item `pi_y_hat`: The point estimate for the quantile in group Y.
#'   \item `eq_region`: The defined equivalence region $[\pi_x - \delta, \pi_x + \delta]$.
#'   \item `alpha`: The nominal significance level.
#'   \item `corrected_alpha`: The adjusted alpha level used (for `alpha-qTOST` method).
#' }
#'
#' @export
#' @examples
#' # https://www.accessdata.fda.gov/drugsatfda_docs/label/2024/021814s030lbl.pdf
#' # C_trough: table 5
#' # reference male group
#' x_bar_orig = 35.6
#' x_sd_orig = 16.7
#' n_x = 106
#' # target female group
#' y_bar_orig = 41.6
#' y_sd_orig = 24.3
#' n_y = 14
#' # data transformation
#' x_bar = log(x_bar_orig^2 / sqrt(x_bar_orig^2 + x_sd_orig^2))
#' x_sd = sqrt(log(1 + (x_sd_orig^2 / x_bar_orig^2)))
#' y_bar = log(y_bar_orig^2 / sqrt(y_bar_orig^2 + y_sd_orig^2))
#' y_sd = sqrt(log(1 + (y_sd_orig^2 / y_bar_orig^2)))
#' x = list(mean=x_bar, sd=x_sd, n=n_x)
#' y = list(mean=y_bar, sd=y_sd, n=n_y)
#' # qTOST
#' sqtost <- qtost(x, y, pi_x = 0.8, delta = 0.15, method = "unadjusted")
#' print(sqtost)
#' aqtost <- qtost(x, y, pi_x = 0.8, delta = 0.15, method = "alpha")
#' print(aqtost)
#'
#' # Example 2: Using raw data to test a single quantile
#' set.seed(12345)
#' x_data <- rnorm(n = 30, mean = 0, sd = 0.1)
#' y_data <- rnorm(n = 10, mean = 0, sd = 0.1)
#'
#' # Test the 90th percentile with a margin of delta = 0.05
#' sqtost <- qtost(x = x_data, y = y_data, pi_x = 0.8, delta = 0.15, method = "unadjusted")
#' print(sqtost)
#' aqtost <- qtost(x_data, y_data, pi_x = 0.8, delta = 0.15, method = "alpha")
#' print(aqtost)
#'
#' # Example 2: Using summary statistics with the unadjusted method
#' x_stats <- list(mean = 100, sd = 15, n = 50)
#' y_stats <- list(mean = 102, sd = 16, n = 50)
#'
#' result_summary <- qtost(x_stats, y_stats, pi_x = 0.95, delta = 0.03,
#'                         method = "unadjusted")
#' print(result_summary)
#'
qtost <- function(x, y, pi_x, delta, alpha = 0.05, method = "alpha", B = NULL, seed = 101010, tol = 1e-6, max_iter=10, tolpower = 1e-3, MC_sup = T, ...) {
  if (is.null(B)) {
   B = if (length(pi_x)==1) 1e5 else 1e4
  }
  # Pre-processing
  x <- .extract_stats(x, "x")
  x_bar = x$mean
  x_sd = x$sd
  n_x = x$n
  y <- .extract_stats(y, "y")
  y_bar = y$mean
  y_sd = y$sd
  n_y = y$n
  # Check inputs
  if (n_x <= 1 || n_y <= 1) {
    stop("Sample sizes (n) for both groups must be greater than 1.")
  }
  if (x_sd < tol || y_sd < tol) {
    stop("Variances must be greater than 0.")
  }
  if ((alpha < tol) || (alpha > (0.5-tol))) {
    stop("alpha must be in (0, 0.5).")
  }
  if (!(method %in% c("unadjusted", "alpha"))) {
    stop("method must be 'unadjusted' (for the standard qTOST) or 'alpha' (for the corrected alpha-qTOST).")
  }
  if (length(delta) > 1 || delta <= 0) {
    stop("Equivalence margin 'delta' is assumed to be a positive scalar (a generalization will be coming soon).")
  } else {
    delta_l = pi_x - delta
    delta_u = pi_x + delta
  }
  if (length(pi_x) > 1) {
    setting = "multiple"
    if (length(pi_x) > 2) {
      stop("No method is currently available to jointly assess more than two quantiles (coming soon).")
    }
  } else {
    setting = "single"
  }
  # Estimates
  l = n_y / n_x
  gamma = y_sd^2 / x_sd^2
  theta = (x_bar + x_sd * qnorm(pi_x) - y_bar) / y_sd
  sigma = sqrt(1 / n_y * (1 + theta^2 / 2 + l / gamma * (1 + (qnorm(pi_x))^2 / 2)))
  if (setting == "single") {
    if (method == "unadjusted") {
      out = qtost_core(theta, sigma, pi_x, delta_l, delta_u, alpha = alpha, corrected_alpha = NULL)
      return(out)
    }
    else if (method == "alpha") {
      alpha_star = get_alpha_qTOST_MC(gamma = gamma, n_x = n_x, n_y = n_y,
                                      pi_x = pi_x, delta_l = delta_l, delta_u = delta_u,
                                      alpha = alpha, B = B,
                                      seed = seed, tol = tol)
      out = qtost_core(theta, sigma, pi_x, delta_l, delta_u, alpha = alpha, corrected_alpha = alpha_star$root)
      return(out)
    }
  } else if (setting == "multiple") {
    if (method == "unadjusted") {
      out = qtost_core(theta, sigma, pi_x, delta_l, delta_u, alpha = alpha, corrected_alpha = NULL)
      return(out)
    } else if (method == "alpha") {
      alpha_star = get_mult_alpha_qTOST_MC(gamma = gamma, n_x = n_x, n_y = n_y,
                                           pi_x = pi_x, delta_l = delta_l, delta_u = delta_u,
                                           alpha = alpha, B = B,
                                           tol = tol, seed = seed, max_iter = max_iter,
                                           tolpower = tolpower, MC_sup = MC_sup)
      out = qtost_core(theta, sigma, pi_x, delta_l, delta_u, alpha = alpha, corrected_alpha = alpha_star$min)
      return(out)
    }
  }
  return(out)
}
