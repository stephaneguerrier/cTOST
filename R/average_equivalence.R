#' @title Finite Sample Adjustment for Average (Bio)Equivalence Assessment
#'
#' @description Computes finite sample corrected versions of the standard (univariate or multivariate) TOST, as developed in Boulaguiem et al. (2024, https://doi.org/10.1002/sim.9993), Boulaguiem et al. (2024, https://doi.org/10.1002/sim.10258), and Insolia et al. (2025).
#'
#' @param theta A \code{numeric} value or vector representing the estimated difference(s) (e.g., between a generic and reference product).
#' @param sigma A \code{numeric} value (univariate) or \code{matrix} (multivariate) corresponding to the estimated variance of \code{theta}.
#' @param nu A \code{numeric} value specifying the degrees of freedom. In the multivariate case, it is assumed to be the same across all dimensions.
#' @param delta A \code{numeric} value or vector defining the (bio)equivalence margin(s). The procedure assumes symmetry, i.e., the (bio)equivalence region is \eqn{(-\delta, \delta)}. In the multivariate case, it is assumed to be the same across all dimensions.
#' @param method A \code{character} string specifying the finite sample adjustment method. Available methods are: \code{"unadjusted"} (standard unadjusted TOST), \code{"alpha"} (alpha-TOST), \code{"delta"} (delta-TOST, not implemented for multivariate settings), and \code{"optimal"} (cTOST, default). See Details.
#' @param alpha A \code{numeric} value specifying the significance level, which must be between 0 and 0.5 (default: \code{alpha = 0.05}).
#' @param B A \code{numeric} value specifying the number of Monte Carlo replications, required for some methods (default: \code{B = 10^4}).
#' @param correction A \code{character} string specifying the correction method. Available options are: \code{"none"} (no correction), \code{"offline"} (offline adjustment), and \code{"bootstrap"} (bootstrap adjustment). In univariate settings, the default is \code{"offline"}; in multivariate settings, the default is \code{"bootstrap"} if \code{nu} < 100, otherwise \code{"none"}.
#' @param seed A \code{numeric} value specifying a seed for reproducibility (default: \code{seed = 101010}).
#' @param ... Additional parameters.
#'
#' @details
#' In univariate settings, three adjustment methods are available: optimal (cTOST, \code{method = "optimal"}) as proposed in Insolia et al. (2025), alpha-TOST (\code{method = "alpha"}), and delta-TOST (\code{method = "delta"}), both proposed in Boulaguiem et al. (2024, https://doi.org/10.1002/sim.9993, https://doi.org/10.1002/sim.10258). For multivariate settings, only cTOST and alpha-TOST are implemented.
#'
#' The cTOST, alpha-TOST, and delta-TOST methods apply different finite sample adjustments. Alpha-TOST corrects the significance level, while delta-TOST adjusts the equivalence limits. The cTOST method is based on a more complex approach, and in small samples (typically less than 30), additional corrections may be beneficial. The \code{correction} argument further adjusts the test level to prevent liberal inference; see Insolia et al. (2025) for details.
#'
#' Generally, cTOST outperforms other methods, with alpha-TOST performing better than delta-TOST. For this reason, delta-TOST is not implemented for multivariate settings and is not recommended.
#'
#' @return An object of class \code{tost} with the following components:
#' \itemize{
#'   \item \code{decision}: Logical; indicates whether (bio)equivalence is accepted.
#'   \item \code{ci}: Confidence region at the \eqn{1 - 2\alpha} level.
#'   \item \code{theta}: The estimated difference(s) used in the test.
#'   \item \code{sigma}: The estimated variance of \code{theta}; a \code{numeric} value (univariate) or \code{matrix} (multivariate).
#'   \item \code{nu}: The degrees of freedom used in the test.
#'   \item \code{alpha}: The significance level used in the test.
#'   \item \code{corrected_alpha}: The significance level after adjustment (if \code{method = "alpha"}).
#'   \item \code{corrected_delta}: The (bio)equivalence limits after adjustment (if \code{method = "delta"}).
#'   \item \code{delta}: The (bio)equivalence limits used in the test.
#'   \item \code{method}: The adjustment method used (optimal, alpha-TOST, or delta-TOST).
#'   \item \code{setting}: The setting used ("univariate" or "multivariate").
#' }
#'
#' @export
#' @examples
#' data(skin)
#' theta_hat = diff(apply(skin, 2, mean))
#' nu = nrow(skin) - 1
#' sig_hat = var(apply(skin, 1, diff)) / nu
#'
#' # alpha-TOST
#' atost = ctost(theta = theta_hat, sigma = sig_hat, nu = nu,
#'               alpha = 0.05, delta = log(1.25), method = "alpha")
#' atost
#' compare_to_tost(atost)
#'
#' # delta-TOST
#' dtost = ctost(theta = theta_hat, sigma = sig_hat, nu = nu,
#'               alpha = 0.05, delta = log(1.25), method = "delta")
#' dtost
#' compare_to_tost(dtost)
ctost = function(theta, sigma, nu, delta, alpha = 0.05, method, B = 10^4, seed = 101010, correction = NULL, ...){

  # Check inputs
  if (alpha < 0.0000001 || alpha > 0.5){
    stop("alpha must be in (0, 0.5).")
  }

  n_theta = length(theta)

  if (n_theta == 1){
    # Univariate setting
    setting = "univariate"
    if (length(sigma) > 1 || length(delta) > 1){
      stop("sigma and delta must be scalars in univariate settings.")
    }
  }else{
    setting = "multivariate"
    if (!is.matrix(sigma) || ncol(sigma) != nrow(sigma)){
      stop("sigma must be a square matrix.")
    }

    if (length(delta) > 1){
      stop("delta is assumed to be a scalar, implying that we consider the alternative theta in (-delta, delta) in each dimension.")
    }
  }

  if (!(method %in% c("unadjusted", "alpha", "delta", "optimal"))) {
    stop("Available methods are 'unadjusted' (standard TOST), 'alpha' (alpha-TOST), 'delta' (delta-TOST), and 'optimal' (cTOST).")
  }

  if (method == "delta" && setting == "multivariate") {
    stop("The delta-TOST method is not implemented for multivariate settings.")
  }

  if (method == "optimal" && setting == "multivariate") {
    stop("The cTOST ('optimal') method is not implemented for multivariate settings.")
  }

  if (!(correction %in% c("none", "bootstrap", "offline"))) {
    stop("Available correction methods are 'none', 'bootstrap', and 'offline'.")
  }

  if (setting == "univariate"){
    # alpha-TOST
    # Transform variance into standard deviation
    sigma = sqrt(sigma)
    if (method == "optimal"){
      return(xtost(theta_hat = theta, sig_hat = sigma, nu = nu, alpha = alpha, delta = delta, correction = correction,
            B = B, seed = seed))
    }

    if (method == "alpha"){
      corrected_alpha = get_alpha_TOST(alpha = alpha, sigma = sigma, nu = nu, delta = delta)$root
      decision = abs(theta) < (delta - qt(1 - corrected_alpha, df = nu) * sigma)
      ci = theta + c(-1, 1) * qt(1 - corrected_alpha, df = nu) * sigma
      out = list(decision = decision, ci = ci, theta = theta,
                 sigma = sigma^2, nu = nu, alpha = alpha,
                 corrected_alpha = corrected_alpha,
                 delta = delta, method = "alpha-TOST",
                 setting = setting)
      class(out) = "tost"
      return(out)
    }

    # delta-TOST
    if (method == "delta"){
      corrected_delta = get_delta_TOST(sigma = sigma, alpha = alpha, delta = delta, nu = nu)$root
      decision = abs(theta) < (corrected_delta - qt(1 - alpha, df = nu) * sigma)
      ci = theta + c(-1, 1) * qt(1 - alpha, df = nu) * sigma
      out = list(decision = decision, ci = ci, theta = theta,
                 sigma = sigma^2, nu = nu, alpha = alpha,
                 corrected_delta = corrected_delta,
                 delta = delta, method = "delta-TOST",
                 setting = setting)
      class(out) = "tost"
      return(out)
    }

    # Unadjusted TOST
    if (method == "unadjusted"){
      return(tost(theta = theta, sigma = sigma, nu = nu, delta = delta, alpha = alpha))
    }

  }else{
    # Multivariate
    if (method == "delta"){
      stop("Multivariate delta-TOST is not implemented.")
    }
    if (method == "alpha"){


      corrected_alpha = get_alpha_TOST_MC_mv(alpha = alpha, Sigma = sigma,
                                           nu = nu, delta = delta, B = B)

      delta_vec = rep(delta, ncol(sigma))
      t_alpha = qt(1 - corrected_alpha$min, df=nu)
      lower = theta - t_alpha * sqrt(diag(sigma))
      upper = theta + t_alpha * sqrt(diag(sigma))
      decision = (lower > -delta_vec) & (upper < delta_vec)

      out = list(decision = decision, ci = cbind(lower, upper), theta = theta,
                 sigma = sigma, nu = nu, alpha = alpha,
                 corrected_alpha = corrected_alpha$min,
                 delta = delta, method = "alpha-TOST", setting = setting)
      class(out) = "mtost"
      return(out)

    }
  }
}
