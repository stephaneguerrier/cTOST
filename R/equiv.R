#' @title Finite Sample Adjusted (Bio)Equivalence Testing
#'
#' @description This function is used to compute finite sample corrected version of the standard (univariate or multivariate) TOST.
#'
#' @param theta                 A \code{numeric} value corresponding to the difference of means.
#' @param sigma                 A \code{numeric} value corresponding to the standard error.
#' @param nu                    A \code{numeric} value corresponding to the number of degrees of freedom.
#' @param delta                 A \code{numeric} value corresponding to (bio)equivalence limit. We assume symmetry, i.e, the (bio)equivalence interval corresponds to (-delta,delta)
#' @param method                A \code{character} value corresponding to the considered finite sample adjustment method, see Details below for more information.
#' @param alpha                 A \code{numeric} value specifying the significance level (default: alpha = 0.05).
#' @param B                     TO DOCUMENT ???
#' @param ...                   Additional parameters.
#'
#'
#' @details
#' In univariate settings, two methods are available: alpha-TOST (using method = 'alpha') and delta-TOST (using method = 'delta').
#' The alpha-TOST and delta-TOST are introduced in Boulaguiem et al. (2024, <https://doi.org/10.1002/sim.9993>). The former is a corrective procedure of the significance level applied to the TOST while the latter
#' adjusts the equivalence limits. In general, the alpha-TOST appears to outperform the delta-TOST.
#'
#' In multivariate setting, the only available method is the (multivariate) alpha-TOST (using method = 'alpha') introduced in Boulaguiem et al. (2024, bioarxiv).
#'
#' @author Younes Boulaguiem, Stéphane Guerrier, Dominique-Laurent Couturier, Luca Insolia
#'
#' @return A \code{tost} object with the structure:
#' \itemize{
#'  \item decision:    A boolean variable indicating whether (bio)equivalence is accepted or not.
#'  \item ci:          Confidence interval at the 1 - 2*alpha level.
#'  \item theta:       The difference of means used for the test.
#'  \item sigma:       The standard error used for the test.
#'  \item nu:          The number of degrees of freedom used for the test.
#'  \item alpha:       The significance level used for the test.
#'  \item delta:       The (bio)equivalence limits used for the test.
#'  \item method:      The method used for the test (alpha-TOST and delta-TOST).
#'  \item setting:     The setting used (univariate or multivariate)
#' }
#' @export
#' @examples
#' data(skin)
#'
#' theta_hat = diff(apply(skin,2,mean))
#' nu = nrow(skin) - 1
#' sig_hat = sd(apply(skin,1,diff))/sqrt(nu)
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
ctost = function(theta, sigma, nu, delta, alpha = 0.05, method = "alpha", B = 10^4, ...){

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

  if (!(method %in% c("alpha", "delta"))){
    stop("Available methods are 'alpha' (for the alpha-TOST) and 'delta' (for the delta-TOST).")
  }

  if (method == "delta" && setting == "multivariate"){
    stop("The delta-TOST is currently not implemented in multivariate settings.")
  }

  if (setting == "univariate"){
    # alpha-TOST
    if (method == "alpha"){
      corrected_alpha = get_alpha_TOST(alpha = alpha, sigma = sigma, nu = nu, delta = delta)$root
      decision = abs(theta) < (delta - qt(1 - corrected_alpha, df = nu) * sigma)
      ci = theta + c(-1, 1) * qt(1 - corrected_alpha, df = nu) * sigma
      out = list(decision = decision, ci = ci, theta = theta,
                 sigma = sigma, nu = nu, alpha = alpha,
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
                 sigma = sigma, nu = nu, alpha = alpha,
                 corrected_delta = corrected_delta,
                 delta = delta, method = "delta-TOST",
                 setting = setting)
      class(out) = "tost"
      return(out)
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


#' @title Power function
#' @author Younes Boulaguiem, Stéphane Guerrier, Dominique-Laurent Couturier, Luca Insolia
#' @param alpha                 A \code{numeric} value specifying the significance level (default = \code{0.05}).
#' @param theta                 A \code{numeric} value corresponding to the estimated parameter of interest (such as a difference of means).
#' @param sigma_nu              A \code{numeric} value corresponding to the estimated standard error.
#' @param nu                    A \code{numeric} value corresponding to the number of degrees of freedom.
#' @param delta                 A \code{numeric} value corresponding to (bio)equivalence limit. We assume symmetry, i.e, the (bio)equivalence interval corresponds to (-delta,delta)
#' @param ...                   Additional parameters.
#' @keywords internal
#' @return The function returns a \code{numeric} value that corresponds to a probability.
power_TOST = function(alpha, theta, sigma_nu, nu, delta, ...){
  tval = qt(1 - alpha, df = nu)
  mu1 = (theta + delta)/sigma_nu
  mu2 = (theta - delta)/sigma_nu
  R = (delta*sqrt(nu))/(tval*sigma_nu)
  p1 = OwensQ(nu, tval, mu1, 0, R)
  p2 = OwensQ(nu, -tval, mu2, 0, R)
  pw = p2-p1
  pw[pw < 0] = 0
  pw
}

#' @title The size
#' @author Younes Boulaguiem, Stéphane Guerrier, Dominique-Laurent Couturier
#' @param alpha                 A \code{numeric} value specifying the significance level (default = \code{0.05}).
#' @param sigma_nu              A \code{numeric} value corresponding to the estimated standard error.
#' @param nu                    A \code{numeric} value corresponding to the number of degrees of freedom.
#' @param delta                 A \code{numeric} value corresponding to (bio)equivalence limit. We assume symmetry, i.e, the (bio)equivalence interval corresponds to (-delta,delta)
#' @param ...                   Additional parameters.
#' @keywords internal
#' @return The function returns a \code{numeric} value that corresponds to a probability.
size_TOST = function(alpha, sigma_nu, nu, delta, ...){
  power_TOST(alpha = alpha, theta = delta, sigma_nu = sigma_nu,
             nu = nu, delta = delta)
}

#' @title Objective function to optimise
#' @author Younes Boulaguiem, Stéphane Guerrier, Dominique-Laurent Couturier
#' @param test                  A \code{numeric} value specifying the significance level to optimise.
#' @param alpha                 A \code{numeric} value specifying the significance level (default = \code{0.05}).
#' @param sigma_nu              A \code{numeric} value corresponding to the estimated standard error.
#' @param nu                    A \code{numeric} value corresponding to the number of degrees of freedom.
#' @param delta                 A \code{numeric} value corresponding to (bio)equivalence limit. We assume symmetry, i.e, the (bio)equivalence interval corresponds to (-delta,delta)
#' @param ...                   Additional parameters.
#' @keywords internal
#' @return The function returns a \code{numeric} value for the objective function.
obj_fun_alpha_star = function(test, alpha, sigma_nu, nu, delta, ...){
  size = size_TOST(alpha = test, sigma_nu = sigma_nu, nu = nu, delta = delta)
  (size - alpha)^2
}

#' @title Get alpha star
#' @author Younes Boulaguiem, Stéphane Guerrier, Dominique-Laurent Couturier
#' @param alpha                 A \code{numeric} value specifying the significance level (default = \code{0.05}).
#' @param sigma_nu              A \code{numeric} value corresponding to the estimated standard error.
#' @param nu                    A \code{numeric} value corresponding to the number of degrees of freedom.
#' @param delta                 A \code{numeric} value corresponding to (bio)equivalence limit. We assume symmetry, i.e, the (bio)equivalence interval corresponds to (-delta,delta)
#' @param ...                   Additional parameters.
#' @keywords internal
#' @return The function returns a \code{numeric} value that corresponds to the solution of the optimisation.
get_alpha_star = function(alpha, sigma_nu, nu, delta, ...){
  out = optimize(obj_fun_alpha_star, c(alpha, 0.5),
                 alpha = alpha, sigma_nu = sigma_nu,
                 nu = nu, delta = delta)

  # size_out = size_TOST(alpha = out$minimum, sigma_nu = sigma_nu, nu = nu, delta = delta)

  out$minimum
}

#' @title Confidence Intervals
#' @author Younes Boulaguiem, Stéphane Guerrier, Dominique-Laurent Couturier
#' @param alpha                 A \code{numeric} value specifying the significance level (default = \code{0.05}).
#' @param theta                 A \code{numeric} value corresponding to the estimated parameter of interest (such as a difference of means).
#' @param sigma_nu              A \code{numeric} value corresponding to the estimated standard error.
#' @param nu                    A \code{numeric} value corresponding to the number of degrees of freedom.
#' @param ...                   Additional parameters.
#' @keywords internal
#' @return The function returns a numerical \code{vector} with the lower and upper bound of the confidence intervals.
ci = function(alpha, theta, sigma_nu, nu, ...){
  tval  = qt(p=alpha,df=nu)
  lower = theta+tval*sigma_nu
  upper = theta-tval*sigma_nu
  cbind(lower,upper)
}


#' @title Two One-Sided Tests (TOST) for (Bio)Equivalence Assessment
#'
#' @description
#' Performs the Two One-Sided Tests (TOST) procedure for (bio)equivalence assessment in both univariate and multivariate settings.
#'
#' @param theta      A \code{numeric} value or vector representing the estimated difference(s) (e.g., between a generic and reference product).
#' @param sigma      A \code{numeric} value (univariate) or \code{matrix} (multivariate) corresponding to the estimated variance of the estimate \code{theta}.
#' @param nu         A \code{numeric} value specifying the degrees of freedom. In the multivariate case, it is assumed to be the same across all dimensions.
#' @param delta      A \code{numeric} value or vector defining the (bio)equivalence margin(s). The procedure assumes symmetry, i.e., the (bio)equivalence region is \eqn{(-\delta, \delta)}. In the multivariate case, it is assumed to be the same across all dimensions.
#' @param alpha      A \code{numeric} value specifying the significance level (default is 0.05).
#' @param ...        Additional arguments.
#'
#' @author Younes Boulaguiem, Stéphane Guerrier, Dominique-Laurent Couturier, Luca Insolia
#'
#' @return A \code{tost} object with the following elements:
#' \itemize{
#'   \item \code{decision}:    A logical indicating whether (bio)equivalence is accepted.
#'   \item \code{ci}:          Confidence region at the \eqn{1 - 2\alpha} level.
#'   \item \code{theta}:       The estimated difference(s) used in the test.
#'   \item \code{sigma}:       The estimated variance of \code{theta}, a \code{numeric} in univariate settings or \code{matrix} in multivariate settings.
#'   \item \code{nu}:          Degrees of freedom.
#'   \item \code{alpha}:       Significance level.
#'   \item \code{delta}:       (Bio)equivalence margin(s).
#'   \item \code{method}:      A character string describing the method used ("TOST").
#'   \item \code{setting}:     A character string indicating the setting ("univariate" or "multivariate").
#' }
#'
#' @examples
#' # Univariate case
#' data(skin)
#' theta_hat <- diff(colMeans(skin))
#' nu <- nrow(skin) - 1
#' sig_hat <- sd(apply(skin, 1, diff)) / sqrt(nu)
#' tost(theta = theta_hat, sigma = sig_hat, nu = nu,
#'      alpha = 0.05, delta = log(1.25))
#'
#' # Multivariate case
#' data(ticlopidine)
#' n = nrow(ticlopidine)
#' nu = n-1
#' theta_hat = colMeans(ticlopidine)
#' Sigma_hat = cov(ticlopidine)/n
#' tost(theta = theta_hat, sigma = Sigma_hat, nu = nu, delta = log(1.25))
#'
#' @export
tost = function(theta, sigma, nu, delta, alpha = 0.05,...){

  n_theta = length(theta)

  if (n_theta == 1){
    # Univariate setting
    setting = "univariate"
    if (length(sigma) > 1 || length(delta) > 1){
      stop("sigma and delta must be scalars in univariate settings.")
    }

  }else{
    # Multivariate setting
    setting = "multivariate"

    if (!is.matrix(sigma) || ncol(sigma) != nrow(sigma)){
      stop("sigma must be a square matrix.")
    }

    if (length(delta) > 1){
      stop("delta is assumed to be a scalar, implying that we consider the alternative theta in (-delta, delta) in each dimension.")
    }
  }

  if (setting == "univariate"){
    decision = abs(theta) < (delta - qt(1 - alpha, df = nu) * sigma)
    ci = theta + c(-1, 1) * qt(1 - alpha, df = nu) * sigma
    out = list(decision = as.vector(decision), ci = ci, theta = theta,
               sigma = sigma, nu = nu, alpha = alpha,
               delta = delta, method = "TOST", setting = setting)
    class(out) = "tost"
    return(out)
  }else{
    delta_vec = rep(delta, ncol(sigma))
    t_alpha = qt(1 - alpha, df=nu)
    lower = theta - t_alpha * sqrt(diag(sigma))
    upper = theta + t_alpha * sqrt(diag(sigma))
    decision = (lower > -delta_vec) & (upper < delta_vec)

    out = list(decision = decision, ci = cbind(lower, upper), theta = theta,
               sigma = sigma, nu = nu, alpha = alpha,
               delta = delta, method = "TOST", setting = setting)
    class(out) = "mtost"
    return(out)
  }
}


#' @title The alpha-TOST Corrective Procedure for (Bio)Equivalence Testing
#'
#' @description This functions is used to compute the alpha-TOST, a corrective procedure of the significance level applied to the Two One-Sided Test (TOST) for (bio)equivalence testing in the univariate framework.
#'
#' @param theta                 A \code{numeric} value corresponding to the difference of means.
#' @param sigma                 A \code{numeric} value corresponding to the standard error.
#' @param nu                    A \code{numeric} value corresponding to the number of degrees of freedom.
#' @param alpha                 A \code{numeric} value specifying the significance level.
#' @param delta                 A \code{numeric} value corresponding to (bio)equivalence limit. We assume symmetry, i.e, the (bio)equivalence interval corresponds to (-delta,delta)
#'
#' @author Younes Boulaguiem, Stéphane Guerrier, Dominique-Laurent Couturier
#'
#' @return A \code{tost} object with the structure:
#' \itemize{
#'  \item decision:    A boolean variable indicating whether (bio)equivalence is accepted or not.
#'  \item ci:          Confidence interval at the 1 - 2*alpha level.
#'  \item theta:       The difference of means used for the test.
#'  \item sigma:       The standard error used for the test.
#'  \item nu:          The number of degrees of freedom used for the test.
#'  \item alpha:       The significance level used for the test.
#'  \item delta:       The (bio)equivalence limits used for the test.
#'  \item method:      The method used for the test (here the "alpha-TOST").
#' }
#' @keywords internal
#'
#' @examples
#' data(skin)
#'
#' theta_hat = diff(apply(skin,2,mean))
#' nu = nrow(skin) - 1
#' sig_hat = sd(apply(skin,1,diff))/sqrt(nu)
#' res_atost = cTOST:::atost(theta = theta_hat, sigma = sig_hat, nu = nu,
#'               alpha = 0.05, delta = log(1.25))
#' compare_to_tost(res_atost)
atost = function(theta, sigma, nu, alpha, delta){
  corrected_alpha = alphahat.fun(sigma = sigma, nu = nu, alpha = alpha, delta = delta)
  decision = abs(theta) < (delta - qt(1 - corrected_alpha, df = nu) * sigma)
  ci = theta + c(-1, 1) * qt(1 - corrected_alpha, df = nu) * sigma
  out = list(decision = decision, ci = ci, theta = theta,
             sigma = sigma, nu = nu, alpha = alpha,
             corrected_alpha = corrected_alpha,
             delta = delta, method = "alpha-TOST")
  class(out) = "tost"
  out
}

#' @title The delta-TOST Corrective Procedure for (Bio)Equivalence Testing
#'
#' @description This functions is used to compute the delta-TOST, a corrective procedure of the (bio)equivalence bounds applied to the Two One-Sided Test (TOST) for (bio)equivalence testing in the univariate framework.
#'
#' @param theta                 A \code{numeric} value corresponding to the difference of means.
#' @param sigma                 A \code{numeric} value corresponding to the standard error.
#' @param nu                    A \code{numeric} value corresponding to the number of degrees of freedom.
#' @param alpha                 A \code{numeric} value specifying the significance level.
#' @param delta                 A \code{numeric} value corresponding to (bio)equivalence limit. We assume symmetry, i.e, the (bio)equivalence interval corresponds to (-delta,delta)
#'
#' @author Younes Boulaguiem, Stéphane Guerrier, Dominique-Laurent Couturier
#'
#' @return A \code{tost} object with the structure:
#' \itemize{
#'  \item decision:    A boolean variable indicating whether (bio)equivalence is accepted or not.
#'  \item ci:          Confidence interval at the 1 - 2*alpha level.
#'  \item theta:       The difference of means used for the test.
#'  \item sigma:       The standard error used for the test.
#'  \item nu:          The number of degrees of freedom used for the test.
#'  \item alpha:       The significance level used for the test.
#'  \item delta:       The (bio)equivalence limits used for the test.
#'  \item method:      The method used for the test (here the "delta-TOST").
#' }
#'
#' @keywords internal
#'
#' @examples
#' data(skin)
#'
#' theta_hat = diff(apply(skin,2,mean))
#' nu = nrow(skin) - 1
#' sig_hat = sd(apply(skin,1,diff))/sqrt(nu)
#' res_dtost = cTOST:::dtost(theta = theta_hat, sigma = sig_hat, nu = nu,
#'               alpha = 0.05, delta = log(1.25))
#' compare_to_tost(res_dtost)
#'
dtost = function(theta, sigma, nu, alpha, delta){
  corrected_delta = deltahat.fun(sigma = sigma, alpha = alpha, delta = delta, nu = nu)
  decision = abs(theta) < (corrected_delta - qt(1 - alpha, df = nu) * sigma)
  ci = theta + c(-1, 1) * qt(1 - alpha, df = nu) * sigma
  out = list(decision = decision, ci = ci, theta = theta,
             sigma = sigma, nu = nu, alpha = alpha,
             corrected_delta = corrected_delta,
             delta = delta, method = "delta-TOST")
  class(out) = "tost"
  out
}

