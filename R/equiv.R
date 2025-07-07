#' @title Finite Sample Adjusted (Bio)Equivalence Testing
#'
#' @description This function is used to compute finite sample corrected version of the standard (univariate or multivariate) TOST.

#' @author Younes Boulaguiem, Stéphane Guerrier, Dominique-Laurent Couturier
#'
#' @param theta                 A \code{numeric} value or vector representing the estimated difference(s) (e.g., between a generic and reference product).
#' @param sigma                 A \code{numeric} value (univariate) or \code{matrix} (multivariate) corresponding to the estimated variance of estimated \code{theta}.
#' @param nu                    A \code{numeric} value specifying the degrees of freedom. In the multivariate case, it is assumed to be the same across all dimensions.
#' @param delta                 A \code{numeric} value or vector defining the (bio)equivalence margin(s). The procedure assumes symmetry, i.e., the (bio)equivalence region is \eqn{(-\delta, \delta)}. In the multivariate case, it is assumed to be the same across all dimensions.
#' @param method                A \code{character} value corresponding to the considered finite sample adjustment method, see Details below for more information.
#' @param alpha                 A \code{numeric} value specifying the significance level (default: alpha = \code{0.05}).
#' @param B                     A \code{numeric} value specifying the number of Monte Carlo replication (default: B = \code{10^4}).
#' @param correction            A \code{character} value corresponding to the considered correction method, see Details below for more information (default: correction = \code{"offline"}).
#' @param seed                  A \code{numeric} value specifying a seed for reproducibility (default: seed = \code{101010}).
#' @param ...                   Additional parameters.
#'
#' @Details
#' In univariate settings, three adjustment methods are available: optimal (using method = 'optimal'), alpha-TOST (using method = 'alpha') and delta-TOST (using method = 'delta').
#' The optimal method is introduced in (add ref).
#' The alpha-TOST and delta-TOST are introduced in Boulaguiem et al. (2024, <https://doi.org/10.1002/sim.9993>). The former is a corrective procedure of the significance level applied to the TOST while the latter
#' adjusts the equivalence limits. In general, the alpha-TOST appears to outperform the delta-TOST.
#'
#' In multivariate setting, the only available method is the (multivariate) alpha-TOST (using method = 'alpha') introduced in Boulaguiem et al. (2024, <https://doi.org/10.48550/arXiv.2411.16429>).
#'
#' The correction method = "offline" refers to ....
#'
#' @return An object of class \code{tost} with the structure:
#' \itemize{
#'  \item \code{decision}:    A boolean variable indicating whether (bio)equivalence is accepted or not.
#'  \item \code{ci}:          Confidence region at the \eqn{1 - 2\alpha} level.
#'  \item \code{theta}:       The estimated difference(s) used in the test.
#'  \item \code{sigma}:       The estimated variance of \code{theta}, a \code{numeric} in univariate settings or \code{matrix} in multivariate settings.
#'  \item \code{nu}:          The number of degrees of freedom used in the test.
#'  \item \code{alpha}:       The significance level used in the test.
#'  \item \code{corrected_alpha}: The significance level corrected by the adjustment ("alpha-TOST" used).
#'  \item \code{corrected_delta}: The (bio)equivalence limits corrected by the adjustment ("delta-TOST" used).
#'  \item \code{delta}:       The (bio)equivalence limits used in the test.
#'  \item \code{method}:      The method used in the test (optimal, alpha-TOST and delta-TOST).
#'  \item \code{setting}:     The setting used (univariate or multivariate).
#' }
#' @export
#' @examples
#' data(skin)
#'
#' theta_hat = diff(apply(skin,2,mean))
#' nu = nrow(skin) - 1
#' sig_hat = var(apply(skin,1,diff))/nu
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
ctost = function(theta, sigma, nu, delta, alpha = 0.05, method, B = 10^4, seed = 101010, correction = "offline", ...){

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

  if (!(method %in% c("alpha", "delta", "optimal"))){
    stop("Available methods are 'alpha' (for the alpha-TOST) and 'delta' (for the delta-TOST).")
  }

  if (method == "delta" && setting == "multivariate"){
    stop("The delta-TOST is currently not implemented in multivariate settings.")
  }

  if (method == "optimal" && setting == "multivariate"){
    stop("The cTOST (optimal) is currently not implemented in multivariate settings.")
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
#'
#' @author Younes Boulaguiem, Stéphane Guerrier, Dominique-Laurent Couturier
#'
#' @param alpha                 A \code{numeric} value specifying the significance level.
#' @param theta                 A \code{numeric} value corresponding to the estimated parameter of interest (such as a difference of means).
#' @param sigma                 A \code{numeric} value corresponding to the estimated standard error of estimated \code{theta}.
#' @param nu                    A \code{numeric} value specifying the degrees of freedom.
#' @param delta                 A \code{numeric} value or vector defining the (bio)equivalence margin(s). The procedure assumes symmetry, i.e., the (bio)equivalence region is \eqn{(-\delta, \delta)}.
#' @param ...                   Additional parameters.
#'
#' @keywords internal
#'
#' @return The function returns a \code{numeric} value that corresponds to a probability.
power_TOST = function(alpha, theta, sigma, nu, delta, ...){
  tval = qt(1 - alpha, df = nu)
  mu1 = (theta + delta)/sigma
  mu2 = (theta - delta)/sigma
  R = (delta*sqrt(nu))/(tval*sigma)
  p1 = OwensQ(nu, tval, mu1, 0, R)
  p2 = OwensQ(nu, -tval, mu2, 0, R)
  pw = p2-p1
  pw[pw < 0] = 0
  pw
}

#' @title The size
#'
#' @author Younes Boulaguiem, Stéphane Guerrier, Dominique-Laurent Couturier
#'
#' @param alpha                 A \code{numeric} value specifying the significance level.
#' @param sigma                 A \code{numeric} value corresponding to the estimated standard error of estimated \code{theta}.
#' @param nu                    A \code{numeric} value specifying the degrees of freedom.
#' @param delta                 A \code{numeric} value or vector defining the (bio)equivalence margin(s). The procedure assumes symmetry, i.e., the (bio)equivalence region is \eqn{(-\delta, \delta)}.
#' @param ...                   Additional parameters.
#' @keywords internal
#'
#' @return The function returns a \code{numeric} value that corresponds to a probability.
size_TOST = function(alpha, sigma, nu, delta, ...){
  power_TOST(alpha = alpha, theta = delta, sigma = sigma,
             nu = nu, delta = delta)
}

#' @title Objective function to optimize
#'
#' @author Younes Boulaguiem, Stéphane Guerrier, Dominique-Laurent Couturier
#'
#' @param test                  A \code{numeric} value specifying the significance level to optimize.
#' @param alpha                 A \code{numeric} value specifying the significance level (default: alpha = \code{0.05}).
#' @param sigma                 A \code{numeric} value corresponding to the estimated standard error of estimated \code{theta}.
#' @param nu                    A \code{numeric} value specifying the degrees of freedom.
#' @param delta                 A \code{numeric} value or vector defining the (bio)equivalence margin(s). The procedure assumes symmetry, i.e., the (bio)equivalence region is \eqn{(-\delta, \delta)}.
#' @param ...                   Additional parameters.
#'
#' @keywords internal
#'
#' @return The function returns a \code{numeric} value for the objective function.
obj_fun_alpha_star = function(test, alpha = 0.05, sigma, nu, delta, ...){
  size = size_TOST(alpha = test, sigma = sigma, nu = nu, delta = delta)
  (size - alpha)^2
}

#' @title Get alpha star
#'
#' @author Younes Boulaguiem, Stéphane Guerrier, Dominique-Laurent Couturier
#'
#' @param alpha                 A \code{numeric} value specifying the significance level (default: alpha = \code{0.05}).
#' @param sigma                 A \code{numeric} value corresponding to the estimated standard error of estimated \code{theta}.
#' @param nu                    A \code{numeric} value specifying the degrees of freedom.
#' @param delta                 A \code{numeric} value or vector defining the (bio)equivalence margin(s). The procedure assumes symmetry, i.e., the (bio)equivalence region is \eqn{(-\delta, \delta)}.
#' @param ...                   Additional parameters.
#'
#' @keywords internal
#'
#' @return The function returns a \code{numeric} value that corresponds to the solution of the optimization.
get_alpha_star = function(alpha=0.05, sigma, nu, delta, ...){
  out = optimize(obj_fun_alpha_star, c(alpha, 0.5),
                 alpha = alpha, sigma = sigma,
                 nu = nu, delta = delta)

  # size_out = size_TOST(alpha = out$minimum, sigma_nu = sigma_nu, nu = nu, delta = delta)

  out$minimum
}

#' @title Confidence Intervals
#'
#' @author Younes Boulaguiem, Stéphane Guerrier, Dominique-Laurent Couturier
#'
#' @param alpha                 A \code{numeric} value specifying the significance level.
#' @param theta                 A \code{numeric} value corresponding to the estimated parameter of interest (such as a difference of means).
#' @param sigma                 A \code{numeric} value corresponding to the estimated standard error of estimated \code{theta}.
#' @param nu                    A \code{numeric} value specifying the degrees of freedom.
#' @param ...                   Additional parameters.
#'
#' @keywords internal
#'
#' @return The function returns a numerical \code{vector} with the lower and upper bound of the confidence interval.
ci = function(alpha, theta, sigma, nu, ...){
  tval  = qt(p=alpha,df=nu)
  lower = theta+tval*sigma
  upper = theta-tval*sigma
  cbind(lower,upper)
}


#' @title Two One-Sided Tests (TOST) for (Bio)Equivalence Assessment
#'
#' @author Younes Boulaguiem, Stéphane Guerrier, Dominique-Laurent Couturier
#'
#' @description
#' This function performs the Two One-Sided Tests (TOST) procedure for (bio)equivalence assessment in both univariate and multivariate settings.
#'
#' @param theta                 A \code{numeric} value or vector representing the estimated difference(s) (e.g., between a generic and reference product).
#' @param sigma                 A \code{numeric} value (univariate) or \code{matrix} (multivariate) corresponding to the estimated variance of estimated \code{theta}.
#' @param nu                    A \code{numeric} value specifying the degrees of freedom. In the multivariate case, it is assumed to be the same across all dimensions.
#' @param delta                 A \code{numeric} value or vector defining the (bio)equivalence margin(s). The procedure assumes symmetry, i.e., the (bio)equivalence region is \eqn{(-\delta, \delta)}. In the multivariate case, it is assumed to be the same across all dimensions.
#' @param alpha                 A \code{numeric} value specifying the significance level (default: alpha = \code{0.05}).
#' @param ...        Additional arguments.
#'
#' @return A \code{tost} object with the following elements:
#' \itemize{
#'   \item \code{decision}:    A boolean variable indicating whether (bio)equivalence is accepted or not.
#'   \item \code{ci}:          Confidence region at the \eqn{1 - 2\alpha} level.
#'   \item \code{theta}:       The estimated difference(s) used in the test.
#'   \item \code{sigma}:       The estimated variance of \code{theta}, a \code{numeric} in univariate settings or \code{matrix} in multivariate settings.
#'   \item \code{nu}:          The number of degrees of freedom used in the test.
#'   \item \code{alpha}:       The significance level used in the test.
#'   \item \code{delta}:       The (bio)equivalence limits used in the test.
#'   \item \code{method}:      A character string describing the method used ("TOST").
#'   \item \code{setting}:     The setting used (univariate or multivariate).
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
    # Variance -> standard dev
    sigma = sqrt(sigma)

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
               sigma = sigma^2, nu = nu, alpha = alpha,
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
#' @author Younes Boulaguiem, Stéphane Guerrier, Dominique-Laurent Couturier
#'
#' @description This function is used to compute the alpha-TOST, a corrective procedure of the significance level applied to the Two One-Sided Test (TOST) for (bio)equivalence testing in the univariate framework.

#' @param theta                 A \code{numeric} value corresponding to the estimated parameter of interest (such as a difference of means).
#' @param sigma                 A \code{numeric} value corresponding to the estimated standard error of estimated \code{theta}.
#' @param nu                    A \code{numeric} value specifying the degrees of freedom.
#' @param alpha                 A \code{numeric} value specifying the significance level.
#' @param delta                 A \code{numeric} value or vector defining the (bio)equivalence margin(s). The procedure assumes symmetry, i.e., the (bio)equivalence region is \eqn{(-\delta, \delta)}.
#'
#' @return An object of class \code{tost} with the structure:
#' \itemize{
#'  \item \code{decision}:    A boolean variable indicating whether (bio)equivalence is accepted or not.
#'  \item \code{ci}:          Confidence region at the \eqn{1 - 2\alpha} level.
#'  \item \code{theta}:       The estimated difference(s) used in the test.
#'  \item \code{sigma}:       The estimated standard error used in the test.
#'  \item \code{nu}:          The number of degrees of freedom used in the test.
#'  \item \code{alpha}:       The significance level used in the test.
#'  \item \code{corrected_alpha}:       The significance level corrected by the adjustment.
#'  \item \code{delta}:       The (bio)equivalence limits used in the test.
#'  \item \code{method}:      The method used in the test (here the "alpha-TOST").
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
#' @author Younes Boulaguiem, Stéphane Guerrier, Dominique-Laurent Couturier
#'
#' @description This function is used to compute the delta-TOST, a corrective procedure of the (bio)equivalence bounds applied to the Two One-Sided Test (TOST) for (bio)equivalence testing in the univariate framework.
#'
#' @param theta                 A \code{numeric} value corresponding to the estimated parameter of interest (such as a difference of means).
#' @param sigma                 A \code{numeric} value corresponding to the estimated standard error of estimated \code{theta}.
#' @param nu                    A \code{numeric} value specifying the degrees of freedom.
#' @param alpha                 A \code{numeric} value specifying the significance level (default: alpha = \code{0.05}).
#' @param delta                 A \code{numeric} value or vector defining the (bio)equivalence margin(s). The procedure assumes symmetry, i.e., the (bio)equivalence region is \eqn{(-\delta, \delta)}.
#'
#' @author Younes Boulaguiem, Stéphane Guerrier, Dominique-Laurent Couturier
#'
#' @return An object of class \code{tost} with the structure:
#' \itemize{
#'  \item \code{decision}:    A boolean variable indicating whether (bio)equivalence is accepted or not.
#'  \item \code{ci}:          Confidence region at the \eqn{1 - 2\alpha} level.
#'  \item \code{theta}:       The estimated difference(s) used in the test.
#'  \item \code{sigma}:       The estimated standard error used in the test.
#'  \item \code{nu}:          The number of degrees of freedom used in the test.
#'  \item \code{alpha}:       The significance level used in the test.
#'  \item \code{delta}:       The (bio)equivalence limits used in the test.
#'  \item \code{corrected_delta}: The (bio)equivalence limits corrected by the adjustment.
#'  \item \code{method}:      The method used in the test (here the "delta-TOST").
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

