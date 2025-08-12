#' @title Power Function of Univariate TOST
#'
#' @description
#' Computes the power function for the univariate Two One-Sided Tests (TOST) procedure.
#'
#' @param alpha A \code{numeric} value specifying the significance level.
#' @param theta A \code{numeric} value representing the parameter of interest (e.g., a difference of means).
#' @param sigma A \code{numeric} value representing the estimated standard error of \code{theta}.
#' @param nu A \code{numeric} value specifying the degrees of freedom.
#' @param delta A \code{numeric} value defining the (bio)equivalence margin. The procedure assumes symmetry, i.e., the (bio)equivalence region is \eqn{(-\delta, \delta)}.
#' @param ... Additional parameters.
#'
#' @keywords internal
#'
#' @return A \code{numeric} value corresponding to the probability (power) of the TOST procedure.
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

#' @title Size of Univariate TOST
#'
#' @description
#' Computes the size (type I error rate) of the univariate Two One-Sided Tests (TOST) procedure.
#'
#' @param alpha A \code{numeric} value specifying the significance level.
#' @param sigma A \code{numeric} value representing the estimated standard error of \code{theta}.
#' @param nu A \code{numeric} value specifying the degrees of freedom.
#' @param delta A \code{numeric} value defining the (bio)equivalence margin. The procedure assumes symmetry, i.e., the (bio)equivalence region is \eqn{(-\delta, \delta)}.
#' @param ... Additional parameters.
#'
#' @keywords internal
#'
#' @return A \code{numeric} value corresponding to the probability (size) of the TOST procedure.
size_TOST = function(alpha, sigma, nu, delta, ...){
  power_TOST(alpha = alpha, theta = delta, sigma = sigma,
             nu = nu, delta = delta)
}

#' @title Objective Function for Optimization in Univariate TOST
#'
#' @description
#' Computes the objective function used to optimize the significance level in the univariate Two One-Sided Tests (TOST) procedure.
#'
#' @param test A \code{numeric} value specifying the significance level to be optimized.
#' @param alpha A \code{numeric} value specifying the target significance level (default: \code{alpha = 0.05}).
#' @param sigma A \code{numeric} value representing the estimated standard error of \code{theta}.
#' @param nu A \code{numeric} value specifying the degrees of freedom.
#' @param delta A \code{numeric} value defining the (bio)equivalence margin. The procedure assumes symmetry, i.e., the (bio)equivalence region is \eqn{(-\delta, \delta)}.
#' @param ... Additional parameters.
#'
#' @keywords internal
#'
#' @return A \code{numeric} value representing the objective function to be minimized.
obj_fun_alpha_star = function(test, alpha = 0.05, sigma, nu, delta, ...){
  size = size_TOST(alpha = test, sigma = sigma, nu = nu, delta = delta)
  (size - alpha)^2
}

#' @title Compute Adjusted Significance Level (Alpha Star) for Univariate TOST
#'
#' @description
#' Calculates the adjusted significance level (\eqn{\alpha^*}) for the univariate Two One-Sided Tests (TOST) procedure.
#'
#' @param alpha A \code{numeric} value specifying the target significance level (default: \code{alpha = 0.05}).
#' @param sigma A \code{numeric} value representing the estimated standard error of \code{theta}.
#' @param nu A \code{numeric} value specifying the degrees of freedom.
#' @param delta A \code{numeric} value defining the (bio)equivalence margin. The procedure assumes symmetry, i.e., the (bio)equivalence region is \eqn{(-\delta, \delta)}.
#' @param ... Additional parameters.
#'
#' @keywords internal
#'
#' @return A \code{numeric} value corresponding to the optimized (adjusted) significance level.
get_alpha_star = function(alpha=0.05, sigma, nu, delta, ...){
  out = optimize(obj_fun_alpha_star, c(alpha, 0.5),
                 alpha = alpha, sigma = sigma,
                 nu = nu, delta = delta)

  # size_out = size_TOST(alpha = out$minimum, sigma_nu = sigma_nu, nu = nu, delta = delta)
  out$minimum
}

#' @title Confidence Interval for Univariate TOST
#'
#' @description
#' Computes the confidence interval for the univariate Two One-Sided Tests (TOST) procedure.
#'
#' @param alpha A \code{numeric} value specifying the significance level.
#' @param theta A \code{numeric} value representing the estimated parameter of interest (e.g., a difference of means).
#' @param sigma A \code{numeric} value representing the estimated standard error of \code{theta}.
#' @param nu A \code{numeric} value specifying the degrees of freedom.
#' @param ... Additional parameters.
#'
#' @keywords internal
#'
#' @return A numeric \code{vector} containing the lower and upper bounds of the confidence interval.
ci = function(alpha, theta, sigma, nu, ...){
  tval  = qt(p=alpha,df=nu)
  lower = theta+tval*sigma
  upper = theta-tval*sigma
  cbind(lower,upper)
}


#' @title Two One-Sided Tests (TOST) for (Bio)Equivalence Assessment
#'
#' @description
#' Performs the Two One-Sided Tests (TOST) procedure for (bio)equivalence assessment in both univariate and multivariate settings.
#'
#' @param theta A \code{numeric} value or vector representing the estimated difference(s) (e.g., between a generic and reference product).
#' @param sigma A \code{numeric} value (univariate) or \code{matrix} (multivariate) corresponding to the estimated variance of \code{theta}.
#' @param nu A \code{numeric} value specifying the degrees of freedom. In the multivariate case, it is assumed to be the same across all dimensions.
#' @param delta A \code{numeric} value or vector defining the (bio)equivalence margin(s). The procedure assumes symmetry, i.e., the (bio)equivalence region is \eqn{(-\delta, \delta)}. In the multivariate case, it is assumed to be the same across all dimensions.
#' @param alpha A \code{numeric} value specifying the significance level (default: \code{alpha = 0.05}).
#' @param ... Additional arguments.
#'
#' @return A \code{tost} object with the following elements:
#' \itemize{
#'   \item \code{decision}: Logical; indicates whether (bio)equivalence is accepted.
#'   \item \code{ci}: Confidence region at the \eqn{1 - 2\alpha} level.
#'   \item \code{theta}: The estimated difference(s) used in the test.
#'   \item \code{sigma}: The estimated variance of \code{theta}; a \code{numeric} value (univariate) or \code{matrix} (multivariate).
#'   \item \code{nu}: The degrees of freedom used in the test.
#'   \item \code{alpha}: The significance level used in the test.
#'   \item \code{delta}: The (bio)equivalence limits used in the test.
#'   \item \code{method}: A character string describing the method used ("TOST").
#'   \item \code{setting}: The setting used ("univariate" or "multivariate").
#' }
#'
#' @keywords internal
#'
#' @examples
#' # Univariate case
#' data(skin)
#' theta_hat = diff(colMeans(skin))
#' nu = nrow(skin) - 1
#' sig_hat = sd(apply(skin, 1, diff)) / sqrt(nu)
#' tost(theta = theta_hat, sigma = sig_hat, nu = nu,
#'      alpha = 0.05, delta = log(1.25))
#'
#' # Multivariate case
#' data(ticlopidine)
#' n = nrow(ticlopidine)
#' nu = n - 1
#' theta_hat = colMeans(ticlopidine)
#' Sigma_hat = cov(ticlopidine) / n
#' tost(theta = theta_hat, sigma = Sigma_hat, nu = nu, delta = log(1.25))
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
#' @description
#' Computes the alpha-TOST, a corrective procedure for the significance level applied to the Two One-Sided Test (TOST) for (bio)equivalence testing in the univariate framework.
#'
#' @param theta A \code{numeric} value corresponding to the estimated parameter of interest (such as a difference of means).
#' @param sigma A \code{numeric} value corresponding to the estimated standard error of \code{theta}.
#' @param nu A \code{numeric} value specifying the degrees of freedom.
#' @param alpha A \code{numeric} value specifying the significance level.
#' @param delta A \code{numeric} value defining the (bio)equivalence margin. The procedure assumes symmetry, i.e., the (bio)equivalence region is \eqn{(-\delta, \delta)}.
#'
#' @return An object of class \code{tost} with the following elements:
#' \itemize{
#'   \item \code{decision}: Logical; indicates whether (bio)equivalence is accepted.
#'   \item \code{ci}: Confidence region at the \eqn{1 - 2\alpha} level.
#'   \item \code{theta}: The estimated difference(s) used in the test.
#'   \item \code{sigma}: The estimated standard error used in the test.
#'   \item \code{nu}: The degrees of freedom used in the test.
#'   \item \code{alpha}: The significance level used in the test.
#'   \item \code{corrected_alpha}: The significance level after adjustment.
#'   \item \code{delta}: The (bio)equivalence limits used in the test.
#'   \item \code{method}: The method used in the test ("alpha-TOST").
#' }
#'
#' @keywords internal
#'
#' @examples
#' data(skin)
#'
#' theta_hat = diff(apply(skin, 2, mean))
#' nu = nrow(skin) - 1
#' sig_hat = sd(apply(skin, 1, diff)) / sqrt(nu)
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
#' @description
#' Computes the delta-TOST, a corrective procedure that adjusts the (bio)equivalence bounds applied to the Two One-Sided Test (TOST) for (bio)equivalence testing in the univariate framework.
#'
#' @param theta A \code{numeric} value corresponding to the estimated parameter of interest (such as a difference of means).
#' @param sigma A \code{numeric} value corresponding to the estimated standard error of \code{theta}.
#' @param nu A \code{numeric} value specifying the degrees of freedom.
#' @param alpha A \code{numeric} value specifying the significance level (default: \code{alpha = 0.05}).
#' @param delta A \code{numeric} value defining the (bio)equivalence margin. The procedure assumes symmetry, i.e., the (bio)equivalence region is \eqn{(-\delta, \delta)}.
#'
#' @return An object of class \code{tost} with the following elements:
#' \itemize{
#'   \item \code{decision}: Logical; indicates whether (bio)equivalence is accepted.
#'   \item \code{ci}: Confidence region at the \eqn{1 - 2\alpha} level.
#'   \item \code{theta}: The estimated difference(s) used in the test.
#'   \item \code{sigma}: The estimated standard error used in the test.
#'   \item \code{nu}: The degrees of freedom used in the test.
#'   \item \code{alpha}: The significance level used in the test.
#'   \item \code{delta}: The (bio)equivalence limits used in the test.
#'   \item \code{corrected_delta}: The (bio)equivalence limits after adjustment.
#'   \item \code{method}: The method used in the test ("delta-TOST").
#' }
#'
#' @keywords internal
#'
#' @examples
#' data(skin)
#'
#' theta_hat = diff(apply(skin, 2, mean))
#' nu = nrow(skin) - 1
#' sig_hat = sd(apply(skin, 1, diff)) / sqrt(nu)
#' res_dtost = cTOST:::dtost(theta = theta_hat, sigma = sig_hat, nu = nu,
#'               alpha = 0.05, delta = log(1.25))
#' compare_to_tost(res_dtost)
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




#' Get Corrected (Bio)Equivalence Bounds
#'
#' This function applies the delta-TOST corrective procedure to obtain the corrected (bio)equivalence bounds
#'
#' @param alpha                 A \code{numeric} value specifying the significance level.
#' @param sigma                 A \code{numeric} value corresponding to the standard error.
#' @param nu                    A \code{numeric} value corresponding to the number of degrees of freedom.
#' @param delta                 A \code{numeric} value corresponding to (bio)equivalence limit. We assume symmetry, i.e, the (bio)equivalence interval corresponds to (-delta,delta).
#' @param l                     A \code{numeric} value corresponding to the upper limit of (bio)equivalence margin to optimize (default: l = \code{100}).
#' @param tol                   A \code{numeric} value specifying a tolerance level (default: tol = \code{.Machine$double.eps}).
#'
#' @keywords internal
#'
#' @importFrom stats uniroot
#'
#' @return The function returns a \code{numeric} vector that minimized the objective function.
get_delta_TOST = function(alpha, sigma, nu, delta, l=100, tol = .Machine$double.eps, ...){
  obj_func_a1=obj_fun_delta_TOST(test=delta,alpha=alpha,sigma=sigma, nu=nu, delta=delta)
  obj_func_al=obj_fun_delta_TOST(test=l,alpha=alpha,sigma=sigma, nu=nu, delta=delta)
  if(obj_func_a1>-tol) return(list(root=delta,f.root=obj_func_a1))
  if(obj_func_al<tol) return(list(root=l,f.root=obj_func_al))

  out = uniroot(obj_fun_delta_TOST, interval=c(delta,l),
                alpha = alpha, sigma=sigma, nu=nu, delta = delta)
  out
}

#' Objective Function of the delta-TOST Corrective Procedure
#'
#' @param test  A \code{numeric} value specifying the significance level to optimize.
#' @param alpha A \code{numeric} value specifying the significance level.
#' @param sigma A \code{numeric} value corresponding to the estimated standard error of estimated \code{theta}.
#' @param nu    A \code{numeric} value specifying the degrees of freedom.
#' @param delta A \code{numeric} value or vector defining the (bio)equivalence margin(s). The procedure assumes symmetry, i.e., the (bio)equivalence region is \eqn{(-\delta, \delta)}.
#' @param theta A \code{numeric} value representing the estimated difference(s) (e.g., between a generic and reference product) or a \code{character} value representing the use of equivalence margin(s) for \eqn{\theta} under \code{NULL} (e.g., \eqn{\theta} = \eqn{(-\delta, \delta)}).
#' @param ...      Additional parameters.
#'
#' @keywords internal
#'
#' @return The function returns a \code{numeric} value for the objective function.
#'
obj_fun_delta_TOST = function(test, alpha, sigma, nu, delta, theta=NULL, ...){
  if(is.null(theta)) theta = delta
  size = power_TOST(alpha = alpha, theta = theta, sigma = sigma,
                    nu = nu, delta = test)
  (size - alpha)
}


#' @title Get alpha star of the alpha-TOST Corrective Procedure
#'
#' @description This function applies the alpha-TOST corrective procedure to obtain the corrected level.
#'
#' @param alpha                 A \code{numeric} value specifying the significance level.
#' @param sigma                 A \code{numeric} value corresponding to the estimated standard error of estimated \code{theta}.
#' @param nu                    A \code{numeric} value corresponding to the number of degrees of freedom.
#' @param delta                 A \code{numeric} value corresponding to (bio)equivalence limit. We assume symmetry, i.e, the (bio)equivalence interval corresponds to (-delta,delta).
#' @param l                     A \code{numeric} value corresponding to the upper limit of the significance level to optimize.
#' @param tol                   A \code{numeric} value specifying a tolerance level (default: tol = \code{.Machine$double.eps}).
#' @param ...                   Additional parameters.
#'
#' @keywords internal
#' @importFrom stats qt uniroot
#'
#' @return A list with at least two components:
#' \itemize{
#'  \item \code{root}:       A \code{numeric} value corresponding to the location of the root.
#'  \item \code{f.root}:     A \code{numeric} value corresponding to the function evaluated at \code{root}.
#' }
get_alpha_TOST = function(alpha, sigma, nu, delta, l=0.5, tol = .Machine$double.eps, ...){
  obj_func_a1=obj_fun_alpha_TOST(test=alpha,alpha=alpha,sigma=sigma, nu=nu, delta=delta)
  obj_func_al=obj_fun_alpha_TOST(test=l,alpha=alpha,sigma=sigma,nu=nu,delta=delta)
  if(obj_func_a1>-tol) return(list(root=alpha,f.root=obj_func_a1)) #size at a=1 close to alpha, so we consider them equal
  if(obj_func_al<tol) return(list(root=l,f.root=obj_func_al)) #size at a=l close to alpha, so we consider them equal

  out = uniroot(obj_fun_alpha_TOST, interval=c(alpha,l),
                alpha = alpha, sigma=sigma,
                nu = nu, delta = delta)
  out
}


#' @title Objective Function of the delta-TOST Corrective Procedure
#'
#' @param test  A \code{numeric} value specifying the significance level to optimize.
#' @param alpha A \code{numeric} value specifying the significance level.
#' @param sigma A \code{numeric} value corresponding to the estimated standard error of estimated \code{theta}.
#' @param nu    A \code{numeric} value corresponding to the number of degrees of freedom.
#' @param delta A \code{numeric} value or vector defining the (bio)equivalence margin(s). The procedure assumes symmetry, i.e., the (bio)equivalence region is \eqn{(-\delta, \delta)}.
#' @param theta A \code{numeric} value representing the estimated difference(s) (e.g., between a generic and reference product) or a \code{character} value representing the use of equivalence margin(s) for \eqn{\theta} under \code{NULL} (e.g., \eqn{\theta} = \eqn{\delta}).
#' @param ...   Additional parameters.
#'
#' @keywords internal
#'
#' @return The function returns a \code{numeric} value for the objective function.

obj_fun_alpha_TOST = function(test, alpha, sigma, nu, delta, theta = NULL, ...){
  if(is.null(theta)) theta = delta
  size = power_TOST(alpha = test, theta = theta, sigma = sigma,
                    nu = nu, delta = delta)
  (size - alpha)
}


#' @title Power function of the xTOST Corrective Procedure
#'
#' @description This function is used to calculate the power by xTOST.
#'
#' @param theta   A \code{character} value specifying the value or vector representing the estimated difference(s).
#' @param sig_hat A \code{numeric} value (univariate) or \code{matrix} (multivariate) corresponding to the estimated variance of estimated \code{theta}.
#' @param delta   A \code{numeric} value or vector defining the (bio)equivalence margin(s). The procedure assumes symmetry, i.e., the (bio)equivalence region is \eqn{(-\delta, \delta)}.
#' @param ...     Additional parameters.
#'
#' @keywords internal
#' @importFrom stats pnorm
#'
#' @return The function returns a \code{numeric} value that corresponds to a probability.
#'
power_xTOST = function(theta, sig_hat, delta, ...){
  pnorm((theta + delta)/sig_hat) - pnorm((theta - delta)/sig_hat)
}

#' @title Size function of the xTOST Corrective Procedure
#'
#' @description
#' This function is used to calculate the size in the xTOST corrective procedure.
#'
#' @param sig_hat A \code{numeric} value (univariate) or \code{matrix} (multivariate) corresponding to the estimated variance of estimated \code{theta}.
#' @param delta   A \code{numeric} value or vector defining the (bio)equivalence margin(s). The procedure assumes symmetry, i.e., the (bio)equivalence region is \eqn{(-\delta, \delta)}.
#' @param delta_star A \code{numeric} value specifying the corrected (bio)equivalence margin(s).
#' @param ... description
#'
#' @keywords internal
#' @importFrom stats qnorm dnorm
#'
#' @return The function returns a \code{numeric} value that corresponds to a probability.
size_xTOST = function(sig_hat, delta, delta_star, ...){
  power_xTOST(theta = delta, sig_hat = sig_hat, delta = delta_star)
}

#' @title Critical value search in the xTOST Corrective Procedure
#'
#' @description This function is used to calculate the critical value such that the size of the xTOST procedure
#' matches the nominal significance level.
#'
#' @param delta   A \code{numeric} value defining the (bio)equivalence margin(s). The procedure assumes symmetry, i.e., the (bio)equivalence is \eqn{(-\delta, \delta)}.
#' @param sigma   A \code{numeric} value corresponding to the estimated variance of estimated \code{theta}.
#' @param alpha   A \code{numeric} value specifying the significance level.
#' @param B       A \code{numeric} value specifying the number of iterations using Newton-Raphson method (default: B = \code{1000}).
#' @param tol     A \code{numeric} value specifying a tolerance level (default: tol = \code{10^(-8)}).
#' @param l       A \code{numeric} value corresponding to the upper limit of the significance level to optimize.
#' @param optim   A \code{character} value representing the method to use (default: optim = \code{"NR"}), see Details below for more information.
#'
#' @details
#' There are two methods available for optimization: Newton-Raphson (using \code{optim} = "NR") and minimization without
#' derivatives (using \code{optim} = "uniroot").
#'
#' @keywords internal
#' @importFrom stats pnorm
#'
#' @return The function returns a \code{list} value with the structure:
#' \itemize{
#'  \item \code{c}:        A numerical variable that corresponds to the estimated critical value.
#'  \item \code{size}:     A numerical variable that corresponds to the size when using the estimated critical value.
#'  \item \code{coverged}: A boolean variable that corresponds to whether Newton-Raphson coverged (only returned if \code{optim = "NR"}).
#'  \item \code{iter}:     A numerical variable that corresponds to the actual iterations used (only returned if \code{optim = "NR"}).
#' }
#'
#'
get_c_of_0 = function(delta, sigma, alpha, B = 1000, tol = 10^(-8), l=1, optim = "NR"){
  c0 = delta
  alpha0 = alpha
  # argzero
  if(optim=="uniroot"){
    c_uniroot_ = uniroot(obj_fun_c_of_0, interval=c(10^-8,l),
                         alpha=alpha0, sigma=sigma, delta=c0,tol=.Machine$double.eps)
    size_uniroot = size_xTOST(sig_hat=sigma,delta=c0,delta_star=c_uniroot_$root)
    out = list(c = c_uniroot_$root, size = size_uniroot)
  }else if(optim == "NR"){
    # Sequence of c's
    m=30
    cte_vect = rep(NA, m)

    # Initial approximation
    c_init = c0 - sigma*qnorm(1 - alpha0)
    if (c_init < 0){
      c_init = c0/2
    }
    cte_vect[1] = c_init

    # Start newton raphson
    for (i in 1:(B-1)){
      delta = (pnorm((c0 + cte_vect[i])/sigma) - pnorm((c0 - cte_vect[i])/sigma) - alpha0)*sigma /
        (dnorm((c0 + cte_vect[i])/sigma) + dnorm((c0 - cte_vect[i])/sigma))
      cte_vect[i+1] = cte_vect[i] - delta

      if (cte_vect[i+1] < tol || cte_vect[i+1] > c0){
        cte_vect[i+1] = c0 - tol
      }

      if (abs(delta) < tol){
        out = list(c = cte_vect[i+1], converged = TRUE, iter = i+1)
        break
      }
    }

    if(i < B-1){
      size_NR = size_xTOST(sig_hat=sigma, delta=c0, delta_star=out$c)
      out = list(c = out$c, size=size_NR, converged = out$converged, iter = out$iter)
    }else{
      size_NR = size_xTOST(sig_hat=sigma, delta=c0, delta_star=cte_vect[B])
      out = list(c = cte_vect[B], size=size_NR, converged = F, iter = B)
    }
  }

  out
}

#' @title Finite Sample Adjusted (Bio)Equivalence Testing of the xTOST Corrective Procedure
#'
#' @description This function is used to compute finite sample corrected version of the multivariate xTOST.
#'
#' @param theta_hat  A \code{numeric} value or vector representing the estimated difference(s) (e.g., between a generic and reference product).
#' @param sig_hat  A \code{matrix} (multivariate) corresponding to the estimated variance of estimated \code{theta}.
#' @param nu         A \code{numeric} value specifying the degrees of freedom. In the multivariate case, it is assumed to be the same across all dimensions.
#' @param alpha      A \code{numeric} value specifying the significance level.
#' @param delta      A \code{numeric} value or vector defining the (bio)equivalence margin(s). The procedure assumes symmetry, i.e., the (bio)equivalence region is \eqn{(-\delta, \delta)}.
#' @param correction A \code{character} value corresponding to the considered correction method, see Details below for more information (default: correction = \code{"no"}).
#' @param B                     A \code{numeric} value specifying the number of Monte Carlo replication (default: B = \code{10^4}).
#' @param seed                  A \code{numeric} value specifying a seed for reproducibility (default: seed = \code{85}).
#'
#' @details
#' #' The correction method = "no" refers to ...., the "bootstrap" refers to ... and the "offline" refers to ...
#'
#' @return An object of class \code{tost} with the structure:
#' \itemize{
#'  \item \code{decision}:    A boolean variable indicating whether (bio)equivalence is accepted or not.
#'  \item \code{ci}:          Confidence region at the \eqn{1 - 2\alpha} level.
#'  \item \code{theta_hat}:   The estimated difference(s) used in the test.
#'  \item \code{sigma}:       The estimated variance of \code{theta}, a \code{numeric} in univariate settings or \code{matrix} in multivariate settings.
#'  \item \code{nu}:          The number of degrees of freedom used in the test.
#'  \item \code{alpha}:       The significance level used in the test.
#'  \item \code{c0}:          The estimated critical value.
#'  \item \code{correction}:  The correction used in the test.
#'  \item \code{corrected_alpha}: The significance level corrected by the adjustment.
#'  \item \code{delta}:       The (bio)equivalence limits used in the test.
#'  \item \code{method}:      The method used in the test (default: method = "cTOST").
#' }
#'
#' @seealso \code{\link{ctost}}
#' @export
#' @examples
#' data(skin)
#'
#' theta_hat = diff(apply(skin,2,mean))
#' nu = nrow(skin) - 1
#' sig_hat = var(apply(skin,1,diff))/nu
#'
#' # x-TOST
#' x_tost = xtost(theta_hat = theta_hat, sig_hat = sig_hat, nu = nu,
#'               alpha = 0.05, delta = log(1.25))
#' x_tost
xtost = function(theta_hat, sig_hat, nu, alpha, delta, correction = "no", B = 10^4, seed = 85){

  if (correction == "bootstrap"){
    res = rep(NA, B)
    for (i in 1:B){
      dat = simulate_data(mu = delta, sigma = sig_hat^2,
                          nu = nu, seed = seed + i)
      c_0_hat = get_c_of_0(delta = delta, sigma = dat$sig_hat, alpha = alpha)
      res[i] = abs(dat$theta_hat) < c_0_hat$c
    }
    correct_alpha = 2*alpha - mean(res)
  }

  if (correction == "offline"){

    index_alpha = which.min(abs(ctost_offline_adj$alphas - alpha))
    index_sigma = which.min(abs(ctost_offline_adj$sigmas - sig_hat))
    index_nu = which.min(abs(ctost_offline_adj$nus - nu))
    correct_alpha = 2*alpha - ctost_offline_adj$tier[index_nu, index_sigma, index_alpha]
    correct_alpha = max(correct_alpha, 1e-6) # CHECK-ME: to avoid close to zero or negatives
    #  if (plot){
    #    correct_alpha_all = 2*alpha - ctost_offline_adj$tier[,,index_alpha]
    #    library(pheatmap)
    #    pheatmap(correct_alpha_all, cluster_rows = FALSE, cluster_cols = FALSE)
    #  }
    correct_alpha
  }

  if (correction == "none"){
    correct_alpha = alpha # i.e. no correction
  }

  c_0_hat = get_c_of_0(delta = delta, sigma = sig_hat, alpha = correct_alpha)
  decision = abs(theta_hat) < c_0_hat$c
  ci_half_length = delta - c_0_hat$c
  ci = theta_hat + c(-1, 1) * ci_half_length
  out = list(decision = decision, ci = ci, theta_hat = theta_hat,
             sigma = sig_hat^2, nu = nu, alpha = alpha,
             c0 = c_0_hat$c,
             correction = correction,
             correct_alpha = correct_alpha,
             delta = delta, method = "cTOST")
  class(out) = "tost"
  out
}

