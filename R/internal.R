#' Get Corrected (Bio)Equivalence Bounds
#'
#' This function applies the  delta-TOST corrective procedure to obtain the corrected (bio)equivalence bounds
#'
#' @param alpha                 A \code{numeric} value specifying the significance level.
#' @param sigma                 A \code{numeric} value corresponding to the standard error.
#' @param nu                    A \code{numeric} value corresponding to the number of degrees of freedom.
#' @param delta                 A \code{numeric} value corresponding to (bio)equivalence limit. We assume symmetry, i.e, the (bio)equivalence interval corresponds to (-delta,delta)
#'
#'
#' @keywords internal
#' @importFrom stats uniroot
#'
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
#' @param test The (bio)equivalence bound parameter to be optimized.
#' @param alpha The nominal level for the test.
#' @param sigma The considered standard error.
#' @param nu The degrees of freedom parameter.
#' @param delta The (bio)equivalence bound used for the TOST decision.
#' @param theta    Parameter value (typically equivalence bounds)
#' @param ...      Additional parameters.
#'
#' @keywords internal
#'
obj_fun_delta_TOST = function(test, alpha, sigma, nu, delta, theta=NULL, ...){
  if(is.null(theta)) theta = delta
  size = power_TOST(alpha = alpha, theta = theta, sigma = sigma,
                    nu = nu, delta = test)
  (size - alpha)
}


#' @title Get Corrected Level
#'
#' @description This function applies the alpha-TOST corrective procedure to obtain the corrected level.
#'
#' @param alpha                 A \code{numeric} value specifying the significance level.
#' @param sigma                 A \code{numeric} value corresponding to the standard error.
#' @param nu                    A \code{numeric} value corresponding to the number of degrees of freedom.
#' @param delta                 A \code{numeric} value corresponding to (bio)equivalence limit. We assume symmetry, i.e, the (bio)equivalence interval corresponds to (-delta,delta)
#' @param tol                   A \code{numeric} value corresponding to the tolerance to be applied during the optimization (see `optim`)
#' @param ...                   Additional parameters.
#'
#' @keywords internal
#' @importFrom stats qt uniroot
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


#' Objective Function of the delta-TOST Corrective Procedure
#'
#' @param test     The (test) level to be evaluated.
#' @param alpha    The nominal level for the test.
#' @param sigma    The considered standard error.
#' @param nu       The degrees of freedom parameter.
#' @param delta    The (bio)equivalence bound used for the TOST decision.
#' @param theta    Parameter value (typically equivalence bound)
#' @param ...      Additional parameters.
#'
#' @keywords internal
#'
obj_fun_alpha_TOST = function(test, alpha, sigma, nu, delta, theta = NULL, ...){
  if(is.null(theta)) theta = delta
  size = power_TOST(alpha = test, theta = theta, sigma = sigma,
                    nu = nu, delta = delta)
  (size - alpha)
}

