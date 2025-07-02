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


#' TITLE
#'
#' What this function does
#'
#' @param theta ...
#' @param sig_hat ...
#' @param delta ...
#' @param ...      Additional parameters.
#'
#' @keywords internal
#' @importFrom stats pnorm
#'
power_xTOST = function(theta, sig_hat, delta, ...){
  pnorm((theta + delta)/sig_hat) - pnorm((theta - delta)/sig_hat)
}

#' TITLE
#'
#' What this function does
#'
#' @param sig_hat ...
#' @param delta ...
#' @param delta_star ...
#' @param ... description
#'
#' @keywords internal
#' @importFrom stats qnorm dnorm
#'
size_xTOST = function(sig_hat, delta, delta_star, ...){
  power_xTOST(theta = delta, sig_hat = sig_hat, delta = delta_star)
}

#' TITLE
#'
#' What this function does
#'
#' @param delta ...
#' @param sig_hat ...
#' @param alpha ...
#' @param B description
#' @param tol description
#'
#' @keywords internal
#' @importFrom stats pnorm
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
  }else if(optim=="NR"){
    # Sequence of c's
    m=30
    cte_vect = rep(NA, m)

    # Initial approximaiton
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
#' TO BE DOCUMENTED
#'
#' @param theta_hat The estimated mean.
#' @param sig_hat The estimated standard deviation.
#' @param nu The degrees of freedom parameter.
#' @param alpha The significance level for the test.
#' @param delta The equivalence bound used for the TOST decision.
#' @param correction TO BE DOCUMENTED
#' @param B TO BE DOCUMENTED (numb of bootstrap for correction based on bootstrap)
#' @param seed TO BE DOCUMENTED (also for bootstrap)

#' @return TO BE DOCUMENTED
#'
#' @keywords internal
#'
xtost = function(theta_hat, sig_hat, nu, alpha, delta, correction = "no", B = 10^4, seed = 85){

  if (!(correction %in% c("no", "bootstrap", "offline"))){
    stop("Finite sample correction method not implemented.")
  }

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

  if (correction == "no"){
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

