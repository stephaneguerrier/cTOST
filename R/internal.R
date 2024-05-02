#' Get Corrected (Bio)Equivalence Bounds
#'
#' This function applies the  delta-TOST corrective procedure to obtain the corrected (bio)equivalence bounds
#'
#' @param sigma                 A \code{numeric} value corresponding to the standard error.
#' @param alpha                 A \code{numeric} value specifying the significance level.
#' @param delta                 A \code{numeric} value corresponding to (bio)equivalence limit. We assume symmetry, i.e, the (bio)equivalence interval corresponds to (-delta,delta)
#' @param nu                    A \code{numeric} value corresponding to the number of degrees of freedom.
#'
#' @keywords internal
#' @importFrom stats optimize
#'
deltahat.fun = function(sigma, alpha, delta, nu){
  #ipowen4 = utils::getFromNamespace("ipowen4", "OwenQ")
  # sigma = dfw$sigma.hat[546]; delta=log(1.25);tol=1e-8

  # Check estimated Type I error, i.e. |estimated size - alpha|
  size_error = sqrt(obj_fun_delta_hat(delta_star = delta, sigma = sigma, alpha = alpha,
                                      delta = delta, nu = nu)/10^8)

  if (size_error < 10^(-4)){
    return(delta)
  }else{
    m = 1000
    delta.m = seq(from = delta, to = 15*delta, length.out = m)
    omega.m = rep(NA,m)
    for (i in 1:m){
      delta_star = delta.m[i]
      tval       = qt(1 - alpha, df = nu)
      delta1     = (delta + delta_star)/sigma
      delta2     = (delta - delta_star)/sigma
      R          = (delta_star*sqrt(nu))/(tval*sigma)
      omega.m[i] = ipowen4(nu, tval, -tval, delta1, delta2)
      if(omega.m[i]>(alpha+0.01)){
        break
      }
    }
    above     = delta.m[omega.m>(alpha+0.01)&!is.na(omega.m)]
    below     = delta.m[omega.m<alpha&!is.na(omega.m)]
    #if(length(above)==0&length(below)==0){
    #    return(delta)
    #}else{
    max_delta = ifelse(length(above)==0,1.1*delta,max(1.1*delta, min(above)))
    min_delta = ifelse(length(above)==0,delta,max(delta, max(below)))
    res = optimize(obj_fun_delta_hat, c(min_delta, max_delta), sigma = sigma, alpha = alpha,
                   delta = delta, nu = nu)

    if (sqrt(res$objective/10^8) > 10^(-4)){

      m = 10000
      delta.m = seq(from = delta, to = 15*delta, length.out = m)
      omega.m = rep(NA,m)
      for (i in 1:m){
        delta_star = delta.m[i]
        tval       = qt(1 - alpha, df = nu)
        delta1     = (delta + delta_star)/sigma
        delta2     = (delta - delta_star)/sigma
        R          = (delta_star*sqrt(nu))/(tval*sigma)
        omega.m[i] = ipowen4(nu, tval, -tval, delta1, delta2)
        if(omega.m[i]>(alpha+0.01)){
          break
        }
      }
      above     = delta.m[omega.m>(alpha+0.01)&!is.na(omega.m)]
      below     = delta.m[omega.m<alpha&!is.na(omega.m)]
      max_delta = ifelse(length(above)==0,1.1*delta,max(1.1*delta, min(above)))
      min_delta = ifelse(length(above)==0,delta,max(delta, max(below)))
      res = optimize(obj_fun_delta_hat, c(min_delta, max_delta), sigma = sigma, alpha = alpha,
                     delta = delta, nu = nu)
      if (sqrt(res$objective/10^8) > 10^(-4)){
        return(NA)
      }else{
        return(res$minimum)
      }
    }else{
      return(res$minimum)
    }
    #}
  }
}

#' Objective Function of the delta-TOST Corrective Procedure
#'
#' @param delta_star The (bio)equivalence bound parameter to be optimized.
#' @param sigma The considered standard error.
#' @param alpha The nominal level for the test.
#' @param delta The (bio)equivalence bound used for the TOST decision.
#' @param nu The degrees of freedom parameter.
#' @keywords internal
#'
obj_fun_delta_hat = function(delta_star, sigma, alpha, delta, nu){
  #ipowen4    = utils::getFromNamespace("ipowen4", "OwenQ")
  # delta_star = delta
  tval       = qt(1 - alpha, df = nu)
  delta1     = (delta + delta_star)/sigma
  delta2     = (delta - delta_star)/sigma
  R          = (delta_star*sqrt(nu))/(tval*sigma)
  omega      = ipowen4(nu, tval, -tval, delta1, delta2)
  10^8*(omega - alpha)^2
}

#' @title Get Corrected Level
#'
#' @description This function applies the  alpha-TOST corrective procedure to obtain the corrected level.
#'
#' @param sigma                 A \code{numeric} value corresponding to the standard error.
#' @param nu                    A \code{numeric} value corresponding to the number of degrees of freedom.
#' @param alpha                 A \code{numeric} value specifying the significance level.
#' @param delta                 A \code{numeric} value corresponding to (bio)equivalence limit. We assume symmetry, i.e, the (bio)equivalence interval corresponds to (-delta,delta)
#' @param tol                   A \code{numeric} value corresponding to the tolerance to be applied during the optimization (see `optim`)
#'
#' @author Younes Boulaguiem, Stéphane Guerrier, Dominique-Laurent Couturier
#' @keywords internal
#' @importFrom stats qt
alphahat.fun = function(sigma, nu, alpha, delta, tol=1e-7){
  K = 10000
  alpha.k = c(alpha,rep(NA,K-1))
  for(k in 2:K){
    tval       = qt(1 - alpha.k[k-1], df = nu)
    delta1     = (2*delta)/sigma
    delta2     = 0
    R          = (delta*sqrt(nu))/(tval*sigma)
    #ipowen4    =  utils::getFromNamespace("ipowen4", "OwenQ")
    # NOTE: OwenQ:::powen4 unreliable
    #       OwenQ:::ipowen4 gets very close results to PowerTOST but faster
    omega      = ipowen4(nu, tval, -tval, delta1, delta2)
    alpha.k[k] = min(c(alpha + alpha.k[k-1] - omega,0.5))
    # alpha.k[k] = alpha + alpha.k[k-1] - omega
    if(abs(alpha.k[k]-alpha.k[k-1])<tol){break}
  }
  # out
  ifelse(k==K,NA,alpha.k[k])
}

#' @title ipowen4 function from the OwenQ package
#' @author Stéphane Laurent
#' @keywords internal
ipowen4 <- function(...) {
  asNamespace("OwenQ")$ipowen4(...)
}
