#' @title ipowen4 function from the OwenQ package
#' @author Stéphane Laurent
#' @keywords internal
ipowen4 <- function(...) {
  asNamespace("OwenQ")$ipowen4(...)
}

#' @title Power function
#' @author Younes Boulaguiem, Stéphane Guerrier, Dominique-Laurent Couturier
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
  tval=qt(p=alpha,df=nu)
  lower <- theta+tval*sigma_nu
  upper <- theta-tval*sigma_nu
  cbind(lower,upper)
}

#' @title Two One-Sided Test (TOST) for (Bio)Equivalence Testing
#'
#' @description This function performs a Two One-Sided Test (TOST) for (bio)equivalence testing.
#'
#' @param theta                 A \code{numeric} value corresponding to the difference of means (e.g. between a generic and reference drug).
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
#'  \item method:      The method used for the test (here the "TOST").
#' }
#' @examples
#' data(skin)
#'
#' theta_hat = diff(apply(skin,2,mean))
#' nu = nrow(skin) - 1
#' sig_hat = sd(apply(skin,1,diff))/sqrt(nu)
#' tost(theta = theta_hat, sigma = sig_hat, nu = nu,
#'      alpha = 0.05, delta = log(1.25))
#'
#' @export
tost = function(theta, sigma, nu, alpha, delta){
  decision = abs(theta) < (delta - qt(1 - alpha, df = nu) * sigma)
  ci = theta + c(-1, 1) * qt(1 - alpha, df = nu) * sigma
  out = list(decision = as.vector(decision), ci = ci, theta = theta,
             sigma = sigma, nu = nu, alpha = alpha,
             delta = delta, method = "TOST")
  class(out) = "tost"
  out
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
#' @examples
#' data(skin)
#'
#' theta_hat = diff(apply(skin,2,mean))
#' nu = nrow(skin) - 1
#' sig_hat = sd(apply(skin,1,diff))/sqrt(nu)
#' res_atost = atost(theta = theta_hat, sigma = sig_hat, nu = nu,
#'               alpha = 0.05, delta = log(1.25))
#' res_atost
#' compare_to_tost(res_atost)
#'
#' @export
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
#' @examples
#' data(skin)
#'
#' theta_hat = diff(apply(skin,2,mean))
#' nu = nrow(skin) - 1
#' sig_hat = sd(apply(skin,1,diff))/sqrt(nu)
#' res_dtost = dtost(theta = theta_hat, sigma = sig_hat, nu = nu,
#'               alpha = 0.05, delta = log(1.25))
#' res_dtost
#' compare_to_tost(res_dtost)
#'
#' @export
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

#' Print Results of (Bio)Equivalence Assessment
#'
#' @param x      A \code{tost} object, which is the output of one of the following functions `tost`, `atost` or `dtost`.
#' @param ticks  Number of ticks to print the confidence interval in the console.
#' @param rn     Number of digits to consider when printing the results.
#' @param ...    Further arguments to be passed to or from methods.
#' @return       Prints object.
#' @importFrom   cli cli_text col_green col_red
#'
#' @rdname print.tost
#'
#' @export
print.tost = function(x, ticks = 30, rn = 5, ...){

  if (x$decision){
    cli_text(col_green("{symbol$tick} Accept (bio)equivalence"))
  }else{
    cli_text(col_red("{symbol$cross} Can't accept (bio)equivalence"))
  }

  if (x$method == "delta-TOST"){
    lower_be = x$ci[1] > -x$corrected_delta
    upper_be = x$ci[2] < x$corrected_delta
    rg = range(c(x$ci, x$corrected_delta, -x$corrected_delta))
    rg_delta = rg[2] - rg[1]
    std_be_interval = round(ticks*(c(-x$corrected_delta, x$corrected_delta) - rg[1])/rg_delta) + 1
    std_zero = round(-ticks*rg[1]/rg_delta) + 1
    std_fit_interval = round(ticks*(x$ci - rg[1])/rg_delta) + 1
    std_fit_interval_center = round(ticks*(sum(x$ci)/2 - rg[1])/rg_delta) + 1
  }else{
    lower_be = x$ci[1] > -x$delta
    upper_be = x$ci[2] < x$delta
    rg = range(c(x$ci, x$delta, -x$delta))
    rg_delta = rg[2] - rg[1]
    std_be_interval = round(ticks*(c(-x$delta, x$delta) - rg[1])/rg_delta) + 1
    std_zero = round(-ticks*rg[1]/rg_delta) + 1
    std_fit_interval = round(ticks*(x$ci - rg[1])/rg_delta) + 1
    std_fit_interval_center = round(ticks*(sum(x$ci)/2 - rg[1])/rg_delta) + 1
  }

  if (x$method == "delta-TOST"){
    cat("Corr. Equiv. Region:  ")
  }else{
    cat("Equiv. Region:  ")
  }

  for (i in 1:(ticks+1)){
    if (i >= std_be_interval[1] && i <= std_be_interval[2]){
      if (i == std_be_interval[1]){
        cat(("|-"))
      }else{
        if (i ==  std_be_interval[2]){
          cat(("-|"))
        }else{
          if (i == std_zero){
            cat(("-0-"))
          }else{
            cat(("-"))
          }
        }
      }
    }else{
      cat(" ")
    }
  }

  cat("\n")

  if (x$method == "delta-TOST"){
    cat("      Estim. Inter.:  ")
  }else{
    cat("Estim. Inter.:  ")
  }

  for (i in 1:(ticks+1)){
    if (i >= std_fit_interval[1] && i <= std_fit_interval[2]){
      if (i == std_fit_interval[1]){
        if (i > std_be_interval[1] && i < std_be_interval[2]){
          cat(col_green("(-"))
        }else{
          if (lower_be){
            cat(col_green("(-"))
          }else{
            cat(col_red("(-"))
          }
        }

      }else{
        if (i ==  std_fit_interval[2]){
          if (i > std_be_interval[1] && i < std_be_interval[2]){
            cat(col_green("-)"))
          }else{
            if (upper_be){
              cat(col_green("-)"))
            }else{
              cat(col_red("-)"))
            }
          }

        }else{
          if (i == std_fit_interval_center){
            if (i >= std_be_interval[1] && i <= std_be_interval[2]){
              cat(col_green("-x-"))
            }else{
              cat(col_red("-x-"))
            }
          }else{
            if (i >= std_be_interval[1] && i <= std_be_interval[2]){
              cat(col_green("-"))
            }else{
              cat(col_red("-"))
            }
          }
        }
      }
    }else{
      cat(" ")
    }
  }
  cat("\n")
  cat("CI =  (")
  cat(format(round(x$ci[1], rn), nsmall = rn))
  cat(" ; ")
  cat(format(round(x$ci[2], rn), nsmall = rn))
  cat(")\n\n")

  cat("Method: ")
  cat(x$method)
  cat("\n")
  cat("alpha = ")
  cat(x$alpha)
  cat("; ")
  cat("Equiv. lim. = +/- ")
  cat(format(round(x$delta, rn), nsmall = rn))
  cat("\n")
  if (x$method == "alpha-TOST"){
    cat("Corrected alpha = ")
    cat(format(round(x$corrected_alpha, rn), nsmall = rn))
    cat("\n")
  }

  if (x$method == "delta-TOST"){
    cat("Corrected Equiv. lim. = +/- ")
    cat(format(round(x$corrected_delta, rn), nsmall = rn))
    cat("\n")
  }

  if (x$method == "x-TOST"){
    cat("Estimated c(0) = ")
    cat(format(round(x$c0, rn), nsmall = rn))
    cat("\n")

    cat("Finite sample correction: ")
    cat(x$correction)
    cat("\n")
    if (x$correction %in% c("bootstrap", "offline")){
      cat("Corrected alpha = ")
      cat(format(round(x$correct_alpha, rn), nsmall = rn))
      cat("\n")
    }
  }
  cat("Mean = ")
  cat(format(round(x$theta, rn), nsmall = rn))
  cat("; ")
  cat("Stand. dev. = ")
  cat(format(round(x$sigma, rn), nsmall = rn))
  cat("; ")
  cat("df = ")
  cat(x$nu)
  cat("\n")
}

#' @title Comparison of a Corrective Procedure to the results of the Two One-Sided Tests (TOST)
#'
#' @description This function renders a comparison of the alpha-TOST or the delta-TOST outputs obtained with the functions `atost` or `dtost`, respectively, to the TOST output obtained with `tost`.
#'
#' @param x A \code{tost} object, which is the output of one of the following functions: `atost` or `dtost`.
#' @param ticks an integer indicating the number of segments that will be printed to represent the confidence intervals.
#' @param rn integer indicating the number of decimals places to be used (see function `round`) for the printed results.
#' @return Pints a comparison between the TOST results (i.e., output of `tost`) and either the alpha-TOST or the delta-TOST results (i.e., outputs of `atost` or `dtost`, respectively).
#'
#' @importFrom cli cli_text col_green col_red
#'
#' @export
compare_to_tost = function(x, ticks = 30, rn = 5){
  result_tost = tost(theta = x$theta, sigma = x$sigma, nu = x$nu,
                     alpha = x$alpha, delta = x$delta)

  if (!(x$method %in% c("alpha-TOST", "delta-TOST", "x-TOST"))){
    stop("This method is not compatible")
  }

  if (x$method == "delta-TOST"){
    cat("TOST:       ")
    if (result_tost$decision){
      cli_text(col_green("{symbol$tick} Accept (bio)equivalence"))
    }else{
      cli_text(col_red("{symbol$cross} Can't accept (bio)equivalence"))
    }

    cat("delta-TOST: ")
    if (x$decision){
      cli_text(col_green("{symbol$tick} Accept (bio)equivalence"))
    }else{
      cli_text(col_red("{symbol$cross} Can't accept (bio)equivalence"))
    }
    cat("\n")

    lower_be_tost = x$ci[1] > -x$delta
    upper_be_tost = x$ci[2] < x$delta

    lower_be_dtost = x$ci[1] > -x$corrected_delta
    upper_be_dtost = x$ci[2] < x$corrected_delta

    rg = range(c(x$ci, result_tost$ci, x$delta, -x$delta, x$corrected_delta, -x$corrected_delta))
    rg_delta = rg[2] - rg[1]
    std_be_interval = round(ticks*(c(-x$delta, x$delta) - rg[1])/rg_delta) + 1
    std_zero = round(-ticks*rg[1]/rg_delta) + 1

    std_be_interval_cor = round(ticks*(c(-x$corrected_delta, x$corrected_delta) - rg[1])/rg_delta) + 1
    std_zero_cor = round(-ticks*rg[1]/rg_delta) + 1

    std_fit_interval_tost = round(ticks*(result_tost$ci - rg[1])/rg_delta) + 1
    std_fit_interval_center_tost = round(ticks*(sum(result_tost$ci)/2 - rg[1])/rg_delta) + 1

    cat("Stand. Equiv. Region:  ")

    for (i in 1:(ticks+1)){
      if (i >= std_be_interval[1] && i <= std_be_interval[2]){
        if (i == std_be_interval[1]){
          cat(("|-"))
        }else{
          if (i ==  std_be_interval[2]){
            cat(("-|"))
          }else{
            if (i == std_zero){
              cat(("-0-"))
            }else{
              cat(("-"))
            }
          }
        }
      }else{
        cat(" ")
      }
    }
    cat("\n")
    cat(" Corr. Equiv. Region:  ")

    for (i in 1:(ticks+1)){
      if (i >= std_be_interval_cor[1] && i <= std_be_interval_cor[2]){
        if (i == std_be_interval_cor[1]){
          cat(("|-"))
        }else{
          if (i ==  std_be_interval_cor[2]){
            cat(("-|"))
          }else{
            if (i == std_zero_cor){
              cat(("-0-"))
            }else{
              cat(("-"))
            }
          }
        }
      }else{
        cat(" ")
      }
    }
    cat("\n")
    cat("                TOST:  ")
    for (i in 1:(ticks+1)){
      if (i >= std_fit_interval_tost[1] && i <= std_fit_interval_tost[2]){
        if (i == std_fit_interval_tost[1]){
          if (i > std_be_interval[1] && i < std_be_interval[2]){
            cat(col_green("(-"))
          }else{
            if (lower_be_tost){
              cat(col_green("(-"))
            }else{
              cat(col_red("(-"))
            }
          }

        }else{
          if (i ==  std_fit_interval_tost[2]){
            if (i > std_be_interval[1] && i < std_be_interval[2]){
              cat(col_green("-)"))
            }else{
              if(upper_be_tost){
                cat(col_green("-)"))
              }else{
                cat(col_red("-)"))
              }
            }

          }else{
            if (i == std_fit_interval_center_tost){
              if (i >= std_be_interval[1] && i <= std_be_interval[2]){
                cat(col_green("-x-"))
              }else{
                cat(col_red("-x-"))
              }
            }else{
              if (i >= std_be_interval[1] && i <= std_be_interval[2]){
                cat(col_green("-"))
              }else{
                cat(col_red("-"))
              }
            }
          }
        }
      }else{
        cat(" ")
      }
    }
    cat("\n")
    cat("          delta-TOST:  ")
    for (i in 1:(ticks+1)){
      if (i >= std_fit_interval_tost[1] && i <= std_fit_interval_tost[2]){
        if (i == std_fit_interval_tost[1]){
          if (i > std_be_interval_cor[1] && i < std_be_interval_cor[2]){
            cat(col_green("(-"))
          }else{
            if(lower_be_dtost){
              cat(col_green("-)"))
            }else{
              cat(col_red("-)"))
            }
          }

        }else{
          if (i ==  std_fit_interval_tost[2]){
            if (i > std_be_interval_cor[1] && i < std_be_interval_cor[2]){
              cat(col_green("-)"))
            }else{
              if(upper_be_dtost){
                cat(col_green("-)"))
              }else{
                cat(col_red("-)"))
              }
            }

          }else{
            if (i == std_fit_interval_center_tost){
              if (i >= std_be_interval_cor[1] && i <= std_be_interval_cor[2]){
                cat(col_green("-x-"))
              }else{
                cat(col_red("-x-"))
              }
            }else{
              if (i >= std_be_interval_cor[1] && i <= std_be_interval_cor[2]){
                cat(col_green("-"))
              }else{
                cat(col_red("-"))
              }
            }
          }
        }
      }else{
        cat(" ")
      }
    }

    cat("\n")
    cat("\n")
    cat(" Standard Equiv. lim. = +/- ")
    cat(format(round(x$delta, rn), nsmall = rn))
    cat("\n")

    cat("Corrected Equiv. lim. = +/- ")
    cat(format(round(x$corrected_delta, rn), nsmall = rn))
    cat("\n")
  }

  if (x$method %in% c("alpha-TOST", "x-TOST")){

    lower_be_atost = x$ci[1] > -x$delta
    upper_be_atost = x$ci[2] < x$delta
    lower_be_tost = result_tost$ci[1] > -x$delta
    upper_be_tost = result_tost$ci[2] < x$delta

    if (x$method == "alpha-TOST"){
      cat("TOST:       ")
    }else{
      cat("TOST:   ")
    }
    if (result_tost$decision){
      cli_text(col_green("{symbol$tick} Accept (bio)equivalence"))
    }else{
      cli_text(col_red("{symbol$cross} Can't accept (bio)equivalence"))
    }


    if (x$method == "alpha-TOST"){
      cat("alpha-TOST: ")
    }else{
      cat("x-TOST: ")
    }

    if (x$decision){
      cli_text(col_green("{symbol$tick} Accept (bio)equivalence"))
    }else{
      cli_text(col_red("{symbol$cross} Can't accept (bio)equivalence"))
    }
    cat("\n")

    rg = range(c(x$ci, result_tost$ci, x$delta, -x$delta))
    rg_delta = rg[2] - rg[1]
    std_be_interval = round(ticks*(c(-x$delta, x$delta) - rg[1])/rg_delta) + 1
    std_zero = round(-ticks*rg[1]/rg_delta) + 1
    std_fit_interval_atost = round(ticks*(x$ci - rg[1])/rg_delta) + 1
    std_fit_interval_center_atost = round(ticks*(sum(x$ci)/2 - rg[1])/rg_delta) + 1
    std_fit_interval_tost = round(ticks*(result_tost$ci - rg[1])/rg_delta) + 1
    std_fit_interval_center_tost = round(ticks*(sum(result_tost$ci)/2 - rg[1])/rg_delta) + 1

    cat("Equiv. Region:  ")

    for (i in 1:(ticks+1)){
      if (i >= std_be_interval[1] && i <= std_be_interval[2]){
        if (i == std_be_interval[1]){
          cat(("|-"))
        }else{
          if (i ==  std_be_interval[2]){
            cat(("-|"))
          }else{
            if (i == std_zero){
              cat(("-0-"))
            }else{
              cat(("-"))
            }
          }
        }
      }else{
        cat(" ")
      }
    }

    cat("\n")
    cat("TOST:           ")

    for (i in 1:(ticks+1)){
      if (i >= std_fit_interval_tost[1] && i <= std_fit_interval_tost[2]){
        if (i == std_fit_interval_tost[1]){
          if (i > std_be_interval[1] && i < std_be_interval[2]){
            cat(col_green("(-"))
          }else{
            if (lower_be_tost){
              cat(col_green("(-"))
            }else{
              cat(col_red("(-"))
            }
          }

        }else{
          if (i ==  std_fit_interval_tost[2]){
            if (i > std_be_interval[1] && i < std_be_interval[2]){
              cat(col_green("-)"))
            }else{
              if (upper_be_tost){
                cat(col_green("-)"))
              }else{
                cat(col_red("-)"))
              }
            }

          }else{
            if (i == std_fit_interval_center_tost){
              if (i >= std_be_interval[1] && i <= std_be_interval[2]){
                cat(col_green("-x-"))
              }else{
                cat(col_red("-x-"))
              }
            }else{
              if (i >= std_be_interval[1] && i <= std_be_interval[2]){
                cat(col_green("-"))
              }else{
                cat(col_red("-"))
              }
            }
          }
        }
      }else{
        cat(" ")
      }
    }

    cat("\n")

    if (x$method == "alpha-TOST"){
      cat("alpha-TOST:     ")
    }else{
      cat("x-TOST:         ")
    }

    for (i in 1:(ticks+1)){
      if (i >= std_fit_interval_atost[1] && i <= std_fit_interval_atost[2]){
        if (i == std_fit_interval_atost[1]){
          if (i > std_be_interval[1] && i < std_be_interval[2]){
            cat(col_green("(-"))
          }else{
            if (lower_be_atost){
              cat(col_green("(-"))
            }else{
              cat(col_red("(-"))
            }
          }
        }else{
          if (i ==  std_fit_interval_atost[2]){
            if (i > std_be_interval[1] && i < std_be_interval[2]){
              cat(col_green("-)"))
            }else{
              if (upper_be_atost){
                cat(col_green("-)"))
              }else{
                cat(col_red("-)"))
              }
            }
          }else{
            if (i == std_fit_interval_center_atost){
              if (i >= std_be_interval[1] && i <= std_be_interval[2]){
                cat(col_green("-x-"))
              }else{
                cat(col_red("-x-"))
              }
            }else{
              if (i >= std_be_interval[1] && i <= std_be_interval[2]){
                cat(col_green("-"))
              }else{
                cat(col_red("-"))
              }
            }
          }
        }
      }else{
        cat(" ")
      }
    }
    cat("\n")
    cat("\n")
    cat("                 CI - low      ")
    cat("CI - high")
    cat("\n")

    cat("TOST:            ")
    cat(format(round(result_tost$ci[1], rn), nsmall = rn))
    cat("       ")
    cat(format(round(result_tost$ci[2], rn), nsmall = rn))
    cat("\n")

    if (x$method == "alpha-TOST"){
      cat("alpha-TOST:      ")
    }else{
      cat("x-TOST:          ")
    }

    cat(format(round(x$ci[1], rn), nsmall = rn))
    cat("       ")
    cat(format(round(x$ci[2], rn), nsmall = rn))
    cat("\n")

    cat("\n")
    cat("Equiv. lim. = +/- ")
    cat(format(round(x$delta, rn), nsmall = rn))
    cat("\n")
  }
}



