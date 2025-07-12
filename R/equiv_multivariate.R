require(mvtnorm)

#' @title Power function of univariate or multivariate TOST using Monte Carlo integration
#'
#' @author Younes Boulaguiem, Luca Insolia, Stéphane Guerrier, Dominique-Laurent Couturier
#'
#' @param alpha                 A \code{numeric} value specifying the significance level.
#' @param theta                 A \code{numeric} value or vector representing the estimated difference(s) (e.g., between a generic and reference product).
#' @param Sigma                 A \code{numeric} value (univariate) or \code{matrix} (multivariate) corresponding to the estimated variance of estimated \code{theta}.
#' @param nu                    A \code{numeric} value specifying the degrees of freedom. In the multivariate case, it is assumed to be the same across all dimensions.
#' @param delta                 A \code{numeric} value or vector defining the (bio)equivalence margin(s). The procedure assumes symmetry, i.e., the (bio)equivalence region is \eqn{(-\delta, \delta)}.
#' @param B                     A \code{numeric} value specifying the number of Monte Carlo replication (default: B = \code{10^5}).
#' @param seed                  A \code{numeric} value specifying a seed for reproducibility (default: seed = \code{10^8}).
#'
#' @keywords internal
#'
#' The function returns a \code{list} value with the structure:
#' \itemize{
#'  \item \code{power_univ}:    A numerical vector that corresponds to a probability in univariate setting.
#'  \item \code{power_mult}:    A numerical vector that corresponds to a probability in multivariate setting.
#' }
#'
power_TOST_MC_mv = function(alpha, theta, Sigma, nu, delta,
                            B = 10^5, seed = 10^8){
  p=ncol(Sigma)
  #
  set.seed(seed)
  X = rWishart(B, nu, Sigma)/nu
  #
  set.seed(seed*10)
  Z = rmvnorm(n = B, mean = rep(0, p), sigma = Sigma)
  #
  Theta_star=t(as.vector(theta)+t(Z))
  diags=apply(X,3,diag)
  t_val = qt(1 - alpha, df = nu)
  tmp=sqrt(diags)*t_val
  ubs = Theta_star+t(tmp)
  lbs = abs(Theta_star-t(tmp))
  ubs_eval0=(ubs<delta)
  lbs_eval0=(lbs<delta)
  eval=ubs_eval0*lbs_eval0
  #
  list(power_univ=apply(eval,2,mean),power_mult=mean(apply(eval,1,prod)))
}

# alphaTOST
#' @title Objective function to optimize for multivariate TOST
#'
#' @author Younes Boulaguiem, Luca Insolia, Stéphane Guerrier, Dominique-Laurent Couturier
#'
#' @param test                  A \code{numeric} value specifying the significance level to optimize.
#' @param alpha                 A \code{numeric} value specifying the significance level.
#' @param Sigma                 A \code{numeric} value (univariate) or \code{matrix} (multivariate) corresponding to the estimated variance of estimated \code{theta}.
#' @param nu                    A \code{numeric} value specifying the degrees of freedom. In the multivariate case, it is assumed to be the same across all dimensions.
#' @param delta                 A \code{numeric} value or vector defining the (bio)equivalence margin(s). The procedure assumes symmetry, i.e., the (bio)equivalence region is \eqn{(-\delta, \delta)}.
#' @param theta                 A \code{numeric} value representing the estimated difference(s) (e.g., between a generic and reference product) or a \code{character} value representing the use of equivalence margin(s) for \eqn{\theta} under \code{NULL} (e.g., \eqn{\theta} = \eqn{(-\delta, \delta)}).
#' @param B                     A \code{numeric} value specifying the number of Monte Carlo replication (default: B = \code{10^5}).
#' @param seed                  A \code{numeric} value specifying a seed for reproducibility representing the use of multivariate power or a \code{character} value representing the use of univariate power under \code{NULL}.
#' @param ...                   Additional parameters.
#'
#' @keywords internal
#'
#' @return The function returns a \code{numeric} value for the objective function.
#'
obj_fun_alpha_TOST_MC_mv = function(test, alpha, Sigma, nu, delta, theta=NULL, B=10^5, seed=NULL, ...){
  if(is.null(theta)) theta = delta
  if(is.null(seed)){
    size = power_TOST_MC_mv(alpha = test, theta = theta, Sigma = Sigma,
                            nu = nu, delta = delta, B=B)
  }else{
    size = power_TOST_MC_mv(alpha = test, theta = theta, Sigma = Sigma,
                            nu = nu, delta = delta, B=B, seed=seed)
  }

  res = 10000*(size$power_mult - alpha)^2
  # cat(res,"\t", test ,"\n")
  res
}

#' @title Get alpha star for multivariate TOST or xTOST
#'
#' @author Younes Boulaguiem, Luca Insolia, Stéphane Guerrier, Dominique-Laurent Couturier
#'
#' @param alpha                 A \code{numeric} value specifying the significance level.
#' @param Sigma                 A \code{numeric} value (univariate) or \code{matrix} (multivariate) corresponding to the estimated variance of estimated \code{theta}.
#' @param nu                    A \code{numeric} value specifying the degrees of freedom. In the multivariate case, it is assumed to be the same across all dimensions.
#' @param delta                 A \code{numeric} value or vector defining the (bio)equivalence margin(s). The procedure assumes symmetry, i.e., the (bio)equivalence region is \eqn{(-\delta, \delta)}.
#' @param theta                 A \code{character} value specifying the value or vector representing the estimated difference(s) obtained by \code{argsup_meth} method under \code{NULL}.
#' @param B                     A \code{numeric} value specifying the number of Monte Carlo replication (default: B = \code{10^5}).
#' @param tol                   A \code{numeric} value specifying a tolerance level (default: tol = \code{.Machine$double.eps^0.5}).
#' @param seed                  A \code{character} value representing the use of default seed under \code{NULL}.
#' @param argsup_meth           A \code{character} value representing the method to use (default: argsup_meth = \code{"x"}), see Details below for more information.
#' @param ...                   Additional parameters.
#'
#' @details
#' In multivariate setting, two methods to compute supremum are available: Monte Carlo (using \code{argsup_meth} = "MC") or ... (using \code{argsup_meth} = "x").
#' The former is introduced in Boulaguiem et al. (2024, <https://doi.org/10.48550/arXiv.2411.16429>) and the latter is introduced in ...
#'
#' @keywords internal
#'
#' @return The function returns a \code{numeric} value that corresponds to the solution of the optimization.
get_alpha_TOST_MC_mv_core = function(alpha, Sigma, nu, delta, theta=NULL, B=10^5, tol = .Machine$double.eps^0.5, seed=NULL, argsup_meth="x", ...){
  if(is.null(seed)) seed=10^8
  # theta = find_sup_MC(alpha, Sigma, nu, delta, B=10^4,seed=seed)
  if(is.null(theta)){
    if(argsup_meth=="MC"){
      theta = find_sup_MC(alpha,Sigma,nu,delta,B)
    }else if(argsup_meth=="x"){
      theta = find_sup_x(alpha,Sigma,delta)
    }
  }
  obj_func_alpha=obj_fun_alpha_TOST_MC_mv(test = alpha, alpha=alpha, Sigma=Sigma, nu=nu, delta=delta,theta=theta, B=B, seed=seed)
  obj_func_onehalf=obj_fun_alpha_TOST_MC_mv(test = 0.5, alpha=alpha, Sigma=Sigma, nu=nu, delta=delta,theta=theta, B=B, seed=seed)
  if(obj_func_alpha<tol) return(list(minimum=alpha,objective=obj_func_alpha,theta_sup = theta)) #size at a=1 close to alpha, so we consider them equal
  if(obj_func_onehalf<tol) return(list(minimum=alpha,objective=obj_func_onehalf,theta_sup = theta)) #size at a=l close to alpha, so we consider them equal
  out = optimize(obj_fun_alpha_TOST_MC_mv, interval=c(alpha, 0.5), tol=tol,
                 alpha = alpha, Sigma=Sigma, nu = nu,
                 delta = delta, theta = theta, B=B, seed=seed)
  out$theta_sup = theta
  out
}

#' @title Get alpha star for multivariate TOST
#'
#' @author Younes Boulaguiem, Luca Insolia, Stéphane Guerrier, Dominique-Laurent Couturier
#'
#' @param alpha                 A \code{numeric} value specifying the significance level.
#' @param Sigma                 A \code{numeric} value (univariate) or \code{matrix} (multivariate) corresponding to the estimated variance of estimated \code{theta}.
#' @param nu                    A \code{numeric} value specifying the degrees of freedom. In the multivariate case, it is assumed to be the same across all dimensions.
#' @param delta                 A \code{numeric} value or vector defining the (bio)equivalence margin(s). The procedure assumes symmetry, i.e., the (bio)equivalence region is \eqn{(-\delta, \delta)}.
#' @param theta                 A \code{numeric} value specifying the initial value for the worst-case parameter configuration. If not provided, it's computed via find_sup_x().
#' @param B                     A \code{numeric} value specifying the number of Monte Carlo replication (default: B = \code{10^5}).
#' @param tol                   A \code{numeric} value specifying a tolerance level (default: tol = \code{.Machine$double.eps^0.5}).
#' @param seed                  A \code{character} value representing the use of default seed under \code{NULL}.
#' @param max_iter              A \code{numeric} value specifying a maximum number of iteration.
#' @param tolpower              A \code{numeric} value specifying the power level to optimize, see Details below for more information.
#' @param ...                   Additional parameters.
#'
#' @details Tolerance of power allowed between simulated power and nominal alpha. If not specified, it's automatically computed from the 1st and 99th percentiles of a binomial distribution.
#'
#' @keywords internal
#'
#' The function returns a \code{list} value with the structure:
#' \itemize{
#'  \item \code{min}:        A numerical vector that corresponds to the corrected significance level.
#'  \item \code{theta_sup}:  A numerical vector that corresponds to the evaluated value of supremum.
#'  \item \code{alphas}:     A numerical vector that corresponds to the significance level used in each iteration.
#'  \item \code{theta_sups}: A numerical matrix that corresponds to the \code{theta} values evaluated.
#'  \item \code{powers}:     A numerical vector that corresponds to the estimated power .
#'  \item \code{iter}:       A numerical variable that corresponds to the number of iterations performed.
#'  \item \code{err_power}:  A numerical variable that corresponds to the final absolute error between simulated power and nominal significance level \eqn{\alpha}.
#' }
get_alpha_TOST_MC_mv = function(alpha, Sigma, nu, delta, theta=NULL, B=10^5, tol = .Machine$double.eps^0.5, seed=NULL, max_iter=10, tolpower=NULL, ...){
  if(is.null(tolpower)) tolpower=max(abs(qbinom(c(0.01,0.99),B,alpha)/B-alpha))
  theta_sups = matrix(NA, (max_iter+1), ncol(Sigma))
  if(is.null(theta)) theta_sups[1,]=find_sup_x(alpha,Sigma,delta) else theta_sups[1,]=theta
  alphas = powers = rep(NA, max_iter)
  alphas[1] = alpha
  i=1
  #
  while (i <= max_iter) {
    powers[i] = power_TOST_MC_mv(alpha = alphas[i],
                                 theta = theta_sups[i,], Sigma = Sigma,
                                 nu = nu, delta = delta, B=B)$power_mult
    err_power = abs(powers[i]-alpha)
    if (err_power>tolpower) {
      theta  = theta_sups[i,]
      sol_aTOST = get_alpha_TOST_MC_mv_core(alpha=alpha,
                                            Sigma=Sigma,
                                            nu=nu, delta=delta, theta=theta,
                                            B=B, tol=tol, seed=seed)
      alphas[i+1] = sol_aTOST$min
      theta_sups[i+1,] = find_sup_x(alphas[i+1],Sigma,delta)
      i = i+1
    } else {
      break
    }
  }
  alphas = alphas[!is.na(alphas)]
  powers = powers[!is.na(powers)]
  colnames(theta_sups)=paste0("argsup", 1:ncol(Sigma))
  theta_sups = theta_sups[complete.cases(theta_sups), ]
  if (is.null(nrow(theta_sups))) {
    theta_sup=theta_sups
  } else {
    theta_sup=theta_sups[nrow(theta_sups), ]
  }
  out = list(min=alphas[length(alphas)],
             theta_sup=theta_sup,
             alphas=alphas,
             theta_sups=theta_sups,
             powers=powers,
             iter=i-1,
             err_power=err_power)
  out
}


#' @title Power function for multivariate xTOST procedure
#'
#' @author Younes Boulaguiem, Luca Insolia, Stéphane Guerrier, Dominique-Laurent Couturier
#'
#' @description This function is used to approximate the power by a standardized multivariate normal vector lies within the (bio)equivalence margins.
#'
#' @param theta                 A \code{numeric} value or vector representing the estimated difference(s) (e.g., between a generic and reference product).
#' @param delta                 A \code{numeric} value or vector defining the (bio)equivalence margin(s). The procedure assumes symmetry, i.e., the (bio)equivalence region is \eqn{(-\delta, \delta)}. In the multivariate case, it is assumed to be the same across all dimensions.
#' @param Sigma                 A \code{numeric} value (univariate) or \code{matrix} (multivariate) corresponding to the estimated variance of estimated \code{theta}.
#' @param alpha                 A \code{numeric} value specifying the significance level (default: alpha = \code{0.5}).
#' @param seed                  A \code{numeric} value specifying a seed for reproducibility (default: seed = \code{10^5}).
#'
#' @keywords internal
#'
#' @return The function returns a \code{numeric} value that corresponds to a probability.
#'
power_xTOST_mv = function(theta, delta, Sigma, alpha=1/2, seed=10^5){
  #delta is what we optimise over
  set.seed(seed)
  mu = rep(0, ncol(Sigma))
  Sig_diag = sqrt(diag(Sigma))
  lb = -delta/Sig_diag-theta/Sig_diag + qnorm(1-alpha)
  ub = delta/Sig_diag-theta/Sig_diag - qnorm(1-alpha)
  if (sum(lb>ub)>0){
    return(0)
  } else {
    tmp_cdf_mvtnorm = mvtnorm::pmvnorm(lower=lb,
                                       upper=ub,
                                       mean=mu,corr=cov2cor(Sigma))[1]
    # tmp_cdf_TN = TruncatedNormal::pmvnorm(mu=mu,sigma=cov2cor(Sigma),lb=lb,
    #                                       ub=ub)[1]
    return(tmp_cdf_mvtnorm)
    # return(list(mvtnorm=tmp_cdf_mvtnorm, TruncatedN=tmp_cdf_TN) )
  }
}

#' @title Objective function to optimize for multivariate xTOST
#'
#' @author Younes Boulaguiem, Luca Insolia, Stéphane Guerrier, Dominique-Laurent Couturier
#'
#' @param test                  A \code{numeric} value specifying the significance level to optimize.
#' @param Sigma                 A \code{numeric} value (univariate) or \code{matrix} (multivariate) corresponding to the estimated variance of estimated \code{theta}.
#' @param cte                   A \code{numeric} value specifying ...(? the constraint?)
#' @param delta                 A \code{numeric} value or vector defining the (bio)equivalence margin(s). The procedure assumes symmetry, i.e., the (bio)equivalence region is \eqn{(-\delta, \delta)}.
#' @param alpha                 A \code{numeric} value specifying the significance level.
#'
#' @keywords internal
#'
#' @return The function returns a \code{numeric} value for the objective function.
#'
obj_fun_xTOST_constr_size_mv = function(test, Sigma, cte, delta, alpha){
  tmp_cdf = power_xTOST_mv(cte, test*delta, Sigma, alpha=1/2)
  10^16*(tmp_cdf - alpha)^2
}

#' @title Confidence interval search in multivariate xTOST
#'
#' @description This function performs a grid search over a specified interval to
#' find the interval where the size is maximized.
#'
#' @author Younes Boulaguiem, Luca Insolia, Stéphane Guerrier, Dominique-Laurent Couturier
#'
#' @param delta                 A \code{numeric} value or vector defining the (bio)equivalence margin(s). The procedure assumes symmetry, i.e., the (bio)equivalence region is \eqn{(-\delta, \delta)}.
#' @param Sigma                 A \code{numeric} value (univariate) or \code{matrix} (multivariate) corresponding to the estimated variance of estimated \code{theta}.
#' @param cte                   A \code{numeric} value specifying ... (the constraint?).
#' @param alpha                 A \code{numeric} value specifying the significance level.
#' @param interval              A \code{numeric} vector specifying range of interest to search.
#' @param B                     A \code{numeric} value specifying the number of point to search.
#'
#' @keywords internal
#'
#' @return The function returns a \code{numeric} vector minimized the objective function.

get_interval = function(delta, Sigma, cte, alpha, interval, B){
  # runif(B, interval[1], interval[2])
  rint = seq(interval[1], interval[2],l=B)
  res_b = rep(NA, B)
  for (b in 1:B) {
    res_b[b] = obj_fun_xTOST_constr_size_mv(rint[b], Sigma, cte, delta, alpha)
  }
  rint[which.min(res_b)]
}


#' @title Size in multivariate xTOST
#'
#' @description This function is used to get the solution that maximizes the supremum for a constrained xTOST procedure.
#'
#' @author Younes Boulaguiem, Luca Insolia, Stéphane Guerrier, Dominique-Laurent Couturier
#'
#' @param alpha                 A \code{numeric} value specifying the significance level.
#' @param Sigma                 A \code{numeric} value (univariate) or \code{matrix} (multivariate) corresponding to the estimated variance of estimated \code{theta}.
#' @param delta                 A \code{numeric} value or vector defining the (bio)equivalence margin(s). The procedure assumes symmetry, i.e., the (bio)equivalence region is \eqn{(-\delta, \delta)}.
#' @param seed                  A \code{numeric} value specifying a seed for reproducibility (default: seed = \code{10^5}).
#'
#' @details
#' The optimization method "Brent" is used for 2-dimensional or less problems only. The optimization method "Nelder-Mead" is used for higher dimensions.
#'
#' @keywords internal
#'
#' @return The function returns a \code{numeric} vector that maximized the objective function.
find_sup_x = function(alpha,Sigma,delta,seed=10^5){
  # set.seed(seed)
  p=ncol(Sigma)
  if(p<=2){
    method="Brent"
    lower=0;upper=delta
  }else{
    method="Nelder-Mead"
    lower=-Inf;upper=Inf
  }  # We evaluate the size by putting the boundary at every dim
  inds_ = combn(1:p, p-1)
  par0=rep(log(1.25),p-1)
  ## The solution stays sometimes too long at the same place unnecessarily, this is because during optimization
  # the difference between numerical fluctuation of the objective function at the same point
  # are not smaller than tol. --> solution= reltol.
  tmp_=t(apply(inds_,2, function(inds){
    tmp = optim(par=par0,fn=argsup_x,
                inds=inds,alpha=alpha,Sigma=Sigma,delta=delta,seed=seed,
                method=method, lower=lower,upper=upper,
                control=list(reltol=.Machine$double.eps^0.5))
    c(unlist(tmp[1]),unlist(tmp[2]))

  }))
  tmp_[,ncol(tmp_)]=round(tmp_[,ncol(tmp_)],4)
  ind_max_size=which.min(tmp_[,ncol(tmp_)])
  inds=inds_[,ind_max_size] #these are the indexes of the dimensions we will let vary, we fix to c the one missing
  thetas=rep(log(1.25),p)
  thetas[inds]=as.vector(tmp_[ind_max_size,-ncol(tmp_)])
  return(thetas)
}


#' @title Objective function in multivariate xTOST
#'
#' @description This function is used to get the maximum power for a constrained xTOST procedure.
#'
#' @author Younes Boulaguiem, Luca Insolia, Stéphane Guerrier, Dominique-Laurent Couturier
#'
#' @param theta                 A \code{numeric} value or vector representing the estimated difference(s) (e.g., between a generic and reference product).
#' @param inds                  A \code{numeric} value or vector representing the interested \code{theta}.
#' @param alpha                 A \code{numeric} value specifying the significance level.
#' @param Sigma                 A \code{numeric} value (univariate) or \code{matrix} (multivariate) corresponding to the estimated variance of estimated \code{theta}.
#' @param delta                 A \code{numeric} value or vector defining the (bio)equivalence margin(s). The procedure assumes symmetry, i.e., the (bio)equivalence region is \eqn{(-\delta, \delta)}.
#' @param seed                  A \code{numeric} value specifying a seed for reproducibility (default: seed = \code{10^5}).
#'
#' @details
#' The optimization method "Brent" is used for 2-dimensional or less problems only. The optimization method "Nelder-Mead" is used for higher dimensions.
#'
#' @keywords internal
#'
#' @return The function returns a \code{numeric} value that corresponds to a probability.
argsup_x = function(theta,inds,alpha,Sigma,delta,seed=10^5){
  set.seed(seed)
  p=ncol(Sigma)
  thetas=rep(log(1.25),p)
  thetas[inds]=theta
  # cat(thetas,"\n")
  -power_xTOST_mv(theta = thetas, Sigma = Sigma,
                  delta = delta, alpha=alpha, seed=seed)[1]
}
