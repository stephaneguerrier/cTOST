#' Extract Summary Statistics
#'
#' @description
#' An internal helper function that processes input data to ensure it is in the
#' correct summary statistics format (mean, standard deviation, and sample size)
#' for equivalence testing on multiple quantiles (see the `qtost` function).
#' If a numeric vector is provided, it calculates these statistics. If a list
#' is provided, it validates that the required statistics are present.
#'
#' @param data A numeric vector, or a list containing the named elements `mean`, `sd`, and `n`.
#' @param group_name A character string representing the name of the group
#'   (e.g., "x" or "y"), used for constructing informative error messages.
#'
#' @return
#' A list containing three named elements:
#' \itemize{
#'   \item `mean`: The mean of the data.
#'   \item `sd`: The standard deviation of the data.
#'   \item `n`: The number of non-missing observations.
#' }
#'
#' @keywords internal
.extract_stats <- function(data, group_name) {
  if (is.numeric(data)) {
    return(list(
      mean = mean(data, na.rm = TRUE),
      sd   = sd(data, na.rm = TRUE),
      n    = length(na.omit(data))
    ))
  } else if (is.list(data) && all(c("mean", "sd", "n") %in% names(data))) {
    return(data)
  } else {
    stop(paste("Input for group '", group_name, "' is not a numeric vector or a valid list.", sep = ""))
  }
}


power_qTOST_MC_biv = function(sol_1, sol_2, delta_l, delta_u, alpha){
  delta_l_1 = delta_l[1]
  delta_u_1 = delta_u[1]
  delta_l_2 = delta_l[2]
  delta_u_2 = delta_u[2]
  theta_hat_1 = sol_1$theta_hat
  theta_hat_2 = sol_2$theta_hat
  sigma_theta_hat_1 = sol_1$sigma_theta_hat
  sigma_theta_hat_2 = sol_2$sigma_theta_hat
  B = length(theta_hat_1)
  b_theta_1 = theta_hat_1 + qnorm(1-alpha)*sigma_theta_hat_1
  a_theta_1 = theta_hat_1 - qnorm(1-alpha)*sigma_theta_hat_1
  b_theta_2 = theta_hat_2 + qnorm(1-alpha)*sigma_theta_hat_2
  a_theta_2 = theta_hat_2 - qnorm(1-alpha)*sigma_theta_hat_2
  tmp = (b_theta_1<qnorm(delta_u_1) & a_theta_1>qnorm(delta_l_1)) &
    (b_theta_2<qnorm(delta_u_2) & a_theta_2>qnorm(delta_l_2))
  res = sum(tmp)/B
  res
}

find_sup_qTOST = function(gamma,n_x,n_y,pi_x,delta_l,delta_u,alpha,B=10^5,seed=12345,side="upper",MC_sup=T){
  p=length(delta_l)
  if(p<=2){
    method="Brent"
    lower=min(qnorm(delta_l))
    upper=max(qnorm(delta_u))
  }else{
    method="Nelder-Mead"
    lower=-Inf;upper=Inf
  }  # We evaluate the size by putting the boundary at every dim
  inds_ = combn(1:p, p-1)
  ## The solution stays sometimes too long at the same place unnecessarily, this is because during optimization
  # the difference between numerical fluctuation of the objective function at the same point
  # are not smaller than tol. --> solution= reltol.
  tmp_=t(apply(inds_,2, function(inds){
    tmp = optim(par=setdiff(1:p,inds),fn=argsup_qTOST,
                inds=inds,gamma=gamma,n_x=n_x,n_y=n_y,
                pi_x=pi_x,delta_l=delta_l,delta_u=delta_u,alpha=alpha,
                B=B,seed=seed, side=side, MC_sup=MC_sup,
                method=method, lower=lower,upper=upper,
                control=list(reltol=.Machine$double.eps^0.5))
    c(unlist(tmp[1]),unlist(tmp[2]))

  }))
  tmp_[,ncol(tmp_)]=round(tmp_[,ncol(tmp_)],4)
  ind_max_size=which.min(tmp_[,ncol(tmp_)])
  inds=inds_[,ind_max_size] #these are the indexes of the dimensions we will let vary, we fix to c the one missing
  if (side=="upper"){ # we try the argsup on both delta_l and delta_u for the element at the boundary
    thetas = qnorm(delta_u)
  } else if (side=="lower") {
    thetas = qnorm(delta_l)
  }
  thetas[inds]=as.vector(tmp_[ind_max_size,-ncol(tmp_)])
  return(thetas)
}

argsup_qTOST = function(theta,inds,gamma,n_x,n_y,pi_x,delta_l,delta_u,alpha,B=1e4,seed=12345,side="upper",MC_sup=T){
  if (side=="upper"){ # we try the argsup on both delta_l and delta_u for the element at the boundary
    thetas = qnorm(delta_u)
  } else if (side=="lower") {
    thetas = qnorm(delta_l)
  } else {
    error("side must be 'upper' or 'lower'")
  }
  thetas[inds] = theta
  l = n_y/n_x
  Sigma = (1 + (thetas%*%t(thetas))/2 + l*(qnorm(pi_x)%*%t(qnorm(pi_x))/2 + 1)/gamma) / n_y
  if (MC_sup){
    sol_1 = get_qTOST_rvs(thetas[1], gamma, n_x, n_y, pi_x[1], B=B, seed=seed)
    sol_2 = get_qTOST_rvs(thetas[2], gamma, n_x, n_y, pi_x[2], B=B, seed=seed)
    res = -power_qTOST_MC_biv(sol_1, sol_2, delta_l, delta_u, alpha)
  } else {
    res = -power_cTOST_mv(theta = thetas, Sigma = Sigma,
                          delta_l=delta_l, delta_u=delta_u, alpha=alpha, seed=seed)[1]
  }
  res
}

tost_sup = function(gamma, n_x, n_y, pi_x, delta_l, delta_u, alpha, B=1e5, seed=12345,MC_sup=T){
  tost_lambda_1 = find_sup_qTOST(gamma,n_x,n_y,pi_x,delta_l,delta_u,alpha,B,seed,side="upper",MC_sup)
  sol_1 = get_qTOST_rvs(tost_lambda_1[1], gamma, n_x, n_y, pi_x[1], B=B, seed=seed)
  sol_2 = get_qTOST_rvs(tost_lambda_1[2], gamma, n_x, n_y, pi_x[2], B=B, seed=seed)
  pw_1 = power_qTOST_MC_biv(sol_1, sol_2, delta_l, delta_u, alpha)
  tost_lambda_2 = find_sup_qTOST(gamma,n_x,n_y,pi_x,delta_l,delta_u,alpha,B,seed,side="lower",MC_sup)
  sol_1 = get_qTOST_rvs(tost_lambda_2[1], gamma, n_x, n_y, pi_x[1], B=B, seed=seed)
  sol_2 = get_qTOST_rvs(tost_lambda_2[2], gamma, n_x, n_y, pi_x[2], B=B, seed=seed)
  pw_2 = power_qTOST_MC_biv(sol_1, sol_2, delta_l, delta_u, alpha)
  if (pw_1 >= pw_2){
    return(tost_lambda_1)
  } else {
    return(tost_lambda_2)
  }
}

power_cTOST_mv = function(theta, Sigma, delta_l, delta_u, alpha=1/2, seed=12345){
  set.seed(seed)
  Sig_diag = sqrt(diag(Sigma))
  lb = qnorm(delta_l)/Sig_diag-theta/Sig_diag + qnorm(1-alpha)
  ub = qnorm(delta_u)/Sig_diag-theta/Sig_diag - qnorm(1-alpha)
  if (sum(lb>ub)>0){
    return(0)
  } else {
    Rho = cov2cor(Sigma)
    tmp_cdf_mvtnorm = mvtnorm::pmvnorm(lower=lb,
                                       upper=ub,
                                       mean=0,
                                       corr=Rho)[1]
    return(tmp_cdf_mvtnorm)
  }
}

get_mult_alpha_qTOST_MC = function(gamma, n_x, n_y, pi_x, delta_l, delta_u, alpha=0.05, B=10^3,
                          tol = .Machine$double.eps^0.5, seed=12345, max_iter=10, tolpower=1e-3, MC_sup=T, ...){
  l = n_y/n_x
  out = list()
  for (side in c("upper", "lower")){ # see argsup_qTOST: we try the argsup on both delta_l and delta_u
    theta_sups = matrix(NA, (max_iter+1), length(pi_x))
    theta_sups[1,]=find_sup_qTOST(gamma,n_x,n_y,pi_x,delta_l,delta_u,alpha,B,seed,side,MC_sup)
    alphas = powers = rep(NA, max_iter)
    alphas[1] = alpha
    i=1
    while (i <= max_iter) {
      sol_1 = get_qTOST_rvs(theta_sups[i,1], gamma, n_x, n_y, pi_x[1], B=B, seed=seed)
      sol_2 = get_qTOST_rvs(theta_sups[i,2], gamma, n_x, n_y, pi_x[2], B=B, seed=seed)
      powers[i] = power_qTOST_MC_biv(sol_1, sol_2, delta_l, delta_u, alphas[i])
      err_power = abs(powers[i]-alpha)
      if (err_power>tolpower) {
        theta  = theta_sups[i,]
        sol_aqTOST = get_qTOST_mv_core(theta=theta, gamma=gamma, alpha=alpha, pi_x=pi_x, n_x=n_x, n_y=n_y,
                                       delta_l=delta_l, delta_u=delta_u,
                                       B=B, tol=tol, seed=seed)
        alphas[i+1] = sol_aqTOST$min
        theta_sups[i+1,] = find_sup_qTOST(gamma,n_x,n_y,pi_x,delta_l,delta_u,alphas[i+1],B,seed,side=side,MC_sup)
        i = i+1
      } else {
        break
      }
    }
    alphas = alphas[!is.na(alphas)]
    powers = powers[!is.na(powers)]
    colnames(theta_sups)=paste0("argsup", 1:length(pi_x))
    theta_sups = theta_sups[complete.cases(theta_sups), ]
    if (is.null(nrow(theta_sups))) {
      theta_sup=theta_sups
    } else {
      theta_sup=theta_sups[nrow(theta_sups), ]
    }
    ind_side = which(c("upper", "lower") %in% side)
    out[[ind_side]] = list(min=alphas[length(alphas)],
                           theta_sup=theta_sup,
                           alphas=alphas,
                           theta_sups=matrix(theta_sups, ncol=length(pi_x)),
                           powers=powers,
                           iter=i-1,
                           err_power=err_power)
  }
  ind_out = which.min(c(out[[1]]$min, out[[2]]$min))
  out = out[[ind_out]]
  out
}

get_qTOST_mv_core = function(theta, gamma, alpha, pi_x, n_x, n_y, delta_l, delta_u, B=10^5, tol=.Machine$double.eps^0.5, seed=12345, ...){
  obj_func_alpha=obj_fun_qTOST_mv(test=alpha, theta=theta, gamma=gamma, alpha=alpha, pi_x=pi_x,
                                  n_x=n_x, n_y=n_y, delta_l=delta_l, delta_u=delta_u, B=B, seed=seed)
  obj_func_onehalf=obj_fun_qTOST_mv(test=0.5, theta=theta, gamma=gamma, alpha=alpha, pi_x=pi_x,
                                    n_x=n_x, n_y=n_y, delta_l=delta_l, delta_u=delta_u, B=B, seed=seed)
  if(obj_func_alpha<tol) return(list(minimum=alpha,objective=obj_func_alpha,theta_sup=theta))   #size at a=1 close to alpha
  if(obj_func_onehalf<tol) return(list(minimum=0.5,objective=obj_func_onehalf,theta_sup=theta)) #size at a=l close to 1/2
  out = optimize(obj_fun_qTOST_mv, interval=c(alpha, 0.5), tol=tol,
                 theta=theta, gamma=gamma, alpha = alpha, pi_x=pi_x,
                 n_x=n_x, n_y=n_y, delta_l=delta_l, delta_u=delta_u,
                 B=B, seed=seed)
  out$theta_sup = theta
  out
}

obj_fun_qTOST_mv = function(test, theta, gamma, alpha, pi_x, n_x, n_y, delta_l, delta_u, B=10^5, seed=12345, ...){
  sol_1 = get_qTOST_rvs(theta[1], gamma, n_x, n_y, pi_x[1], B=B, seed=seed)
  sol_2 = get_qTOST_rvs(theta[2], gamma, n_x, n_y, pi_x[2], B=B, seed=seed)
  size = power_qTOST_MC_biv(sol_1, sol_2, delta_l, delta_u, test)
  1e3*(size - alpha)^2
}
