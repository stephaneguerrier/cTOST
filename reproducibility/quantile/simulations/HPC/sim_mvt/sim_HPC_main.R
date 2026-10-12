rm(list=ls())
options(echo=TRUE)
library(tictoc)

sim_dir = Sys.getenv("SIM2_DIR")
setwd(sim_dir)

##########################################################
# number of simulations and arrays (inputs from bash)
args <- commandArgs(trailingOnly = TRUE)
print(args)
sim_setting = as.numeric(args[1])
n_simu = as.numeric(args[2])
n_array = as.numeric(args[3])
run_id = Sys.getenv("SLURM_ARRAY_JOB_ID")
id_slurm <- Sys.getenv("SLURM_ARRAY_TASK_ID")
id_slurm <- as.numeric(id_slurm)
print("id_slurm: ")
print(id_slurm)

# id of simu to be run on array
ind_mat <- matrix(1:n_simu, nr = n_array, byr = T)
Rw <- ind_mat[id_slurm, ]
destinationDir <- paste0("sim_HPC_out/sim_setting_", sim_setting, "/tmp")
if(!dir.exists(destinationDir)) dir.create(destinationDir,recursive=T)

source("sim_options.R")

############################################################

tic()
cat("Start -", date(),"\n\n\n")

pow_emp_TOST_ = pow_emp_piTOST_ = array(NA,dim=c(M, M, length(Rw)),
                                       dimnames=list(paste0("theta_2=",theta_2),
                                                     paste0("theta_1=",theta_1),
                                                     paste0("b=",Rw)))
for(rw in Rw){
  tic()
  cat("start Rw:",which(Rw==rw),"/",length(Rw),"\n")
  set.seed(1000000+rw)
  Z = rnorm(1)
  W1 = sqrt(rchisq(n=1,df=nu_x))
  W2 = sqrt(rchisq(n=1,df=nu_y))
  gamma_hat =  W2^2/W1^2*gamma * nu_x/nu_y
  outer_rng = .Random.seed
  alpha_qtost_calibration = tryCatch(cTOST:::get_mult_alpha_qTOST_MC(gamma_hat, n_x, n_y, pi_x, delta_l, delta_u,
      alpha=alpha, B=B, seed=seed, max_iter=max_iter, tol=tol, tolpower=tolpower, MC_sup=MC_sup),
      finally={.Random.seed = outer_rng})
  alpha_star = alpha_qtost_calibration$min
  sgamma = sqrt(gamma)
  snu_y = sqrt(nu_y)
  snu_x = sqrt(nu_x)
  q_x_1 = qnorm(pi_x[1])
  q_x_2 = qnorm(pi_x[2])
  b2 = snu_y/sgamma*sqrt(1/n_x + gamma/n_y)
  for (i in 1:M){
    for (j in 1:M){
      b1 = snu_y/sgamma*(theta_1[j]*sgamma-q_x_1)
      b3 = snu_y*q_x_1/(snu_x*sgamma)
      theta_hat_1 = b1/W2 + b2*(Z/W2) + b3*W1/W2
      b1 = snu_y/sgamma*(theta_2[i]*sgamma-q_x_2)
      b3 = snu_y*q_x_2/(snu_x*sgamma)
      theta_hat_2 = b1/W2 + b2*(Z/W2) + b3*W1/W2
      sigma_hat_theta_1 = sqrt((1/n_y) * (1+theta_hat_1^2/2) + (1/n_x)/gamma_hat*(1+qnorm(pi_x[1])^2/2))
      sigma_hat_theta_2 = sqrt((1/n_y) * (1+theta_hat_2^2/2) + (1/n_x)/gamma_hat*(1+qnorm(pi_x[2])^2/2)) 
      # qTOST
      qtost_sol = cTOST:::qtost_core(c(theta_hat_1, theta_hat_2),
          c(sigma_hat_theta_1, sigma_hat_theta_2), pi_x, delta_l, delta_u, alpha=alpha)
      pow_emp_TOST_[i,j,which(Rw %in% rw)] = all(qtost_sol$decision)
      # alpha-qTOST
      alpha_qtost_sol = cTOST:::qtost_core(c(theta_hat_1, theta_hat_2),
          c(sigma_hat_theta_1, sigma_hat_theta_2), pi_x, delta_l, delta_u,
          alpha=alpha, corrected_alpha=alpha_star)
      pow_emp_piTOST_[i,j,which(Rw %in% rw)] = all(alpha_qtost_sol$decision)
    }
  } # end theta           
  tt_Rw=toc(quiet=T)
  cat("end Rw:",which(Rw==rw),"/",length(Rw),"- ",tt_Rw$callback_msg,"\n\n")
} # end rw

save(
  pow_emp_TOST_,
  pow_emp_piTOST_,
  Rw, run_id,
  file=paste0(destinationDir,"/",Rw[1]))