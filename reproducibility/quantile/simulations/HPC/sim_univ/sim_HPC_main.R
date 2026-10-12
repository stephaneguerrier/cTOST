rm(list=ls())
options(echo=TRUE)
library(tictoc)

sim_sigma.hat = function(sigma, nu) sqrt(rchisq(n=1,df=nu)*sigma^2/nu)

##########################################################
# number of simulations and arrays (inputs from bash)
sim_dir = Sys.getenv("SIM1_DIR")
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
destinationDir <- file.path(sim_dir, paste0("sim_HPC_out/sim_setting_", sim_setting, "/tmp"))
if(!dir.exists(destinationDir)) dir.create(destinationDir,recursive=T)

source(file.path(sim_dir, "sim_options.R"))

############################################################

tic()
cat("Start -", date(),"\n\n\n")

pow_emp_TOST_ =
  pow_emp_alphaTOST_ = array(NA,dim=c(length(n_),length(l_),length(pi_x_),length(gamma_),length(theta_), length(Rw)),
                             dimnames=list(paste0("n=",n_),
                                           paste0("l=",l_),
                                           paste0("pi_x=",pi_x_,"; cte=", delta_),
                                           paste0("gamma=",round(gamma_,2)),
                                           paste0("theta=",theta_),
                                           paste0("b=",Rw)))

for(rw in Rw){
    tic()
    cat("start Rw:",which(Rw==rw),"/",length(Rw),"\n")
    set.seed(rw)
    for(u in 1:length(n_)){
      tic()
      if(u>1){
        cat("Rw:",which(Rw==rw),"/",length(Rw),"\n")
        cat("\tstart n_:",u,"/",length(n_),"\n")
      }
      counter_pb = 0
      pb = txtProgressBar(min = 0, max = length(l_)*length(pi_x_)*length(gamma_), style = 3)
      n_y = n_[u]
      nu_y = n_y-1
      for (i_l in 1:length(l_)){
        l = l_[i_l]
        n_x = n_y * l
        nu_x = n_x-1
        for(v in 1:length(pi_x_)){
          pi_x = pi_x_[v]
          delta = delta_[v]
          delta_u = pi_x + delta
          delta_l = pi_x - delta
          if (sim_size){
            theta_ = c(qnorm(pi_x-delta), qnorm(pi_x+delta)) # size
          } else {
            theta_ = seq(qnorm(pi_x-delta*1.2), qnorm(pi_x+delta*1.2), len=50) # power curves  
          }
          for(r in 1:length(gamma_)){ 
            gamma = gamma_[r]
            sigma_x = sigma
            sigma_y = sqrt(gamma)*sigma
            if (!sim_size){
              set.seed(1000000+rw)
            }
            sigma_x_hat = sim_sigma.hat(sigma_x, nu_x)
            sigma_y_hat = sim_sigma.hat(sigma_y, nu_y)
            gamma_hat = sigma_y_hat^2/sigma_x_hat^2
            cal_seed = seed
            #####################
            data_state = .Random.seed
            calibration = tryCatch(cTOST:::get_alpha_qTOST_MC(
              gamma=gamma_hat, n_x=n_x, n_y=n_y, pi_x=pi_x,
              delta_l=delta_l, delta_u=delta_u, alpha=alpha, B=B, seed=cal_seed, tol=1e-6),
              error=function(e) list(root=NA_real_, error=conditionMessage(e)))
            .Random.seed = data_state
            alpha_star = calibration$root
            stopifnot(length(alpha_star)==1L,
                      is.na(alpha_star) || (is.finite(alpha_star) && alpha_star>=alpha && alpha_star<=.5))
            if (!is.null(calibration$error)) cat("Calibration failed at rw=",rw,
              ", case=",u,i_l,v,r,": ",calibration$error,"\n",file=stderr())
            for(k in 1:length(theta_)){
              if (!sim_size){
                set.seed(1000000+n_simu+rw)
              }
              diff_mu = theta_[k]*sigma_y - qnorm(pi_x)*sigma_x
              diff_sig = sqrt(sigma_x^2/n_x + sigma_y^2/n_y)
              zbar = rnorm(n=1, mean=diff_mu, sd=diff_sig)
              qtost_sol = cTOST::qtost(x=list(mean=zbar,sd=sigma_x_hat,n=n_x),
                y=list(mean=0,sd=sigma_y_hat,n=n_y), pi_x=pi_x,delta=delta,
                alpha=alpha,method="unadjusted",B=B,seed=cal_seed,tol=1e-6)
              theta_hat = qtost_sol$theta
              sigma_hat_theta = qtost_sol$sigma
              # qTOST
              pow_emp_TOST_[u,i_l,v,r,k,which(Rw %in% rw)] = as.numeric(qtost_sol$decision)
              # alpha-qTOST
              if (!is.na(alpha_star)) {
                alpha_qtost_sol = cTOST:::qtost_core(theta=theta_hat,sigma=sigma_hat_theta,
                  pi_x=pi_x,delta_l=delta_l,delta_u=delta_u,alpha=alpha,corrected_alpha=alpha_star)
                pow_emp_alphaTOST_[u,i_l,v,r,k,which(Rw %in% rw)] = as.numeric(alpha_qtost_sol$decision)
              }
            } # end theta
            counter_pb = counter_pb + 1
            setTxtProgressBar(pb, counter_pb)
          } # end sigma
        } # end pi_x
      } # end l
      tt_n=toc(quiet=T)
      if(u==length(n_)){
        cat("\n\tend n_:",u,"/",length(n_),"- ",tt_n$callback_msg,"\n")
      }else{
        cat("\n\tend n_:",u,"/",length(n_),"- ",tt_n$callback_msg,"\n\n")
      }
    } # end n
    tt_Rw=toc(quiet=T)
    cat("end Rw:",which(Rw==rw),"/",length(Rw),"- ",tt_Rw$callback_msg,"\n\n")
} # end rw

save(
  pow_emp_TOST_,
  pow_emp_alphaTOST_,
  Rw, run_id,
  file=paste0(destinationDir,"/",Rw[1]))