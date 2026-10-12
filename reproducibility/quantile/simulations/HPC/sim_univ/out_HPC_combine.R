rm(list=ls())

sim_dir = Sys.getenv("SIM1_DIR")

# simulation inputs from bash
args <- commandArgs(trailingOnly = TRUE)
sim_setting = as.numeric(args[1])
n_simu = as.numeric(args[2])
n_array = as.numeric(args[3])
run_id = Sys.getenv("SIM_RUN_ID")

# create matrix of indices
ind_mat <- matrix(1:n_simu, nr = n_array, byr = T)

indexes_i = ind_mat[,1]

load_dir       = file.path(sim_dir, paste0("sim_HPC_out/sim_setting_", sim_setting, "/tmp"))
destinationDir = file.path(sim_dir, paste0("sim_HPC_out/sim_setting_", sim_setting))
cleanup_marker = file.path(destinationDir, ".cleanup_ok")
if (file.exists(cleanup_marker)) stop("Completion marker already exists")

source(file.path(sim_dir, "sim_options.R"))

if (sim_size){
  name_sim = paste0("results_sim_setting_size_", sim_setting, "_n_simu_", n_simu, "_n_array_", n_array)
} else {
  name_sim = paste0("results_sim_setting_", sim_setting, "_n_simu_", n_simu, "_n_array_", n_array)
}
result_file = file.path(destinationDir,paste0(name_sim,".rdata"))

####################################

sol_TOST_ =
  sol_alphaTOST_ = array(NA,dim=c(length(n_),length(l_),length(pi_x_),length(gamma_),length(theta_), n_simu),
                             dimnames=list(paste0("n=",n_),
                                           paste0("l=",l_),
                                           paste0("pi_x=",pi_x_,"; cte=", delta_),
                                           paste0("gamma=",gamma_),
                                           paste0("theta=",theta_),
                                           paste0("b=",1:n_simu)))


for (i in indexes_i) {
  chunk = new.env(parent = emptyenv())
  load(file.path(load_dir, i), envir = chunk)
  required = c("pow_emp_TOST_", "pow_emp_alphaTOST_", "Rw", "run_id")
  expected_rw = as.vector(ind_mat[which(indexes_i == i), ])
  expected_names = dimnames(sol_TOST_)
  expected_names[[6]] = paste0("b=", expected_rw)
  if (!setequal(required,ls(chunk)) ||
      !identical(chunk$run_id,run_id) ||
      !identical(chunk$Rw,expected_rw)) stop("Incompatible or incorrectly indexed chunk: ", i)
  for (object in c("pow_emp_TOST_", "pow_emp_alphaTOST_")) {
    decisions = chunk[[object]]
    if (!identical(dim(decisions), as.integer(c(dim(sol_TOST_)[1:5], length(expected_rw)))) ||
        !identical(dimnames(decisions), expected_names) ||
        !(is.logical(decisions) || is.numeric(decisions)) ||
        anyNA(decisions) || any(!decisions %in% c(0, 1))) stop("Invalid decisions in chunk: ", i)
  }
  index_end = ind_mat[which(ind_mat[,1] == i), ncol(ind_mat)]
  sol_TOST_[, , , , , i:index_end] = chunk$pow_emp_TOST_
  sol_alphaTOST_[, , , , , i:index_end] = chunk$pow_emp_alphaTOST_
  rm(chunk, decisions)
}

sol_TOST_avg = apply(sol_TOST_, 1:5, mean) * 100
sol_alphaTOST_avg = apply(sol_alphaTOST_, 1:5, mean) * 100

if (file.exists(result_file)) stop("Combined result already exists")
save(sol_TOST_avg, sol_alphaTOST_avg, file=result_file)
marker_connection = file(cleanup_marker, open="wb")
writeLines(c(run_id, basename(result_file)), marker_connection)
close(marker_connection)
