rm(list=ls())

args <- commandArgs(trailingOnly = TRUE)
sim_setting = as.numeric(args[1])
n_simu = as.numeric(args[2])
n_array = as.numeric(args[3])
run_id = Sys.getenv("SIM_RUN_ID")

sim_dir = Sys.getenv("SIM2_DIR")
setwd(sim_dir)

ind_mat <- matrix(1:n_simu, nr = n_array, byr = T)
indexes_i = ind_mat[,1]
load_dir = paste0("sim_HPC_out/sim_setting_", sim_setting, "/tmp")
destinationDir = paste0("sim_HPC_out/sim_setting_", sim_setting)
cleanup_marker = file.path(destinationDir, ".cleanup_ok")
if (file.exists(cleanup_marker)) stop("Completion marker already exists")

source("sim_options.R")

name_sim = paste0("results_mvt_sim_setting_", sim_setting, "_n_simu_", n_simu, "_n_array_", n_array, ".rdata")
result_file = file.path(destinationDir, name_sim)

####################################

# combine simulation chunks into rejection percentages.
sol_TOST_sum = sol_piTOST_sum = matrix(0, M, M,
                                       dimnames = list(paste0("theta_2=", theta_2), 
                                                       paste0("theta_1=", theta_1)))
if (!setequal(list.files(load_dir, all.files = TRUE, no.. = TRUE), as.character(indexes_i)))
  stop("Missing or unexpected chunk files")
n_complete = 0
for (i in indexes_i) {
  x = new.env(parent = emptyenv())
  load(file.path(load_dir, i), envir = x)
  required = c("pow_emp_TOST_", "pow_emp_piTOST_", "Rw", "run_id")
  expected_rw = as.vector(ind_mat[which(indexes_i == i), ])
  expected_names = c(dimnames(sol_TOST_sum), list(paste0("b=", expected_rw)))
  if (!setequal(required,ls(x)) ||
      !identical(x$run_id,run_id) ||
      !identical(x$Rw, expected_rw)) stop("Incompatible or incorrectly indexed chunk: ", i)
  for (object in c("pow_emp_TOST_", "pow_emp_piTOST_")) {
    decisions = x[[object]]
    if (!identical(dim(decisions), as.integer(c(M, M, length(expected_rw)))) ||
        !identical(dimnames(decisions), expected_names) ||
        !(is.logical(decisions) || is.numeric(decisions)) ||
        anyNA(decisions) || any(!decisions %in% c(0, 1))) stop("Invalid decisions in chunk: ", i)
  }
  sol_TOST_sum = sol_TOST_sum + apply(x$pow_emp_TOST_, 1:2, sum)
  sol_piTOST_sum = sol_piTOST_sum + apply(x$pow_emp_piTOST_, 1:2, sum)
  n_complete = n_complete + length(expected_rw)
  rm(x, decisions)
}
stopifnot(n_complete == n_simu)
sol_TOST_avg = sol_TOST_sum / n_complete * 100
sol_piTOST_avg = sol_piTOST_sum / n_complete * 100
if (file.exists(result_file)) stop("Combined result already exists")
save(sol_TOST_avg, sol_piTOST_avg, file = result_file)
marker_connection = file(cleanup_marker, open="wb")
writeLines(c(run_id, basename(result_file)), marker_connection)
close(marker_connection)
