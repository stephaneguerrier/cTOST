stopifnot(requireNamespace("cTOST",quietly=TRUE),
          as.character(packageVersion("cTOST"))=="1.1.0")

pi_x_ = c(0.05, 0.1, 0.25, 0.5)
delta_ = c(0.025, 0.05, 0.1, 0.2)
gamma_ = c(1, 2)
l_ = c(1, 2, 3)
n_ = c(30, 60, 90, 120)
sigma = 1

stopifnot(sim_setting %in% 1:2)
sim_size = if (sim_setting==1) 1 else 0
if (sim_size){
  theta_ = 1:2
} else {
  theta_ = 1:50  
}

B = 50000
seed = 2026
alpha = 0.05