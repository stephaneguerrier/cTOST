stopifnot(requireNamespace("cTOST",quietly=TRUE),
          as.character(packageVersion("cTOST"))=="1.1.0")

if (sim_setting==1){
  delta = 0.15
  pi_x = c(0.2, 0.8)
  delta_l = pi_x - delta
  delta_u = pi_x + delta
  gamma = 1
  n_y = 30
  n_x = 30
} else if (sim_setting==2){
  delta = 0.15
  pi_x = c(0.2, 0.8)
  delta_l = pi_x - delta
  delta_u = pi_x + delta
  gamma = 1
  n_y = 10
  n_x = 10
} else if (sim_setting==3){
  delta = 0.15
  pi_x = c(0.2, 0.8)
  delta_l = pi_x - delta
  delta_u = pi_x + delta
  gamma = 1
  n_y = 50
  n_x = 50
} 

nu_y = n_y-1
nu_x = n_x-1
alpha = 0.05
B = 50000
seed = 2026
tol = .Machine$double.eps^0.5
max_iter = 10
tolpower = 1e-3
MC_sup = TRUE
M = 50

theta_1 = (seq(qnorm(delta_l[1]), qnorm(delta_u[1]), len=M-4))
theta_1 = c(theta_1[1]-diff(theta_1)[1]*2, 
            theta_1[1]-diff(theta_1)[1], 
            theta_1,
            theta_1[length(theta_1)]+diff(theta_1)[1],
            theta_1[length(theta_1)]+diff(theta_1)[1]*2)
theta_2 = rev(seq(qnorm(delta_l[2]), qnorm(delta_u[2]), len=M-4))
theta_2 = c(theta_2[1]-diff(theta_2)[1]*2, 
            theta_2[1]-diff(theta_2)[1], 
            theta_2,
            theta_2[length(theta_2)]+diff(theta_2)[1],
            theta_2[length(theta_2)]+diff(theta_2)[1]*2)
