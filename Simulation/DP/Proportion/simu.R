library(cTOST)

# Get simulation parameters from environment variables
epsilon <- as.numeric(Sys.getenv('epsilon_val'))
n_val <- as.numeric(Sys.getenv('n_val'))
p1 <- as.numeric(Sys.getenv('p1_val'))
scenario <- as.numeric(Sys.getenv('scenario_val'))

# Fixed parameters
alpha <- 0.05
lower <- -0.1
upper <- 0.1
n1 <- n_val
n2 <- n_val  # Balanced design

# Set p2 based on scenario
if (scenario == 1) {
  p2 <- p1 + 0.1  # True diff = -0.1 (at lower bound)
} else {
  p2 <- p1 - 0.1  # True diff = +0.1 (at upper bound)
}

# Get SLURM array task ID
id_slurm <- Sys.getenv("SLURM_ARRAY_TASK_ID")
id_slurm <- as.numeric(id_slurm)

# Monte Carlo iterations per job
iters_per_job <- 3334

# Initialize results matrix
# Columns: decision_dp, ci_lower_dp, ci_upper_dp,
#          decision_classical, ci_lower_classical, ci_upper_classical,
#          p1_hat, p2_hat, p1_obs, p2_obs,
#          p1_true, p2_true, n, epsilon, scenario, id_slurm, seed
res <- matrix(NA, nrow = iters_per_job, ncol = 17)

# Critical value for 90% CI (1 - 2*alpha = 0.90)
z_crit <- qnorm(1 - alpha)

# Create unique parameter IDs for seed generation
# Ensures all seeds are unique across all parameter combinations
eps_id <- match(epsilon, c(0.1, 0.25, 0.5, 1, 2, 5, 10))
n_id <- match(n_val, seq(200, 1000, 100))
p1_id <- match(p1, c(0.5, 0.65, 0.8))

# Run Monte Carlo iterations
for (i in 1:iters_per_job) {
  # Unique seed incorporating all parameters
  # Base seed encodes: epsilon, n, p1, scenario
  base_seed <- eps_id * 1e7 + n_id * 1e6 + p1_id * 1e5 + scenario * 1e4
  seed_i <- base_seed + id_slurm * iters_per_job + i
  set.seed(seed_i)

  # Generate Bernoulli data
  x1 <- rbinom(n1, size = 1, prob = p1)
  x2 <- rbinom(n2, size = 1, prob = p2)

  # Observed proportions (non-private)
  p1_obs <- mean(x1)
  p2_obs <- mean(x2)

  # Classical non-DP test using CLT
  diff_obs <- p1_obs - p2_obs
  se_classical <- sqrt(p1_obs * (1 - p1_obs) / n1 + p2_obs * (1 - p2_obs) / n2)
  ci_classical_lower <- diff_obs - z_crit * se_classical
  ci_classical_upper <- diff_obs + z_crit * se_classical
  decision_classical <- (ci_classical_lower > lower) && (ci_classical_upper < upper)

  # Add Laplace noise (differential privacy)
  scale1 <- 1 / (n1 * epsilon)
  scale2 <- 1 / (n2 * epsilon)

  u1 <- runif(1)
  u2 <- runif(1)

  laplace_noise1 <- -scale1 * sign(u1 - 0.5) * log(1 - 2 * abs(u1 - 0.5))
  laplace_noise2 <- -scale2 * sign(u2 - 0.5) * log(1 - 2 * abs(u2 - 0.5))

  p1_hat <- p1_obs + laplace_noise1
  p2_hat <- p2_obs + laplace_noise2

  # Run DP-TOST (using B = 10^3)
  test_result_dp <- prop_test_equiv_dp(
    p_hat = p1_hat,
    n = n1,
    p_hat2 = p2_hat,
    n2 = n2,
    lower = lower,
    upper = upper,
    epsilon = c(epsilon, epsilon),
    alpha = alpha,
    B = 10^3,
    seed = seed_i + 1000000  # Different seed for the test
  )

  # Store results
  res[i, ] <- c(
    as.numeric(test_result_dp$decision),
    test_result_dp$conf.int[1],
    test_result_dp$conf.int[2],
    as.numeric(decision_classical),
    ci_classical_lower,
    ci_classical_upper,
    p1_hat,
    p2_hat,
    p1_obs,
    p2_obs,
    p1,
    p2,
    n1,
    epsilon,
    scenario,
    id_slurm,
    seed_i
  )
}

# Define filename
name_file <- paste0("ctost/proportion/data_temp/results",
                   "_eps_", epsilon,
                   "_n_", n_val,
                   "_p1_", p1,
                   "_scen_", scenario,
                   "_id_", id_slurm,
                   ".rda")

# Save
save(res, file = name_file)
