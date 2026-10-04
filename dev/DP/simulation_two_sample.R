# Simulation study for DP-TOST two-sample proportion test
# Power curve analysis across different true proportion differences

# Load functions
source("R/dp_proportion_test.R")

# Simulation parameters
n1 <- 300
n2 <- 300
p1_true <- 0.5  # Fixed proportion for group 1
diff_true_grid <- seq(-0.15, 0.15, length.out = 13)  # Grid of true differences (p1 - p2)
delta <- 0.1  # Equivalence margin
epsilon1 <- 1
epsilon2 <- 1
alpha <- 0.05
n_sim <- 100

# Storage for results
results <- matrix(NA, nrow = n_sim, ncol = length(diff_true_grid))
colnames(results) <- paste0("diff_", round(diff_true_grid, 4))

# Set overall seed for reproducibility
set.seed(1337)

cat("Running DP-TOST Two-Sample Simulation\n")
cat("======================================\n")
cat("Parameters:\n")
cat("  n1 =", n1, ", n2 =", n2, "\n")
cat("  Fixed p1 =", p1_true, "\n")
cat("  True diff (p1-p2) grid: [", min(diff_true_grid), ",", max(diff_true_grid), "] with",
    length(diff_true_grid), "values\n")
cat("  Delta =", delta, "\n")
cat("  Epsilon1 =", epsilon1, ", Epsilon2 =", epsilon2, "\n")
cat("  Alpha =", alpha, "\n")
cat("  Number of simulations per diff =", n_sim, "\n\n")

# Progress tracking
total_sims <- length(diff_true_grid) * n_sim
pb <- txtProgressBar(min = 0, max = total_sims, style = 3)
counter <- 0

# Loop over true differences
for (j in 1:length(diff_true_grid)) {

  true_diff <- diff_true_grid[j]
  p2_true <- p1_true - true_diff  # Compute p2 from difference

  # Skip if p2 is outside [0, 1]
  if (p2_true < 0 || p2_true > 1) {
    results[, j] <- NA
    counter <- counter + n_sim
    setTxtProgressBar(pb, counter)
    next
  }

  # Run simulations for this true difference
  for (i in 1:n_sim) {

    # Generate data for both groups
    x1 <- rbinom(n1, 1, p1_true)
    x2 <- rbinom(n2, 1, p2_true)

    p1_obs <- mean(x1)
    p2_obs <- mean(x2)

    # Add privacy noise to each proportion
    p1_hat <- p1_obs + rlaplace_custom(1, scale = 1 / (n1 * epsilon1))
    p2_hat <- p2_obs + rlaplace_custom(1, scale = 1 / (n2 * epsilon2))

    # Compute equivalence bounds for difference
    lower <- -delta
    upper <- delta

    # Run DP-TOST two-sample test with unique seed (internal function)
    test_result <- prop_test_dp_two_sample(
      p1_hat = p1_hat,
      p2_hat = p2_hat,
      n1 = n1,
      n2 = n2,
      lower = lower,
      upper = upper,
      epsilon1 = epsilon1,
      epsilon2 = epsilon2,
      alpha = alpha,
      B = 10^3,
      seed = 1337 + j * 1000 + i
    )

    # Test wrapper function with same inputs
    test_result_wrapper <- prop_test_equiv_dp(
      p_hat = p1_hat,
      n = n1,
      p_hat2 = p2_hat,
      n2 = n2,
      lower = lower,
      upper = upper,
      epsilon = c(epsilon1, epsilon2),
      alpha = alpha,
      B = 10^3,
      seed = 1337 + j * 1000 + i
    )

    # Verify both give same decision
    if (test_result$decision != test_result_wrapper$decision) {
      stop("Mismatch between internal and wrapper function!")
    }

    # Store decision
    results[i, j] <- test_result$decision

    # Update progress
    counter <- counter + 1
    setTxtProgressBar(pb, counter)
  }
}

close(pb)

cat("\n\nSimulation Results\n")
cat("==================\n\n")

# Compute empirical probabilities
empirical_prob <- colMeans(results, na.rm = TRUE)

# Print summary table
cat("True diff | Empirical Prob H1 | In Equiv Region\n")
cat("----------|-------------------|----------------\n")
for (j in 1:length(diff_true_grid)) {
  in_region <- (diff_true_grid[j] >= -delta) && (diff_true_grid[j] <= delta)
  cat(sprintf("%8.4f  | %.4f            | %s\n",
              diff_true_grid[j], empirical_prob[j], in_region))
}
cat("\n")

# Plot power curve
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))

plot(diff_true_grid, empirical_prob, type = "b",
     pch = 19, col = "steelblue", lwd = 2,
     xlab = "True difference (p1 - p2)",
     ylab = "Empirical probability of accepting H1",
     main = "DP-TOST Two-Sample Power Curve",
     ylim = c(0, 1),
     panel.first = grid())

# Add equivalence region
abline(v = -delta, col = "red", lty = 2, lwd = 2)
abline(v = delta, col = "red", lty = 2, lwd = 2)
abline(v = 0, col = "gray50", lty = 3)

# Add reference line at alpha
abline(h = alpha, col = "darkgreen", lty = 2)

# Add shaded equivalence region
rect(-delta, -1, delta, 2,
     col = rgb(1, 0, 0, 0.1), border = NA)

# Legend
legend("topleft",
       legend = c("Empirical probability",
                  "Equivalence bounds",
                  "Type I error (α)"),
       col = c("steelblue", "red", "darkgreen"),
       lty = c(1, 2, 2), lwd = 2, pch = c(19, NA, NA),
       bty = "n")

# Save results
save_data <- list(
  diff_true_grid = diff_true_grid,
  empirical_prob = empirical_prob,
  results = results,
  params = list(n1 = n1, n2 = n2, p1_true = p1_true, delta = delta,
                epsilon1 = epsilon1, epsilon2 = epsilon2, alpha = alpha, n_sim = n_sim)
)

save_path <- "tests/DP/simulation_two_sample_results.rds"
saveRDS(save_data, save_path)
cat("Results saved to:", save_path, "\n")

cat("\n✓ Simulation completed!\n")
