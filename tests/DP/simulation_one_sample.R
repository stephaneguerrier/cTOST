# Simulation study for DP-TOST one-sample proportion test
# Power curve analysis across different true proportions

# Load functions
source("R/dp_proportion_test.R")

# Simulation parameters
n <- 500
p_true_grid <- seq(0.4, 0.52, length.out = 13)
p0 <- 0.4  # Reference proportion
delta <- 0.1
epsilon <- 0.1
alpha <- 0.05
n_sim <- 100

# Storage for results
results <- matrix(NA, nrow = n_sim, ncol = length(p_true_grid))
colnames(results) <- paste0("p_", round(p_true_grid, 4))

# Set overall seed for reproducibility
set.seed(1337)

cat("Running DP-TOST Simulation\n")
cat("==========================\n")
cat("Parameters:\n")
cat("  n =", n, "\n")
cat("  True p grid: [", min(p_true_grid), ",", max(p_true_grid), "] with",
    length(p_true_grid), "values\n")
cat("  Reference p0 =", p0, "\n")
cat("  Delta =", delta, "\n")
cat("  Epsilon =", epsilon, "\n")
cat("  Alpha =", alpha, "\n")
cat("  Number of simulations per p_true =", n_sim, "\n\n")

# Progress tracking
total_sims <- length(p_true_grid) * n_sim
pb <- txtProgressBar(min = 0, max = total_sims, style = 3)
counter <- 0

# Loop over true proportions
for (j in 1:length(p_true_grid)) {

  true_p <- p_true_grid[j]

  # Run simulations for this true_p
  for (i in 1:n_sim) {

    # Generate data
    x <- rbinom(n, 1, true_p)
    p_obs <- mean(x)

    # Add privacy noise
    p_hat <- p_obs + rlaplace_custom(1, scale = 1 / (n * epsilon))

    # Compute equivalence bounds
    lower <- p0 - delta
    upper <- p0 + delta

    # Run DP-TOST with unique seed (internal function)
    test_result <- prop_test_dp_one_sample(
      p_hat = p_hat,
      n = n,
      lower = lower,
      upper = upper,
      epsilon = epsilon,
      alpha = alpha,
      B = 10^3,
      seed = 1337 + j * 1000 + i
    )

    # Test wrapper function with same inputs
    test_result_wrapper <- prop_test_equiv_dp(
      p_hat = p_hat,
      n = n,
      lower = lower,
      upper = upper,
      epsilon = epsilon,
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
empirical_prob <- colMeans(results)

# Print summary table
cat("True p    | Empirical Prob H1 | In Equiv Region\n")
cat("----------|-------------------|----------------\n")
for (j in 1:length(p_true_grid)) {
  in_region <- (p_true_grid[j] >= (p0 - delta)) && (p_true_grid[j] <= (p0 + delta))
  cat(sprintf("%.4f    | %.4f            | %s\n",
              p_true_grid[j], empirical_prob[j], in_region))
}
cat("\n")

# Plot power curve
par(mfrow = c(1, 1), mar = c(5, 4, 4, 2))

plot(p_true_grid, empirical_prob, type = "b",
     pch = 19, col = "steelblue", lwd = 2,
     xlab = "True proportion (p)",
     ylab = "Empirical probability of accepting H1",
     main = "DP-TOST Power Curve",
     ylim = c(0, 1),
     panel.first = grid())

# Add equivalence region
abline(v = p0 - delta, col = "red", lty = 2, lwd = 2)
abline(v = p0 + delta, col = "red", lty = 2, lwd = 2)
abline(v = p0, col = "gray50", lty = 3)

# Add reference line at alpha
abline(h = alpha, col = "darkgreen", lty = 2)

# Add shaded equivalence region
rect(p0 - delta, -1, p0 + delta, 2,
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
  p_true_grid = p_true_grid,
  empirical_prob = empirical_prob,
  results = results,
  params = list(n = n, p0 = p0, delta = delta,
                epsilon = epsilon, alpha = alpha, n_sim = n_sim)
)

save_path <- "tests/DP/simulation_results.rds"
saveRDS(save_data, save_path)
cat("Results saved to:", save_path, "\n")

cat("\n✓ Simulation completed!\n")
