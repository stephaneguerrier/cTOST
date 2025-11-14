# Test script for DP-TOST one-sample proportion test
# Load the function
source("R/dp_proportion_test.R")

# Test 1: Basic functionality
cat("Test 1: Basic one-sample test\n")
cat("==============================\n")
set.seed(123)
n <- 200
true_p <- 0.55
x <- rbinom(n, 1, true_p)
p_obs <- mean(x)

# Add privacy noise
epsilon <- 1
p_hat <- p_obs + rlaplace_custom(1, scale = 1/(n * epsilon))

cat("True proportion:", true_p, "\n")
cat("Observed proportion:", p_obs, "\n")
cat("Privatized proportion:", p_hat, "\n\n")

# Run DP-TOST
lower <- 0.4
upper <- 0.6
result <- prop_test_dp_one_sample(
  p_hat = p_hat,
  n = n,
  lower = lower,
  upper = upper,
  epsilon = epsilon,
  alpha = 0.05,
  B = 1000,  # Using smaller B for quick test
  seed = 2024
)

cat("Confidence Interval:", result$conf.int, "\n")
cat("Equivalence decision:", result$decision, "\n")
cat("Equivalence bounds: [", result$lower, ",", result$upper, "]\n\n")

# Test 2: Test with different epsilon values
cat("\nTest 2: Different privacy budgets\n")
cat("==================================\n")
for (eps in c(0.5, 1, 2)) {
  p_hat_test <- p_obs + rlaplace_custom(1, scale = 1/(n * eps))
  result <- prop_test_dp_one_sample(
    p_hat = p_hat_test,
    n = n,
    lower = 0.4,
    upper = 0.6,
    epsilon = eps,
    B = 1000,
    seed = 2024
  )
  cat("epsilon =", eps, "| CI:", sprintf("[%.4f, %.4f]", result$conf.int[1], result$conf.int[2]),
      "| Decision:", result$decision, "\n")
}

# Test 3: Reproducibility
cat("\nTest 3: Reproducibility check\n")
cat("==============================\n")
result1 <- prop_test_dp_one_sample(p_hat = 0.6, n = 100, lower = 0.4, upper = 0.6, seed = 1337, B = 100)
result2 <- prop_test_dp_one_sample(p_hat = 0.6, n = 100, lower = 0.4, upper = 0.6, seed = 1337, B = 100)
cat("Same seed gives same results:",
    all.equal(result1$conf.int, result2$conf.int), "\n")

result3 <- prop_test_dp_one_sample(p_hat = 0.6, n = 100, lower = 0.4, upper = 0.6, seed = 9999, B = 100)
cat("Different seed gives different results:",
    !isTRUE(all.equal(result1$conf.int, result3$conf.int)), "\n")

# Test 4: Edge cases
cat("\nTest 4: Edge cases\n")
cat("==================\n")

# Proportion close to 0
result_low <- prop_test_dp_one_sample(p_hat = 0.05, n = 100, lower = 0, upper = 0.2, B = 500, seed = 1337)
cat("Low proportion (0.05): CI =", sprintf("[%.4f, %.4f]",
    result_low$conf.int[1], result_low$conf.int[2]), "\n")

# Proportion close to 1
result_high <- prop_test_dp_one_sample(p_hat = 0.95, n = 100, lower = 0.8, upper = 1, B = 500, seed = 1337)
cat("High proportion (0.95): CI =", sprintf("[%.4f, %.4f]",
    result_high$conf.int[1], result_high$conf.int[2]), "\n")

cat("\n✓ All tests completed!\n")
