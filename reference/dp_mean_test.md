# Differentially Private Mean Equivalence Test (Internal)

Internal unified interface for DP-TOST that automatically selects the
best implementation. Uses ultra-fast C++ by default (15-250x speedup),
with option to use pure R.

## Usage

``` r
dp_mean_test(
  a,
  b,
  n,
  epsilon,
  mean_private_obs,
  sd_private_obs,
  lower,
  upper,
  B = 1000,
  alpha = 0.05,
  seed = NULL,
  method = c("cpp", "r")
)
```

## Arguments

- a:

  Lower bound for data truncation

- b:

  Upper bound for data truncation

- n:

  Sample size

- epsilon:

  Privacy budget

- mean_private_obs:

  Observed privatized mean

- sd_private_obs:

  Observed privatized standard deviation

- lower:

  Lower equivalence bound for TOST

- upper:

  Upper equivalence bound for TOST

- B:

  Number of bootstrap iterations (default: 1000)

- alpha:

  Significance level (default: 0.05)

- seed:

  Random seed for reproducibility (default: NULL)

- method:

  Implementation to use: "cpp" (default, fast) or "r" (exact R)

## Value

List with components:

- decision:

  1 if equivalence established, 0 otherwise

- ci_lower:

  Lower bound of confidence interval

- ci_upper:

  Upper bound of confidence interval

- mean_private_obs:

  Observed privatized mean (input)

- sd_private_obs:

  Observed privatized SD (input)

- mu_estimates:

  Bootstrap estimates of mu

- sigma_estimates:

  Bootstrap estimates of sigma

- method:

  Method used ("cpp" or "r")

## Details

This function performs equivalence testing under differential privacy.

\*\*Method Selection:\*\* - \`"cpp"\` (default): Ultra-fast C++
implementation using Nelder-Mead - 15-250x faster depending on problem
size - Uses R's random number generator through Rcpp, so results are
reproducible with \`seed\` - Recommended for production and large-scale
simulations

\- \`"r"\`: Pure R implementation using L-BFGS-B - Uses R's optim() with
gradient-based optimization - Same random draws as the C++ path, paired
differently, so results differ slightly - Useful for validation and
debugging

\*\*Performance Guide:\*\* - Small problems (n~100, B~1000): C++ is
60-250x faster - Large problems (n~1000, B~10000): C++ is 15-20x faster
