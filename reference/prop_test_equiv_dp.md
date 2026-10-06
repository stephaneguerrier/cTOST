# Differentially Private TOST for Proportion Equivalence Testing

Performs equivalence testing for proportions under differential privacy
using the DP-TOST procedure. Supports both one-sample and two-sample
tests.

## Usage

``` r
prop_test_equiv_dp(
  p_hat,
  n,
  p_hat2 = NULL,
  n2 = NULL,
  lower,
  upper,
  epsilon,
  alpha = 0.05,
  B = 10^4,
  max_resample = 100,
  seed = 1337,
  ...
)
```

## Arguments

- p_hat:

  Privatized sample proportion for the first (or only) group. This
  should be the observed proportion with differential privacy noise
  already added.

- n:

  Sample size for the first (or only) group.

- p_hat2:

  Privatized sample proportion for the second group (optional). If NULL,
  a one-sample test is performed. Default is NULL.

- n2:

  Sample size for the second group (optional). Required if p_hat2 is
  provided. Default is NULL.

- lower:

  Lower equivalence bound. For one-sample tests, this is the lower bound
  for the proportion. For two-sample tests, this is the lower bound for
  the difference (p1 - p2). Typical value: -delta.

- upper:

  Upper equivalence bound. For one-sample tests, this is the upper bound
  for the proportion. For two-sample tests, this is the upper bound for
  the difference (p1 - p2). Typical value: delta.

- epsilon:

  Privacy budget. Can be a single value (used for all samples) or a
  vector of length 2 (epsilon\[1\] for first group, epsilon\[2\] for
  second group). Larger values provide less privacy but more statistical
  power.

- alpha:

  Significance level for the test. The confidence interval will have
  level (1 - 2\*alpha). Default is 0.05.

- B:

  Number of Monte Carlo replications for reconstructing the sampling
  distribution. Larger values provide more accurate results but take
  longer. Default is 10000.

- max_resample:

  Maximum number of resampling attempts when solving the moment matching
  equation. Default is 100.

- seed:

  Random seed for reproducibility. Default is 1337.

- ...:

  Additional arguments (currently unused).

## Value

A list of class "prop_test_dp" containing:

- decision:

  Logical. TRUE if equivalence is established, FALSE otherwise.

- conf.int:

  Confidence interval for the parameter of interest (proportion for
  one-sample, difference for two-sample).

- lower:

  Lower equivalence bound.

- upper:

  Upper equivalence bound.

- alpha:

  Significance level used.

- epsilon:

  Privacy budget(s) used.

- estimate:

  Point estimate(s) of the privatized proportion(s).

- sample_size:

  Sample size(s).

- test_type:

  Character string indicating "one-sample" or "two-sample".

- seed:

  Random seed used.

- B:

  Number of Monte Carlo replications.

## Details

This function implements the DP-TOST (Differentially Private Two
One-Sided Tests) procedure for testing equivalence of proportions while
maintaining differential privacy guarantees. The test uses Monte Carlo
simulation to reconstruct the sampling distribution from privatized
statistics and constructs confidence intervals using the percentile
method.

For the one-sample case, the test evaluates whether a proportion is
equivalent to a reference range \[lower, upper\]. For the two-sample
case, it tests whether the difference between two proportions (p1 - p2)
falls within the equivalence bounds \[lower, upper\].

The privacy mechanism adds Laplace noise to the observed proportions.
The privacy budget epsilon controls the amount of noise: larger epsilon
values provide less privacy but more accurate inference.

## Examples

``` r
# One-sample test: Is proportion equivalent to range [0.3, 0.5]?
set.seed(123)
n <- 500
p_true <- 0.42
epsilon <- 1

# Generate Bernoulli data
x <- rbinom(n, 1, p_true)
p_obs <- mean(x)

# Add Laplace noise using inverse CDF method
u <- runif(1)
scale <- 1 / (n * epsilon)
laplace_noise <- -scale * sign(u - 0.5) * log(1 - 2 * abs(u - 0.5))
p_hat <- p_obs + laplace_noise

result <- prop_test_equiv_dp(
  p_hat = p_hat,
  n = n,
  lower = 0.3,
  upper = 0.5,
  epsilon = 1,
  B = 1000
)
result
#> ✔ Accept equivalence
#> Equiv. Region:  |-------------------------------|
#> Estim. Inter.:             (-------x------)        
#> CI = (0.37602 ; 0.44468)
#> 
#> Method: DP-TOST (proportion)
#> alpha = 0.05; Equiv. bounds = [0.3, 0.5]
#> epsilon = 1; B = 1000
#> n = 500; p_hat = 0.40931

# Two-sample test: Is difference (p1-p2) equivalent to [-0.1, 0.1]?
set.seed(456)
n1 <- 300
n2 <- 300
p1_true <- 0.52
p2_true <- 0.48
epsilon1 <- 1
epsilon2 <- 1

# Generate Bernoulli data for both groups
x1 <- rbinom(n1, 1, p1_true)
x2 <- rbinom(n2, 1, p2_true)
p1_obs <- mean(x1)
p2_obs <- mean(x2)

# Add Laplace noise to each proportion
u1 <- runif(1)
u2 <- runif(1)
scale1 <- 1 / (n1 * epsilon1)
scale2 <- 1 / (n2 * epsilon2)
p1_hat <- p1_obs - scale1 * sign(u1 - 0.5) * log(1 - 2 * abs(u1 - 0.5))
p2_hat <- p2_obs - scale2 * sign(u2 - 0.5) * log(1 - 2 * abs(u2 - 0.5))

result <- prop_test_equiv_dp(
  p_hat = p1_hat,
  n = n1,
  p_hat2 = p2_hat,
  n2 = n2,
  lower = -0.1,
  upper = 0.1,
  epsilon = c(1, 1),
  B = 1000
)
result
#> ✔ Accept equivalence
#> Equiv. Region:  |----------------0----------------|
#> Estim. Inter.:      (-----------x-----------)      
#> CI = (-0.07409 ; 0.06052)
#> 
#> Method: DP-TOST (two-sample proportion)
#> alpha = 0.05; Equiv. bounds = [-0.1, 0.1]
#> epsilon = (1, 1); B = 1000
#> n = (300, 300); p_hat = (0.49216, 0.49930)

```
