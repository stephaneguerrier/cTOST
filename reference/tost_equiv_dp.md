# Differentially Private TOST for Mean Equivalence Testing

Performs equivalence testing for bounded means under differential
privacy using the DP-TOST procedure. Supports both one-sample and
two-sample tests.

## Usage

``` r
tost_equiv_dp(
  mean_private_obs,
  sd_private_obs,
  a,
  b,
  n,
  mean_private_obs2 = NULL,
  sd_private_obs2 = NULL,
  a2 = NULL,
  b2 = NULL,
  n2 = NULL,
  lower,
  upper,
  epsilon,
  alpha = 0.05,
  B = 10^4,
  seed = 1337,
  method = "cpp",
  ...
)
```

## Arguments

- mean_private_obs:

  Privatized sample mean for the first (or only) group. This should be
  the observed mean with differential privacy noise already added.

- sd_private_obs:

  Privatized sample standard deviation for the first (or only) group
  with differential privacy noise already added.

- a:

  Lower bound for data truncation in the first (or only) group.

- b:

  Upper bound for data truncation in the first (or only) group.

- n:

  Sample size for the first (or only) group.

- mean_private_obs2:

  Privatized sample mean for the second group (optional). If NULL, a
  one-sample test is performed. Default is NULL.

- sd_private_obs2:

  Privatized sample standard deviation for the second group (optional).
  Required if mean_private_obs2 is provided. Default is NULL.

- a2:

  Lower bound for data truncation in the second group (optional). If
  NULL but mean_private_obs2 is provided, defaults to a. Default is
  NULL.

- b2:

  Upper bound for data truncation in the second group (optional). If
  NULL but mean_private_obs2 is provided, defaults to b. Default is
  NULL.

- n2:

  Sample size for the second group (optional). Required if
  mean_private_obs2 is provided. Default is NULL.

- lower:

  Lower equivalence bound. For one-sample tests, this is the lower bound
  for the mean. For two-sample tests, this is the lower bound for the
  difference (\\\mu_1 - \mu_2\\). Typical value: -delta.

- upper:

  Upper equivalence bound. For one-sample tests, this is the upper bound
  for the mean. For two-sample tests, this is the upper bound for the
  difference (\\\mu_1 - \mu_2\\). Typical value: delta.

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

- seed:

  Random seed for reproducibility. Default is 1337.

- method:

  Implementation to use: "cpp" (default, ultra-fast) or "r" (pure R).
  The C++ implementation is 15-250x faster depending on problem size.

- ...:

  Additional arguments (currently unused).

## Value

A list of class "tost_dp" containing:

- decision:

  Logical. TRUE if equivalence is established, FALSE otherwise.

- conf.int:

  Confidence interval for the parameter of interest (mean for
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

  Point estimate(s) of the privatized mean(s).

- sample_size:

  Sample size(s).

- test_type:

  Character string indicating "one-sample" or "two-sample".

- method:

  Method used ("cpp" or "r").

- seed:

  Random seed used.

- B:

  Number of Monte Carlo replications.

## Details

This function implements the DP-TOST (Differentially Private Two
One-Sided Tests) procedure for testing equivalence of means while
maintaining differential privacy guarantees. The test uses Monte Carlo
simulation with moment-matching optimization to reconstruct the sampling
distribution from privatized statistics.

For the one-sample case, the test evaluates whether a mean is equivalent
to a reference range \[lower, upper\]. For the two-sample case, it tests
whether the difference between two means (\\\mu_1 - \mu_2\\) falls
within the equivalence bounds \[lower, upper\].

The privacy mechanism adds Laplace noise to both the observed mean and
standard deviation, with scales \\(b-a)/(n\\\epsilon/2)\\ and
\\(b-a)/(\sqrt{n-1}\\\epsilon/2)\\ respectively: the privacy budget
epsilon is split equally between the two statistics. Larger epsilon
values provide less privacy but more accurate inference. When the noise
on the standard deviation is comparable to the standard deviation itself
(small epsilon, small n or a wide interval \[a, b\]), the test is valid
but conservative: its type I error is below the nominal level and its
power is reduced.

The function fixes the random seed internally and restores the caller's
random number stream on exit; with `seed = NULL` a seed is drawn from
the caller's stream instead.

## Examples

``` r
# One-sample test: Is mean equivalent to range [1.5, 3.5]?
set.seed(123)
n <- 500
a <- 0
b <- 5
mu_true <- 2.5
sigma_true <- 1.0
epsilon <- 2

# Generate truncated normal data
z <- rnorm(n)
x <- pmin(pmax(mu_true + sigma_true * z, a), b)
mean_obs <- mean(x)
sd_obs <- sd(x)

# Add Laplace noise
delta_f_mean <- (b - a) / n
delta_f_sd <- (b - a) / sqrt(n - 1)
scale_mean <- delta_f_mean / (epsilon / 2)
scale_sd <- delta_f_sd / (epsilon / 2)

u <- runif(2)
mean_private <- mean_obs - scale_mean * sign(u[1] - 0.5) * log(1 - 2 * abs(u[1] - 0.5))
sd_private <- sd_obs - scale_sd * sign(u[2] - 0.5) * log(1 - 2 * abs(u[2] - 0.5))

result <- tost_equiv_dp(
  mean_private_obs = mean_private,
  sd_private_obs = sd_private,
  a = a, b = b, n = n,
  lower = 1.5,
  upper = 3.5,
  epsilon = 2,
  B = 1000
)
result
#> ✔ Accept equivalence
#> Equiv. Region:  |-------------------------------|
#> Estim. Inter.:                (---x--)             
#> CI = (2.45368 ; 2.61987)
#> 
#> Method: DP-TOST (mean)
#> alpha = 0.05; Equiv. bounds = [1.5, 3.5]
#> epsilon = 2; B = 1000
#> Truncation: (a, b) = (0, 5)

# Two-sample test: Is difference (mu1 - mu2) equivalent to [-1, 1]?
set.seed(456)
n1 <- 100
n2 <- 100
a1 <- 0; b1 <- 5
a2 <- 0; b2 <- 5
mu1_true <- 2.5
mu2_true <- 2.3
sigma1_true <- 1.0
sigma2_true <- 1.2
epsilon <- 1

# Generate data for both groups
z1 <- rnorm(n1)
z2 <- rnorm(n2)
x1 <- pmin(pmax(mu1_true + sigma1_true * z1, a1), b1)
x2 <- pmin(pmax(mu2_true + sigma2_true * z2, a2), b2)

# Add Laplace noise to both groups
scale_mean1 <- ((b1 - a1) / n1) / (epsilon / 2)
scale_sd1 <- ((b1 - a1) / sqrt(n1 - 1)) / (epsilon / 2)
scale_mean2 <- ((b2 - a2) / n2) / (epsilon / 2)
scale_sd2 <- ((b2 - a2) / sqrt(n2 - 1)) / (epsilon / 2)

u1 <- runif(2)
u2 <- runif(2)
mean_private1 <- mean(x1) - scale_mean1 * sign(u1[1] - 0.5) * log(1 - 2 * abs(u1[1] - 0.5))
sd_private1 <- sd(x1) - scale_sd1 * sign(u1[2] - 0.5) * log(1 - 2 * abs(u1[2] - 0.5))
mean_private2 <- mean(x2) - scale_mean2 * sign(u2[1] - 0.5) * log(1 - 2 * abs(u2[1] - 0.5))
sd_private2 <- sd(x2) - scale_sd2 * sign(u2[2] - 0.5) * log(1 - 2 * abs(u2[2] - 0.5))

result <- tost_equiv_dp(
  mean_private_obs = mean_private1,
  sd_private_obs = sd_private1,
  a = a1, b = b1, n = n1,
  mean_private_obs2 = mean_private2,
  sd_private_obs2 = sd_private2,
  a2 = a2, b2 = b2, n2 = n2,
  lower = -1,
  upper = 1,
  epsilon = 1,
  B = 1000
)
result
#> ✖ Can't accept equivalence
#> Equiv. Region:  |---------0---------|              
#> Estim. Inter.:    (---------------x---------------)
#> CI = (-0.80773 ; 2.65630)
#> 
#> Method: DP-TOST (two-sample mean)
#> alpha = 0.05; Equiv. bounds = [-1, 1]
#> epsilon = (1, 1); B = 1000
#> Truncation: (a1, a2) = (0, 0); (b1, b2) = (5, 5)

```
