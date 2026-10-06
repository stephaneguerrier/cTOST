# Quantile equivalence testing procedures

Performs a Two One-Sided Test (TOST) to assess the equivalence of a
quantile from a test population (\\Y\\) with the corresponding quantile
from a reference population (\\X\\), assuming the data are normally
distributed. The test evaluates if the true quantile \\\pi_y\\ is within
a pre-specified equivalence margin \\\delta\\ around the reference
quantile \\\pi_x\\.

The null hypotheses for the two one-sided tests are: \\H\_{01}: \pi_y
\ge \pi_x - \delta\\ and \\H\_{02}: \pi_y \le \pi_x + \delta\\.
Equivalence is concluded if both null hypotheses are rejected.

## Usage

``` r
qtost(
  x,
  y,
  pi_x,
  delta,
  alpha = 0.05,
  method = "alpha",
  B = NULL,
  seed = 101010,
  tol = 1e-06,
  max_iter = 10,
  tolpower = 0.001,
  MC_sup = TRUE,
  ...
)
```

## Arguments

- x:

  A `numeric` vector of data for the reference group (\\X\\), or a
  `list` containing \`mean\`, \`sd\`, and \`n\`.

- y:

  A `numeric` vector of data for the test group (\\Y\\), or a `list`
  containing \`mean\`, \`sd\`, and \`n\`.

- pi_x:

  A `numeric` scalar or vector specifying the quantile(s) of interest in
  the reference group \\X\\ (e.g., 0.9 for the 90th percentile).

- delta:

  A `numeric` value or vector defining the (bio)equivalence margin(s).
  The procedure assumes symmetry, i.e., the (bio)equivalence region is
  \\(\pi_x-\delta, \pi_x+\delta)\\. In the multivariate case, it is
  assumed to be the same across all dimensions.

- alpha:

  A `numeric` value specifying the significance level, which must be
  between 0 and 0.5 (default: `alpha = 0.05`).

- method:

  A `character` string specifying the finite sample adjustment method.
  Available methods are: `"unadjusted"` (standard unadjusted qTOST),
  `"alpha"` (alpha-qTOST). Default: `method = "alpha"`.

- B:

  A `numeric` value specifying the number of Monte Carlo replications,
  required for the \`"alpha"\` method (default: `B = NULL`, which uses
  `10^5` replications for a single quantile and `10^4` for two
  quantiles).

- seed:

  A `numeric` value specifying a seed for reproducibility of the Monte
  Carlo step (default: `seed = 101010`). The caller's random number
  stream is restored on exit.

- tol:

  A `numeric` value specifying a tolerance level (default:
  `tol = 1e-6`).

- max_iter:

  A `numeric` value specifying a maximum number of iteration to compute
  the supremum at which the size is assessed (default: `max_iter = 10`).

- tolpower:

  A `numeric` value specifying the tolerance for power convergence when
  computing corrected alpha (default: `tolpower = 1e-3`).

- MC_sup:

  A `logical` value indicating whether to use Monte Carlo simulation to
  find the supremum when computing corrected alpha (default:
  `MC_sup = TRUE`).

- ...:

  Additional arguments (not currently used).

## Value

An object of class \`qtost\` (for a single quantile) or \`mqtost\` (for
multiple quantiles) with the following components:

- \`decision\`: The (component-wise) equivalence decision (\`TRUE\` or
  \`FALSE\`).

- \`method\`: The method used (\`"qTOST"\` or \`"alpha-qTOST"\`).

- \`ci\`: The \\(1 - 2\alpha)\\ confidence interval for the estimated
  quantile \\\hat{\pi}\_y\\.

- \`pi_y_hat\`: The point estimate for the quantile in group Y.

- \`eq_region\`: The defined equivalence region \\\[\pi_x - \delta,
  \pi_x + \delta\]\\.

- \`alpha\`: The nominal significance level.

- \`corrected_alpha\`: The adjusted alpha level used (for
  \`alpha-qTOST\` method).

## Examples

``` r
# https://www.accessdata.fda.gov/drugsatfda_docs/label/2024/021814s030lbl.pdf
# C_trough: table 5
# reference male group
x_bar_orig = 35.6
x_sd_orig = 16.7
n_x = 106
# target female group
y_bar_orig = 41.6
y_sd_orig = 24.3
n_y = 14
# data transformation
x_bar = log(x_bar_orig^2 / sqrt(x_bar_orig^2 + x_sd_orig^2))
x_sd = sqrt(log(1 + (x_sd_orig^2 / x_bar_orig^2)))
y_bar = log(y_bar_orig^2 / sqrt(y_bar_orig^2 + y_sd_orig^2))
y_sd = sqrt(log(1 + (y_sd_orig^2 / y_bar_orig^2)))
x = list(mean=x_bar, sd=x_sd, n=n_x)
y = list(mean=y_bar, sd=y_sd, n=n_y)
# qTOST
sqtost <- qtost(x, y, pi_x = 0.8, delta = 0.15, method = "unadjusted")
print(sqtost)
#> ✖ Can't accept quantile (bio)equivalence
#> Equiv. Region:            |---------------------|
#> Estim. Inter.: (------------x------------)        
#> CI =  (0.50105 ; 0.83712)
#> 
#> Method: qTOST
#> alpha = 0.05; Equiv. lim. = (0.65000 ; 0.95000)
#> theta_hat = 0.49265; Stand. dev. = 0.29791
aqtost <- qtost(x, y, pi_x = 0.8, delta = 0.15, method = "alpha")
print(aqtost)
#> ✖ Can't accept quantile (bio)equivalence
#> Equiv. Region:           |----------------------|
#> Estim. Inter.: (------------x------------)        
#> CI =  (0.51340 ; 0.82937)
#> 
#> Method: alpha-qTOST
#> alpha = 0.05; Equiv. lim. = (0.65000 ; 0.95000)
#> Corrected alpha = 0.06167
#> theta_hat = 0.49265; Stand. dev. = 0.29791

# Example 2: Using raw data to test a single quantile
set.seed(12345)
x_data <- rnorm(n = 30, mean = 0, sd = 0.1)
y_data <- rnorm(n = 10, mean = 0, sd = 0.1)

# Test the 90th percentile with a margin of delta = 0.05
sqtost <- qtost(x = x_data, y = y_data, pi_x = 0.8, delta = 0.15, method = "unadjusted")
print(sqtost)
#> ✖ Can't accept quantile (bio)equivalence
#> Equiv. Region:                  |---------------|
#> Estim. Inter.: (-----------x------------)         
#> CI =  (0.31955 ; 0.75964)
#> 
#> Method: qTOST
#> alpha = 0.05; Equiv. lim. = (0.65000 ; 0.95000)
#> theta_hat = 0.11810; Stand. dev. = 0.35690
aqtost <- qtost(x_data, y_data, pi_x = 0.8, delta = 0.15, method = "alpha")
print(aqtost)
#> ✖ Can't accept quantile (bio)equivalence
#> Equiv. Region:                |-----------------|
#> Estim. Inter.: (----------x----------)            
#> CI =  (0.37345 ; 0.71189)
#> 
#> Method: alpha-qTOST
#> alpha = 0.05; Equiv. lim. = (0.65000 ; 0.95000)
#> Corrected alpha = 0.10839
#> theta_hat = 0.11810; Stand. dev. = 0.35690

# Example 2: Using summary statistics with the unadjusted method
x_stats <- list(mean = 100, sd = 15, n = 50)
y_stats <- list(mean = 102, sd = 16, n = 50)

result_summary <- qtost(x_stats, y_stats, pi_x = 0.95, delta = 0.03,
                        method = "unadjusted")
print(result_summary)
#> ✖ Can't accept quantile (bio)equivalence
#> Equiv. Region:                    |-------------|
#> Estim. Inter.: (---------------x---------------)  
#> CI =  (0.82835 ; 0.97038)
#> 
#> Method: qTOST
#> alpha = 0.05; Equiv. lim. = (0.92000 ; 0.98000)
#> theta_hat = 1.41705; Stand. dev. = 0.28537
```
