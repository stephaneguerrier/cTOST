# Generate distribution of estimated parameter of interest

This function is used to generate estimated parameter of interest and
its standard deviation using Monte Carlo simulation.

## Usage

``` r
get_qTOST_rvs(theta, gamma, n_x, n_y, pi_x, B = 10^4, seed = 12345)
```

## Arguments

- theta:

  A `numeric` value corresponding to the parameter of interest.

- gamma:

  A `numeric` value corresponding to the ratio of sample variance
  between the target and reference group.

- n_x:

  An `integer` value representing the sample size of the reference
  group.

- n_y:

  An `integer` value representing the sample size of the target group.

- pi_x:

  A `numeric` value defining the quantile of interest.

- B:

  An optional `integer` value specifying the number of Monte Carlo
  replications (default: `B = 10^5`).

- seed:

  An optional `integer` value specifying the random seed for
  reproducibility (default: `seed = 12345`).

## Value

A list with the structure:

- theta_hat: B `numeric` values corresponding to the \\\hat{\theta}\\
  used in the test.

- sigma_theta_hat: B `numeric` values corresponding to the
  \\\hat{\sigma}\\ of \\\hat{\theta}\\ used in the test.
