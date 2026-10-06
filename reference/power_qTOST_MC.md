# Power calculation for quantile TOST using Monte Carlo simulation

This function estimates the statistical power of the quantile-based Two
One-Sided Tests (qTOST) procedure using Monte Carlo simulations.

## Usage

``` r
power_qTOST_MC(
  theta,
  gamma = gamma,
  n_x,
  n_y,
  pi_x,
  delta_l,
  delta_u,
  alpha,
  B,
  seed,
  ...
)
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

- delta_l:

  A `numeric` value defining the equivalence lower limit.

- delta_u:

  A `numeric` value defining the equivalence upper limit.

- alpha:

  A `numeric` value specifying the significance level.

- B:

  An `integer` value specifying the number of Monte Carlo replications.

- seed:

  An `integer` value specifying the random seed for reproducibility.

- ...:

  Additional parameters

## Value

The function returns a `numeric` value that corresponds to a
probability.

## Details

The estimated power is the proportion of Monte Carlo replicates in which
equivalence is concluded.
