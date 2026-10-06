# Objective function to optimize for quantile TOST

Objective function to optimize for quantile TOST

## Usage

``` r
obj_fun_qTOST_MC(
  test,
  theta,
  gamma,
  n_x,
  n_y,
  pi_x,
  delta_l,
  delta_u,
  alpha,
  B,
  seed
)
```

## Arguments

- test:

  A `numeric` value corresponding to the significance level to optimize.

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

## Value

The function returns a `numeric` value for the objective function.
