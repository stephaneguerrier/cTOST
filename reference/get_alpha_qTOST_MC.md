# Get alpha star for quantile TOST

Get alpha star for quantile TOST

## Usage

``` r
get_alpha_qTOST_MC(
  gamma,
  n_x,
  n_y,
  pi_x,
  delta_l,
  delta_u,
  alpha,
  B,
  seed,
  tol = .Machine$double.eps,
  ...
)
```

## Arguments

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

- tol:

  An optional `numeric` value defining the tolerance for root-finding
  (default: `tol = .Machine$double.eps^0.5`).

- ...:

  Additional parameters

## Value

A list with at least two components:

- `root`: A `numeric` value corresponding to the location of the root.

- `f.root`: A `numeric` value corresponding to the function evaluated at
  `root`.
