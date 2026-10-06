# Objective function to optimize for multivariate TOST

Objective function to optimize for multivariate TOST

## Usage

``` r
obj_fun_alpha_TOST_MC_mv(
  test,
  alpha,
  Sigma,
  nu,
  delta,
  theta = NULL,
  B = 10^5,
  seed = NULL,
  ...
)
```

## Arguments

- test:

  A `numeric` value specifying the significance level to optimize.

- alpha:

  A `numeric` value specifying the significance level.

- Sigma:

  A `numeric` value (univariate) or `matrix` (multivariate)
  corresponding to the estimated variance of estimated `theta`.

- nu:

  A `numeric` value specifying the degrees of freedom. In the
  multivariate case, it is assumed to be the same across all dimensions.

- delta:

  A `numeric` value or vector defining the (bio)equivalence margin(s).
  The procedure assumes symmetry, i.e., the (bio)equivalence region is
  \\(-\delta, \delta)\\.

- theta:

  A `numeric` value representing the estimated difference(s) (e.g.,
  between a generic and reference product) or a `character` value
  representing the use of equivalence margin(s) for \\\theta\\ under
  `NULL` (e.g., \\\theta\\ = \\(-\delta, \delta)\\).

- B:

  A `numeric` value specifying the number of Monte Carlo replication
  (default: B = `10^5`).

- seed:

  A `numeric` value specifying a seed for reproducibility representing
  the use of multivariate power or a `character` value representing the
  use of univariate power under `NULL`.

- ...:

  Additional parameters.

## Value

The function returns a `numeric` value for the objective function.

## Author

Younes Boulaguiem, Luca Insolia, Stéphane Guerrier, Dominique-Laurent
Couturier
