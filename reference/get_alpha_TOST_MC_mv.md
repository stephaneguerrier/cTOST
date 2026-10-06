# Get alpha star for multivariate TOST

Get alpha star for multivariate TOST

## Usage

``` r
get_alpha_TOST_MC_mv(
  alpha,
  Sigma,
  nu,
  delta,
  theta = NULL,
  B = 10^5,
  tol = .Machine$double.eps^0.5,
  seed = NULL,
  max_iter = 10,
  tolpower = NULL,
  ...
)
```

## Arguments

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

  A `numeric` value specifying the initial value for the worst-case
  parameter configuration. If not provided, it's computed via
  find_sup_x().

- B:

  A `numeric` value specifying the number of Monte Carlo replication
  (default: B = `10^5`).

- tol:

  A `numeric` value specifying a tolerance level (default: tol =
  `.Machine$double.eps^0.5`).

- seed:

  A `character` value representing the use of default seed under `NULL`.

- max_iter:

  A `numeric` value specifying a maximum number of iteration to compute
  the supremum at which the size is assessed (default: `max_iter = 10`).

- tolpower:

  A `numeric` value specifying the power level to optimize, see Details
  below for more information.

- ...:

  Additional parameters.

## Value

The function returns a `list` value with the structure:

- `min`: A numerical vector that corresponds to the corrected
  significance level.

- `theta_sup`: A numerical vector that corresponds to the evaluated
  value of supremum.

- `alphas`: A numerical vector that corresponds to the significance
  level used in each iteration.

- `theta_sups`: A numerical matrix that corresponds to the `theta`
  values evaluated.

- `powers`: A numerical vector that corresponds to the estimated power .

- `iter`: A numerical variable that corresponds to the number of
  iterations performed.

- `err_power`: A numerical variable that corresponds to the final
  absolute error between simulated power and nominal significance level
  \\\alpha\\.

## Details

Tolerance of power allowed between simulated power and nominal alpha. If
not specified, it's automatically computed from the 1st and 99th
percentiles of a binomial distribution.

## Author

Younes Boulaguiem, Luca Insolia, Stéphane Guerrier, Dominique-Laurent
Couturier
