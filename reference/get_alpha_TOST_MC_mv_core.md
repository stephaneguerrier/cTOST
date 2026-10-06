# Get alpha star for multivariate TOST or xTOST

Get alpha star for multivariate TOST or xTOST

## Usage

``` r
get_alpha_TOST_MC_mv_core(
  alpha,
  Sigma,
  nu,
  delta,
  theta = NULL,
  B = 10^5,
  tol = .Machine$double.eps^0.5,
  seed = NULL,
  argsup_meth = "x",
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

  A `character` value specifying the value or vector representing the
  estimated difference(s) obtained by `argsup_meth` method under `NULL`.

- B:

  A `numeric` value specifying the number of Monte Carlo replication
  (default: B = `10^5`).

- tol:

  A `numeric` value specifying a tolerance level (default: tol =
  `.Machine$double.eps^0.5`).

- seed:

  A `character` value representing the use of default seed under `NULL`.

- argsup_meth:

  A `character` value representing the method to use (default:
  argsup_meth = `"x"`), see Details below for more information.

- ...:

  Additional parameters.

## Value

The function returns a `numeric` value that corresponds to the solution
of the optimization.

## Details

The supremum of the size over the boundary of the equivalence region is
found by direct optimisation (`argsup_meth = "x"`, the only implemented
option). The former is introduced in Boulaguiem et al. (2024,
\<https://doi.org/10.48550/arXiv.2411.16429\>) and the latter is
introduced in ...

## Author

Younes Boulaguiem, Luca Insolia, Stéphane Guerrier, Dominique-Laurent
Couturier
