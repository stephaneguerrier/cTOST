# Objective function in multivariate xTOST

This function is used to get the maximum power for a constrained xTOST
procedure.

## Usage

``` r
argsup_x(theta, inds, alpha, Sigma, delta, seed = 10^5)
```

## Arguments

- theta:

  A `numeric` value or vector representing the estimated difference(s)
  (e.g., between a generic and reference product).

- inds:

  A `numeric` value or vector representing the interested `theta`.

- alpha:

  A `numeric` value specifying the significance level.

- Sigma:

  A `numeric` value (univariate) or `matrix` (multivariate)
  corresponding to the estimated variance of estimated `theta`.

- delta:

  A `numeric` value or vector defining the (bio)equivalence margin(s).
  The procedure assumes symmetry, i.e., the (bio)equivalence region is
  \\(-\delta, \delta)\\.

- seed:

  A `numeric` value specifying a seed for reproducibility (default: seed
  = `10^5`).

## Value

The function returns a `numeric` value that corresponds to a
probability.

## Details

The optimization method "Brent" is used for 2-dimensional or less
problems only. The optimization method "Nelder-Mead" is used for higher
dimensions.

## Author

Younes Boulaguiem, Luca Insolia, Stéphane Guerrier, Dominique-Laurent
Couturier
