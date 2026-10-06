# Size in multivariate xTOST

This function is used to get the solution that maximizes the supremum
for a constrained xTOST procedure.

## Usage

``` r
find_sup_x(alpha, Sigma, delta, seed = 10^5)
```

## Arguments

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

The function returns a `numeric` vector that maximized the objective
function.

## Details

The optimization method "Brent" is used for 2-dimensional or less
problems only. The optimization method "Nelder-Mead" is used for higher
dimensions.

## Author

Younes Boulaguiem, Luca Insolia, Stéphane Guerrier, Dominique-Laurent
Couturier
