# Power function for multivariate xTOST procedure

This function is used to approximate the power by a standardized
multivariate normal vector lies within the (bio)equivalence margins.

## Usage

``` r
power_xTOST_mv(theta, delta, Sigma, alpha = 1/2, seed = 10^5)
```

## Arguments

- theta:

  A `numeric` value or vector representing the estimated difference(s)
  (e.g., between a generic and reference product).

- delta:

  A `numeric` value or vector defining the (bio)equivalence margin(s).
  The procedure assumes symmetry, i.e., the (bio)equivalence region is
  \\(-\delta, \delta)\\. In the multivariate case, it is assumed to be
  the same across all dimensions.

- Sigma:

  A `numeric` value (univariate) or `matrix` (multivariate)
  corresponding to the estimated variance of estimated `theta`.

- alpha:

  A `numeric` value specifying the significance level (default: alpha =
  `0.5`).

- seed:

  A `numeric` value specifying a seed for reproducibility (default: seed
  = `10^5`).

## Value

The function returns a `numeric` value that corresponds to a
probability.

## Author

Younes Boulaguiem, Luca Insolia, Stéphane Guerrier, Dominique-Laurent
Couturier
