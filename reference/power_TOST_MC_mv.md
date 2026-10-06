# Power function of univariate or multivariate TOST using Monte Carlo integration

Power function of univariate or multivariate TOST using Monte Carlo
integration

## Usage

``` r
power_TOST_MC_mv(alpha, theta, Sigma, nu, delta, B = 10^5, seed = 10^8)
```

## Arguments

- alpha:

  A `numeric` value specifying the significance level.

- theta:

  A `numeric` value or vector representing the estimated difference(s)
  (e.g., between a generic and reference product).

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

- B:

  A `numeric` value specifying the number of Monte Carlo replication
  (default: B = `10^5`).

- seed:

  A `numeric` value specifying a seed for reproducibility (default: seed
  = `10^8`).

## Value

The function returns a `list` value with the structure:

- `power_univ`: A numerical vector that corresponds to a probability in
  univariate setting.

- `power_mult`: A numerical vector that corresponds to a probability in
  multivariate setting.

## Author

Younes Boulaguiem, Luca Insolia, Stéphane Guerrier, Dominique-Laurent
Couturier
