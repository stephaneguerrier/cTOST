# Simulate data from the canonical form of the univariate average equivalence problem

Simulate data from the canonical form of the univariate average
equivalence problem

## Usage

``` r
simulate_data(mu, sigma, nu, seed = 18)
```

## Arguments

- mu:

  Numeric value indicating the true mean of the data.

- sigma:

  Numeric value indicating the standard deviation of the data.

- nu:

  Numeric value indicating the degrees of freedom parameter of the
  considered setting.

- seed:

  Numeric value (optional) used to set the seed for reproducibility of
  random number generation. Default value is 18.

## Value

A list containing the estimated mean and estimated standard deviation of
the simulated data.

## Examples

``` r
cTOST:::simulate_data(mu = 10, sigma = 2, nu = 10)
#> $theta_hat
#> [1] 11.85292
#> 
#> $sig_hat
#> [1] 2.712557
#> 
```
