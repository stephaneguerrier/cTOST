# Confidence Interval for Univariate TOST

Computes the confidence interval for the univariate Two One-Sided Tests
(TOST) procedure.

## Usage

``` r
ci(alpha, theta, sigma, nu, ...)
```

## Arguments

- alpha:

  A `numeric` value specifying the significance level.

- theta:

  A `numeric` value representing the estimated parameter of interest
  (e.g., a difference of means).

- sigma:

  A `numeric` value representing the estimated standard error of
  `theta`.

- nu:

  A `numeric` value specifying the degrees of freedom.

- ...:

  Additional parameters.

## Value

A numeric `vector` containing the lower and upper bounds of the
confidence interval.
