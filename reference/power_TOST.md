# Power Function of Univariate TOST

Computes the power function for the univariate Two One-Sided Tests
(TOST) procedure.

## Usage

``` r
power_TOST(alpha, theta, sigma, nu, delta, ...)
```

## Arguments

- alpha:

  A `numeric` value specifying the significance level.

- theta:

  A `numeric` value representing the parameter of interest (e.g., a
  difference of means).

- sigma:

  A `numeric` value representing the estimated standard error of
  `theta`.

- nu:

  A `numeric` value specifying the degrees of freedom.

- delta:

  A `numeric` value defining the (bio)equivalence margin. The procedure
  assumes symmetry, i.e., the (bio)equivalence region is \\(-\delta,
  \delta)\\.

- ...:

  Additional parameters.

## Value

A `numeric` value corresponding to the probability (power) of the TOST
procedure.
