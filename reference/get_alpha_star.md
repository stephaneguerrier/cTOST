# Compute Adjusted Significance Level (Alpha Star) for Univariate TOST

Calculates the adjusted significance level (\\\alpha^\*\\) for the
univariate Two One-Sided Tests (TOST) procedure.

## Usage

``` r
get_alpha_star(alpha = 0.05, sigma, nu, delta, ...)
```

## Arguments

- alpha:

  A `numeric` value specifying the target significance level (default:
  `alpha = 0.05`).

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

A `numeric` value corresponding to the optimized (adjusted) significance
level.
