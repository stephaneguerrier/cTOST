# Size of Univariate TOST

Computes the size (type I error rate) of the univariate Two One-Sided
Tests (TOST) procedure.

## Usage

``` r
size_TOST(alpha, sigma, nu, delta, ...)
```

## Arguments

- alpha:

  A `numeric` value specifying the significance level.

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

A `numeric` value corresponding to the probability (size) of the TOST
procedure.
