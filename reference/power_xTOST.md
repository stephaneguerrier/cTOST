# Power function of the xTOST Corrective Procedure

This function is used to calculate the power by xTOST.

## Usage

``` r
power_xTOST(theta, sig_hat, delta, ...)
```

## Arguments

- theta:

  A `character` value specifying the value or vector representing the
  estimated difference(s).

- sig_hat:

  A `numeric` value (univariate) or `matrix` (multivariate)
  corresponding to the estimated variance of estimated `theta`.

- delta:

  A `numeric` value or vector defining the (bio)equivalence margin(s).
  The procedure assumes symmetry, i.e., the (bio)equivalence region is
  \\(-\delta, \delta)\\.

- ...:

  Additional parameters.

## Value

The function returns a `numeric` value that corresponds to a
probability.
