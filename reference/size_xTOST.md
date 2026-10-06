# Size function of the xTOST Corrective Procedure

This function is used to calculate the size in the xTOST corrective
procedure.

## Usage

``` r
size_xTOST(sig_hat, delta, delta_star, ...)
```

## Arguments

- sig_hat:

  A `numeric` value (univariate) or `matrix` (multivariate)
  corresponding to the estimated variance of estimated `theta`.

- delta:

  A `numeric` value or vector defining the (bio)equivalence margin(s).
  The procedure assumes symmetry, i.e., the (bio)equivalence region is
  \\(-\delta, \delta)\\.

- delta_star:

  A `numeric` value specifying the corrected (bio)equivalence margin(s).

- ...:

  Additional arguments (currently unused).

## Value

The function returns a `numeric` value that corresponds to a
probability.
