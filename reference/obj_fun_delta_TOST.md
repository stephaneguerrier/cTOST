# Objective Function of the delta-TOST Corrective Procedure

Objective Function of the delta-TOST Corrective Procedure

## Usage

``` r
obj_fun_delta_TOST(test, alpha, sigma, nu, delta, theta = NULL, ...)
```

## Arguments

- test:

  A `numeric` value specifying the significance level to optimize.

- alpha:

  A `numeric` value specifying the significance level.

- sigma:

  A `numeric` value corresponding to the estimated standard error of
  estimated `theta`.

- nu:

  A `numeric` value specifying the degrees of freedom.

- delta:

  A `numeric` value or vector defining the (bio)equivalence margin(s).
  The procedure assumes symmetry, i.e., the (bio)equivalence region is
  \\(-\delta, \delta)\\.

- theta:

  A `numeric` value representing the estimated difference(s) (e.g.,
  between a generic and reference product) or a `character` value
  representing the use of equivalence margin(s) for \\\theta\\ under
  `NULL` (e.g., \\\theta\\ = \\(-\delta, \delta)\\).

- ...:

  Additional parameters.

## Value

The function returns a `numeric` value for the objective function.
