# Objective Function for Optimization in Univariate TOST

Computes the objective function used to optimize the significance level
in the univariate Two One-Sided Tests (TOST) procedure.

## Usage

``` r
obj_fun_alpha_star(test, alpha = 0.05, sigma, nu, delta, ...)
```

## Arguments

- test:

  A `numeric` value specifying the significance level to be optimized.

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

A `numeric` value representing the objective function to be minimized.
