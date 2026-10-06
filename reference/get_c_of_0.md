# Critical value search in the xTOST Corrective Procedure

This function is used to calculate the critical value such that the size
of the xTOST procedure matches the nominal significance level.

## Usage

``` r
get_c_of_0(delta, sigma, alpha, B = 1000, tol = 10^(-8), l = 1, optim = "NR")
```

## Arguments

- delta:

  A `numeric` value defining the (bio)equivalence margin(s). The
  procedure assumes symmetry, i.e., the (bio)equivalence is \\(-\delta,
  \delta)\\.

- sigma:

  A `numeric` value corresponding to the estimated variance of estimated
  `theta`.

- alpha:

  A `numeric` value specifying the significance level.

- B:

  A `numeric` value specifying the number of iterations using
  Newton-Raphson method (default: B = `1000`).

- tol:

  A `numeric` value specifying a tolerance level (default: tol =
  `10^(-8)`).

- l:

  A `numeric` value corresponding to the upper limit of the significance
  level to optimize.

- optim:

  A `character` value representing the method to use (default: optim =
  `"NR"`), see Details below for more information.

## Value

The function returns a `list` value with the structure:

- `c`: A numerical variable that corresponds to the estimated critical
  value.

- `size`: A numerical variable that corresponds to the size when using
  the estimated critical value.

- `converged`: A boolean variable that corresponds to whether
  Newton-Raphson converged (only returned if `optim = "NR"`).

- `iter`: A numerical variable that corresponds to the actual iterations
  used (only returned if `optim = "NR"`).

## Details

Only the Newton-Raphson method (`optim = "NR"`) is implemented.
