# Get alpha star of the alpha-TOST Corrective Procedure

This function applies the alpha-TOST corrective procedure to obtain the
corrected level.

## Usage

``` r
get_alpha_TOST(
  alpha,
  sigma,
  nu,
  delta,
  l = 0.5,
  tol = .Machine$double.eps,
  ...
)
```

## Arguments

- alpha:

  A `numeric` value specifying the significance level.

- sigma:

  A `numeric` value corresponding to the estimated standard error of
  estimated `theta`.

- nu:

  A `numeric` value corresponding to the number of degrees of freedom.

- delta:

  A `numeric` value corresponding to (bio)equivalence limit. We assume
  symmetry, i.e, the (bio)equivalence interval corresponds to
  (-delta,delta).

- l:

  A `numeric` value corresponding to the upper limit of the significance
  level to optimize.

- tol:

  A `numeric` value specifying a tolerance level (default:
  `tol = .Machine$double.eps`).

- ...:

  Additional parameters.

## Value

A list with at least two components:

- `root`: A `numeric` value corresponding to the location of the root.

- `f.root`: A `numeric` value corresponding to the function evaluated at
  `root`.
