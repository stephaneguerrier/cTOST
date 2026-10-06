# Get Corrected (Bio)Equivalence Bounds

This function applies the delta-TOST corrective procedure to obtain the
corrected (bio)equivalence bounds

## Usage

``` r
get_delta_TOST(
  alpha,
  sigma,
  nu,
  delta,
  l = 100,
  tol = .Machine$double.eps,
  ...
)
```

## Arguments

- alpha:

  A `numeric` value specifying the significance level.

- sigma:

  A `numeric` value corresponding to the standard error.

- nu:

  A `numeric` value corresponding to the number of degrees of freedom.

- delta:

  A `numeric` value corresponding to (bio)equivalence limit. We assume
  symmetry, i.e, the (bio)equivalence interval corresponds to
  (-delta,delta).

- l:

  A `numeric` value corresponding to the upper limit of (bio)equivalence
  margin to optimize (default: l = `100`).

- tol:

  A `numeric` value specifying a tolerance level (default:
  `tol = .Machine$double.eps`).

## Value

The function returns a `numeric` vector that minimized the objective
function.
