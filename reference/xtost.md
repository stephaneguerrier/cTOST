# Finite Sample Adjusted (Bio)Equivalence Testing of the xTOST Corrective Procedure

This function is used to compute finite sample corrected version of the
multivariate xTOST.

## Usage

``` r
xtost(
  theta_hat,
  sig_hat,
  nu,
  alpha,
  delta,
  correction = "none",
  B = 10^4,
  seed = 85
)
```

## Arguments

- theta_hat:

  A `numeric` value or vector representing the estimated difference(s)
  (e.g., between a generic and reference product).

- sig_hat:

  A `numeric` value corresponding to the estimated standard error of
  `theta_hat`.

- nu:

  A `numeric` value specifying the degrees of freedom. In the
  multivariate case, it is assumed to be the same across all dimensions.

- alpha:

  A `numeric` value specifying the significance level.

- delta:

  A `numeric` value or vector defining the (bio)equivalence margin(s).
  The procedure assumes symmetry, i.e., the (bio)equivalence region is
  \\(-\delta, \delta)\\.

- correction:

  A `character` value corresponding to the considered correction method,
  see Details below for more information (default: correction =
  `"none"`).

- B:

  A `numeric` value specifying the number of Monte Carlo replication
  (default: B = `10^4`).

- seed:

  A `numeric` value specifying a seed for reproducibility (default: seed
  = `85`).

## Value

An object of class `tost` with the structure:

- `decision`: A boolean variable indicating whether (bio)equivalence is
  accepted or not.

- `ci`: Confidence region at the \\1 - 2\alpha\\ level.

- `theta`: The estimated difference(s) used in the test.

- `sigma`: The estimated variance of `theta`, a `numeric` in univariate
  settings or `matrix` in multivariate settings.

- `nu`: The number of degrees of freedom used in the test.

- `alpha`: The significance level used in the test.

- `c0`: The estimated critical value.

- `correction`: The correction used in the test.

- `corrected_alpha`: The significance level corrected by the adjustment.

- `delta`: The (bio)equivalence limits used in the test.

- `method`: The method used in the test (default: method = "cTOST").

## Details

The significance level used to compute the critical value is corrected
additively, \\2\alpha - \mathrm{TIER}\\, where TIER is the type I error
rate of the uncorrected procedure at the equivalence boundary. With
`correction = "none"` no correction is applied; with
`correction = "bootstrap"` the TIER is estimated by parametric bootstrap
with `B` replications; with `correction = "offline"` it is read from a
precomputed table shipped with the package.

## See also

[`ctost`](https://stephaneguerrier.github.io/cTOST/reference/ctost.md)

## Examples

``` r
data(skin)

theta_hat = diff(apply(skin,2,mean))
nu = nrow(skin) - 1
sig_hat = sd(apply(skin,1,diff))/sqrt(nu)

# x-TOST
x_tost = cTOST:::xtost(theta_hat = theta_hat, sig_hat = sig_hat, nu = nu,
              alpha = 0.05, delta = log(1.25))
x_tost
#> ✔ Accept (bio)equivalence
#> Equiv. Region:  |----------------0----------------|
#> Estim. Inter.:      (--------------x-------------) 
#> CI =  (-0.16754 ; 0.21295)
#> 
#> Method: cTOST
#> alpha = 0.05; Equiv. lim. = +/- 0.22314
#> Estimated c(0) = 0.03290
#> Finite sample correction: none
#> Mean = 0.02270; Stand. dev. = 0.13428; df = 16
```
