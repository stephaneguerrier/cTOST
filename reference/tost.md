# Two One-Sided Tests (TOST) for (Bio)Equivalence Assessment

Performs the Two One-Sided Tests (TOST) procedure for (bio)equivalence
assessment in both univariate and multivariate settings.

## Usage

``` r
tost(theta, sigma, nu, delta, alpha = 0.05, ...)
```

## Arguments

- theta:

  A `numeric` value or vector representing the estimated difference(s)
  (e.g., between a generic and reference product).

- sigma:

  A `numeric` value (univariate) or `matrix` (multivariate)
  corresponding to the estimated variance of `theta`.

- nu:

  A `numeric` value specifying the degrees of freedom. In the
  multivariate case, it is assumed to be the same across all dimensions.

- delta:

  A `numeric` value or vector defining the (bio)equivalence margin(s).
  The procedure assumes symmetry, i.e., the (bio)equivalence region is
  \\(-\delta, \delta)\\. In the multivariate case, it is assumed to be
  the same across all dimensions.

- alpha:

  A `numeric` value specifying the significance level (default:
  `alpha = 0.05`).

- ...:

  Additional arguments.

## Value

A `tost` object with the following elements:

- `decision`: Logical; indicates whether (bio)equivalence is accepted.

- `ci`: Confidence region at the \\1 - 2\alpha\\ level.

- `theta`: The estimated difference(s) used in the test.

- `sigma`: The estimated variance of `theta`; a `numeric` value
  (univariate) or `matrix` (multivariate).

- `nu`: The degrees of freedom used in the test.

- `alpha`: The significance level used in the test.

- `delta`: The (bio)equivalence limits used in the test.

- `method`: A character string describing the method used ("TOST").

- `setting`: The setting used ("univariate" or "multivariate").

## Examples

``` r
# Univariate case
data(skin)
theta_hat = diff(colMeans(skin))
nu = nrow(skin) - 1
sig_hat = sd(apply(skin, 1, diff)) / sqrt(nu)
cTOST:::tost(theta = theta_hat, sigma = sig_hat, nu = nu,
     alpha = 0.05, delta = log(1.25))
#> ✖ Can't accept (bio)equivalence
#> Equiv. Region:  |---------------0---------------|  
#> Estim. Inter.:   (---------------x----------------)
#> CI =  (-0.21174 ; 0.25715)
#> 
#> Method: TOST
#> alpha = 0.05; Equiv. lim. = +/- 0.22314
#> Mean = 0.02270; Stand. dev. = 0.13428; df = 16

# Multivariate case
data(ticlopidine)
n = nrow(ticlopidine)
nu = n - 1
theta_hat = colMeans(ticlopidine)
Sigma_hat = cov(ticlopidine) / n
cTOST:::tost(theta = theta_hat, sigma = Sigma_hat, nu = nu, delta = log(1.25))
#> ✖ Can't accept (bio)equivalence
#> Equiv. Region:   |----------------0----------------|
#> t_half               (-----------x----------)       
#> AUC                 (-------x--------)              
#> AUC_inf             (--------x-------)              
#> C_max            (---------x---------)              
#> 
#> CIs:
#> t_half   (-0.15767 ; 0.12503)     
#> ✔
#> AUC      (-0.18553 ; 0.00992)     
#> ✔
#> AUC_inf  (-0.17914 ; 0.01620)     
#> ✔
#> C_max    (-0.22379 ; 0.02154)     
#> ✖
#> 
#> Method: TOST
#> alpha = 0.05; Equiv. lim. = +/- 0.22314
```
