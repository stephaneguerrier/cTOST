# Finite Sample Adjustment for Average (Bio)Equivalence Assessment

Computes finite sample corrected versions of the standard (univariate or
multivariate) TOST, as developed in Boulaguiem et al. (2024,
https://doi.org/10.1002/sim.9993), Boulaguiem et al. (2025,
https://doi.org/10.1002/sim.10258), and Insolia et al. (2025,
https://doi.org/10.48550/arXiv.2507.22756).

## Usage

``` r
ctost(
  theta,
  sigma,
  nu,
  delta,
  alpha = 0.05,
  method = "optimal",
  B = 10^4,
  seed = 101010,
  correction = NULL,
  ...
)
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

  A `numeric` value specifying the significance level, which must be
  between 0 and 0.5 (default: `alpha = 0.05`).

- method:

  A `character` string specifying the finite sample adjustment method.
  Available methods are: `"unadjusted"` (standard unadjusted TOST),
  `"alpha"` (alpha-TOST), `"delta"` (delta-TOST, not implemented for
  multivariate settings), and `"optimal"` (cTOST, default). See Details.

- B:

  A `numeric` value specifying the number of Monte Carlo replications,
  required for some methods (default: `B = 10^4`).

- seed:

  A `numeric` value specifying a seed for reproducibility of the Monte
  Carlo and bootstrap steps (default: `seed = 101010`). The caller's
  random number stream is restored on exit.

- correction:

  A `character` string specifying the finite sample correction of the
  significance level used by the cTOST (`method = "optimal"`) in
  univariate settings: `"offline"` (default; precomputed table),
  `"bootstrap"` or `"none"`. No correction is available in multivariate
  settings, where `"none"` is used.

- ...:

  Additional parameters.

## Value

An object of class `tost` (univariate) or `mtost` (multivariate) with
the following components:

- `decision`: Logical; indicates whether (bio)equivalence is accepted.

- `ci`: Confidence region at the \\1 - 2\alpha\\ level.

- `theta`: The estimated difference(s) used in the test.

- `sigma`: The estimated variance of `theta`; a `numeric` value
  (univariate) or `matrix` (multivariate).

- `nu`: The degrees of freedom used in the test.

- `alpha`: The significance level used in the test.

- `corrected_alpha`: The significance level after adjustment (if
  `method = "alpha"` or `"optimal"`).

- `corrected_c`: The corrected critical value \\c(0)\\ (univariate
  cTOST); `c0` holds the per-coordinate values in multivariate settings.

- `correction`: The finite sample correction applied (cTOST only).

- `corrected_ci`: The TOST interval re-expressed on the nominal
  \\\pm\delta\\ scale (if `method = "delta"`).

- `corrected_delta`: The (bio)equivalence limits after adjustment (if
  `method = "delta"`).

- `delta`: The (bio)equivalence limits used in the test.

- `method`: The adjustment method used (optimal, alpha-TOST, or
  delta-TOST).

- `setting`: The setting used ("univariate" or "multivariate").

## Details

In univariate settings, three adjustment methods are available: optimal
(cTOST, `method = "optimal"`) as proposed in Insolia et al. (2025),
alpha-TOST (`method = "alpha"`), and delta-TOST (`method = "delta"`),
both proposed in Boulaguiem et al. (2024,
https://doi.org/10.1002/sim.9993, https://doi.org/10.1002/sim.10258).
For multivariate settings, only cTOST and alpha-TOST are implemented.

The cTOST, alpha-TOST, and delta-TOST methods apply different finite
sample adjustments. Alpha-TOST corrects the significance level, while
delta-TOST adjusts the equivalence limits. The cTOST method is based on
a more complex approach, and in small samples (typically less than 30),
additional corrections may be beneficial. The `correction` argument
further adjusts the test level to prevent liberal inference; see Insolia
et al. (2025) for details.

Generally, cTOST outperforms other methods, with alpha-TOST performing
better than delta-TOST. For this reason, delta-TOST is not implemented
for multivariate settings and is not recommended.

In multivariate settings, the alpha-TOST calibrates the corrected level
by Monte Carlo with `B` replications and a fixed internal seed; the
achieved size is accurate to roughly \\\pm
2.3\sqrt{\alpha(1-\alpha)/B}\\ (about 0.005 at `B = 10^4`) and can
differ in the fourth decimal between platforms. The multivariate cTOST
computes its critical values under a normal approximation that does not
depend on `nu`, and no finite sample correction of the level is
available for it yet.

## Examples

``` r
data(skin)
theta_hat = diff(apply(skin, 2, mean))
nu = nrow(skin) - 1
sig_hat = var(apply(skin, 1, diff)) / nu

# alpha-TOST
atost = ctost(theta = theta_hat, sigma = sig_hat, nu = nu,
              alpha = 0.05, delta = log(1.25), method = "alpha")
atost
#> ✔ Accept (bio)equivalence
#> Equiv. Region:  |----------------0----------------|
#> Estim. Inter.:     (---------------x--------------)
#> CI =  (-0.17654 ; 0.22194)
#> 
#> Method: alpha-TOST
#> alpha = 0.05; Equiv. lim. = +/- 0.22314
#> Corrected alpha = 0.07866
#> Mean = 0.02270; Stand. dev. = 0.13428; df = 16
compare_to_tost(atost)
#> TOST:       
#> ✖ Can't accept (bio)equivalence
#> alpha-TOST: 
#> ✔ Accept (bio)equivalence
#> 
#> Equiv. Region:  |---------------0---------------|  
#> TOST:            (---------------x----------------)
#> alpha-TOST:        (-------------x--------------)  
#> 
#>                  CI - low      CI - high
#> TOST:            -0.21174       0.25715
#> alpha-TOST:      -0.17654       0.22194
#> 
#> Equiv. lim. = +/- 0.22314

# delta-TOST
dtost = ctost(theta = theta_hat, sigma = sig_hat, nu = nu,
              alpha = 0.05, delta = log(1.25), method = "delta")
dtost
#> ✖ Can't accept (bio)equivalence
#> Corr. Equiv. Region:  |----------------0----------------|
#>       Estim. Inter.:     (--------------x---------------)
#> CI =  (-0.21174 ; 0.25715)
#> 
#> Method: delta-TOST
#> alpha = 0.05; Equiv. lim. = +/- 0.22314
#> Corrected Equiv. lim. = +/- 0.25473
#> Mean = 0.02270; Stand. dev. = 0.13428; df = 16
compare_to_tost(dtost)
#> TOST:       
#> ✖ Can't accept (bio)equivalence
#> delta-TOST: 
#> ✖ Can't accept (bio)equivalence
#> 
#> Stand. Equiv. Region:    |--------------0--------------|  
#>  Corr. Equiv. Region:  |----------------0----------------|
#>                 TOST:     (--------------x---------------)
#>           delta-TOST:     (--------------x---------------)
#> 
#>  Standard Equiv. lim. = +/- 0.22314
#> Corrected Equiv. lim. = +/- 0.25473
```
