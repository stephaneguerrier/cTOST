
<!-- README.md is generated from README.Rmd. Please edit this file -->

# `cTOST` Overview <a href="https://yboulag.github.io/cTOST/"><img src="man/figures/hex-cTOST.png" alt="" align="right" height="138" width="125" /></a>

<!-- badges: start -->

[![Licence](https://img.shields.io/badge/licence-AGPL--3.0-blue.svg)](https://opensource.org/licenses/AGPL-3.0)
[![Last-changedate](https://img.shields.io/badge/last%20change-2025--07--10-green.svg)](https://github.com/yboulag/cTOST)
[![R-CMD-check](https://github.com/yboulag/cTOST/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/stephaneguerrier/cTOST/actions/workflows/R-CMD-check.yaml)
[![CRAN RStudio mirror
downloads](http://cranlogs.r-pkg.org/badges/cTOST)](https://www.r-pkg.org/pkg/cTOST)
[![CRAN RStudio mirror
downloads](https://cranlogs.r-pkg.org/badges/grand-total/cTOST)](https://www.r-pkg.org/pkg/cTOST)
<!-- badges: end -->

## 1. `cTOST` Overview

This R package contains functions for testing equivalence in both
univariate and multivariate settings, based on the Two One-Sided Tests
(TOST). The `cTOST` package implements the $\alpha$-TOST and
$\delta$-TOST methods proposed by Boulaguiem et al. (2024a, 2024b).
These two corrective procedures that can be applied to the standard TOST
to adjust the size of the procedure to the desired nominal level,
resulting in more powerful procedures.

## 2. Install Instructions

The `cTOST` package is available on cran and on GitHub. The version on
GitHub is subject to ongoing updates that may lead to stability issues.

The package can from cran as follows:

``` r
install.packages("cTOST")
```

The package can also be install from GitHub using the `devtools`
package:

``` r
install.packages("devtools")
devtools::install_github("stephaneguerrier/cTOST")
```

Note that Windows users are assumed that have Rtools installed (if this
is not the case, please visit this
[link](https://cran.r-project.org/bin/windows/Rtools/).

## 3. How to use

We provide here a few examples on the usage of the `cTOST` package.

### 3.1. Univariate settings

To illustrate the use of the proposed method in the univariate settings,
we consider the `skin` dataset analyzed in Boulaguiem et al. (2024a),
which can be loaded as follows:

``` r
data(skin)
theta_hat = diff(apply(skin,2,mean))
nu = nrow(skin) - 1
sig_hat = var(apply(skin,1,diff))/nu
```

### 3.1.1. Standard TOST

The standard TOST can be used as follows:

``` r
stost = tost(theta = theta_hat, sigma = sig_hat, nu = nu, delta = log(1.25))
stost
```

<img src="README_files/figure-gfm//unnamed-chunk-6.svg" width="100%" />

#### 3.1.2. cTOST

The cTOST can be used through the function `ctost` as follows:

``` r
opt_tost = ctost(theta = theta_hat, sigma = sig_hat, nu = nu, 
              delta = log(1.25), method = "optimal")
opt_tost
```

<img src="README_files/figure-gfm//unnamed-chunk-8.svg" width="100%" />

#### 3.1.3. $\alpha$-TOST

The $\alpha$-TOST can be used through the function `ctost` as follows:

``` r
atost = ctost(theta = theta_hat, sigma = sig_hat, nu = nu, 
              delta = log(1.25), method = "alpha")
atost
```

<img src="README_files/figure-gfm//unnamed-chunk-10.svg" width="100%" />

It is possible to compare the results of the $\alpha$-TOST (or
$\delta$-TOST, see below) with the standard TOST as follows:

``` r
compare_to_tost(atost)
```

<img src="README_files/figure-gfm//unnamed-chunk-12.svg" width="100%" />

### 3.1.4. $\delta$-TOST

The $\delta$-TOST can be used through the function `ctost` as follows:

``` r
dtost = ctost(theta = theta_hat, sigma = sig_hat, nu = nu, 
              delta = log(1.25), method = "delta")
dtost
#> ✖ Can't accept (bio)equivalence
#> Corr. Equiv. Region:  |----------------0----------------|
#>       Estim. Inter.:     (--------------x---------------)
#> CI =  (-0.21174 ; 0.25715)
#> 
#> Method: delta-TOST
#> alpha = 0.05; Equiv. lim. = +/- 0.22314
#> Corrected Equiv. lim. = +/- 0.25470
#> Mean = 0.02270; Stand. dev. = 0.13428; df = 16
```

<img src="README_files/figure-gfm//unnamed-chunk-14.svg" width="100%" />

### 3.2. Multivariate settings

``` r
data(ticlopidine)
n = nrow(ticlopidine)
p = ncol(ticlopidine)
nu = n-1
theta_hat = colMeans(ticlopidine)
Sigma_hat = cov(ticlopidine)/n
```

#### 3.2.1. Multivariate TOST

``` r
mtost = tost(theta = theta_hat, sigma = Sigma_hat, nu = nu, delta = log(1.25))
mtost
```

<img src="README_files/figure-gfm//unnamed-chunk-17.svg" width="100%" />

#### 3.2.2. Multivariate $\alpha$-TOST

``` r
matost = ctost(theta = theta_hat, sigma = Sigma_hat, nu = nu, delta = log(1.25))
matost
```

<img src="README_files/figure-gfm//unnamed-chunk-19.svg" width="100%" />

## 4. How to cite

    @Manual{boulaguiem2024ctost,
      title = {cTOST: Finite Sample Correction of The TOST in The Univariate Framework},
      author = {Boulaguiem, Y. and Insolia, L. and Couturier, D.-L. and Guerrier, S.},
      year = {2024},
      note = {R package version 1.1.0},
      url = {https://github.com/stephaneguerrier},
    }

## 5. License

The license this source code is released under is the GNU AFFERO GENERAL
PUBLIC LICENSE (AGPL) v3.0. Please see the LICENSE file for full text.
Otherwise, please consult
[GNU](https://www.gnu.org/licenses/agpl-3.0.en.html) which will provide
a synopsis of the restrictions placed upon the code.

## 6. References

Boulaguiem, Y., Quartier, J., Lapteva, M., Kalia, Y. N., Victoria-Feser,
M. P., Guerrier, S. & Couturier, D. L., “*Finite Sample Adjustments for
Average Equivalence Testing*”, Statistics in Medicine, 2024a,
<https://doi.org/10.1002/sim.9993>.

Boulaguiem, Y., Insolia, L., Victoria-Feser, M. P., Couturier, D. L. &
Guerrier, S., “*Multivariate Adjustments for Average Equivalence
Testing*”, submitted manuscript, 2024b.
