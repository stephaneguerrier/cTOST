
<!-- README.md is generated from README.Rmd. Please edit this file -->

# `cTOST` Overview <a href="https://yboulag.github.io/cTOST/"><img src="man/figures/hex-cTOST.png" alt="" align="right" height="138" width="125" /></a>

<!-- badges: start -->

[![Licence](https://img.shields.io/badge/licence-AGPL--3.0-blue.svg)](https://opensource.org/licenses/AGPL-3.0)
[![Last-changedate](https://img.shields.io/badge/last%20change-2024--05--02-green.svg)](https://github.com/yboulag/cTOST)
[![R-CMD-check](https://github.com/yboulag/cTOST/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/yboulag/cTOST/actions/workflows/R-CMD-check.yaml)
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

<svg xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" width="840" height="300.52">
<rect width="840" height="300.52" rx="0" ry="0" class="a"></rect><svg height="260.52" viewBox="0 0 80 26.052" width="800" x="20" xmlns="http://www.w3.org/2000/svg" xmlns:xlink="http://www.w3.org/1999/xlink" y="20">
<style>.a{fill:rgb(40,45,53)}.b{font-family:'Fira Code',Monaco,Consolas,Menlo,'Bitstream Vera Sans Mono','Powerline Symbols',monospace}.c{fill:transparent}.d{fill:rgb(185,192,203);white-space:pre}.e{fill:rgb(232,131,136);white-space:pre}.f{fill:rgb(168,204,140);white-space:pre}</style>
<g font-family="&#x27;Fira Code&#x27;,Monaco,Consolas,Menlo,&#x27;Bitstream Vera Sans Mono&#x27;,&#x27;Powerline Symbols&#x27;,monospace" font-size="1.67" class="b"><defs><symbol id="a"><rect height="12" width="80" x="0" y="0" class="c"></rect></symbol></defs><rect height="26.052" width="80" class="a"></rect><svg x="0" y="0" width="80"><svg x="0"><use xlink:href="#a"></use><text font-size="1.67" x="0" y="1.67" class="d">TOST:</text><text font-size="1.67" x="12.024" y="1.67" class="e">✖</text><text font-size="1.67" x="14.027999999999999" y="1.67" class="e">Can't</text><text font-size="1.67" x="20.04" y="1.67" class="e">accept</text><text font-size="1.67" x="27.054" y="1.67" class="e">(bio)equivalence</text><text font-size="1.67" x="0" y="3.8409999999999997" class="d">alpha-TOST:</text><text font-size="1.67" x="12.024" y="3.8409999999999997" class="f">✔</text><text font-size="1.67" x="14.027999999999999" y="3.8409999999999997" class="f">Accept</text><text font-size="1.67" x="21.041999999999998" y="3.8409999999999997" class="f">(bio)equivalence</text><text font-size="1.67" x="0" y="8.183" class="d">Equiv.</text><text font-size="1.67" x="7.013999999999999" y="8.183" class="d">Region:</text><text font-size="1.67" x="16.032" y="8.183" class="d">\|—————0—————\|</text><text font-size="1.67" x="0" y="10.354" class="d">TOST:</text><text font-size="1.67" x="17.034" y="10.354" class="f">(—————x————–</text><text font-size="1.67" x="48.096" y="10.354" class="e">–)</text><text font-size="1.67" x="0" y="12.525" class="d">alpha-TOST:</text><text font-size="1.67" x="19.037999999999997" y="12.525" class="f">(————-x————–)</text><text font-size="1.67" x="17.034" y="16.866999999999997" class="d">CI</text><text font-size="1.67" x="20.04" y="16.866999999999997" class="d">-</text><text font-size="1.67" x="22.043999999999997" y="16.866999999999997" class="d">low</text><text font-size="1.67" x="31.061999999999998" y="16.866999999999997" class="d">CI</text><text font-size="1.67" x="34.068" y="16.866999999999997" class="d">-</text><text font-size="1.67" x="36.071999999999996" y="16.866999999999997" class="d">high</text><text font-size="1.67" x="0" y="19.037999999999997" class="d">TOST:</text><text font-size="1.67" x="17.034" y="19.037999999999997" class="d">-0.21174</text><text font-size="1.67" x="32.064" y="19.037999999999997" class="d">0.25715</text><text font-size="1.67" x="0" y="21.209000000000003" class="d">alpha-TOST:</text><text font-size="1.67" x="17.034" y="21.209000000000003" class="d">-0.17654</text><text font-size="1.67" x="32.064" y="21.209000000000003" class="d">0.22194</text><text font-size="1.67" x="0" y="25.550999999999995" class="d">Equiv.</text><text font-size="1.67" x="7.013999999999999" y="25.550999999999995" class="d">lim.</text><text font-size="1.67" x="12.024" y="25.550999999999995" class="d">=</text><text font-size="1.67" x="14.027999999999999" y="25.550999999999995" class="d">+/-</text><text font-size="1.67" x="18.035999999999998" y="25.550999999999995" class="d">0.22314</text>
</svg>
</svg>

</g></svg></svg>

We provide here a few examples on the usage of the `cTOST` package. More
information can be found
[here](https://stephaneguerrier.github.io/cTOST/)

### 3.1. Univariate settings

To illustrate the use of the proposed method in the univariate settings,
we consider the `skin` dataset analyzed in Boulaguiem et al. (2024a),
which can be loaded as follows:

``` r
data(skin)
theta_hat = diff(apply(skin,2,mean))
nu = nrow(skin) - 1
sig_hat = sd(apply(skin,1,diff))/sqrt(nu)
```

### 3.1.1. Standard TOST

The standard TOST can be used as follows:

``` r
stost = tost(theta = theta_hat, sigma = sig_hat, nu = nu, delta = log(1.25))
stost
```

<img src="README_files/figure-gfm//unnamed-chunk-6.svg" width="100%" />

#### 3.1.2. $\alpha$-TOST

The $\alpha$-TOST can be used through the function `ctost` as follows:

``` r
atost = ctost(theta = theta_hat, sigma = sig_hat, nu = nu, 
              delta = log(1.25), method = "alpha")
atost
```

<img src="README_files/figure-gfm//unnamed-chunk-8.svg" width="100%" />

It is possible to compare the results of the $\alpha$-TOST (or
$\delta$-TOST, see below) with the standard TOST as follows:

``` r
compare_to_tost(atost)
```

<img src="README_files/figure-gfm//unnamed-chunk-10.svg" width="100%" />

### 3.1.3. $\delta$-TOST

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
#> Corrected Equiv. lim. = +/- 0.25473
#> Mean = 0.02270; Stand. dev. = 0.13428; df = 16
```

<img src="README_files/figure-gfm//unnamed-chunk-12.svg" width="100%" />

### 3.2. Multivariate settings

COMING SOON

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
