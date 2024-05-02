
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

## `cTOST` Overview

This R package contains the functions to test equivalence in univariate
and multivariate settings based on the Two One-Sided Tests (TOST). In
addition, the package contains different corrective procedures applied
to the standard TOST in order to adjust the size of the procedure at the
desired nominal level leading to more powerful procedures. This package
implements the $\alpha$-TOST and $\delta$-TOST methods proposed in
Boulaguiem et al. (2024a) and in Boulaguiem et al. (2024b).

## Install Instructions

The `cTOST` package is available on GitHub at the moment. It is subject
to ongoing updates that may lead to stability issues.

In order to install the package, it is required to pre-install the
`devtools` dependency. Run the following command if you do not have it
already installed:

``` r
install.packages("devtools")
```

The package is then installed with the following command:

``` r
devtools::install_github("stephaneguerrier/cTOST")
```

Note that Windows users are assumed that have Rtools installed (if this
is not the case, please visit this
[link](https://cran.r-project.org/bin/windows/Rtools/).

## How to use

TO DO

## How to cite

    @Manual{boulaguiem2024ctost,
      title = {cTOST: Finite Sample Correction of The TOST in The Univariate Framework},
      author = {Boulaguiem, Y. and Insolia, L. and Couturier, D.-L. and Guerrier, S.},
      year = {2024},
      note = {R package version 1.1.0},
      url = {https://github.com/stephaneguerrier},
    }

## License

The license this source code is released under is the GNU AFFERO GENERAL
PUBLIC LICENSE (AGPL) v3.0. Please see the LICENSE file for full text.
Otherwise, please consult
[GNU](https://www.gnu.org/licenses/agpl-3.0.en.html) which will provide
a synopsis of the restrictions placed upon the code.

## References

Boulaguiem, Y., Quartier, J., Lapteva, M., Kalia, Y. N., Victoria-Feser,
M. P., Guerrier, S. & Couturier, D. L., “*Finite Sample Adjustments for
Average Equivalence Testing*”, Statistics in Medicine, 2024a,
<https://doi.org/10.1002/sim.9993>.

Boulaguiem, Y., Insolia, L., Victoria-Feser, M. P., Couturier, D. L. &
Guerrier, S., “*Multivariate Adjustments for Average Equivalence
Testing*”, submitted manuscript, 2024b.
