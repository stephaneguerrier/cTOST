# Changelog

## cTOST 1.1.0

### Maintainer

- Stéphane Guerrier is the new maintainer (previously Younes
  Boulaguiem). Development has moved to
  <https://github.com/stephaneguerrier/cTOST>, with the documentation
  site at <https://stephaneguerrier.github.io/cTOST/>.

### New features

- Multivariate average equivalence:
  [`ctost()`](https://stephaneguerrier.github.io/cTOST/reference/ctost.md)
  accepts a vector `theta` and a covariance matrix `sigma` and
  implements the multivariate alpha-TOST and cTOST of Boulaguiem et
  al. (2025, <doi:10.1002/sim.10258>).
- Quantile equivalence:
  [`qtost()`](https://stephaneguerrier.github.io/cTOST/reference/qtost.md)
  and
  [`compare_to_qtost()`](https://stephaneguerrier.github.io/cTOST/reference/compare_to_qtost.md)
  implement the qTOST and alpha-qTOST of Wu et al. (2025,
  <doi:10.48550/arXiv.2510.17514>) for one or two quantiles, from raw
  data or summary statistics.
- Differentially private equivalence tests for bounded means
  ([`tost_equiv_dp()`](https://stephaneguerrier.github.io/cTOST/reference/tost_equiv_dp.md),
  with a C++ implementation) and for proportions
  ([`prop_test_equiv_dp()`](https://stephaneguerrier.github.io/cTOST/reference/prop_test_equiv_dp.md)),
  following Pareek et al. (2026, <doi:10.48550/arXiv.2604.06499>).
- A [`plot()`](https://rdrr.io/r/graphics/plot.default.html) method for
  `tost`, `mtost`, `qtost` and `mqtost` objects.
- New datasets `ticlopidine`, `skin_mvt` and `biodistribution`.
- New vignettes: getting started, univariate and multivariate average
  equivalence, quantile equivalence, differential privacy, mathematical
  background and a reproduction of Boulaguiem et al. (2024).

### Breaking changes

- [`tost()`](https://stephaneguerrier.github.io/cTOST/reference/tost.md),
  `atost()` and `dtost()` are no longer exported. Use
  `ctost(method = "unadjusted")`, `ctost(method = "alpha")` and
  `ctost(method = "delta")`, which return the same results. Note that
  [`ctost()`](https://stephaneguerrier.github.io/cTOST/reference/ctost.md)
  takes the **variance** of `theta` as `sigma`, whereas `atost()` and
  `dtost()` took its standard error.
- `PowerTOST`, `mvtnorm` and `cli` moved from Depends to Imports:
  [`library(cTOST)`](https://stephaneguerrier.github.io/cTOST/) no
  longer attaches them.

### Numerical change

- The corrected level of the alpha-TOST and the corrected margin of the
  delta-TOST are again solved to `.Machine$double.eps^0.5`. Development
  versions after 1.0.1 solved them only to about 1e-4, so the fifth
  significant digit could differ from CRAN 1.0.1: on the `skin` example
  the corrected alpha is 0.078658 (CRAN 1.0.1: 0.078658; development:
  0.078649) and the corrected margin 0.254728 (0.254731; 0.254700).
  Decisions are unaffected.

### Bug fixes

- [`print()`](https://rdrr.io/r/base/print.html) works for every `qtost`
  object; two-quantile objects (class `mqtost`) have a print method;
  [`compare_to_qtost()`](https://stephaneguerrier.github.io/cTOST/reference/compare_to_qtost.md)
  accepts an alpha-qTOST object (it could not accept any input).
- [`plot()`](https://rdrr.io/r/graphics/plot.default.html) dispatches on
  multivariate and quantile objects; the decision marks fall back to
  `v`/`x` on devices and locales that cannot render the Unicode glyphs.
- `qtost(..., MC_sup = FALSE)` for two quantiles runs (a shadowed
  internal function always made it fail).
- [`ctost()`](https://stephaneguerrier.github.io/cTOST/reference/ctost.md)
  with the default offline correction no longer fails on combinations of
  `nu`, `sigma` and `alpha` that are missing from the precomputed table:
  it warns and uses the bootstrap correction. A non-positive corrected
  level (very small `nu` or `alpha`) is floored at 1e-6 with a warning
  instead of an error.
- `tost_equiv_dp(B = 0)` crashed the R session; all four exported
  functions now validate their inputs and fail with informative
  messages.
- [`compare_to_tost()`](https://stephaneguerrier.github.io/cTOST/reference/compare_to_tost.md)
  refuses multivariate objects with a message;
  [`print()`](https://rdrr.io/r/base/print.html) of multivariate results
  works for unnamed inputs; the multivariate cTOST no longer warns on
  every default call.
- [`ctost()`](https://stephaneguerrier.github.io/cTOST/reference/ctost.md),
  [`qtost()`](https://stephaneguerrier.github.io/cTOST/reference/qtost.md),
  [`tost_equiv_dp()`](https://stephaneguerrier.github.io/cTOST/reference/tost_equiv_dp.md)
  and
  [`prop_test_equiv_dp()`](https://stephaneguerrier.github.io/cTOST/reference/prop_test_equiv_dp.md)
  fix their seeds internally and now restore the caller’s random number
  stream on exit, so simulation loops around them no longer repeat
  identical data. `seed = NULL` is accepted by the DP tests.
- 1 x 1 matrix inputs to
  [`ctost()`](https://stephaneguerrier.github.io/cTOST/reference/ctost.md)
  no longer trigger recycling warnings on R \>= 4.6.
- [`prop_test_equiv_dp()`](https://stephaneguerrier.github.io/cTOST/reference/prop_test_equiv_dp.md)
  stops when no Monte Carlo replicate has a moment-matching solution
  instead of returning `NA`.

### Documentation and infrastructure

- R CMD check is clean (0 errors, 0 warnings) and runs on GitHub Actions
  for five platforms; a testthat suite with 283 expectations locks the
  results of every exported function; the pkgdown site is deployed from
  GitHub Actions.
- Non-portable compiler flags removed; documentation math, references
  and typos fixed; DOIs in CRAN form.

## cTOST 1.0.1

CRAN release: 2025-02-10

- CRAN release (2025-02-10).
