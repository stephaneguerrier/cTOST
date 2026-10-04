# Shared inputs and tolerances for the characterization tests.
#
# Provenance of the reference values quoted in the test files:
# - "cdd680d": computed from commit cdd680d (main, 2025-11-27), built without
#   src/Makevars (PR 1), R 4.6.0, macOS arm64, 2026-10-04.
# - "CRAN 1.0.1": computed from the CRAN 1.0.1 source tarball installed on the
#   same machine. In 1.0.1 the sigma argument of tost()/atost()/dtost() is the
#   standard error; in 1.1.0 ctost()'s sigma is the variance.
# - "README": numbers printed in README.md, which was knitted from cdd680d.

# Univariate: skin data (econazole), as in ?ctost and the README
skin_stats = function() {
  skin = cTOST::skin
  theta_hat = diff(apply(skin, 2, mean))
  nu = nrow(skin) - 1
  sig_hat = var(apply(skin, 1, diff)) / nu   # variance of the mean difference
  list(theta = theta_hat, sigma = sig_hat, nu = nu)
}

# Multivariate: ticlopidine data (log-differences), as in the README
ticlopidine_stats = function() {
  ticlopidine = cTOST::ticlopidine
  n = nrow(ticlopidine)
  list(theta = colMeans(ticlopidine), sigma = cov(ticlopidine) / n, nu = n - 1)
}

# Quantile: FDA C_trough example from ?qtost (summary statistics)
qtost_fda_inputs = function() {
  x_bar_orig = 35.6; x_sd_orig = 16.7; n_x = 106
  y_bar_orig = 41.6; y_sd_orig = 24.3; n_y = 14
  x_bar = log(x_bar_orig^2 / sqrt(x_bar_orig^2 + x_sd_orig^2))
  x_sd = sqrt(log(1 + (x_sd_orig^2 / x_bar_orig^2)))
  y_bar = log(y_bar_orig^2 / sqrt(y_bar_orig^2 + y_sd_orig^2))
  y_sd = sqrt(log(1 + (y_sd_orig^2 / y_bar_orig^2)))
  list(x = list(mean = x_bar, sd = x_sd, n = n_x),
       y = list(mean = y_bar, sd = y_sd, n = n_y))
}

# Tolerances (relative, as interpreted by testthat's expect_equal)
tol_exact = 1e-10   # deterministic values locked from cdd680d (quoted to 12 significant digits)
tol_cran101 = 1e-3  # alpha-/delta-TOST versus CRAN 1.0.1: 1.1.0 solves the corrected
                    # level/margin only to uniroot's default tolerance (audit F27).
                    # Tighten to 1e-7 once the tolerance is passed to uniroot again.
