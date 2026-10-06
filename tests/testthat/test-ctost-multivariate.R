# Characterization tests: multivariate average equivalence, ctost() on ticlopidine.
# See helper-data.R for the provenance of the reference values.

test_that("multivariate TOST on ticlopidine is locked and matches the closed form", {
  t = ticlopidine_stats()
  res = ctost(theta = t$theta, sigma = t$sigma, nu = t$nu, delta = log(1.25),
              method = "unadjusted")
  expect_s3_class(res, "mtost")
  expect_equal(res$setting, "multivariate")
  expect_equal(res$method, "TOST")
  expect_equal(unname(res$decision), c(TRUE, TRUE, TRUE, FALSE))
  ci_ref = matrix(c(-0.15767110413, -0.185532468496, -0.179142662887, -0.223791545503,
                    0.125026438423, 0.00991821734256, 0.0161961122682, 0.0215381798994),
                  ncol = 2)
  expect_equal(unname(res$ci), ci_ref, tolerance = tol_exact)
  se = sqrt(diag(t$sigma))
  expect_equal(unname(res$ci),
               unname(cbind(t$theta - qt(0.95, t$nu) * se, t$theta + qt(0.95, t$nu) * se)),
               tolerance = tol_exact)
  expect_equal(unname(res$decision), unname(apply(abs(res$ci) < log(1.25), 1, all)))
})

test_that("multivariate alpha-TOST on ticlopidine is locked (Monte Carlo, fixed internal seed)", {
  t = ticlopidine_stats()
  res = ctost(theta = t$theta, sigma = t$sigma, nu = t$nu, delta = log(1.25),
              method = "alpha")
  expect_s3_class(res, "mtost")
  expect_equal(res$method, "alpha-TOST")
  expect_true(all(res$decision))
  # cdd680d, default B, macOS value. alpha* is a Monte Carlo estimate (accuracy about
  # +/- 0.005, audit F34) and depends on the BLAS/LAPACK used by rmvnorm (0.0590800 on
  # Linux/Windows CI), hence tol_mc. The multivariate vignette prints 0.05909. The
  # value moves if the Monte Carlo size, the seed handling (F14) or the sup search
  # (F12) change.
  expect_equal(unname(res$corrected_alpha), 0.0590946010753, tolerance = tol_mc)
  ci_ref = matrix(c(-0.150099920315, -0.180297923278, -0.173911114853, -0.217221143689,
                    0.117455254607, 0.00468367212509, 0.010964564234, 0.0149677780853),
                  ncol = 2)
  expect_equal(unname(res$ci), ci_ref, tolerance = tol_mc)
  # the intervals use the corrected level per coordinate
  se = sqrt(diag(t$sigma))
  q = qt(1 - res$corrected_alpha, t$nu)
  expect_equal(unname(res$ci), unname(cbind(t$theta - q * se, t$theta + q * se)),
               tolerance = tol_exact)
})

test_that("multivariate cTOST on ticlopidine is locked and warns about the default correction", {
  t = ticlopidine_stats()
  # audit F09 (fixed in PR 4a): the default correction is 'none' in multivariate settings
  expect_no_warning(res <- ctost(theta = t$theta, sigma = t$sigma, nu = t$nu, delta = log(1.25)))
  expect_equal(res$correction, "none")
  expect_warning(ctost(theta = t$theta, sigma = t$sigma, nu = t$nu, delta = log(1.25),
                       correction = "bootstrap"), "none")
  expect_s3_class(res, "mtost")
  expect_equal(res$method, "cTOST")
  expect_true(all(res$decision))
  # cdd680d. c0 is reported to three decimals; audit F30 would move it (stopping tolerance).
  expect_equal(unname(res$c0), c(0.094, 0.134, 0.134, 0.111), tolerance = tol_exact)
  ci_ref = matrix(c(-0.145202660226, -0.176937417696, -0.170552533423, -0.213002287799,
                    0.112557994518, 0.00132316654268, 0.00760598280364, 0.0107489221954),
                  ncol = 2)
  expect_equal(unname(res$ci), ci_ref, tolerance = tol_exact)
})

test_that("multivariate ctost() rejects a vector delta and a non-square sigma", {
  t = ticlopidine_stats()
  expect_error(ctost(theta = t$theta, sigma = t$sigma, nu = t$nu, delta = rep(log(1.25), 4),
                     method = "unadjusted"))
  expect_error(ctost(theta = t$theta, sigma = t$sigma[, 1:3], nu = t$nu, delta = log(1.25),
                     method = "unadjusted"))
})
