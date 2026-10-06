# Characterization tests: univariate average equivalence, ctost().
# See helper-data.R for the provenance of the reference values.

test_that("standard TOST on skin is identical to CRAN 1.0.1 and to the closed form", {
  s = skin_stats()
  res = ctost(theta = s$theta, sigma = s$sigma, nu = s$nu, delta = log(1.25),
              method = "unadjusted")
  expect_s3_class(res, "tost")
  expect_equal(res$method, "TOST")
  expect_equal(res$setting, "univariate")
  # CRAN 1.0.1: tost(theta_hat, sig_hat, nu, 0.05, log(1.25))$ci
  expect_equal(unname(res$ci), c(-0.211741483133, 0.257145789481), tolerance = tol_exact)
  expect_false(unname(res$decision))
  # README: CI = (-0.21174 ; 0.25715)
  expect_equal(round(unname(res$ci), 5), c(-0.21174, 0.25715))
  # closed form: theta +/- t_{1 - alpha, nu} * SE, decision = CI inside +/- delta
  se = sqrt(s$sigma)
  expect_equal(unname(res$ci), unname(s$theta) + c(-1, 1) * qt(0.95, s$nu) * se,
               tolerance = tol_exact)
  expect_equal(unname(res$decision), all(abs(res$ci) < log(1.25)))
})

test_that("alpha-TOST on skin matches CRAN 1.0.1 (audit F27 tolerance) and is locked", {
  s = skin_stats()
  res = ctost(theta = s$theta, sigma = s$sigma, nu = s$nu, delta = log(1.25),
              method = "alpha")
  expect_s3_class(res, "tost")
  expect_equal(res$method, "alpha-TOST")
  expect_true(unname(res$decision))
  # CRAN 1.0.1: atost(theta_hat, sig_hat, nu, 0.05, log(1.25))
  expect_equal(unname(res$corrected_alpha), 0.078657919501, tolerance = tol_cran101)
  expect_equal(unname(res$ci), c(-0.176536587848, 0.221940894197), tolerance = tol_cran101)
  # README: Corrected alpha = 0.07865; CI = (-0.17655 ; 0.22195)
  expect_equal(round(unname(res$corrected_alpha), 5), 0.07865)
  expect_equal(round(unname(res$ci), 5), c(-0.17655, 0.22195))
  # cdd680d (differs from 1.0.1 in the 6th significant digit, audit F27)
  expect_equal(unname(res$corrected_alpha), 0.0786485339274427, tolerance = tol_exact)
  expect_equal(unname(res$ci), c(-0.176546187381565, 0.221950493730324), tolerance = tol_exact)
  # the interval uses the corrected level
  expect_equal(unname(res$ci),
               unname(s$theta) + c(-1, 1) * qt(1 - res$corrected_alpha, s$nu) * sqrt(s$sigma),
               tolerance = tol_exact)
})

test_that("delta-TOST on skin matches CRAN 1.0.1, keeps the TOST interval and is locked", {
  s = skin_stats()
  res = ctost(theta = s$theta, sigma = s$sigma, nu = s$nu, delta = log(1.25),
              method = "delta")
  expect_s3_class(res, "tost")
  expect_equal(res$method, "delta-TOST")
  expect_false(unname(res$decision))
  # CRAN 1.0.1: dtost(theta_hat, sig_hat, nu, 0.05, log(1.25))
  expect_equal(unname(res$corrected_delta), 0.254730957152, tolerance = tol_cran101)
  # the interval is the plain TOST interval (identical in 1.0.1)
  expect_equal(unname(res$ci), c(-0.211741483133, 0.257145789481), tolerance = tol_exact)
  # README: Corrected Equiv. lim. = +/- 0.25470
  expect_equal(round(unname(res$corrected_delta), 5), 0.25470)
  # cdd680d
  expect_equal(unname(res$corrected_delta), 0.254700490986057, tolerance = tol_exact)
  expect_equal(unname(res$corrected_ci), c(-0.180184543460833, 0.225588849809592),
               tolerance = tol_exact)
  # corrected_ci is the TOST interval re-expressed on the nominal +/- delta scale
  expect_equal(unname(res$corrected_ci),
               unname(res$ci) + c(1, -1) * (res$corrected_delta - res$delta),
               tolerance = tol_exact)
  expect_equal(unname(res$decision), all(abs(res$ci) < res$corrected_delta))
})

test_that("alpha-TOST and delta-TOST match CRAN 1.0.1 on a synthetic cell", {
  # CRAN 1.0.1: atost(theta = 0.05, sigma = 0.08, nu = 10, alpha = 0.05, delta = log(1.25))
  # (sigma = SE there; variance 0.08^2 here)
  a = ctost(theta = 0.05, sigma = 0.08^2, nu = 10, delta = log(1.25), method = "alpha")
  expect_equal(unname(a$corrected_alpha), 0.050193259564, tolerance = tol_cran101)
  expect_equal(unname(a$ci), c(-0.094807707292, 0.194807707292), tolerance = tol_cran101)
  d = ctost(theta = 0.05, sigma = 0.08^2, nu = 10, delta = log(1.25), method = "delta")
  expect_equal(unname(d$corrected_delta), 0.223314094889, tolerance = tol_cran101)
  expect_equal(unname(d$ci), c(-0.094996889825, 0.194996889825), tolerance = tol_cran101)
})

test_that("cTOST on skin (offline default, bootstrap, none) is locked to cdd680d", {
  s = skin_stats()

  off = ctost(theta = s$theta, sigma = s$sigma, nu = s$nu, delta = log(1.25))
  expect_s3_class(off, "tost")
  expect_equal(off$method, "cTOST")
  expect_equal(off$correction, "offline")
  expect_true(unname(off$decision))
  expect_equal(unname(off$corrected_alpha), 0.038532, tolerance = tol_exact)
  expect_equal(unname(off$corrected_c), 0.025525302221127, tolerance = tol_exact)
  expect_equal(unname(off$ci), c(-0.174916095918703, 0.220320402267462), tolerance = tol_exact)
  # README: Estimated c(0) = 0.02553; Corrected alpha = 0.03853; CI = (-0.17492 ; 0.22032)
  expect_equal(round(unname(off$corrected_c), 5), 0.02553)
  expect_equal(round(unname(off$corrected_alpha), 5), 0.03853)
  expect_equal(round(unname(off$ci), 5), c(-0.17492, 0.22032))
  # cdd680d stores the estimate as 'theta_hat' for this method only (audit F79)
  expect_true("theta_hat" %in% names(off))
  # interval = theta +/- (delta - c), decision = interval inside +/- delta
  expect_equal(unname(off$ci), unname(s$theta) + c(-1, 1) * (off$delta - off$corrected_c),
               tolerance = tol_exact)
  expect_equal(unname(off$decision), all(abs(off$ci) < off$delta))

  boot = ctost(theta = s$theta, sigma = s$sigma, nu = s$nu, delta = log(1.25),
               correction = "bootstrap")
  expect_equal(boot$correction, "bootstrap")
  expect_true(unname(boot$decision))
  expect_equal(unname(boot$corrected_alpha), 0.0404, tolerance = tol_exact)
  expect_equal(unname(boot$corrected_c), 0.0267358462950114, tolerance = tol_exact)
  expect_equal(unname(boot$ci), c(-0.173705551844819, 0.219109858193578), tolerance = tol_exact)

  none = ctost(theta = s$theta, sigma = s$sigma, nu = s$nu, delta = log(1.25),
               correction = "none")
  expect_equal(none$correction, "none")
  expect_equal(unname(none$corrected_alpha), 0.05)
  expect_equal(unname(none$corrected_c), 0.032897629062446, tolerance = tol_exact)
  expect_equal(unname(none$ci), c(-0.167543769077384, 0.212948075426143), tolerance = tol_exact)
})

test_that("ctost() rejects an invalid method or correction", {
  s = skin_stats()
  expect_error(ctost(theta = s$theta, sigma = s$sigma, nu = s$nu, delta = log(1.25),
                     method = "nope"))
  expect_error(ctost(theta = s$theta, sigma = s$sigma, nu = s$nu, delta = log(1.25),
                     correction = "nope"))
  expect_error(ctost(theta = s$theta, sigma = s$sigma, nu = s$nu, delta = log(1.25),
                     alpha = 0.6))
})
