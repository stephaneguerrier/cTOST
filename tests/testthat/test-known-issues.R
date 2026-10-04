# Known defects of cdd680d, documented by the audit of 2026-10-03.
# Each test asserts the CURRENT (broken) behaviour so that a fix is visible:
# when a fix lands (PR 3 / PR 4a), flip the expectation in the same PR
# (expect_error -> expect_no_error, expect_null -> expect_s3_class, ...).

test_that("F01: print() of a single-quantile qtost object errors", {
  q = qtost(list(mean = 100, sd = 15, n = 50), list(mean = 102, sd = 16, n = 50),
            pi_x = 0.95, delta = 0.03, method = "unadjusted")
  expect_error(suppressMessages(capture.output(print(q))))
  i = qtost_fda_inputs()
  a = qtost(i$x, i$y, pi_x = 0.8, delta = 0.15, method = "alpha")
  expect_error(suppressMessages(capture.output(print(a))))
})

test_that("F04: no print method is registered for two-quantile ('mqtost') objects", {
  expect_null(getS3method("print", "mqtost", optional = TRUE))
  # the method exists under a different class name
  expect_type(getS3method("print", "m_qtost", optional = TRUE), "closure")
})

test_that("F05: compare_to_qtost() cannot take a qtost object", {
  i = qtost_fda_inputs()
  a = qtost(i$x, i$y, pi_x = 0.8, delta = 0.15, method = "alpha")
  expect_error(capture.output(compare_to_qtost(a)))
})

test_that("F03: two-quantile alpha-qTOST with MC_sup = FALSE fails (shadowed power_cTOST_mv)", {
  i = qtost_fda_inputs()
  expect_error(qtost(i$x, i$y, pi_x = c(0.25, 0.75), delta = 0.15, method = "alpha",
                     MC_sup = FALSE), "unused argument")
})

test_that("F06: plot() is not registered for multivariate objects", {
  t = ticlopidine_stats()
  m = ctost(theta = t$theta, sigma = t$sigma, nu = t$nu, delta = log(1.25),
            method = "unadjusted")
  pdf(NULL)
  on.exit(dev.off())
  expect_error(plot(m))
})

test_that("F07: compare_to_tost() errors on multivariate objects", {
  t = ticlopidine_stats()
  m = ctost(theta = t$theta, sigma = t$sigma, nu = t$nu, delta = log(1.25),
            method = "alpha")
  expect_error(suppressMessages(suppressWarnings(capture.output(compare_to_tost(m)))), "length")
})

test_that("F08: print.mtost() errors when sigma has no column names", {
  m = ctost(theta = c(0.01, 0.02), sigma = diag(2) * 0.004, nu = 11, delta = log(1.25),
            method = "unadjusted")
  expect_error(suppressMessages(suppressWarnings(capture.output(print(m)))))
})

test_that("F20: seed = NULL is rejected by the C++ path of tost_equiv_dp()", {
  expect_error(tost_equiv_dp(mean_private_obs = 2.5, sd_private_obs = 1, a = 0, b = 5,
                             n = 100, lower = 1.5, upper = 3.5, epsilon = 1, B = 100,
                             seed = NULL))
})

# F02 (tost_equiv_dp(B = 0) crashes the R session) cannot be tested until the
# input guard exists; add expect_error(..., "B") together with the fix.
