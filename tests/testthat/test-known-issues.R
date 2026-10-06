# Known defects of cdd680d, documented by the audit of 2026-10-03.
# F01, F03, F04, F05 and F06 were fixed in PR 3 and are now covered by positive tests in
# test-qtost.R and test-ctost-multivariate.R.
# Each test asserts the CURRENT (broken) behaviour so that a fix is visible:
# when a fix lands (PR 3 / PR 4a), flip the expectation in the same PR
# (expect_error -> expect_no_error, expect_null -> expect_s3_class, ...).

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
