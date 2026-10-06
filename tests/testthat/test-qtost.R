# Characterization tests: quantile equivalence, qtost().
# See helper-data.R for the provenance of the reference values.
# Values are locked from cdd680d. Audit F13 (boundary selection of the single-quantile
# alpha-qTOST) would change the corrected alphas below when applied; update them in
# the same PR and record the before/after values in NEWS.md.

test_that("qTOST on the FDA C_trough example (summary statistics) is locked", {
  i = qtost_fda_inputs()
  res = qtost(i$x, i$y, pi_x = 0.8, delta = 0.15, method = "unadjusted")
  expect_s3_class(res, "qtost")
  expect_equal(res$method, "qTOST")
  expect_false(unname(res$decision))
  expect_equal(unname(res$theta), 0.492648245202, tolerance = tol_exact)
  expect_equal(unname(res$sigma), 0.297912152114, tolerance = tol_exact)
  expect_equal(as.numeric(res$ci), c(0.501047765356, 0.837115091621), tolerance = tol_exact)  # ci is a 1 x 2 matrix
  expect_equal(as.numeric(res$eq_region), c(0.65, 0.8, 0.95))
  expect_null(res$corrected_alpha)
})

test_that("alpha-qTOST on the FDA C_trough example is locked", {
  i = qtost_fda_inputs()
  res = qtost(i$x, i$y, pi_x = 0.8, delta = 0.15, method = "alpha")
  expect_s3_class(res, "qtost")
  expect_equal(res$method, "alpha-qTOST")
  expect_false(unname(res$decision))
  expect_equal(unname(res$corrected_alpha), 0.0616720401485, tolerance = tol_exact)
  expect_equal(as.numeric(res$ci), c(0.513401550455, 0.829374782194), tolerance = tol_exact)
  # the corrected interval is narrower than the unadjusted one
  una = qtost(i$x, i$y, pi_x = 0.8, delta = 0.15, method = "unadjusted")
  expect_gt(res$ci[1], una$ci[1])
  expect_lt(res$ci[2], una$ci[2])
})

test_that("qtost() on raw data (?qtost Example 2) is locked", {
  set.seed(12345)
  x_data = rnorm(n = 30, mean = 0, sd = 0.1)
  y_data = rnorm(n = 10, mean = 0, sd = 0.1)
  una = qtost(x_data, y_data, pi_x = 0.8, delta = 0.15, method = "unadjusted")
  expect_equal(unname(una$theta), 0.118101999078, tolerance = tol_exact)
  expect_equal(unname(una$sigma), 0.356904728063, tolerance = tol_exact)
  expect_equal(as.numeric(una$ci), c(0.319551244918, 0.75964405336), tolerance = tol_exact)
  expect_false(unname(una$decision))
  adj = qtost(x_data, y_data, pi_x = 0.8, delta = 0.15, method = "alpha")
  expect_equal(unname(adj$corrected_alpha), 0.108391096569, tolerance = tol_exact)
  expect_equal(as.numeric(adj$ci), c(0.373453018423, 0.711893775721), tolerance = tol_exact)
  expect_false(unname(adj$decision))
})

test_that("qTOST at the 95th percentile (?qtost Example 3) is locked", {
  res = qtost(list(mean = 100, sd = 15, n = 50), list(mean = 102, sd = 16, n = 50),
              pi_x = 0.95, delta = 0.03, method = "unadjusted")
  expect_equal(unname(res$theta), 1.41705027527, tolerance = tol_exact)
  expect_equal(unname(res$sigma), 0.285372791872, tolerance = tol_exact)
  expect_equal(as.numeric(res$ci), c(0.828347136993, 0.970382610918), tolerance = tol_exact)
  expect_equal(as.numeric(res$eq_region), c(0.92, 0.95, 0.98))
  expect_false(unname(res$decision))
})

test_that("two-quantile qTOST and alpha-qTOST on the FDA example are locked", {
  i = qtost_fda_inputs()
  una = qtost(i$x, i$y, pi_x = c(0.25, 0.75), delta = 0.15, method = "unadjusted")
  expect_s3_class(una, "mqtost")
  expect_equal(unname(una$decision), c(FALSE, FALSE))
  expect_equal(unname(una$theta), c(-0.755268110718, 0.355081725384), tolerance = tol_exact)
  expect_equal(unname(una$sigma), c(0.315668239072, 0.289442404404), tolerance = tol_exact)
  expect_equal(unname(una$ci),
               matrix(c(0.10124381658, 0.451842086654, 0.406700791376, 0.797061797148), ncol = 2),
               tolerance = tol_exact)
  expect_equal(unname(una$eq_region), matrix(c(0.1, 0.6, 0.25, 0.75, 0.4, 0.9), ncol = 3))
  # default MC_sup = TRUE path (the MC_sup = FALSE path is broken, see test-known-issues.R)
  adj = qtost(i$x, i$y, pi_x = c(0.25, 0.75), delta = 0.15, method = "alpha")
  expect_s3_class(adj, "mqtost")
  expect_equal(adj$method, "alpha-qTOST")
  expect_equal(unname(adj$decision), c(TRUE, FALSE))
  expect_equal(unname(adj$corrected_alpha), 0.121657816212, tolerance = tol_exact)
  expect_equal(unname(adj$ci),
               matrix(c(0.130597600639, 0.506932419165, 0.34939085839, 0.755777937985), ncol = 2),
               tolerance = tol_exact)
})

test_that("qtost() rejects a non-positive delta and an unknown method", {
  i = qtost_fda_inputs()
  expect_error(qtost(i$x, i$y, pi_x = 0.8, delta = 0, method = "unadjusted"))
  expect_error(qtost(i$x, i$y, pi_x = 0.8, delta = 0.15, method = "nope"))
})
