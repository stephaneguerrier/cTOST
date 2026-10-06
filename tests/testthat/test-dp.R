# Characterization tests: differential-privacy equivalence tests.
# See helper-data.R for the provenance of the reference values (cdd680d).
# The DP methods stay exactly as published (owner decision, 2026-10-04); these tests
# lock their output on the documented examples and on one informative setting.

dp_one_sample_inputs = function() {
  # ?tost_equiv_dp, one-sample example: privatized mean and sd
  set.seed(123)
  n = 100; a = 0; b = 5
  z = rnorm(n)
  x = pmin(pmax(2.5 + 1 * z, a), b)
  scale_mean = ((b - a) / n) / (1 / 2)
  scale_sd = ((b - a) / sqrt(n - 1)) / (1 / 2)
  u = runif(2)
  list(mean = mean(x) - scale_mean * sign(u[1] - 0.5) * log(1 - 2 * abs(u[1] - 0.5)),
       sd = sd(x) - scale_sd * sign(u[2] - 0.5) * log(1 - 2 * abs(u[2] - 0.5)),
       n = n, a = a, b = b)
}

test_that("tost_equiv_dp(): the documented one-sample example is locked (degenerate interval, audit F71)", {
  i = dp_one_sample_inputs()
  expect_equal(i$mean, 2.51647675511, tolerance = tol_exact)
  expect_equal(i$sd, 3.51235877981, tolerance = tol_exact)
  res = tost_equiv_dp(mean_private_obs = i$mean, sd_private_obs = i$sd, a = i$a, b = i$b,
                      n = i$n, lower = 1.5, upper = 3.5, epsilon = 1, B = 1000)
  expect_s3_class(res, "tost_dp")
  expect_equal(res$test_type, "one-sample")
  expect_equal(res$method, "cpp")
  expect_false(res$decision)
  # the privatized sd (3.51) exceeds (b - a) / 2, so the interval collapses to the bounds
  expect_equal(unname(res$conf.int), c(0, 5))
  res_r = tost_equiv_dp(mean_private_obs = i$mean, sd_private_obs = i$sd, a = i$a, b = i$b,
                        n = i$n, lower = 1.5, upper = 3.5, epsilon = 1, B = 1000, method = "r")
  expect_equal(unname(res_r$conf.int), c(0, 5))
  expect_false(res_r$decision)
})

test_that("tost_equiv_dp(): an informative one-sample setting is locked and reproducible", {
  res = tost_equiv_dp(mean_private_obs = 2.5, sd_private_obs = 1, a = 0, b = 5, n = 500,
                      epsilon = 2, lower = 1.5, upper = 3.5, seed = 42)   # default B = 10^4
  expect_true(res$decision)
  expect_equal(unname(res$conf.int), c(2.41170063019, 2.58292349577), tolerance = tol_exact)
  expect_equal(res$estimate, 2.5)
  expect_equal(res$sd_estimate, 1)
  expect_equal(res$sample_size, 500)
  expect_equal(res$bounds, c(0, 5))
  expect_equal(res$B, 10^4)
  expect_equal(res$seed, 42)
  expect_lt(res$conf.int[1], res$conf.int[2])
  expect_equal(res$decision, res$conf.int[1] >= res$lower && res$conf.int[2] <= res$upper)
  # same seed -> identical result; another seed -> different interval
  again = tost_equiv_dp(mean_private_obs = 2.5, sd_private_obs = 1, a = 0, b = 5, n = 500,
                        epsilon = 2, lower = 1.5, upper = 3.5, seed = 42)
  expect_identical(again$conf.int, res$conf.int)
  other = tost_equiv_dp(mean_private_obs = 2.5, sd_private_obs = 1, a = 0, b = 5, n = 500,
                        epsilon = 2, lower = 1.5, upper = 3.5, seed = 43)
  expect_false(identical(other$conf.int, res$conf.int))
  expect_equal(unname(other$conf.int), c(2.41667607427, 2.58597375154), tolerance = tol_exact)
  # pure-R path (slower; B = 1000). It uses the same draws differently (audit F32), so it
  # is locked separately rather than compared with the C++ path.
  res_r = tost_equiv_dp(mean_private_obs = 2.5, sd_private_obs = 1, a = 0, b = 5, n = 500,
                        epsilon = 2, lower = 1.5, upper = 3.5, seed = 42, B = 1000, method = "r")
  expect_true(res_r$decision)
  expect_equal(unname(res_r$conf.int), c(2.41271767013, 2.59073412689), tolerance = tol_exact)
})

test_that("tost_equiv_dp(): the documented two-sample example is locked", {
  set.seed(456)
  n1 = 100; n2 = 100
  z1 = rnorm(n1); z2 = rnorm(n2)
  x1 = pmin(pmax(2.5 + 1.0 * z1, 0), 5)
  x2 = pmin(pmax(2.3 + 1.2 * z2, 0), 5)
  scale_mean = (5 / 100) / (1 / 2)
  scale_sd = (5 / sqrt(99)) / (1 / 2)
  u1 = runif(2); u2 = runif(2)
  m1 = mean(x1) - scale_mean * sign(u1[1] - 0.5) * log(1 - 2 * abs(u1[1] - 0.5))
  s1 = sd(x1) - scale_sd * sign(u1[2] - 0.5) * log(1 - 2 * abs(u1[2] - 0.5))
  m2 = mean(x2) - scale_mean * sign(u2[1] - 0.5) * log(1 - 2 * abs(u2[1] - 0.5))
  s2 = sd(x2) - scale_sd * sign(u2[2] - 0.5) * log(1 - 2 * abs(u2[2] - 0.5))
  expect_equal(c(m1, s1, m2, s2), c(2.64686142331, 0.699937479251, 2.29428695374, 1.09270613574),
               tolerance = tol_exact)
  res = tost_equiv_dp(mean_private_obs = m1, sd_private_obs = s1, a = 0, b = 5, n = n1,
                      mean_private_obs2 = m2, sd_private_obs2 = s2, a2 = 0, b2 = 5, n2 = n2,
                      lower = -1, upper = 1, epsilon = 1, B = 1000)
  expect_s3_class(res, "tost_dp")
  expect_equal(res$test_type, "two-sample")
  expect_false(res$decision)
  expect_equal(unname(res$conf.int), c(-0.807733487334, 2.6562973178), tolerance = tol_exact)
})

test_that("tost_equiv_dp() rejects an unknown method and an incomplete second sample", {
  expect_error(tost_equiv_dp(mean_private_obs = 2.5, sd_private_obs = 1, a = 0, b = 5, n = 100,
                             lower = 1.5, upper = 3.5, epsilon = 1, B = 100, method = "nope"))
  expect_error(tost_equiv_dp(mean_private_obs = 2.5, sd_private_obs = 1, a = 0, b = 5, n = 100,
                             mean_private_obs2 = 2.4, lower = -1, upper = 1, epsilon = 1, B = 100))
})

test_that("prop_test_equiv_dp(): the documented one-sample example is locked and reproducible", {
  set.seed(123)
  n = 500
  x = rbinom(n, 1, 0.42)
  u = runif(1)
  p_hat = mean(x) + (-(1 / (n * 1)) * sign(u - 0.5) * log(1 - 2 * abs(u - 0.5)))
  expect_equal(p_hat, 0.409307150831, tolerance = tol_exact)
  res = prop_test_equiv_dp(p_hat = p_hat, n = n, lower = 0.3, upper = 0.5, epsilon = 1, B = 1000)
  expect_s3_class(res, "prop_test_dp")
  expect_equal(res$test_type, "one-sample")
  expect_true(res$decision)
  # conf.int carries quantile names in cdd680d (audit F81); compare values only
  expect_equal(unname(res$conf.int), c(0.376021334984, 0.444683331219), tolerance = tol_exact)
  expect_equal(res$estimate, p_hat)
  expect_equal(res$sample_size, n)
  expect_equal(res$decision, res$conf.int[[1]] >= res$lower && res$conf.int[[2]] <= res$upper)
  again = prop_test_equiv_dp(p_hat = p_hat, n = n, lower = 0.3, upper = 0.5, epsilon = 1, B = 1000)
  expect_identical(again$conf.int, res$conf.int)
})

test_that("prop_test_equiv_dp(): the documented two-sample example is locked", {
  set.seed(456)
  x1 = rbinom(300, 1, 0.52); x2 = rbinom(300, 1, 0.48)
  u1 = runif(1); u2 = runif(1)
  p1 = mean(x1) - (1 / 300) * sign(u1 - 0.5) * log(1 - 2 * abs(u1 - 0.5))
  p2 = mean(x2) - (1 / 300) * sign(u2 - 0.5) * log(1 - 2 * abs(u2 - 0.5))
  expect_equal(c(p1, p2), c(0.492157635095, 0.499298117345), tolerance = tol_exact)
  res = prop_test_equiv_dp(p_hat = p1, n = 300, p_hat2 = p2, n2 = 300, lower = -0.1, upper = 0.1,
                           epsilon = c(1, 1), B = 1000)
  expect_equal(res$test_type, "two-sample")
  expect_true(res$decision)
  expect_equal(unname(res$conf.int), c(-0.0740910662797, 0.0605186688686), tolerance = tol_exact)
})
