# Behaviour fixes of PR 4a (no numerical change): input validation, crash guards,
# RNG restoration, object consistency. Audit IDs in the test names.

test_that("F02: tost_equiv_dp() rejects B < 1 and n < 2 instead of crashing the session", {
  base = list(mean_private_obs = 2.5, sd_private_obs = 1, a = 0, b = 5, n = 100,
              lower = 1.5, upper = 3.5, epsilon = 1)
  expect_error(do.call(tost_equiv_dp, c(base, B = 0)), "B")
  expect_error(do.call(tost_equiv_dp, c(base, B = 0.5)), "B")
  expect_error(do.call(tost_equiv_dp, c(base, B = -1)), "B")
  expect_error(tost_equiv_dp(mean_private_obs = 2.5, sd_private_obs = 1, a = 0, b = 5, n = 1,
                             lower = 1.5, upper = 3.5, epsilon = 1, B = 100), "n")
  # the C++ entry point is guarded too
  expect_error(cTOST:::tost_dp_one_sample_ultra_fast(0, 5, 100L, 1, 2.5, 1, 1.5, 3.5, 0L, 0.05, 1L), "B")
})

test_that("F18/F20: tost_equiv_dp() validates its inputs and resolves seed = NULL", {
  base = list(mean_private_obs = 2.5, sd_private_obs = 1, a = 0, b = 5, n = 100,
              lower = 1.5, upper = 3.5, epsilon = 1, B = 100)
  bad = list(list(epsilon = 0), list(epsilon = NA), list(alpha = 0.7), list(alpha = 0),
             list(n = 10.7), list(a = 5, b = 0), list(lower = 3.5, upper = 1.5),
             list(mean_private_obs = NA), list(sd_private_obs = "a"), list(seed = "a"),
             list(seed = 2^31), list(seed = 1.5), list(method = "nope"))
  for (b in bad) {
    args = base; args[names(b)] = b
    expect_error(suppressWarnings(do.call(tost_equiv_dp, args)), info = paste(names(b), collapse = ","))
  }
  expect_warning(do.call(tost_equiv_dp, c(base[names(base) != "mean_private_obs"], mean_private_obs = 7)), "outside")
  # seed = NULL works and gives different results on repeated calls
  r1 = do.call(tost_equiv_dp, c(base, seed = list(NULL)))
  r2 = do.call(tost_equiv_dp, c(base, seed = list(NULL)))
  expect_s3_class(r1, "tost_dp")
  expect_false(identical(r1$conf.int, r2$conf.int))
})

test_that("F19: prop_test_equiv_dp() validates its inputs and reports an unsolvable configuration", {
  base = list(p_hat = 0.42, n = 500, lower = 0.3, upper = 0.5, epsilon = 1, B = 100)
  bad = list(list(epsilon = -1), list(alpha = 0.5), list(n = 0), list(lower = 0.5, upper = 0.3),
             list(B = 0), list(max_resample = -1), list(p_hat = NA), list(seed = "a"))
  for (b in bad) {
    args = base; args[names(b)] = b
    expect_error(suppressWarnings(do.call(prop_test_equiv_dp, args)), info = paste(names(b), collapse = ","))
  }
  # a privatized proportion far outside [0, 1] has no moment-matching solution
  expect_error(prop_test_equiv_dp(p_hat = 5, n = 1000, lower = 0.3, upper = 0.5, epsilon = 10, B = 50),
               "moment-matching")
  r1 = do.call(prop_test_equiv_dp, c(base, seed = list(NULL)))
  r2 = do.call(prop_test_equiv_dp, c(base, seed = list(NULL)))
  expect_false(identical(r1$conf.int, r2$conf.int))
})

test_that("F16/F80: ctost() validates its inputs and accepts 1 x 1 matrices without warnings", {
  s = skin_stats()
  ok = function(...) ctost(theta = s$theta, sigma = s$sigma, nu = s$nu, delta = log(1.25), ...)
  expect_error(ok(alpha = 0.5)); expect_error(ok(alpha = NA)); expect_error(ok(method = c("alpha", "delta")))
  expect_error(ok(B = 0)); expect_error(ok(seed = "a"))
  expect_error(ctost(theta = s$theta, sigma = 0, nu = s$nu, delta = log(1.25)), "sigma")
  expect_error(ctost(theta = s$theta, sigma = s$sigma, nu = 0, delta = log(1.25)), "nu")
  expect_error(ctost(theta = s$theta, sigma = s$sigma, nu = s$nu, delta = 0), "delta")
  expect_error(ctost(theta = NA, sigma = s$sigma, nu = s$nu, delta = log(1.25)), "theta")
  t = ticlopidine_stats()
  expect_error(ctost(theta = t$theta[1:3], sigma = t$sigma, nu = t$nu, delta = log(1.25), method = "unadjusted"), "p x p")
  expect_error(ctost(theta = t$theta, sigma = t$sigma, nu = 2, delta = log(1.25), method = "alpha"), "nu")
  # 1 x 1 matrices (e.g. cov(X)/n of a one-column matrix): same result, no recycling warning
  ref = ok(method = "alpha")
  expect_no_warning(m <- ctost(theta = matrix(s$theta), sigma = matrix(s$sigma), nu = s$nu, delta = log(1.25), method = "alpha"))
  expect_equal(unname(m$ci), unname(ref$ci), tolerance = tol_exact)
  expect_equal(unname(m$corrected_alpha), unname(ref$corrected_alpha), tolerance = tol_exact)
  expect_no_warning(ctost(theta = matrix(s$theta), sigma = matrix(s$sigma), nu = s$nu, delta = log(1.25)))
})

test_that("F17: qtost() validates pi_x, the margins and the Monte Carlo settings", {
  i = qtost_fda_inputs()
  expect_error(qtost(i$x, i$y, pi_x = 0, delta = 0.1), "pi_x")
  expect_error(qtost(i$x, i$y, pi_x = 1.5, delta = 0.1), "pi_x")
  expect_error(qtost(i$x, i$y, pi_x = NA, delta = 0.1), "pi_x")
  expect_error(qtost(i$x, i$y, pi_x = 0.95, delta = 0.1), "pi_x")   # pi_x + delta >= 1
  expect_error(qtost(i$x, i$y, pi_x = 0.8, delta = 0.15, B = 0), "B")
  expect_error(qtost(i$x, i$y, pi_x = 0.8, delta = 0.15, max_iter = 0), "max_iter")
  expect_error(qtost(i$x, i$y, pi_x = 0.8, delta = 0.15, MC_sup = NA), "MC_sup")
  expect_error(qtost(i$x, i$y, pi_x = 0.8, delta = 0.15, seed = "a"), "seed")
  # tiny samples still have a corrected level once the margins are valid
  expect_s3_class(qtost(list(mean = 0, sd = 1, n = 3), list(mean = 0, sd = 1, n = 3), pi_x = 0.8, delta = 0.1, B = 2000), "qtost")
})

test_that("F07/F08: compare_to_tost() refuses multivariate objects; print.mtost() handles unnamed inputs", {
  t = ticlopidine_stats()
  m = ctost(theta = t$theta, sigma = t$sigma, nu = t$nu, delta = log(1.25), method = "alpha")
  expect_error(compare_to_tost(m), "univariate")
  u = ctost(theta = c(0.01, 0.02), sigma = diag(2) * 0.004, nu = 11, delta = log(1.25), method = "unadjusted")
  out = capture.output(suppressMessages(print(u)))
  expect_true(any(grepl("theta1", out, fixed = TRUE)))
  expect_true(any(grepl("theta2", out, fixed = TRUE)))
  u2 = ctost(theta = c(a = 0.01, b = 0.02), sigma = diag(2) * 0.004, nu = 11, delta = log(1.25), method = "unadjusted")
  out = capture.output(suppressMessages(print(u2)))
  expect_true(any(grepl("^a ", out)))
})

test_that("F15: the caller's random number stream is restored by every exported function", {
  s = skin_stats(); t = ticlopidine_stats(); i = qtost_fda_inputs()
  calls = list(
    boot = function() ctost(theta = s$theta, sigma = s$sigma, nu = s$nu, delta = log(1.25), correction = "bootstrap", B = 200),
    mv_alpha = function() ctost(theta = t$theta, sigma = t$sigma, nu = t$nu, delta = log(1.25), method = "alpha", B = 500),
    mv_opt = function() ctost(theta = t$theta, sigma = t$sigma, nu = t$nu, delta = log(1.25)),
    qtost = function() qtost(i$x, i$y, pi_x = 0.8, delta = 0.15, method = "alpha", B = 2000),
    dp = function() tost_equiv_dp(mean_private_obs = 2.5, sd_private_obs = 1, a = 0, b = 5, n = 100, lower = 1.5, upper = 3.5, epsilon = 1, B = 500),
    prop = function() prop_test_equiv_dp(p_hat = 0.42, n = 500, lower = 0.3, upper = 0.5, epsilon = 1, B = 200))
  for (nm in names(calls)) {
    set.seed(1); u_ref = runif(3)
    set.seed(1); res = calls[[nm]](); u_after = runif(3)
    expect_identical(u_after, u_ref, info = nm)
    # and the function's own result is unchanged by the caller's state
    set.seed(999); res2 = calls[[nm]]()
    res$call = NULL; res2$call = NULL
    expect_equal(res, res2, info = nm)
  }
  # with no seed set at all, nothing is left behind either
  if (exists(".Random.seed", envir = globalenv(), inherits = FALSE)) rm(".Random.seed", envir = globalenv())
  calls$qtost()
  expect_false(exists(".Random.seed", envir = globalenv(), inherits = FALSE))
})

test_that("F10/F11: the univariate cTOST warns and falls back instead of crashing on table gaps", {
  # (nu, sigma) cell with no offline table entry -> bootstrap with a warning, same as explicit bootstrap
  expect_warning(r <- ctost(theta = 0, sigma = 0.1, nu = 100, delta = log(1.25)), "bootstrap")
  b = ctost(theta = 0, sigma = 0.1, nu = 100, delta = log(1.25), correction = "bootstrap")
  expect_equal(r$correction, "bootstrap")
  expect_equal(r$corrected_alpha, b$corrected_alpha)
  expect_equal(unname(r$ci), unname(b$ci), tolerance = tol_exact)
  # bootstrap correction non-positive (tiny nu): floored with a warning instead of NaN crash
  expect_warning(r2 <- ctost(theta = 0, sigma = 0.0001, nu = 2, delta = log(1.25), correction = "bootstrap", B = 200, seed = 1), "floored")
  expect_true(is.finite(r2$corrected_c))
  expect_equal(r2$corrected_alpha, 1e-6)
})

test_that("F23/F24/F25: dead optimizer branches fail with a clear message", {
  expect_error(cTOST:::get_alpha_TOST_MC_mv_core(alpha = 0.05, Sigma = diag(2), nu = 10, delta = log(1.25), argsup_meth = "MC"), "not implemented")
  expect_error(cTOST:::get_c_of_0(delta = log(1.25), sigma = 0.1, alpha = 0.05, optim = "uniroot"), "not implemented")
  expect_error(cTOST:::obj_fun_one_sample(c(2, 1), 0, 5, 10, 0.1, 0.1, rnorm(10), 0.5, 0.5, 2, 1, use_cpp = TRUE), "not available")
})

test_that("F78/F79/F81: result objects are consistent", {
  i = qtost_fda_inputs()
  q = qtost(i$x, i$y, pi_x = 0.8, delta = 0.15, method = "unadjusted")
  expect_false("corrected_alpha" %in% names(q))
  a = qtost(i$x, i$y, pi_x = 0.8, delta = 0.15, method = "alpha")
  expect_true("corrected_alpha" %in% names(a))
  s = skin_stats()
  off = ctost(theta = s$theta, sigma = s$sigma, nu = s$nu, delta = log(1.25))
  expect_identical(unname(off[["theta"]]), unname(s$theta))
  p = prop_test_equiv_dp(p_hat = 0.42, n = 500, lower = 0.3, upper = 0.5, epsilon = 1, B = 200)
  expect_null(names(p$conf.int))
})
