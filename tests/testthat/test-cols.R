## Corrected ordinary least squares, sfm(estimator = "cols").

test_that("COLS recovers the DGP for NHN and NE", {
  skip_on_cran()
  d <- data_gen_cs(N = 4000, rand = 1, sig_u = 1, sig_v = 0.3, cons = 0.5,
                   beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.1)

  f <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = d, estimator = "cols")
  p <- coef(f)
  expect_equal(unname(p[["sigv"]]), 0.3, tolerance = 0.15)
  expect_equal(unname(p[["sigu"]]), 1.0, tolerance = 0.15)
  expect_equal(unname(p[["(Intercept)"]]), 0.5, tolerance = 0.15)
  expect_equal(unname(p[["x1"]]), 0.5, tolerance = 0.15)

  g <- sfm(y_pcs_e ~ x1 + x2, model_name = "NE", data = d, estimator = "cols")
  q <- coef(g)
  expect_equal(unname(q[["sigv"]]), 0.3, tolerance = 0.2)
  expect_equal(unname(q[["sigu"]]), 1.0, tolerance = 0.2)
})

test_that("the intercept is actually corrected by E[u]", {
  ## The whole content of the "corrected" in COLS. An earlier version located
  ## the intercept column by a name that data_i_vars does not carry (it holds
  ## make.names()-mangled labels), so match() returned NA and the correction
  ## was silently skipped -- leaving a plain OLS intercept that looked
  ## plausible but was E[u] too low.
  skip_on_cran()
  d <- data_gen_cs(N = 3000, rand = 4, sig_u = 1, sig_v = 0.3, cons = 0.5,
                   beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.1)
  f  <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = d, estimator = "cols")
  b0 <- coef(f)[["(Intercept)"]]
  su <- coef(f)[["sigu"]]
  ols <- stats::lm(y_pcs ~ x1 + x2, data = d)$coefficients[["(Intercept)"]]

  expect_equal(b0 - ols, su*sqrt(2/pi), tolerance = 1e-8)
  expect_gt(b0, ols)                       ## the frontier lies above the mean
  ## and the slopes are left exactly as OLS gave them
  expect_equal(unname(coef(f)[["x1"]]),
               unname(stats::lm(y_pcs ~ x1 + x2, data = d)$coefficients[["x1"]]),
               tolerance = 1e-10)
})

test_that("NE and NHN corrections use their own E[u], not a shared one", {
  ## E[u] = sigma_u*sqrt(2/pi) for the half-normal but sigma_u for the
  ## exponential; using one formula for both would pass the NHN test above and
  ## quietly mis-shift every NE intercept.
  skip_on_cran()
  d <- data_gen_cs(N = 3000, rand = 5, sig_u = 1, sig_v = 0.3, cons = 0.5,
                   beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.1)
  g   <- sfm(y_pcs_e ~ x1 + x2, model_name = "NE", data = d, estimator = "cols")
  ols <- stats::lm(y_pcs_e ~ x1 + x2, data = d)$coefficients[["(Intercept)"]]
  expect_equal(coef(g)[["(Intercept)"]] - ols, coef(g)[["sigu"]], tolerance = 1e-8)
})

test_that("wrong-skew residuals are reported, not silently absorbed", {
  ## Olson, Schmidt and Waldman's Type I failure: with positively skewed
  ## residuals the moment equations have no admissible solution.
  set.seed(3)
  n <- 400
  x1 <- stats::runif(n); x2 <- stats::runif(n)
  ## deliberately POSITIVE skew: an "inefficiency" term with the wrong sign
  y  <- 0.5 + 0.5*x1 + 0.5*x2 + stats::rnorm(n, 0, 0.3) + abs(stats::rnorm(n))
  dd <- data.frame(y, x1, x2)

  expect_warning(f <- sfm(y ~ x1 + x2, model_name = "NHN", data = dd,
                          estimator = "cols"), "WRONG")
  expect_equal(unname(coef(f)[["sigu"]]), 0)
  expect_true(f$wrong_skew)
  expect_gte(f$residual_moments[["m3"]], 0)
  ## no efficiency prediction is possible when sigma_u collapses
  expect_true(all(is.na(f$exp_u_hat)))
  ## and the intercept is then left uncorrected, since E[u] = 0
  ols <- stats::lm(y ~ x1 + x2, data = dd)$coefficients[["(Intercept)"]]
  expect_equal(coef(f)[["(Intercept)"]], ols, tolerance = 1e-10)
})

test_that("COLS is deterministic and needs no optimizer", {
  skip_on_cran()
  d <- data_gen_cs(N = 1000, rand = 6, sig_u = 1, sig_v = 0.3, cons = 0.5,
                   beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.1)
  a <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = d, estimator = "cols")
  b <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = d, estimator = "cols")
  expect_equal(coef(a), coef(b))
  expect_null(a$opt)                       ## no optimizer output to store
  expect_equal(a$estimator, "cols")
})

test_that("bootstrap standard errors are produced and reproducible", {
  skip_on_cran()
  d <- data_gen_cs(N = 1000, rand = 7, sig_u = 1, sig_v = 0.3, cons = 0.5,
                   beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.1)
  a <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = d, estimator = "cols",
           cols_boot = 100, rand.cols = 11)
  b <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = d, estimator = "cols",
           cols_boot = 100, rand.cols = 11)
  expect_equal(a$std.errors, b$std.errors)
  expect_true(all(is.finite(a$std.errors)))
  expect_equal(dim(a$cols_boot_draws), c(100L, length(coef(a))))
  ## without the bootstrap, the moment-based parameters carry no standard error
  ## rather than a misleading one
  c0 <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = d, estimator = "cols")
  expect_true(is.na(c0$std.errors[["sigv"]]))
  expect_true(is.na(c0$std.errors[["sigu"]]))
  expect_true(is.na(c0$std.errors[["(Intercept)"]]))
  expect_true(is.finite(c0$std.errors[["x1"]]))   ## OLS slope SEs are valid
})

test_that("the bootstrap restores the caller's RNG stream", {
  skip_on_cran()
  d <- data_gen_cs(N = 400, rand = 8, sig_u = 1, sig_v = 0.3, cons = 0.5,
                   beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.1)
  set.seed(99); ref <- stats::rnorm(3)
  set.seed(99)
  invisible(sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = d,
                estimator = "cols", cols_boot = 20, rand.cols = 5))
  expect_equal(stats::rnorm(3), ref)
})

test_that("unsupported models and incompatible options error clearly", {
  skip_on_cran()
  d <- data_gen_cs(N = 300, rand = 9, sig_u = 1, sig_v = 0.3, cons = 0.5,
                   beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.1)
  expect_error(sfm(y_pcs ~ x1 + x2, model_name = "NR", data = d,
                   estimator = "cols"), "NHN")
  expect_error(sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = d,
                   estimator = "cols", robust = "mlqe"), "moment estimator")
})

test_that("estimator defaults to mle, leaving existing calls untouched", {
  expect_equal(eval(formals(sfm)$estimator)[1], "mle")
})
