## vhet / uhet / muhet: the heteroskedastic cross-sectional engine.

.het_data <- function(n = 700, seed = 42) {
  set.seed(seed)
  x1 <- runif(n, 1, 5)
  zv <- runif(n, -1, 1)
  zu <- runif(n, -1, 1)
  s_v <- exp(-1.2 + 0.7 * zv)
  s_u <- exp(-0.6 + 0.9 * zu)
  y <- 1 + 0.5 * x1 + rnorm(n, 0, s_v) - abs(rnorm(n, 0, s_u))
  data.frame(y = y, x1 = x1, zv = zv, zu = zu)
}

test_that("vhet and the pipe segment recover both variance blocks", {
  skip_on_cran()
  d <- .het_data(n = 1200)
  f <- sfm(y ~ x1 | zu, model_name = "NHN_Z", data = d, vhet = ~zv)

  expect_s3_class(f, "sfareg")
  expect_identical(
    names(coef(f)),
    c("(Intercept)", "x1", "Zv.(Intercept)", "Zv.zv", "Zu.(Intercept)", "Zu.zu")
  )
  ## Every coefficient within 3 SE of truth.
  truth <- c(1, 0.5, -1.2, 0.7, -0.6, 0.9)
  expect_true(all(abs(coef(f) - truth) < 3 * f$std.errors))
  ## Per-observation scales, not scalars.
  expect_length(f$sigma_u, nrow(d))
  expect_true(all(f$sigma_v > 0))
  expect_true(all(f$exp_u_hat >= 0 & f$exp_u_hat <= 1))
})

test_that("the two z-links are the same model, reparameterized", {
  skip_on_cran()
  d <- .het_data(n = 500)
  a <- sfm(y ~ x1 | zu, "NHN_Z", d, vhet = ~zv, z_link = "sd")
  b <- sfm(y ~ x1 | zu, "NHN_Z", d, vhet = ~zv, z_link = "var")

  expect_equal(as.numeric(logLik(a)), as.numeric(logLik(b)), tolerance = 1e-6)
  ## Frontier identical; every variance coefficient exactly doubled.
  expect_equal(coef(a)[1:2], coef(b)[1:2], tolerance = 1e-5)
  expect_equal(2 * coef(a)[3:6], coef(b)[3:6], tolerance = 1e-4, ignore_attr = TRUE)
})

test_that("muhet fits Battese-Coelli (1995) and needs NTN", {
  skip_on_cran()
  set.seed(11)
  n <- 1500
  x1 <- runif(n, 1, 5)
  zm <- runif(n, -1, 1)
  mu_i <- 0.5 + 1.0 * zm
  ## N+(mu_i, 0.6^2) by inverse-CDF truncation.
  u <- mu_i + 0.6 * qnorm(runif(n, pnorm(-mu_i / 0.6), 1))
  d <- data.frame(y = 1 + 0.5 * x1 + rnorm(n, 0, 0.25) - u, x1 = x1, zm = zm)

  f <- sfm(y ~ x1, model_name = "NTN", data = d, muhet = ~zm)
  expect_true(all(c("Zmu.(Intercept)", "Zmu.zm") %in% names(coef(f))))
  expect_lt(abs(coef(f)[["Zmu.zm"]] - 1.0), 3 * f$std.errors[["Zmu.zm"]])
  expect_lt(abs(coef(f)[["x1"]] - 0.5), 3 * f$std.errors[["x1"]])
  expect_false(is.null(f$het$mu))

  expect_error(sfm(y ~ x1, "NHN", d, muhet = ~zm), "PRE-TRUNCATION MEAN")
})

test_that("unsupported combinations refuse rather than ignore", {
  skip_on_cran()
  d <- .het_data(n = 200)
  expect_error(sfm(y ~ x1, "NGE", d, vhet = ~zv), "implemented for model_name")
  expect_error(sfm(y ~ x1 | zu, "NHN_Z", d, uhet = ~zv), "specified twice")
  expect_error(sfm(y ~ x1, "NHN", d, vhet = ~nosuchvar), "not found in")
  expect_error(sfm(y ~ x1, "NHN", d, vhet = ~zv, estimator = "mols"), "moment estimator")
  expect_error(sfm(y ~ x1, "NHN", d, vhet = ~zv, robust = "mdpd"), "homoskedastic")
})

test_that("uhet is the pipe segment by another name", {
  skip_on_cran()
  d <- .het_data(n = 500)
  a <- sfm(y ~ x1 | zu, "NHN_Z", d, vhet = ~zv)
  b <- sfm(y ~ x1, "NHN", d, uhet = ~zu, vhet = ~zv)
  expect_equal(coef(a), coef(b), tolerance = 1e-5)
})

test_that("marginal_effects reads a heteroskedastic fit unchanged", {
  skip_on_cran()
  d <- .het_data(n = 400)
  f <- sfm(y ~ x1 | zu, "NHN_Z", d, vhet = ~zv)
  me <- marginal_effects(f)
  expect_true(all(c("sigma_u", "E_u", "Var_u", "dE_u.dzu") %in% colnames(me)))
  expect_equal(nrow(me), nrow(d))
  ## delta_zu > 0, so raising zu must raise expected inefficiency.
  expect_true(all(me[, "dE_u.dzu"] > 0))
})

test_that("omitting the het arguments changes nothing", {
  skip_on_cran()
  d <- .het_data(n = 400)
  a <- sfm(y ~ x1, "NHN", d)
  b <- sfm(y ~ x1, "NHN", d, vhet = NULL, uhet = NULL, muhet = NULL)
  expect_identical(coef(a), coef(b))
})

test_that("rows missing a het covariate are dropped, not silently misaligned", {
  skip_on_cran()
  d <- .het_data(n = 400)
  d$zv[c(3, 17, 200)] <- NA
  f <- sfm(y ~ x1, "NHN", d, vhet = ~zv)
  expect_equal(nobs(f), nrow(d) - 3L)
  expect_length(f$sigma_v, nrow(d) - 3L)
})
