## Coelli (1995) Section 3 and Appendix 1.

.ct_data <- function(N, seed, sig_u = 1, sig_v = 0.5) {
  as.data.frame(data_gen_cs(N = N, rand = seed, sig_u = sig_u, sig_v = sig_v,
    cons = 0.5, beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.1
  ))
}

test_that("the chi-bar-square mixture is the documented one", {
  ## One boundary restriction: half the chi2(1) tail, because chi2(0) is a
  ## point mass at zero and contributes nothing above it.
  expect_equal(.chi_bar_p(3.0, 1L), 0.5 * pchisq(3.0, 1, lower.tail = FALSE))
  ## Coelli's own statement of the rule: the size-alpha critical value equals
  ## the chi2(1) critical value for size 2*alpha, i.e. 2.71 at 5%.
  cv <- qchisq(2 * 0.05, df = 1, lower.tail = FALSE)
  expect_equal(cv, 2.705543, tolerance = 1e-5)
  expect_equal(.chi_bar_p(cv, 1L), 0.05, tolerance = 1e-8)
  ## Two restrictions: 0.25/0.5/0.25 over chi2(0,1,2).
  expect_equal(
    .chi_bar_p(3.0, 2L),
    0.5 * pchisq(3, 1, lower.tail = FALSE) + 0.25 * pchisq(3, 2, lower.tail = FALSE)
  )
  expect_equal(.chi_bar_p(0, 1L), 1)
  expect_equal(.chi_bar_p(-1, 1L), 1)
})

test_that("the one-sided LR p-value is exactly half the naive one", {
  skip_on_cran()
  f <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = .ct_data(300, 7))
  tt <- inefficiency_test(f, test = c("lr_1sided", "lr"))
  expect_equal(tt$statistic[1], tt$statistic[2])
  expect_equal(tt$p.value[1], 0.5 * tt$p.value[2], tolerance = 1e-10)
  ## and the naive test can therefore never reject when the one-sided does not
  expect_true(tt$p.value[1] <= tt$p.value[2])
})

test_that("inefficiency_test() returns a well-formed table and rejects on strong data", {
  skip_on_cran()
  f <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = .ct_data(600, 11, sig_u = 2, sig_v = 0.5))
  tt <- inefficiency_test(f)
  expect_s3_class(tt, "data.frame")
  expect_equal(nrow(tt), 4L)
  expect_named(tt, c("test", "statistic", "null", "p.value", "reject"))
  expect_true(all(tt$p.value >= 0 & tt$p.value <= 1))
  ## sigma_u = 4 sigma_v is not a subtle alternative
  expect_true(tt$reject[tt$test == "LR (one-sided)"])
  expect_true(is.finite(attr(tt, "logLik_H1")))
  expect_true(attr(tt, "logLik_H1") >= attr(tt, "logLik_H0") - 1e-6)
})

test_that("inefficiency_test() refuses models Coelli's mixture is not stated for", {
  skip_on_cran()
  f <- sfm(y_pcs_e ~ x1 + x2, model_name = "NE", data = .ct_data(200, 3))
  expect_error(inefficiency_test(f), "half-normal")
})

test_that("Coelli's A13 inversion is the package's own moment inversion", {
  ## A13 written out directly from the paper, against .gtre_two_step() --
  ## the same estimator, not merely a similar one.
  set.seed(4)
  e <- rnorm(50000, 0, 0.6) - abs(rnorm(50000, 0, 1.2))
  e <- e - mean(e)
  m2 <- mean(e^2); m3 <- mean(e^3)
  su2_paper <- (m3 / (sqrt(2 / pi) * (1 - 4 / pi)))^(2 / 3)
  sig2_paper <- m2 + (2 / pi) * su2_paper
  ts <- .gtre_two_step(e, e, 0)
  expect_equal(ts$sigmaSq_uv, sig2_paper, tolerance = 1e-10)
  expect_equal(ts$gamma_uv, su2_paper / sig2_paper, tolerance = 1e-10)
})

test_that(".nhn_central_moments() is shared and production-signed", {
  m <- .nhn_central_moments(g = 0.8, ssq = 1.25)
  expect_named(m, c("mu2", "mu3", "mu4", "mu5", "mu6"))
  ## (1 - 4/pi) < 0, so the odd central moments come back NEGATIVE already:
  ## a production frontier's composed error is left-skewed. Flipping them
  ## downstream would be a double negation.
  expect_lt(m[["mu3"]], 0)
  expect_lt(m[["mu5"]], 0)
  expect_gt(m[["mu2"]], 0)
  expect_gt(m[["mu4"]], 0)
  expect_gt(m[["mu6"]], 0)
  ## against a direct simulation of v - u at the same (g, ssq)
  set.seed(3)
  su <- sqrt(0.8 * 1.25); sv <- sqrt(0.2 * 1.25)
  e <- rnorm(4e6, 0, sv) - abs(rnorm(4e6, 0, su))
  e <- e - mean(e)
  expect_equal(m[["mu2"]], mean(e^2), tolerance = 0.01)
  expect_equal(m[["mu3"]], mean(e^3), tolerance = 0.02)
  expect_equal(m[["mu4"]], mean(e^4), tolerance = 0.02)
})

test_that(".cols_se_nhn() is finite and tracks the bootstrap as n grows", {
  skip_on_cran()
  ## cst = 1/(sqrt(2/pi)(1 - 4/pi)) is NEGATIVE, so any formulation taking its
  ## cube root directly returns NaN. This is a regression test for that.
  se <- .cols_se_nhn(sigma_u = 0.94, sigma_v = 0.54, n = 800)
  expect_true(all(is.finite(se)))
  expect_true(all(se > 0))
  expect_named(se, c("sigma_v", "sigma_u", "eu"))
  ## sqrt(n)-consistent: quadrupling n halves the standard error
  se4 <- .cols_se_nhn(0.94, 0.54, n = 3200)
  expect_equal(unname(se4[["sigma_u"]] / se[["sigma_u"]]), 0.5, tolerance = 1e-8)

  ## and it agrees with a nonparametric bootstrap of the same estimator
  d <- .ct_data(3200, 7)
  an <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = d, estimator = "cols")
  bo <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = d, estimator = "cols",
    cols_boot = 400, rand.cols = 2
  )
  bse <- apply(bo$cols_boot_draws, 2, stats::sd, na.rm = TRUE)
  expect_equal(unname(an$out["sigu", "st_err"]), unname(bse[["sigu"]]), tolerance = 0.12)
  expect_equal(unname(an$out["sigv", "st_err"]), unname(bse[["sigv"]]), tolerance = 0.15)
})

test_that("estimator = \"cols\" reports analytic errors for NHN and NA elsewhere", {
  skip_on_cran()
  d <- .ct_data(600, 7)
  fh <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = d, estimator = "cols")
  expect_true(is.finite(fh$out["sigu", "st_err"]))
  expect_true(is.finite(fh$out["sigv", "st_err"]))
  ## Coelli derives the variance for the half-normal only; the moment
  ## inversion is distribution-specific, so NE keeps NA and the bootstrap.
  fe <- sfm(y_pcs_e ~ x1 + x2, model_name = "NE", data = d, estimator = "cols")
  expect_true(is.na(fe$out["sigu", "st_err"]))
})
