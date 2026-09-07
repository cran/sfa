## TIC() and vuong(): non-nested model selection (L1).
##
## What is asserted here is the THEORY, not a measured number. The Takeuchi
## penalty equals the parameter count exactly when the information matrix
## equality holds, and the equality holds exactly when the model is correctly
## specified -- so a correctly specified fit is a case where the right answer
## is known in advance and does not depend on the platform's optimizer path.

ms_data <- function(seed = 42, n = 1500, u = c("hn", "exp")) {
  u <- match.arg(u)
  set.seed(seed)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  ui <- if (u == "hn") abs(rnorm(n, 0, 1)) else rexp(n, rate = 1)
  data.frame(
    y = 1 + 0.5 * x1 + 0.5 * x2 + rnorm(n, 0, 0.4) - ui,
    x1 = x1, x2 = x2
  )
}

ms_fit <- function(d, nm) {
  sfm(y ~ x1 + x2, model_name = nm, data = d, keep_objective = TRUE)
}

test_that("the per-observation log-likelihood sums to the fitted value", {
  skip_on_cran()
  ## Everything below is built on this. If per_obs did not decompose the
  ## objective, TIC and vuong would both be quietly wrong rather than broken.
  d <- ms_data()
  f <- ms_fit(d, "NHN")
  expect_equal(sum(f$objective(f$opt$par, per_obs = TRUE)), -f$opt$value,
    tolerance = 1e-8
  )
  expect_length(f$objective(f$opt$par, per_obs = TRUE), nrow(d))
})

test_that("the Takeuchi penalty collapses to df under correct specification", {
  skip_on_cran()
  ## H = I when the model is right, so tr[H I^-1] = p. This is the property
  ## that makes TIC correct; a bug in either matrix breaks it immediately.
  ## The tolerance is loose because the penalty is a sampling quantity, not an
  ## identity -- at n = 1500 it lands within a few percent.
  d <- ms_data(u = "hn")
  z <- TIC(ms_fit(d, "NHN"), detail = TRUE)
  expect_equal(z$penalty, z$df, tolerance = 0.15)
  expect_equal(z$ratio, 1, tolerance = 0.15)
  ## And TIC therefore nearly equals AIC in this case, which is the sense in
  ## which TIC "generalizes" AIC rather than disagreeing with it.
  expect_equal(z$TIC, z$AIC, tolerance = 0.01)
})

test_that("TIC's pieces are internally consistent", {
  skip_on_cran()
  d <- ms_data()
  f <- ms_fit(d, "NHN")
  z <- TIC(f, detail = TRUE)
  expect_equal(z$TIC, -2 * z$logLik + 2 * z$penalty)
  expect_equal(z$AIC, unname(AIC(f)))
  expect_equal(z$logLik, as.numeric(logLik(f)))
  expect_equal(z$df, length(coef(f)))
  ## The bare call returns just the number.
  expect_equal(TIC(f), z$TIC)
})

test_that("vuong points at the correctly specified model", {
  skip_on_cran()
  ## The DIRECTION of the statistic is the property: whichever model is true,
  ## the statistic must lean toward it. SIGNIFICANCE is not asserted here --
  ## it is a sampling outcome, and on this DGP at n = 1500 the half-normal
  ## case lands at z = 1.95 against a critical value of 1.96, which would flip
  ## from platform to platform for reasons that have nothing to do with the
  ## code. Two runs with different true distributions, so a routine that
  ## always answered "model 1" could not pass both.
  d_hn <- ms_data(u = "hn")
  v1 <- vuong(ms_fit(d_hn, "NHN"), ms_fit(d_hn, "NE"))
  expect_s3_class(v1, "sfa_vuong")
  expect_gt(v1$statistic, 0)

  d_ex <- ms_data(u = "exp")
  v2 <- vuong(ms_fit(d_ex, "NHN"), ms_fit(d_ex, "NE"))
  expect_lt(v2$statistic, 0)
})

test_that("vuong can actually reject, and names the winner when it does", {
  skip_on_cran()
  ## The exponential case is the well-powered one on this DGP (z is around
  ## -3.7, nowhere near the boundary), so it is the one used to check that a
  ## verdict is reached rather than only that a number comes back.
  d <- ms_data(u = "exp")
  v <- vuong(ms_fit(d, "NHN"), ms_fit(d, "NE"))
  expect_lt(v$statistic, -v$critical)
  expect_identical(v$favoured, "NE")
  expect_lt(v$p.value, 0.05)
})

test_that("vuong is antisymmetric in its arguments", {
  skip_on_cran()
  ## Swapping the models must negate the statistic and leave the p-value
  ## alone. A sign error would otherwise be invisible whenever the favoured
  ## model happened to be passed first.
  d <- ms_data()
  a <- ms_fit(d, "NHN")
  b <- ms_fit(d, "NE")
  v1 <- vuong(a, b)
  v2 <- vuong(b, a)
  expect_equal(v1$statistic, -v2$statistic)
  expect_equal(v1$p.value, v2$p.value)
  expect_equal(v1$lr, -v2$lr)
  expect_identical(v1$favoured, v2$favoured)
})

test_that("the log-likelihood difference is the sum of the per-observation ones", {
  skip_on_cran()
  d <- ms_data()
  a <- ms_fit(d, "NHN")
  b <- ms_fit(d, "NE")
  v <- vuong(a, b)
  expect_equal(v$lr, as.numeric(logLik(a)) - as.numeric(logLik(b)),
    tolerance = 1e-6
  )
  expect_equal(v$n, nobs(a))
})

test_that("the corrections shift the statistic and are flagged as inexact", {
  skip_on_cran()
  ## Lai and Huang derive the N(0,1) limit for the UNCORRECTED statistic only,
  ## and say the corrected one needs further investigation. The object must
  ## carry that distinction so print() can report it.
  d <- ms_data()
  a <- ms_fit(d, "NHN")
  b <- ms_fit(d, "NE")
  plain <- vuong(a, b)
  aic_c <- vuong(a, b, correction = "aic")
  tic_c <- vuong(a, b, correction = "tic")

  expect_true(plain$exact_null)
  expect_false(aic_c$exact_null)
  expect_false(tic_c$exact_null)
  ## Equal parameter counts, so the AIC correction is exactly a no-op here --
  ## which is itself the check that the penalty enters as a difference.
  expect_equal(aic_c$lr_adjusted, plain$lr)
  expect_equal(tic_c$penalty, c(TIC(a, detail = TRUE)$penalty,
    TIC(b, detail = TRUE)$penalty))
  expect_output(print(tic_c), "further")
})

test_that("print reports the models, the verdict and the p-value", {
  skip_on_cran()
  d <- ms_data()
  v <- vuong(ms_fit(d, "NHN"), ms_fit(d, "NE"))
  expect_output(print(v), "NHN")
  expect_output(print(v), "NE")
  expect_output(print(v), "Vuong")
  expect_output(print(v), "p-value")
})

test_that("the refusals name the fix", {
  skip_on_cran()
  d <- ms_data()
  bare <- sfm(y ~ x1 + x2, model_name = "NHN", data = d)   # no keep_objective

  expect_error(vuong(bare, bare), "keep_objective")
  expect_error(TIC(list(a = 1)), "sfareg")
  expect_error(vuong(ms_fit(d, "NHN"), "not a fit"), "sfareg")

  ## Same model twice: the differences are identically zero, so the statistic
  ## is 0/0 and must be refused rather than returned as NaN.
  f <- ms_fit(d, "NHN")
  expect_error(vuong(f, f), "identical")

  ## Different numbers of observations.
  f2 <- ms_fit(ms_data(n = 900), "NE")
  expect_error(vuong(f, f2), "different numbers of observations")
})

test_that("a moment-based fit has no likelihood-based criterion", {
  skip_on_cran()
  ## C2SLS maximises nothing, so logLik() is NA and TIC must say so rather
  ## than return -2*NA + 2*penalty.
  set.seed(3)
  n <- 600
  w1 <- rnorm(n); w2 <- rnorm(n); x1 <- rnorm(n); eta <- rnorm(n)
  v <- 0.5 * (0.6 * eta + sqrt(1 - 0.36) * rnorm(n))
  x2 <- 0.9 * w1 - 0.7 * w2 + 0.5 * x1 + eta
  d <- data.frame(
    y = 0.5 + 0.8 * x1 - 0.6 * x2 + v - abs(rnorm(n)),
    x1 = x1, x2 = x2, w1 = w1, w2 = w2
  )
  f <- ivsfm(y ~ x1 + x2, endogenous = ~x2, instruments = ~ w1 + w2,
    data = d, model_name = "C2SLS"
  )
  expect_error(suppressWarnings(TIC(f)))
})
