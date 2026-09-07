## The sample-selection stochastic frontier: selsfm(), Greene (2010).

## Greene's Equation (10) exactly: the selection disturbance w is correlated
## with the frontier NOISE v, not with the inefficiency u.
sel_data <- function(seed, n = 600, rho = 0.6, su = 1.0, sv = 0.5,
                     b = c(0.5, 0.8, -0.3), al = c(0.2, 0.7, -0.5),
                     cost = FALSE) {
  set.seed(seed)
  z1 <- rnorm(n)
  z2 <- rnorm(n)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  w <- rnorm(n)
  d <- as.numeric(al[1] + al[2] * z1 + al[3] * z2 + w > 0)
  v <- sv * (rho * w + sqrt(1 - rho^2) * rnorm(n))
  u <- su * abs(rnorm(n))
  y <- b[1] + b[2] * x1 + b[3] * x2 + v + if (cost) u else -u
  y[d == 0] <- NA
  data.frame(y = y, x1 = x1, x2 = x2, z1 = z1, z2 = z2, d = d)
}

## Small Nsim keeps these tests fast. That trips the draw-count warning, which
## is asserted on its own below; silence it where it is incidental so it does
## not surface as noise in R CMD check.
sel_fit <- function(...) suppressWarnings(selsfm(...))

test_that("selsfm returns a well-formed sfareg with the documented layout", {
  f <- sel_fit(d ~ z1 + z2, y ~ x1 + x2, sel_data(1), Nsim = 60, seed = 1)

  expect_s3_class(f, "sfareg")
  expect_identical(f$model_name, "SEL")

  ## out is p x 3 -- one ROW per parameter -- and the convenience vectors are
  ## its columns. The source builds it 3 x p and stores t(out), so this is the
  ## orientation a user actually sees.
  expect_equal(ncol(f$out), 3L)
  expect_equal(nrow(f$out), 6L)
  expect_identical(colnames(f$out), c("par", "st_err", "t-val"))
  expect_identical(
    rownames(f$out),
    c("sigma_u", "sigma_v", "rho", "(Intercept)", "x1", "x2")
  )
  expect_equal(unname(f$coefficients), unname(f$out[, "par"]))
  expect_equal(unname(f$std.errors), unname(f$out[, "st_err"]))

  ## Scales are reported as magnitudes; rho stays inside the unit interval.
  expect_gt(f$out["sigma_u", "par"], 0)
  expect_gt(f$out["sigma_v", "par"], 0)
  expect_lt(abs(f$out["rho", "par"]), 1)

  ## One inefficiency score per SELECTED observation, not per row of data.
  expect_equal(length(f$jlms), f$n_selected)
  expect_true(all(f$jlms >= 0))
  expect_equal(f$efficiency, exp(-f$jlms))
  expect_lt(f$n_selected, f$n_total)
})

test_that("the sfareg generics use the SELECTED count, not the rows supplied", {
  f <- sel_fit(d ~ z1 + z2, y ~ x1 + x2, sel_data(1), Nsim = 60, seed = 1)

  ## The second-stage likelihood is a sum over the selected observations only,
  ## so that is the n BIC() has to divide by. Without an explicit `nobs` on the
  ## object, nobs.sfareg() falls back to re-evaluating the call's `data` and
  ## reports the whole sample -- the same defect that was fixed once before for
  ## rows lost to missingness.
  expect_equal(nobs(f), f$n_selected)
  expect_lt(nobs(f), f$n_total)

  k <- nrow(f$out)
  expect_equal(unname(AIC(f)), unname(-2 * as.numeric(logLik(f)) + 2 * k))
  expect_equal(
    unname(BIC(f)),
    unname(-2 * as.numeric(logLik(f)) + log(f$n_selected) * k)
  )

  expect_equal(unname(coef(f)), unname(f$out[, "par"]))
  expect_equal(dim(vcov(f)), c(k, k))
  expect_silent(invisible(capture.output(print(f))))
})

test_that("selsfm recovers the parameters it was simulated from", {
  f <- sel_fit(d ~ z1 + z2, y ~ x1 + x2, sel_data(7, n = 2000),
    Nsim = 200, seed = 7
  )
  p <- f$out[, "par"]

  ## The frontier slopes are the best identified parameters and are held to a
  ## tighter bar than the scales; rho is the loosest, as it is throughout the
  ## selection literature.
  expect_equal(unname(p["x1"]), 0.8, tolerance = 0.12)
  expect_equal(unname(p["x2"]), -0.3, tolerance = 0.12)
  expect_equal(unname(p["sigma_u"]), 1.0, tolerance = 0.25)
  expect_equal(unname(p["sigma_v"]), 0.5, tolerance = 0.25)
  expect_gt(p["rho"], 0)

  ## The first stage is an ordinary probit and should be close.
  expect_equal(unname(f$alpha[2]), 0.7, tolerance = 0.2)
  expect_equal(unname(f$alpha[3]), -0.5, tolerance = 0.2)
})

test_that("a cost frontier flips the sign of the inefficiency term", {
  dd <- sel_data(3, n = 1200, cost = TRUE)
  fc <- sel_fit(d ~ z1 + z2, y ~ x1 + x2, dd, Nsim = 100, seed = 3,
    inefdec = FALSE
  )
  expect_equal(fc$S, -1)
  expect_equal(unname(fc$out["sigma_u", "par"]), 1.0, tolerance = 0.35)
  expect_equal(unname(fc$out["x1", "par"]), 0.8, tolerance = 0.15)

  ## Fitting cost data as a production frontier should not reproduce sigma_u:
  ## the skew is the wrong way round for it.
  fp <- sel_fit(d ~ z1 + z2, y ~ x1 + x2, dd, Nsim = 100, seed = 3)
  expect_equal(fp$S, 1)
  expect_gt(
    abs(fp$out["sigma_u", "par"] - 1.0),
    abs(fc$out["sigma_u", "par"] - 1.0)
  )
})

test_that("the same seed gives the same fit, and the caller's RNG is untouched", {
  dat <- sel_data(11, n = 400)
  a <- sel_fit(d ~ z1 + z2, y ~ x1 + x2, dat, Nsim = 50, seed = 99)
  b <- sel_fit(d ~ z1 + z2, y ~ x1 + x2, dat, Nsim = 50, seed = 99)
  expect_equal(a$out[, "par"], b$out[, "par"])

  ## .sml_draws() snapshots and restores the stream; a fit must not move it.
  set.seed(1234)
  before <- runif(1)
  set.seed(1234)
  invisible(sel_fit(d ~ z1 + z2, y ~ x1 + x2, dat, Nsim = 50, seed = 99))
  expect_equal(runif(1), before)
})

test_that("selsfm rejects malformed calls", {
  dat <- sel_data(2, n = 300)

  expect_error(sel_fit("d ~ z1", y ~ x1, dat), "must be a two-sided formula")
  expect_error(sel_fit(d ~ z1, ~x1, dat), "must have a left-hand side")
  expect_error(sel_fit(d ~ z1 + z2, y ~ x1, as.matrix(dat)), "must be a data.frame")
  expect_error(
    sel_fit(d ~ z1 + z2, y ~ x1 | x2, dat),
    "must not contain a `\\|` segment"
  )
  expect_error(
    sel_fit(d ~ z1 + z2, y ~ x1 + x2, dat, inefdec = NA),
    "must be TRUE or FALSE"
  )
  expect_error(
    sel_fit(d ~ z1 + z2, y ~ x1 + x2, dat, maxit.optim = 0),
    "must be a single positive number"
  )

  ## A non-binary selection response is a modelling error, not something to
  ## coerce silently.
  bad <- dat
  bad$d <- bad$d + 0.5
  expect_error(sel_fit(d ~ z1 + z2, y ~ x1 + x2, bad), "must be binary")

  ## Nothing to correct for if everything is selected.
  allsel <- dat
  allsel$d <- 1
  expect_error(
    sel_fit(d ~ z1 + z2, y ~ x1 + x2, allsel),
    "both selected and unselected"
  )
})

test_that("too few draws warns rather than failing quietly", {
  ## Simulated ML is consistent only if the draw count grows with n, so a count
  ## below the automatic floor is a correctness problem worth saying out loud.
  expect_warning(
    selsfm(d ~ z1 + z2, y ~ x1 + x2, sel_data(5, n = 400), Nsim = 12, seed = 5),
    "below the"
  )
})

test_that("a binary factor or logical selector is accepted", {
  dat <- sel_data(13, n = 400)
  num <- sel_fit(d ~ z1 + z2, y ~ x1 + x2, dat, Nsim = 50, seed = 4)

  lg <- dat
  lg$d <- as.logical(lg$d)
  expect_equal(
    sel_fit(d ~ z1 + z2, y ~ x1 + x2, lg, Nsim = 50, seed = 4)$out[, "par"],
    num$out[, "par"]
  )

  fc <- dat
  fc$d <- factor(fc$d, levels = c(0, 1), labels = c("no", "yes"))
  expect_equal(
    sel_fit(d ~ z1 + z2, y ~ x1 + x2, fc, Nsim = 50, seed = 4)$out[, "par"],
    num$out[, "par"]
  )
})
