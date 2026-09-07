## Tests for the robust-divergence tuning additions.
##
## These exercise the parts that do not need a fitted object: the H-score in
## both evaluations, the weight-matching calibration, and the density-power
## weights.
##
## The second half tests the ADAPTERS, which is where integration actually went
## wrong. The staged versions of these functions assumed sfm() reports NHN as
## (sigma_v, sigma_u) and that start_from takes a parameter vector. Neither is
## true -- sfm() reports (lambda, sigma), and start_from takes an "sfareg"
## object and matches by name -- so every candidate fit failed and
## hscore_select() aborted. Both assumptions are pinned below, because both
## are silent when wrong: the first stops, but the second is swallowed by the
## tryCatch that exists to tolerate a genuinely non-converging candidate.

test_that("the two H-score evaluations agree where both are computable", {
  set.seed(20260903)
  e <- rnorm(500, 0, 0.30) - abs(rnorm(500, 0, 0.60))
  for (cc in c(0, 0.05, 0.10, 0.217, 0.40)) {
    a <- hscore(e, 0.30, 0.60, cc, stable = TRUE)
    b <- hscore(e, 0.30, 0.60, cc, stable = FALSE)
    expect_true(is.finite(a))
    if (is.finite(b)) expect_equal(a, b, tolerance = 1e-8)
  }
})

test_that("the log-space evaluation survives residuals that overflow the natural one", {
  ## An outlier far enough out that its fitted density underflows to zero. The
  ## natural-scale form then contains f^(c-2) = Inf; the log-space form does not.
  ## This is the failure the paper documents: it is silent, and it biases the
  ## selected tuning parameter toward maximum likelihood.
  set.seed(1)
  e <- c(rnorm(200, 0, 0.30) - abs(rnorm(200, 0, 0.60)), -60)
  stable  <- hscore(e, 0.30, 0.60, 0.15, stable = TRUE)
  natural <- hscore(e, 0.30, 0.60, 0.15, stable = FALSE)
  expect_true(is.finite(stable))
  expect_false(is.finite(natural))
})

test_that("hscore rejects invalid arguments", {
  e <- rnorm(50)
  expect_error(hscore(e, -1,  0.6, 0.2), "sigma_v")
  expect_error(hscore(e, 0.3, 0,   0.2), "sigma_u")
  expect_error(hscore(e, 0.3, 0.6, -0.1), "non-negative")
})

test_that("calibrate_c reproduces the paper's fixed calibration", {
  ## Design of the Monte Carlo study: sigma_v = 0.30, sigma_u = 0.60, so
  ## lambda = 2. The paper reports c_MLqE = 0.217 and c_MDPD = 0.540.
  cal <- calibrate_c(sigma_v = 0.30, sigma_u = 0.60, method = "all")
  expect_named(cal, c("mlqe", "psi", "mdpd"))
  expect_equal(unname(cal[["mlqe"]]), 0.217, tolerance = 5e-3)
  expect_equal(unname(cal[["mdpd"]]), 0.540, tolerance = 5e-3)
})

test_that("psi and mdpd calibrate to the same value", {
  ## They are the same estimator: the objectives differ by the factor (1 + c),
  ## which is constant in the parameters and so cancels from the ratio the
  ## calibration matches.
  cal <- calibrate_c(0.30, 0.60, method = "all")
  expect_equal(unname(cal[["psi"]]), unname(cal[["mdpd"]]), tolerance = 1e-6)
})

test_that("the calibrated c does what it says", {
  ## By construction a residual 3 scale units out should retain 10 percent of
  ## its clean-data estimating-equation influence. To solver precision, not to
  ## a grid step: the tolerance here is what caught the grid fallback.
  for (m in c("mlqe", "psi", "mdpd")) {
    cc <- calibrate_c(0.30, 0.60, method = m)
    expect_equal(sfa:::.weight_ratio(m, 0.30, 0.60, unname(cc), k = 3),
                 0.10, tolerance = 1e-8)
  }
})

test_that("the influence ratio is non-monotone, and both roots are reported", {
  ## Not a curiosity: it is why uniroot() over the whole interval is refused,
  ## and why the first implementation silently returned a grid point instead of
  ## a root. Psi/MDPD cross the target twice; MLqE crosses once.
  cal <- calibrate_c(0.30, 0.60, method = "all")
  rt  <- attr(cal, "roots")
  expect_length(rt[["mlqe"]], 1L)
  expect_length(rt[["mdpd"]], 2L)
  ## The larger root is the one returned.
  expect_equal(unname(cal[["mdpd"]]), max(rt[["mdpd"]]))
  expect_lt(rt[["mdpd"]][1], rt[["mdpd"]][2])
  ## Both really are roots, and the ratio really does dip below the target in
  ## between -- which is what makes them two crossings rather than one.
  for (r in rt[["mdpd"]]) {
    expect_equal(sfa:::.weight_ratio("mdpd", 0.30, 0.60, r, k = 3), 0.10,
                 tolerance = 1e-8)
  }
  mid <- mean(rt[["mdpd"]])
  expect_lt(sfa:::.weight_ratio("mdpd", 0.30, 0.60, mid, k = 3), 0.10)
})

test_that("density_weights are in (0, 1] and order with the density", {
  set.seed(2)
  e <- rnorm(300, 0, 0.30) - abs(rnorm(300, 0, 0.60))
  w <- density_weights(e, sigma_v = 0.30, sigma_u = 0.60, c = 0.217)
  expect_length(w, length(e))
  expect_true(all(w > 0 & w <= 1))
  expect_equal(max(w), 1)
  f <- sfa:::.dens_nhn(e, 0.30, 0.60)
  expect_equal(order(w), order(f))
})

test_that("density_weights are all one at the maximum likelihood endpoint", {
  e <- rnorm(50, 0, 0.3) - abs(rnorm(50, 0, 0.6))
  expect_equal(density_weights(e, 0.30, 0.60, c = 0), rep(1, 50))
})

test_that("density_weights requires the scale parameters for a bare vector", {
  expect_error(density_weights(rnorm(10)), "required")
})

## ---------------------------------------------------------------------------
## The adapters, against a real sfm() fit
## ---------------------------------------------------------------------------

nhn_fit <- function(seed = 3, n = 200, ...) {
  set.seed(seed)
  x <- runif(n, 1, 10)
  dat <- data.frame(y = 1 + 0.5 * log(x) + rnorm(n, 0, 0.3) -
                        abs(rnorm(n, 0, 0.6)), x = log(x))
  sfm(y ~ x, data = dat, model_name = "NHN", ...)
}

test_that("the scale accessor inverts sfm()'s (lambda, sigma) reporting", {
  skip_on_cran()
  fit <- nhn_fit()
  ## The layout the adapter has to cope with. If this assertion ever fails the
  ## package changed how it reports NHN, and .robust_scales() must change too.
  expect_true(all(c("lambda", "sigma") %in% names(fit$coefficients)))
  expect_false(any(c("sigma_v", "sigma_u") %in% names(fit$coefficients)))

  sv <- sfa:::.robust_get(fit, "sigma_v")
  su <- sfa:::.robust_get(fit, "sigma_u")
  ## The defining identities, not a restatement of the same two lines.
  expect_equal(su / sv, unname(fit$coefficients[["lambda"]]))
  expect_equal(sqrt(sv^2 + su^2), unname(fit$coefficients[["sigma"]]))
  expect_true(sv > 0 && su > 0)

  ## A direct (sigma_v, sigma_u) layout is still honoured.
  expect_equal(sfa:::.robust_scales(c(sigma_v = 0.3, sigma_u = 0.6)),
               c(sigma_v = 0.3, sigma_u = 0.6))
  expect_null(sfa:::.robust_scales(c(a = 1, b = 2)))
  expect_error(sfa:::.robust_get(structure(list(coefficients = c(a = 1)),
                                           class = "sfareg"), "sigma_v"),
               "neither")
})

test_that("residuals are rebuilt in the v - u orientation", {
  skip_on_cran()
  fit <- nhn_fit()
  e <- sfa:::.robust_residuals(fit)
  expect_length(e, fit$nobs)
  ## Reconstructed independently here, so this is a check and not an echo.
  d <- eval(fit$call$data, envir = environment(fit$formula))
  X <- stats::model.matrix(fit$formula, stats::model.frame(fit$formula, d))
  b <- fit$coefficients[colnames(X)]
  expect_equal(e, as.numeric(d$y - X %*% b))
  ## A production frontier subtracts u, so the residuals must lean negative.
  expect_lt(mean(e), 0)
})

test_that("the refit closure passes a fit, not a parameter vector", {
  skip_on_cran()
  ## sfm(start_from = ) requires an "sfareg" object and matches by parameter
  ## name; handed a numeric vector it stops. The closure's tryCatch would
  ## swallow that and report a failed candidate, so the bug is only visible
  ## from the outside as "no candidate was scorable". Both halves are asserted.
  fit <- nhn_fit()
  expect_error(sfm(y ~ x, data = eval(fit$call$data,
                                      envir = environment(fit$formula)),
                   model_name = "NHN",
                   start_from = as.numeric(fit$coefficients)),
               "sfareg")
  expect_s3_class(sfa:::.robust_ref_start(fit), "sfareg")

  refit <- sfa:::.robust_refit_fn(fit, "mlqe")
  r0 <- refit(0, fit)
  r2 <- refit(0.2, fit)
  expect_true(r0$ok)
  expect_true(r2$ok)
  ## c = 0 is the maximum likelihood endpoint, so the refit must reproduce the
  ## original fit rather than merely converge somewhere.
  expect_equal(r0$sigma_v, sfa:::.robust_get(fit, "sigma_v"), tolerance = 1e-6)
  expect_equal(r0$objective, as.numeric(stats::logLik(fit)), tolerance = 1e-6)
  ## Robustification pulls the scales in on clean data with a few long tails.
  expect_lt(r2$sigma_u, r0$sigma_u)

  expect_error(sfa:::.robust_refit_fn(fit, "mlqe")(0.1, "not a fit")$ok, NA)
  expect_false(refit(0.1, "not a fit")$ok)
})

test_that("density_weights reads the tuning value off a robust fit", {
  skip_on_cran()
  fit <- nhn_fit(robust = "mlqe", c_mlqe = 0.217)
  expect_equal(fit$robust_c, 0.217)
  w <- density_weights(fit)
  expect_length(w, fit$nobs)
  expect_true(all(w > 0 & w <= 1))
  expect_equal(max(w), 1)
  ## The weight orders inversely with distance from the fitted surface.
  e <- sfa:::.robust_residuals(fit)
  expect_equal(order(w), order(sfa:::.dens_nhn(e,
    sfa:::.robust_get(fit, "sigma_v"), sfa:::.robust_get(fit, "sigma_u"))))
  ## An MLE fit carries no tuning value, so every weight is one.
  expect_equal(density_weights(nhn_fit()), rep(1, 200))
})

## ---------------------------------------------------------------------------
## Integration test (skipped unless a fitted NHN model can be produced quickly)
## ---------------------------------------------------------------------------
test_that("hscore_select returns a coherent object", {
  skip_on_cran()
  skip_if_not(exists("sfm"), "sfm() not available")
  set.seed(3)
  n  <- 200
  x  <- runif(n, 1, 10)
  y  <- 1 + 0.5 * log(x) + rnorm(n, 0, 0.3) - abs(rnorm(n, 0, 0.6))
  dat <- data.frame(y = y, x = log(x))
  fit <- try(sfm(y ~ x, data = dat, model_name = "NHN"), silent = TRUE)
  skip_if(inherits(fit, "try-error"), "could not fit the reference model")

  sel <- hscore_select(fit, method = "mlqe", range = c(0, 0.2), coarse = 0.05,
                       fine = 0.01, half_width = 0.05)
  expect_s3_class(sel, "sfa_hscore")
  expect_true(sel$c >= 0 && sel$c <= 0.2)
  expect_true(is.finite(sel$hscore))
  expect_true(all(c("c", "hscore", "stage") %in% names(sel$path)))
  ## the reported minimum really is the minimum over the scorable path
  expect_equal(sel$hscore, min(sel$path$hscore, na.rm = TRUE), tolerance = 1e-10)
})
