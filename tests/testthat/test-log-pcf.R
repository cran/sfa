## log D_nu(z), the parabolic cylinder function behind NNAK (part of G8).

test_that("the integral fallback agrees with gsl::hyperg_U where both work", {
  ## Two independent routes to the same quantity: the hypergeometric identity
  ## and the integral representation. Agreement across the grid is what makes
  ## the fallback trustworthy where only one of them survives.
  worst <- 0
  for (nu in c(-0.5, -1, -2, -4, -8, -12)) {
    for (z in c(0.6, 1, 2, 5, 10, 20)) {
      u <- tryCatch(gsl::hyperg_U(-nu / 2, 0.5, z^2 / 2), error = function(e) NA_real_)
      if (!is.finite(u) || u <= 0) next
      a <- (nu / 2) * log(2) - z^2 / 4 + log(u)
      b <- .log_pcf_integral(nu, z)
      expect_true(is.finite(b))
      worst <- max(worst, abs(a - b) / max(abs(a), 1e-12))
    }
  }
  expect_lt(worst, 1e-6)
})

test_that("a large shape no longer returns NA", {
  ## gsl::hyperg_U THROWS for a/2 >~ 16 at small arguments. The call is
  ## vectorized, so one bad element used to null the whole vector -- which is
  ## why NNAK "failed outright" once its shape wandered past about 8. The
  ## fallback is now applied per element.
  for (nu in c(-32, -64)) {
    v <- .log_pcf(nu, c(1, 5, 20, 40))
    expect_true(all(is.finite(v)), info = paste("nu =", nu))
    ## log D_nu is decreasing in z for these orders.
    expect_true(all(diff(v) < 0), info = paste("nu =", nu))
  }
})

test_that("small shapes, where the integrand has no interior peak, still work", {
  ## For a <= 1 the exponent is decreasing on (0, Inf) with no finite maximum
  ## to factor out, which is a different branch. It used to warn "NaNs produced"
  ## from sqrt() of a negative discriminant.
  expect_silent(v <- .log_pcf_integral(-0.5, c(0.6, 1, 5)))
  expect_true(all(is.finite(v)))
  expect_silent(.log_pcf_integral(-1, 2))
})

test_that("NNAK fits at shapes that used to kill it", {
  skip_on_cran()
  for (m in c(1, 4)) {
    d <- data_gen_cs(N = 400, rand = 7, sig_u = 1, sig_v = 0.3, cons = 0.5,
      beta1 = 0.5, beta2 = 0.5, a = 1, mu = 0.5, m_nak = m
    )
    f <- try(suppressWarnings(
      sfm(y_pcs_nak ~ x1 + x2, model_name = "NNAK", data = d)
    ), silent = TRUE)
    expect_false(inherits(f, "try-error"), info = paste("m_nak =", m))
    expect_true(is.finite(f$opt$value), info = paste("m_nak =", m))
  }
  ## NOTE: this fixes the hard FAILURE, not NNAK's identification. Recovery of
  ## sigma_u at large shape remains poor and is a separate problem.
})

test_that("log D_nu accepts a VECTOR order, one per observation", {
  ## This is what lets the gamma/Nakagami shape depend on covariates (G5).
  z <- c(-2, -0.3, 0.6, 2, 8)

  ## A scalar order must recycle to exactly what it did before.
  expect_equal(.log_pcf(-4, z), .log_pcf(rep(-4, length(z)), z))

  ## And a genuinely varying order must match element-by-element evaluation.
  nu <- c(-1, -2, -4, -8, -16)
  v <- .log_pcf(nu, z)
  ref <- vapply(seq_along(z), function(i) .log_pcf(nu[i], z[i]), numeric(1))
  expect_equal(v, ref)
  expect_true(all(is.finite(v)))
})

test_that("NG's fitted values are unchanged by the log_pcf work", {
  skip_on_cran()
  ## A REGRESSION PIN, and it earned its place. An earlier version of the
  ## large-order fix also replaced the series branch's floor
  ## (`log(pmax(br, xmin))`) with the integral representation, on the reasoning
  ## that a finite-but-wrong value is worse than a correct one. That is true in
  ## isolation and wrong here: NG evaluates on the series branch for 99.8% of
  ## observations at its own optimum, and the floor acts as a barrier keeping
  ## the optimizer out of the sigma_v -> 0 corner where the composed likelihood
  ## is unbounded. Removing it moved NG from (sigv 0.185, x1 0.560) to
  ## (sigv 0.000, x1 0.775) against a true x1 of 0.5 -- at a HIGHER likelihood,
  ## which is exactly why no likelihood-based check would have caught it.
  d <- data_gen_cs(N = 400, rand = 1, sig_u = 1, sig_v = 0.3, cons = 0.5,
    beta1 = 0.5, beta2 = 0.5, a = 1, mu = 0.5
  )
  f <- suppressWarnings(sfm(y_pcs_g ~ x1 + x2, model_name = "NG", data = d))
  p <- f$out[, "par"]

  expect_equal(unname(p[["sigv"]]), 0.185, tolerance = 0.02)
  expect_equal(unname(p[["sigu"]]), 0.373, tolerance = 0.02)
  expect_equal(unname(p[["x1"]]), 0.560, tolerance = 0.02)
  ## The specific failure mode: sigma_v must not be at the boundary.
  expect_gt(p[["sigv"]], 0.05)
})
