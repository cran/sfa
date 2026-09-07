## The NE log-density's second catastrophic cancellation, at sigma_v -> 0.
##
## 1.2.0 fixed the sigma_u -> 0 cancellation with .log_phi_tilt(). That fix
## created a mirror of itself: .log_phi_tilt() returns log Phi(z) + z^2/2, and
## for eps < 0 with sigma_v small the argument goes large POSITIVE, so the
## returned value IS essentially z^2/2 -- and the caller then subtracted
## eps^2/(2 sigma_v^2) of the same magnitude. These pin the shape of the
## objective, because the symptom was a fit that looked entirely ordinary.

ne_obj <- function(seed = 11, n = 1000) {
  d <- as.data.frame(data_gen_cs(N = n, rand = seed, sig_u = 1, sig_v = 0.3,
    cons = 0.5, beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.1,
    m_nak = 1, lam_tsl = 1.5))
  f <- suppressWarnings(sfm(y_pcs_e ~ x1 + x2, model_name = "NE", data = d,
    keep_objective = TRUE))
  list(fit = f, lf = f$objective, data = d)
}

test_that("the objective does not dive as sigma_v goes to zero", {
  skip_on_cran()
  ## The defect, stated as the property that failed. Before the fix this
  ## objective read -5188 at sigma_v = 1e-8 against +4398 at 1e-6: a spuriously
  ## excellent fit in a band a minimiser walks straight into.
  o <- ne_obj()
  p <- function(sv) c(sv, 0.9792, 0.8001, 0.9785, 1.2157)
  ref <- o$lf(p(1e-6))
  expect_true(is.finite(ref))
  for (sv in 10^seq(-7, -12)) {
    v <- o$lf(p(sv))
    ## Never better than the reference, and never negative: the negative
    ## summed log-likelihood of a degenerate scale cannot be a good value.
    expect_gt(v, 0)
    expect_gte(v, ref - 1e-6)
  }
})

test_that("the two algebraic forms agree across the switch", {
  ## The fix replaces the subtraction with its closed form above
  ## NE_TILT_SWITCH. If the two disagree, one of them is wrong.
  su <- 1
  for (sv in c(0.3, 0.05, 0.01)) {
    for (zt in c(20, 30, 40, 60)) {
      eps <- -(zt + sv / su) * sv
      z <- -(eps / sv) - (sv / su)
      tilt <- -log(su) - eps^2 / (2 * sv^2) + sfa:::.log_phi_tilt(z)
      anal <- -log(su) + eps / su + sv^2 / (2 * su^2) + pnorm(z, log.p = TRUE)
      expect_equal(tilt, anal, tolerance = 1e-10)
    }
  }
  ## And the switch sits where log Phi(z) is negligible, which is what makes
  ## the closed form the whole story above it.
  expect_lt(pnorm(sfa:::.SFA_CONSTANTS$NE_TILT_SWITCH, lower.tail = FALSE), 1e-190)
})

test_that("the replication that failed now reaches the right optimum", {
  skip_on_cran()
  ## Seed 11 of the convergence DGP returned log-likelihood -1381.09 with the
  ## frontier coefficients pinned on optim()'s lower bounds, where every other
  ## start reached -1228.74. The gradient there was ~1667 on beta2 and the
  ## reported convergence code was 0, so nothing about the object looked wrong.
  o <- ne_obj(seed = 11)
  expect_gt(as.numeric(logLik(o$fit)), -1229)

  ## Not merely a better number: a genuine interior stationary point.
  g <- numDeriv::grad(o$lf, as.numeric(o$fit$coefficients))
  expect_lt(max(abs(g)), 1e-2)
})

test_that("disabling nlminb no longer changes the answer", {
  skip_on_cran()
  ## nlminb was the stage that found the bad band; with the objective repaired
  ## there is nothing to find, so the two paths must agree. This is the
  ## regression guard: if the cancellation ever comes back, they will diverge
  ## again, and by a lot rather than a little.
  d <- ne_obj(seed = 11)$data
  a <- suppressWarnings(sfm(y_pcs_e ~ x1 + x2, model_name = "NE", data = d))
  b <- suppressWarnings(sfm(y_pcs_e ~ x1 + x2, model_name = "NE", data = d,
    use.nlminb = FALSE))
  expect_equal(as.numeric(logLik(a)), as.numeric(logLik(b)), tolerance = 1e-4)
  expect_equal(unname(coef(a)), unname(coef(b)), tolerance = 1e-3)
})
