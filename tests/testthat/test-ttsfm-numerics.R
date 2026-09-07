## ttsfm() numerical guards, from a code review of the TTHN branch (2026-09-07).
##
## The review's own finding was that the DENSITY is right -- it was checked
## against direct numerical convolution of (v + w - u) and agrees. Everything
## below is about the arithmetic around it: the scale the optimiser starts on,
## what happens to a log-density that cannot be computed, and what the third
## optimiser stage is allowed to return.

test_that(".tt_objective() refuses a draw rather than inventing a number", {
  expect_equal(sfa:::.tt_objective(c(-1, -2, -3)), 6)
  ## Every non-finite route to the same answer: the finite barrier, and
  ## nothing else. What this replaced substituted
  ## -sqrt(.Machine$double.xmax / n) per bad entry -- about -9.5e152 at
  ## n = 200, i.e. 140 orders of magnitude past the barrier it sat beside.
  for (bad in list(c(-1, -Inf), c(-1, Inf), c(-1, NaN), c(-1, NA_real_))) {
    expect_identical(sfa:::.tt_objective(bad), sfa:::.TT_PENALTY)
  }
  expect_identical(sfa:::.tt_objective(NULL), sfa:::.TT_PENALTY)
  expect_identical(sfa:::.tt_objective(numeric(0)), sfa:::.TT_PENALTY)
  ## +Inf used to become a large NEGATIVE log-density, so an unusable draw
  ## could be REWARDED. It cannot now.
  expect_gt(sfa:::.tt_objective(c(-1, Inf)), sfa:::.tt_objective(c(-1, -2)))
})

test_that("the barrier is the finite one, not a near-overflow substitute", {
  expect_equal(sfa:::.TT_PENALTY, 1e12)
  expect_lt(sfa:::.TT_PENALTY, sqrt(.Machine$double.xmax / 200))
  ## rho must stay strictly inside (-1, 1): pmnorm() is handed a singular
  ## varcov at exactly +-1, which sigma_v small enough does reach.
  expect_lt(sfa:::.TT_RHO_MAX, 1)
  expect_gt(sfa:::.TT_RHO_MAX, 1 - 1e-6)
})

test_that("the sigma_v start is on the scale the likelihood reads it on", {
  ## Every branch forms sigv as exp(p[nr + 1]), so the slot holds log sigma_v.
  ## The literal 0.2 came from the reference implementation in
  ## `base code/ttsfm/2TierR.Rnw`, which parameterises sigma_v DIRECTLY -- so
  ## the optimiser used to start at exp(0.2) = 1.22.
  skip_on_cran()
  d <- data_gen_cs(N = 120, rand = 3, sig_u = 1, sig_v = 0.3, cons = 1,
    beta1 = 0.5, beta2 = 0.5, a = 1, mu = 0.5)
  f <- suppressWarnings(ttsfm(y_ttne ~ x1 + x2, model_name = "TTNE", data = d,
    optHessian = FALSE, PSopt = FALSE, maxit.bobyqa = 2))
  ## start_v is returned as handed to stage 1; slot n_x + 1 is sigma_v's.
  expect_equal(exp(f$start_v[4]), 0.2, tolerance = 1e-8)
  expect_lt(f$start_v[4], 0)
})

test_that("stage 3 is allowed to escape a spike at the point it started from", {
  skip_on_cran()
  ## Regression test for a fix that was tried and REVERTED, and for the reason.
  ##
  ## The obvious guard -- stages 1 and 2 refuse to move backwards, so make
  ## stage 3 do the same -- is wrong, because the point stage 3 starts from is
  ## not always trustworthy. On exactly this data the NR likelihood has a
  ## spike: the stage-2 point evaluates to an objective of -3.6e17, optim()
  ## correctly escapes it to 541.8, and such a guard drags the fit back onto
  ## the spike. sfma() then saw logLik = +3.6e17 for NR and gave it all the
  ## weight, over an NE fit that was right.
  ##
  ## L-BFGS-B is a descent method. When it ends above its own start, that is
  ## evidence about the START, not about the search.
  set.seed(3)
  n <- 400
  x1 <- stats::rnorm(n); x2 <- stats::rnorm(n)
  d <- data.frame(y = 1 + 0.5 * x1 + 0.5 * x2 + stats::rnorm(n, 0, 0.4) -
    stats::rexp(n, 1), x1 = x1, x2 = x2)
  f <- suppressWarnings(sfm(y ~ x1 + x2, model_name = "NR", data = d))
  ll <- as.numeric(logLik(f))
  expect_true(is.finite(ll))
  ## A log-likelihood for 400 observations belongs in the hundreds. The bug
  ## produced +3.6e17.
  expect_lt(ll, 0)
  expect_gt(ll, -1e4)
})

test_that("a non-zero convergence code with an IMPROVED value is kept", {
  ## L-BFGS-B reports ABNORMAL_TERMINATION_IN_LNSRCH (code 52) whenever its
  ## line search meets a discontinuity, and it does so having improved the
  ## objective. ttsfm() used to stop() on any non-zero code, which threw away
  ## the whole fit in exactly that case.
  p0 <- c(1, 1)
  fn <- function(p) { v <- sum(p^2); if (v < 1.5) v + 40 else v }
  res <- suppressWarnings(opt.optim(
    fn = fn, start_v = p0, lower.optim = c(-5, -5), upper.optim = c(5, 5),
    maxit.optim = 100, opt.TF = TRUE, method = "L-BFGS-B",
    optHessian = TRUE, trace = 0, verbose = FALSE
  ))
  expect_identical(res$opt$convergence, 52L)
  expect_lt(res$opt$value, fn(p0))
})
