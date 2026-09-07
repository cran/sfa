## efficiency(): the three point predictors, on either scale of y.

## The data frame is kept as a NAMED object rather than built inline in the
## call. fitted.sfareg() re-evaluates the call's `data` argument, so an inline
## data.frame() built from variables local to a helper is gone by the time
## fitted() runs and the level-scale branch cannot form a frontier. Naming it
## is also how a user actually works.
eff_data <- function(seed = 2, n = 400, intercept = 6) {
  set.seed(seed)
  x1 <- rnorm(n); x2 <- rnorm(n)
  ## A large intercept keeps the fitted frontier well away from zero, so the
  ## level-scale branch is exercised where it is actually defined.
  y <- intercept + 0.5 * x1 + 0.5 * x2 + rnorm(n, 0, 0.3) - abs(rnorm(n, 0, 0.8))
  data.frame(y = y, x1 = x1, x2 = x2)
}

test_that("the three predictors are ordered the way the theory says", {
  d <- eff_data()
  f <- sfm(y ~ x1 + x2, model_name = "NHN", data = d)
  bc <- efficiency(f, "bc")
  jl <- efficiency(f, "jlms")
  md <- efficiency(f, "mode")

  expect_length(bc, nobs(f))
  expect_true(all(bc > 0 & bc <= 1))
  expect_true(all(jl > 0 & jl <= 1))

  ## bc reports the STORED number rather than a round trip through log/exp.
  expect_equal(bc, as.numeric(f$exp_u_hat))

  ## E[exp(-u)] >= exp(-E[u]) by Jensen, so bc is never below jlms.
  expect_true(all(jl <= bc + 1e-12))

  ## The mode is the only predictor that reaches exactly 1: u = 0 wherever the
  ## posterior mean is negative, which is most of the efficient tail.
  expect_true(any(md >= 1 - 1e-12))
  expect_true(all(md <= 1 + 1e-12))
  expect_false(any(bc >= 1 - 1e-12))
})

test_that("logDepVar = FALSE gives the level-scale ratio", {
  d <- eff_data()
  f <- sfm(y ~ x1 + x2, model_name = "NHN", data = d)
  lg <- efficiency(f, "jlms", logDepVar = TRUE)
  lv <- efficiency(f, "jlms", logDepVar = FALSE, newdata = d)

  expect_length(lv, nobs(f))
  expect_false(isTRUE(all.equal(lg, lv)))

  ## TE = 1 - u/f exactly.
  u <- -log(lg)
  expect_equal(lv, 1 - u / as.numeric(sfa:::.sfa_xb(f, d)$xb))

  ## On a frontier bounded well away from zero the scores stay in (0, 1].
  expect_true(all(lv > 0 & lv <= 1))
})

test_that("a frontier that crosses zero warns rather than returning nonsense", {
  ## 1 - u/f is a ratio of estimates and goes negative where the fitted
  ## frontier does. Silence there would be the bug.
  d <- eff_data(intercept = 0)
  f <- sfm(y ~ x1 + x2, model_name = "NHN", data = d)
  expect_warning(efficiency(f, "jlms", logDepVar = FALSE, newdata = d),
    "negative or non-finite"
  )
})

test_that("efficiency() refuses what it cannot compute", {
  d <- eff_data()
  f <- sfm(y ~ x1 + x2, model_name = "NHN", data = d)
  expect_error(efficiency(list()), "must be an \"sfareg\" fit")
  expect_error(efficiency(f, "nonsense"))
  expect_error(efficiency(f, logDepVar = NA), "must be TRUE or FALSE")

  ## The mode needs a truncated-normal posterior, which only some models have.
  set.seed(5); n <- 300
  x1 <- rnorm(n); x2 <- rnorm(n)
  y <- 2 + 0.5 * x1 + 0.5 * x2 + rnorm(n, 0, 0.3) - rexp(n)
  dne <- data.frame(y = y, x1 = x1, x2 = x2)
  g <- sfm(y ~ x1 + x2, model_name = "NE", data = dne)
  if (is.null(g$u_posterior)) {
    expect_error(efficiency(g, "mode"), "posterior")
  } else {
    expect_length(efficiency(g, "mode"), nobs(g))
  }
})
