## npsfm(): nonparametric stochastic frontier models.
## Every test needs np, which is only in Suggests.
##
## Deliberately does NOT call library(np). npsfm() must work with np merely
## installed, not attached -- attaching it here once hid a bug where np's own
## generic gradients() was called unqualified and only resolved because the
## test harness had put np on the search path.

np_data <- function(n = 120, seed = 42, su = 0.6, sv = 0.25) {
  set.seed(seed)
  x1 <- runif(n, 1, 4)
  x2 <- runif(n, 1, 4)
  m <- 1 + 0.6 * log(x1) + 0.4 * sqrt(x2)
  list(
    d = data.frame(y = m + rnorm(n, 0, sv) - abs(rnorm(n, 0, su)), x1 = x1, x2 = x2),
    m = m, su = su, sv = sv
  )
}

test_that("npsfm() rejects bad input before doing any work", {
  g <- np_data(30)
  ## a pipe formula is an error, not silently ignored
  expect_error(npsfm(y ~ x1 | x2, data = g$d, method = "FLW"), "single-part formula")
  ## SVKZ is half-normal only
  expect_error(npsfm(y ~ x1 + x2, data = g$d, method = "SVKZ", dist = "exp"),
    "normal-half normal"
  )
  ## since 1.1.4 the message comes from .match_model_name(), not match.arg()
  expect_error(npsfm(y ~ x1 + x2, data = g$d, method = "nope"),
               "not a recognized choice")
  ## and method is case-insensitive like model_name
  expect_equal(sfa:::.match_model_name("flw", c("FLW","SVKZ"), arg = "method"), "FLW")
})

test_that(".flw_neg_cll is finite, positive-lambda only, and minimized near truth", {
  set.seed(3)
  e <- rnorm(400, 0, 0.25) - abs(rnorm(400, 0, 0.6))
  e <- e - mean(e)
  expect_equal(sfa:::.flw_neg_cll(-1, e), sfa:::.SFA_CONSTANTS$MAX_VALUE)
  expect_equal(sfa:::.flw_neg_cll(0, e), sfa:::.SFA_CONSTANTS$MAX_VALUE)
  expect_true(is.finite(sfa:::.flw_neg_cll(2.4, e)))
  best <- optimize(sfa:::.flw_neg_cll, c(1e-4, 50), e = e)$minimum
  ## true lambda is 2.4; the concentrated likelihood should land in the region
  expect_gt(best, 0.8)
  expect_lt(best, 8)
})

test_that("FLW recovers a known normal-half normal frontier", {
  skip_if_not_installed("np")
  skip_on_cran()
  g <- np_data(250, seed = 11)
  f <- npsfm(y ~ x1 + x2, data = g$d, method = "FLW", dist = "hn")

  expect_s3_class(f, "npsfareg")
  expect_equal(f$sigma.u, g$su, tolerance = 0.30)
  expect_equal(f$sigma.v, g$sv, tolerance = 0.40)
  expect_equal(f$lambda, g$su / g$sv, tolerance = 0.50)
  ## the frontier tracks the truth
  expect_gt(cor(f$frontier, g$m), 0.85)
  ## and lies above the conditional mean (production frontier)
  expect_true(all(f$frontier >= f$conditional.mean))
  ## efficiency predictions are on (0, 1]
  expect_true(all(f$exp_u_hat > 0 & f$exp_u_hat <= 1))
  expect_true(all(f$u_hat >= 0))
  expect_equal(nobs(f), nrow(g$d))
  expect_length(fitted(f), nrow(g$d))
  expect_equal(residuals(f), g$d$y - f$frontier, ignore_attr = TRUE)
})

test_that("SVKZ returns per-observation sigmas and flags wrong skew", {
  skip_if_not_installed("np")
  skip_on_cran()
  g <- np_data(250, seed = 11)
  f <- npsfm(y ~ x1 + x2, data = g$d, method = "SVKZ")

  expect_length(f$sigma.u, nrow(g$d))
  expect_length(f$sigma.v, nrow(g$d))
  expect_true(all(f$sigma.u >= 0))
  expect_true(all(f$sigma.v >= 0))
  expect_type(f$wrong.skew, "logical")
  ## wrong-skew points are floored at zero and contribute no gradient
  expect_true(all(f$sigma.u[f$wrong.skew] == 0))
  expect_true(all(f$sigma.u.grad[f$wrong.skew, ] == 0))
  expect_equal(dim(f$frontier.grad), c(nrow(g$d), 2L))
})

test_that("cost frontiers recover a cost DGP and sit below the conditional mean", {
  skip_if_not_installed("np")
  skip_on_cran()
  ## A cost frontier adds inefficiency: y = m(x) + v + u.
  set.seed(5); n <- 250
  x1 <- runif(n, 1, 4); x2 <- runif(n, 1, 4)
  m <- 1 + 0.6 * log(x1) + 0.4 * sqrt(x2)
  dcost <- data.frame(y = m + rnorm(n, 0, 0.25) + abs(rnorm(n, 0, 0.6)),
                      x1 = x1, x2 = x2)
  f <- npsfm(y ~ x1 + x2, data = dcost, method = "FLW", dist = "hn", cost = TRUE)
  expect_equal(f$sigma.u, 0.6, tolerance = 0.40)
  ## the cost frontier is the lower envelope: below the conditional mean
  expect_true(all(f$frontier <= f$conditional.mean))
  expect_gt(cor(f$frontier, m), 0.85)
})

test_that("the wrong orientation warns instead of silently returning sigma_u = 0", {
  skip_if_not_installed("np")
  skip_on_cran()
  g <- np_data(150, seed = 5)          # production data
  expect_warning(
    f <- npsfm(y ~ x1 + x2, data = g$d, method = "FLW", dist = "hn", cost = TRUE),
    "not negatively skewed"
  )
  expect_lt(f$sigma.u, 0.05)
})

test_that("supplying bandwidths skips cross-validation and changes the fit", {
  skip_if_not_installed("np")
  skip_on_cran()
  g <- np_data(120, seed = 8)
  f1 <- npsfm(y ~ x1 + x2, data = g$d, method = "FLW", dist = "hn",
              bw = c(0.5, 0.5))
  f2 <- npsfm(y ~ x1 + x2, data = g$d, method = "FLW", dist = "hn",
              bw = c(5, 5))
  expect_false(isTRUE(all.equal(f1$frontier, f2$frontier)))
  expect_equal(as.numeric(f1$bws$bw), c(0.5, 0.5))
})

test_that("the exponential branch inverts its moments", {
  skip_if_not_installed("np")
  skip_on_cran()
  set.seed(21); n <- 250
  x1 <- runif(n, 1, 4); x2 <- runif(n, 1, 4)
  m <- 1 + 0.6 * log(x1) + 0.4 * sqrt(x2)
  d <- data.frame(y = m + rnorm(n, 0, 0.25) - rexp(n, rate = 1 / 0.6), x1 = x1, x2 = x2)
  f <- npsfm(y ~ x1 + x2, data = d, method = "FLW", dist = "exp")
  expect_equal(f$sigma.u, 0.6, tolerance = 0.45)
  expect_equal(1 / f$theta, f$sigma.u, tolerance = 1e-8)
  expect_true(all(f$exp_u_hat > 0 & f$exp_u_hat <= 1))
})

test_that("PSZ returns per-observation sigmas and a convergence code", {
  skip_if_not_installed("np")
  skip_on_cran()
  g <- np_data(80, seed = 3)
  f <- npsfm(y ~ x1 + x2, data = g$d, method = "PSZ")
  expect_s3_class(f, "npsfareg")
  expect_length(f$sigma.u, nrow(g$d))
  expect_length(f$sigma.v, nrow(g$d))
  expect_true(all(f$sigma.u > 0))
  expect_length(f$convergence, nrow(g$d))
  expect_equal(dim(f$frontier.grad), c(nrow(g$d), 2L))
  ## "KPST" is an accepted alias and must give the identical fit
  f2 <- npsfm(y ~ x1 + x2, data = g$d, method = "KPST")
  expect_equal(f$frontier, f2$frontier)
  expect_equal(f2$method, "PSZ")
})

test_that("MY iterates and reports its convergence state", {
  skip_if_not_installed("np")
  skip_on_cran()
  g <- np_data(80, seed = 3)
  f <- npsfm(y ~ x1 + x2, data = g$d, method = "MY", iter = 5)
  expect_length(f$sigma.u, 1L)
  expect_true(f$iterations >= 1L && f$iterations <= 5L)
  expect_type(f$converged, "logical")
  expect_true(is.finite(f$lambda) && f$lambda > 0)
  expect_equal(f$sigma.u / f$sigma.v, f$lambda, tolerance = 1e-6)
})

test_that("local-likelihood frontiers carry no half-normal mean shift", {
  skip_if_not_installed("np")
  skip_on_cran()
  ## PSZ/MY maximize the composed-error likelihood, so their local intercept
  ## IS the frontier. Adding sqrt(2/pi)*sigma_u on top of it -- as the
  ## least-squares methods require -- biases the frontier up by about E[u].
  ## Guard against that regression: the fitted frontier must sit close to the
  ## truth, not a systematic ~E[u] above it.
  g <- np_data(150, seed = 701)
  for (meth in c("PSZ", "MY")) {
    f <- npsfm(y ~ x1 + x2, data = g$d, method = meth)
    bias <- mean(f$frontier - g$m)
    expect_lt(abs(bias), 0.35)
    ## and the conditional mean must sit BELOW the frontier by E[u]
    expect_true(all(f$conditional.mean <= f$frontier + 1e-8))
  }
})

test_that("non-half-normal distributions are refused by every method but FLW", {
  skip_if_not_installed("np")
  g <- np_data(40)
  for (meth in c("SVKZ", "PSZ", "MY", "SZ")) {
    expect_error(
      npsfm(y ~ x1 + x2, data = g$d, method = meth, dist = "exp"),
      "normal-half normal"
    )
  }
})

test_that("SZ monotonizes a prior fit and refuses cost frontiers", {
  skip_if_not_installed("np")
  skip_if_not_installed("Benchmarking")
  skip_on_cran()
  g <- np_data(60, seed = 9)
  expect_error(
    npsfm(y ~ x1 + x2, data = g$d, method = "SZ", cost = TRUE),
    "production frontiers only"
  )
  f <- npsfm(y ~ x1 + x2, data = g$d, method = "SZ")
  expect_length(f$dea.efficiency, nrow(g$d))
  ## output-oriented DEA efficiencies are >= 1, so the DEA frontier is never
  ## below the smooth fit it monotonizes
  expect_true(all(f$dea.efficiency >= 1 - 1e-8))
  expect_true(all(f$frontier >= f$prior.fit - 1e-8))
  ## a supplied prior fit of the wrong length is rejected
  expect_error(
    npsfm(y ~ x1 + x2, data = g$d, method = "SZ", prior.fit = c(1, 2, 3)),
    "one fitted value per observation"
  )
})
