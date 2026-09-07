## Chen and Wang (2012) centered-residuals moment estimator and test.

.cw_data <- function(n, seed, dgp = c("NHN", "NE"), sv = 1, su = 1) {
  dgp <- match.arg(dgp)
  set.seed(seed)
  X <- matrix(rnorm(n * 3), n, 3)
  u <- if (dgp == "NHN") abs(rnorm(n, 0, su)) else rexp(n, rate = 1 / su)
  y <- as.numeric(X %*% c(1, 1, 1) + rnorm(n, 0, sv) - u)
  as.numeric(residuals(lm(y ~ X)))
}

test_that(".cw_theory() matches a direct simulation", {
  taus <- c(1, 1.5, 2)
  tt <- .cw_theory("NHN", c(lambda = 1, sigma = sqrt(2)), taus = taus, kmax = 3L)
  set.seed(7)
  R <- 6e5
  u <- abs(rnorm(R, 0, 1)); v <- rnorm(R, 0, 1)
  ec <- (v - u) - mean(v - u)
  expect_equal(tt$mu1, sqrt(2 / pi), tolerance = 1e-8)      # E|N(0,1)|
  expect_equal(tt$moments[1], mean(ec^2), tolerance = 0.01)
  expect_equal(tt$moments[2], mean(ec^3), tolerance = 0.02)
  for (j in seq_along(taus)) {
    expect_equal(tt$Ecos[j], mean(cos(taus[j] * ec)), tolerance = 0.005)
    expect_equal(tt$Esin[j], mean(sin(taus[j] * ec)), tolerance = 0.005)
  }
})

test_that("the moment estimator is Chen and Wang's closed form (32)", {
  ## For NHN and NE their (32) is the same moment inversion COLS uses, so the
  ## numerical solve must land on it rather than merely near it.
  for (mn in c("NHN", "NE")) {
    e <- .cw_data(4000, 11, mn)
    ec <- e - mean(e)
    m2 <- mean(ec^2); m3 <- mean(ec^3)
    k3 <- if (mn == "NHN") (1 - 4 / pi) * sqrt(2 / pi) else -2
    k2 <- if (mn == "NHN") 1 - 2 / pi else 1
    su <- (m3 / k3)^(1 / 3)
    sv <- sqrt(m2 - k2 * su^2)
    z <- cw_test(e, mn)
    expect_equal(unname(z$estimate[["sigma_u"]]), su, tolerance = 1e-5)
    expect_equal(unname(z$estimate[["sigma_v"]]), sv, tolerance = 1e-5)
  }
})

test_that("the moment estimator recovers the truth on average", {
  skip_on_cran()
  ## One draw is not the right assertion: sigma_u comes from a cube root of a
  ## third moment and is genuinely high-variance, which is what Chen and Wang's
  ## own Table 3 shows. Average over seeds instead.
  for (mn in c("NHN", "NE")) {
    est <- vapply(1:20, function(r) {
      z <- try(cw_test(.cw_data(2000, 300 + r, mn), mn), silent = TRUE)
      if (inherits(z, "try-error")) c(NA, NA) else unname(z$estimate)
    }, numeric(2))
    ok <- stats::complete.cases(t(est))
    expect_gt(sum(ok), 15)
    expect_equal(mean(est[1, ok]), 1, tolerance = 0.07)  # sigma_v
    expect_equal(mean(est[2, ok]), 1, tolerance = 0.10)  # sigma_u
  }
})

test_that("the test statistic is a proper chi-square quantity", {
  e <- .cw_data(1000, 3, "NHN")
  for (ty in c("cosine", "sine")) {
    z <- cw_test(e, "NHN", tau = 1, type = ty)
    expect_s3_class(z, "cw_test")
    expect_equal(z$df, 1L)
    expect_gte(z$statistic, 0)
    expect_equal(z$p.value, pchisq(z$statistic, 1, lower.tail = FALSE))
  }
  ## q frequencies give q degrees of freedom
  expect_equal(cw_test(e, "NHN", tau = c(1, 2))$df, 2L)
  expect_output(print(cw_test(e, "NHN")), "Chen and Wang")
})

test_that("the centering correction is not optional", {
  ## Dropping the -E[dpsi/deps]*eps_c term shrinks the variance and so inflates
  ## D. This pins that the shipped statistic is the corrected one.
  e <- .cw_data(2000, 5, "NHN")
  ec <- e - mean(e)
  z <- cw_test(e, "NHN", tau = 1)
  th <- c(z$estimate[["sigma_v"]], z$estimate[["sigma_u"]])
  tf <- function(t) {
    .cw_theory("NHN", c(lambda = t[2] / t[1], sigma = sqrt(sum(t^2))),
      taus = 1, kmax = 3L
    )
  }
  TT <- tf(th)
  psi1 <- cbind(ec^2 - TT$moments[1], ec^3 - TT$moments[2])
  psi2 <- matrix(cos(ec) - TT$Ecos, ncol = 1)
  K1 <- -numDeriv::jacobian(function(t) tf(t)$moments, th)
  K2 <- matrix(-numDeriv::jacobian(function(t) tf(t)$Ecos, th), nrow = 1)
  ## UNcorrected: no xi transformation at all
  Zbad <- psi2 - psi1 %*% t(K2 %*% solve(K1))
  Dbad <- length(ec) * mean(psi2)^2 / (crossprod(Zbad) / length(ec))[1, 1]
  expect_false(isTRUE(all.equal(as.numeric(Dbad), z$statistic, tolerance = 1e-6)))
})

test_that("wrong skew is refused rather than answered", {
  set.seed(9)
  n <- 800
  X <- matrix(rnorm(n * 3), n, 3)
  ## cost orientation read as production: third moment comes out positive
  y <- as.numeric(X %*% c(1, 1, 1) + rnorm(n, 0, 1) + abs(rnorm(n, 0, 1)))
  e <- residuals(lm(y ~ X))
  expect_error(cw_test(e, "NHN"), "WRONG WAY|wrong way")
})

test_that("unsupported models and bad input error clearly", {
  e <- .cw_data(400, 2, "NHN")
  expect_error(cw_test(e, "NTN"), "two-parameter normal-noise")
  expect_error(cw_test(e, "THT"), "two-parameter normal-noise")
  expect_error(cw_test(e), "model_name")
  expect_error(cw_test("nope", "NHN"), "numeric")
  expect_error(cw_test(e, "NHN", tau = -1), "tau")
  expect_error(cw_test(e, "NHN", tau = numeric(0)), "tau")
})

test_that("it has power against the wrong inefficiency distribution", {
  skip_on_cran()
  ## Power at tau = 1, n = 2000 is 0.95 in the study behind ?cw_test, so a
  ## single well-separated draw is a fair assertion here.
  e <- .cw_data(4000, 21, "NE")
  expect_lt(cw_test(e, "NHN", tau = 1)$p.value, 0.05)   # wrong pair
  expect_gt(cw_test(e, "NE", tau = 1)$p.value, 0.05)    # right pair
})

test_that("the Gauss-Legendre cache returns identical rules", {
  a <- .gauss_legendre_01(64L)
  b <- .gauss_legendre_01(64L)
  expect_identical(a, b)
  expect_equal(sum(a$weights), 1, tolerance = 1e-12)
  expect_equal(sum(a$weights * a$nodes), 0.5, tolerance = 1e-12)  # mean of U(0,1)
  expect_false(identical(.gauss_legendre_01(32L)$nodes, a$nodes))
})
