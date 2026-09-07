## Normal / truncated skew-Laplace, sfm(model_name = "TSL").
##
## The density is a SIGNED mixture of two exponentials, so two things can go
## wrong that would not show up as an error: the composed likelihood can be
## evaluated by subtracting raw terms that overflow, and the efficiency
## predictor can be built with convex weights instead of the signed ones.
## Both are pinned below, along with the closed form itself against numerical
## convolution.

## The implied u-density, written out independently of the package.
tsl_fu <- function(u, su, lam) {
  ((1 + lam) / (su * (2 * lam + 1))) * (2 * exp(-u / su) - exp(-(1 + lam) * u / su))
}

test_that("the implied u-density is a proper density", {
  for (su in c(0.4, 0.8)) {
    for (lam in c(0.5, 1.5, 4)) {
      expect_equal(integrate(tsl_fu, 0, Inf, su = su, lam = lam)$value, 1,
                   tolerance = 1e-8)
      expect_true(all(tsl_fu(seq(0, 30, length.out = 2000), su, lam) >= -1e-12))
    }
  }
})

test_that("the composed log-density matches numerical convolution", {
  skip_on_cran()
  d <- cs_small(N = 120)
  fit <- suppressWarnings(sfm(y_pcs_tsl ~ x1 + x2, model_name = "TSL", data = d))
  p <- fit$out[, "par"]
  sv <- unname(p[["sigv"]]); su <- unname(p[["sigu"]]); lam <- unname(p[["lambda"]])

  ## Rebuild the closed form from the fitted parameters and compare against
  ## direct integration of f_u against the normal noise density.
  closed <- function(e) {
    A <- sv^2 / (2 * su^2) + e / su
    B <- (1 + lam)^2 * sv^2 / (2 * su^2) + e * (1 + lam) / su
    a <- -sv / su - e / sv
    b <- -sv * (1 + lam) / su - e / sv
    l1 <- log(2) + A + pnorm(a, log.p = TRUE)
    l2 <- B + pnorm(b, log.p = TRUE)
    log1p(lam) - log(2 * lam + 1) - log(su) + l1 + log(-expm1(pmin(l2 - l1, -1e-300)))
  }
  numeric_ll <- function(e) {
    log(integrate(function(u) tsl_fu(u, su, lam) * dnorm((e + u) / sv) / sv,
                  0, Inf, rel.tol = 1e-11, subdivisions = 2000L)$value)
  }
  for (e in c(-4, -1.5, -0.25, 0, 0.75)) {
    expect_equal(closed(e), numeric_ll(e), tolerance = 1e-6, info = paste("eps =", e))
  }
})

test_that("the log-domain form survives where subtracting raw terms does not", {
  ## Both exponentials carry sigma_v^2/(2 sigma_u^2), which overflows once
  ## sigma_u is small relative to sigma_v -- a region the optimizer visits.
  naive <- function(e, su, sv, lam) {
    A <- sv^2 / (2 * su^2) + e / su
    B <- (1 + lam)^2 * sv^2 / (2 * su^2) + e * (1 + lam) / su
    log(2 * exp(A) * pnorm(-sv / su - e / sv) - exp(B) * pnorm(-sv * (1 + lam) / su - e / sv))
  }
  stable <- function(e, su, sv, lam) {
    A <- sv^2 / (2 * su^2) + e / su
    B <- (1 + lam)^2 * sv^2 / (2 * su^2) + e * (1 + lam) / su
    l1 <- log(2) + A + pnorm(-sv / su - e / sv, log.p = TRUE)
    l2 <- B + pnorm(-sv * (1 + lam) / su - e / sv, log.p = TRUE)
    l1 + log(-expm1(pmin(l2 - l1, -1e-300)))
  }
  grid <- expand.grid(su = c(0.05, 0.1, 0.3), sv = c(0.3, 1), e = c(-2, 0, 2, 5, 20))
  st <- mapply(stable, grid$e, grid$su, grid$sv, 1.5)
  nv <- mapply(naive, grid$e, grid$su, grid$sv, 1.5)
  expect_true(all(is.finite(st)))
  ## The naive form must actually fail somewhere, or this guard is pointless.
  expect_true(any(!is.finite(nv)))
  ## Where it does survive, the two must agree.
  ok <- is.finite(nv)
  expect_equal(st[ok], nv[ok], tolerance = 1e-6)
})

test_that("TSL recovers its DGP and predicts efficiency", {
  skip_on_cran()
  d <- data_gen_cs(N = 1500, rand = 3, sig_u = 0.8, sig_v = 0.3, cons = 0.5,
                   beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.5, lam_tsl = 1.5)
  fit <- suppressWarnings(sfm(y_pcs_tsl ~ x1 + x2, model_name = "TSL", data = d))
  p <- fit$out[, "par"]
  expect_named(fit$out[, "par"], c("sigv", "sigu", "lambda",
                                   "(Intercept)", "x1", "x2"))
  expect_equal(unname(p[["sigv"]]), 0.3, tolerance = 0.12)
  expect_equal(unname(p[["sigu"]]), 0.8, tolerance = 0.15)
  expect_equal(unname(p[["x1"]]), 0.5, tolerance = 0.1)
  expect_equal(unname(p[["x2"]]), 0.5, tolerance = 0.1)

  ## Signed-mixture weights: getting these convex instead produces scores
  ## essentially uncorrelated with the truth, so a loose correlation floor is
  ## enough to catch it.
  expect_length(fit$u_hat, nrow(d))
  expect_length(fit$exp_u_hat, nrow(d))
  expect_true(all(fit$u_hat >= 0))
  expect_true(all(fit$exp_u_hat >= 0 & fit$exp_u_hat <= 1))
  expect_gt(cor(fit$u_hat, d$u_tsl), 0.8)
  expect_gt(cor(fit$exp_u_hat, exp(-d$u_tsl)), 0.8)
})

test_that("the TSL generator matches the density it claims to draw from", {
  d <- data_gen_cs(N = 60000, rand = 1, sig_u = 0.8, sig_v = 0.3, cons = 0.5,
                   beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.5, lam_tsl = 1.5)
  u <- d$u_tsl
  su <- 0.8; lam <- 1.5
  w1 <- 2 * (1 + lam) / (2 * lam + 1)
  w2 <- 1 / (2 * lam + 1)
  expect_true(all(u >= 0))
  expect_equal(mean(u), w1 * su - w2 * su / (1 + lam), tolerance = 0.02)
  expect_equal(mean(u^2), w1 * 2 * su^2 - w2 * 2 * (su / (1 + lam))^2, tolerance = 0.05)
})

test_that("the TSL start uses the exponential inversion, not the half-normal one", {
  skip_on_cran()
  ## The half-normal inversion over-attributes variance to u on skew-Laplace
  ## data, driving m2 - su^2(1-2/pi) negative and pinning sigma_v at zero --
  ## worth 215 log-likelihood points on one sample in four before this was
  ## changed. sigma_v must start strictly positive.
  d <- data_gen_cs(N = 800, rand = 3, sig_u = 0.8, sig_v = 0.3, cons = 0.5,
                   beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.5, lam_tsl = 1.5)
  sc <- sfa:::start_cs(y_pcs_tsl ~ x1 + x2, d, c("(Intercept)", "x1", "x2"),
                       1, "TSL", 3, FALSE, 0, NULL)
  expect_length(sc$start_v, 6)
  expect_gt(sc$start_v[1], 0)
  expect_gt(sc$start_v[2], 0)
  expect_equal(sc$start_v[3], 1)
})

test_that("the TSL density is the thesis's, and its moments are too", {
  ## Wang (2012), "A Normal Truncated Skewed-Laplace Model in Stochastic
  ## Frontier Analysis", WKU master's thesis, equations (3.7) and (3.8). The
  ## distribution itself is Aryal and Rao (2005); the thesis applies it to the
  ## frontier model. Checking the package against the SOURCE rather than
  ## against itself, which is what makes the citation in ?sfm a claim rather
  ## than a decoration.
  f_thesis <- function(u, phi, lam) {
    (1 + lam) / (phi * (2 * lam + 1)) *
      (2 * exp(-u / phi) - exp(-(1 + lam) * u / phi))
  }
  ## eq (3.8): E[X^k] = phi^k (1+lam) Gamma(k+1)/(2 lam+1) * {2 - (1+lam)^-(k+1)}
  mom_thesis <- function(k, phi, lam) {
    phi^k * (1 + lam) * gamma(k + 1) / (2 * lam + 1) * (2 - 1 / (1 + lam)^(k + 1))
  }
  for (lam in c(0.1, 1 / sqrt(2), 1.5, 3)) {
    for (phi in c(0.7, 1, 1.6)) {
      ## The density integrates to 1, so it is the truncated version.
      ## rel.tol tightened from integrate()'s default: over (0, Inf) the
      ## default leaves ~4e-8, which is quadrature error rather than any
      ## disagreement with the thesis, and loosening the assertion instead
      ## would weaken the check for no reason.
      expect_equal(stats::integrate(function(u) f_thesis(u, phi, lam),
        0, Inf, rel.tol = 1e-12)$value, 1, tolerance = 1e-9)
      for (k in 1:3) {
        num <- stats::integrate(function(u) u^k * f_thesis(u, phi, lam),
          0, Inf, rel.tol = 1e-12, subdivisions = 400L)$value
        expect_equal(num, mom_thesis(k, phi, lam), tolerance = 1e-7,
          info = paste("lambda", lam, "phi", phi, "k", k))
      }
    }
  }
  ## lambda = 0 must collapse to the exponential, which is what makes TSL a
  ## generalisation rather than a different family.
  expect_equal(mom_thesis(1, 1, 0), 1, tolerance = 1e-12)
})

test_that("the scale-shape multiplier peaks at 4 - 2*sqrt(2)", {
  ## The reason lambda is weakly identified, as an exact statement rather than
  ## a measured one: E[U]/sigma_u equals 1 in BOTH limits and its maximum over
  ## all lambda is 4 - 2 sqrt(2) = 1.1716, at lambda = 1/sqrt(2). The shape
  ## parameter can move the standardised mean by at most 17 percent.
  m <- function(lam) (1 + lam) / (2 * lam + 1) * (2 - 1 / (1 + lam)^2)
  expect_equal(m(1e-10), 1, tolerance = 1e-8)
  expect_equal(m(1e10), 1, tolerance = 1e-6)
  o <- stats::optimize(m, c(1e-3, 50), maximum = TRUE)
  expect_equal(o$maximum, 1 / sqrt(2), tolerance = 1e-4)
  expect_equal(o$objective, 4 - 2 * sqrt(2), tolerance = 1e-8)
  ## And the convergence DGP's lambda = 1.5 sits in that flat region.
  expect_gt(m(1.5), 0.97 * (4 - 2 * sqrt(2)))
})
