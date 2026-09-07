## Hafner, Manner and Simar (2018) extended stochastic frontier.

.hms_gen <- function(n, gam, sv, seed, dist = "halfnormal") {
  set.seed(seed)
  lx1 <- rnorm(n, 1.5, 0.3); lx2 <- rnorm(n, 1.8, 0.3)
  u <- if (dist == "halfnormal") abs(rnorm(n, 0, gam)) else rexp(n, rate = 1 / gam)
  data.frame(y = 0.9 + 0.6 * lx1 + 0.5 * lx2 + rnorm(n, 0, sv) - u, lx1, lx2)
}

test_that("the fixed constants are the paper's", {
  ## a0 solves E[u] under the mirrored density = E[u] under the original.
  f_hn <- function(a) a + (dnorm(a) - dnorm(0)) / (pnorm(a) - pnorm(0)) - sqrt(2 / pi)
  expect_equal(.HMS$halfnormal$a0, uniroot(f_hn, c(0.5, 3), tol = 1e-12)$root,
    tolerance = 1e-9
  )
  expect_equal(.HMS$halfnormal$a0, 1.3892032925, tolerance = 1e-9)  # paper
  ## exponential: a0 = 2(1 - exp(-a0))
  expect_equal(.HMS$exponential$a0,
    uniroot(function(a) 2 * (1 - exp(-a)) - a, c(0.5, 3), tol = 1e-12)$root,
    tolerance = 1e-9
  )
  expect_equal(.HMS$halfnormal$k1, sqrt(2 / pi))
  expect_equal(.HMS$exponential$k1, 1)
})

test_that("every branch of the density integrates to 1 with the right mean", {
  for (dist in c("halfnormal", "exponential")) {
    k1 <- .HMS[[dist]]$k1
    for (g in c(0.8, -0.8, 0.3, -0.3)) {
      f <- function(w) exp(.esfm_logdens(w, g, 0.5, dist))
      I <- integrate(f, -30, 30, rel.tol = 1e-10, subdivisions = 2000)$value
      M <- integrate(function(w) w * f(w), -30, 30, rel.tol = 1e-10,
        subdivisions = 2000)$value
      expect_equal(I, 1, tolerance = 1e-7, info = paste(dist, g))
      ## E[w] = -k1|gamma| whichever sign gamma has: that is what fixing a0 buys
      expect_equal(M, -k1 * abs(g), tolerance = 1e-6, info = paste(dist, g))
    }
    ## and the SIGN of the skewness flips with the sign of gamma
    m3 <- function(g) {
      f <- function(w) exp(.esfm_logdens(w, g, 0.5, dist))
      M <- integrate(function(w) w * f(w), -30, 30, rel.tol = 1e-10)$value
      integrate(function(w) (w - M)^3 * f(w), -30, 30, rel.tol = 1e-9)$value
    }
    expect_lt(m3(0.8), 0)   # classical: negative skew
    expect_gt(m3(-0.8), 0)  # extended: positive skew, the wrong-skew case
  }
})

test_that("gamma > 0 reproduces the package's own classical densities", {
  w <- c(-2, -0.5, 0, 0.5, 2)
  g <- 0.8; sv <- 0.5
  s <- sqrt(g^2 + sv^2)
  nhn <- log(2 / s) + dnorm(w / s, log = TRUE) + pnorm(-(w / s) * (g / sv), log.p = TRUE)
  expect_equal(.esfm_logdens(w, g, sv, "halfnormal"), nhn, tolerance = 1e-10)
  ne <- -log(g) + w / g + sv^2 / (2 * g^2) + pnorm(-w / sv - sv / g, log.p = TRUE)
  expect_equal(.esfm_logdens(w, g, sv, "exponential"), ne, tolerance = 1e-9)
})

test_that("the density is continuous at gamma = 0", {
  for (dist in c("halfnormal", "exponential")) {
    w <- c(-1, -0.2, 0.4, 1.5)
    lim <- dnorm(w, 0, 0.5, log = TRUE)
    expect_equal(.esfm_logdens(w, 0, 0.5, dist), lim)
    expect_equal(.esfm_logdens(w, 1e-6, 0.5, dist), lim, tolerance = 1e-4, info = dist)
    expect_equal(.esfm_logdens(w, -1e-6, 0.5, dist), lim, tolerance = 1e-4, info = dist)
  }
})

test_that("esfm() recovers gamma on a correctly skewed sample", {
  skip_on_cran()
  f <- esfm(y ~ lx1 + lx2, data = .hms_gen(500, 0.5, 0.25, 3))
  expect_s3_class(f, "esfm")
  expect_equal(unname(f$gamma), 0.5, tolerance = 0.12)
  expect_equal(unname(f$sigma_v), 0.25, tolerance = 0.08)
  expect_equal(unname(f$out["lx1", "par"]), 0.6, tolerance = 0.12)
  ## and it agrees with the classical fit, which it nests at gamma > 0
  cf <- sfm(y ~ lx1 + lx2, model_name = "NHN", data = .hms_gen(500, 0.5, 0.25, 3))
  lam <- cf$out["lambda", "par"]; sg <- cf$out["sigma", "par"]
  expect_equal(unname(f$gamma), unname(sg * lam / sqrt(1 + lam^2)), tolerance = 0.02)
})

test_that("a wrongly skewed sample gives a negative gamma where sfm() collapses", {
  skip_on_cran()
  d <- .hms_gen(50, 0.3, 0.25, 1)
  skip_if_not(mean(residuals(lm(y ~ lx1 + lx2, d))^3) > 0)
  cf <- suppressWarnings(sfm(y ~ lx1 + lx2, model_name = "NHN", data = d))
  lam <- cf$out["lambda", "par"]; sg <- cf$out["sigma", "par"]
  expect_lt(sg * lam / sqrt(1 + lam^2), 0.01)   # classical collapses
  f <- esfm(y ~ lx1 + lx2, data = d)
  expect_lt(f$gamma, 0)                         # extended does not
  expect_true(all(is.finite(f$Efficiency)))
  expect_true(all(f$Efficiency > 0 & f$Efficiency <= 1))
  expect_lt(mean(f$Efficiency), 0.999)          # genuine inefficiency reported
})

test_that("efficiency is a proper score on both sides of zero", {
  skip_on_cran()
  for (dist in c("halfnormal", "exponential")) {
    f <- esfm(y ~ lx1 + lx2, data = .hms_gen(200, 0.5, 0.25, 8, dist), dist = dist)
    expect_length(f$Efficiency, 200L)
    expect_true(all(f$Efficiency > 0 & f$Efficiency <= 1))
  }
})

test_that("symmetry_test() is chi-square(1), not a mixture", {
  skip_on_cran()
  f <- esfm(y ~ lx1 + lx2, data = .hms_gen(500, 0.5, 0.25, 3))
  tt <- symmetry_test(f)
  expect_equal(tt$df, 1L)
  expect_equal(tt$p.value, pchisq(tt$statistic, 1, lower.tail = FALSE))
  expect_lt(tt$p.value, 0.05)   # sigma_u = 2 sigma_v is not subtle
  expect_error(symmetry_test(structure(list(), class = "lm")), "esfm")
  ## pure noise: should not reject
  set.seed(21); n <- 400
  lx1 <- rnorm(n, 1.5, 0.3); lx2 <- rnorm(n, 1.8, 0.3)
  d0 <- data.frame(y = 0.9 + 0.6 * lx1 + 0.5 * lx2 + rnorm(n, 0, 1), lx1, lx2)
  expect_gt(symmetry_test(esfm(y ~ lx1 + lx2, data = d0))$p.value, 0.05)
})
