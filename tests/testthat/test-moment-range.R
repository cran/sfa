## moment_range(): the Papadopoulos-Parmeter (2021) Type II failure check (L19).
##
## The ranges are derived from PP2021 equations (5) and (10), so the checks
## here are load-bearing: the paper prints only a handful of numbers, and every
## one of them appears below. Where the paper prints nothing, the closed forms
## are checked against brute force or against simulation.

test_that("the truncated-normal zeta formulas reproduce PP2021's table", {
  ## PP2021 section 2.1 and 2.2 tabulate gamma1(u) and gamma2(u) against
  ## z = mu/sigma_u. z = -30 is deliberately excluded: see the comment in
  ## moment_range.R -- zeta4 loses precision there and quadrature, the zeta
  ## recursion and the paper's own table give three different answers.
  z <- c(-10, -5, 0, 1, 2, 3, 5)
  m <- sfa:::.mr_tn_moments(z)
  expect_equal(round(m$g1, 3), c(1.946, 1.831, 0.995, 0.592, 0.221, 0.036, 0.000))
  expect_equal(round(m$g2, 3), c(5.580, 4.759, 0.869, 0.001, -0.239, -0.083, 0.000))
})

test_that("z = 0 is the half-normal, in closed form", {
  m <- sfa:::.mr_tn_moments(0)
  expect_equal(m$g1, sqrt(2) * (4 - pi) / (pi - 2)^(3 / 2), tolerance = 1e-10)
  expect_equal(m$g2, 8 * (pi - 3) / (pi - 2)^2, tolerance = 1e-10)
})

test_that("the truncated-normal moments agree with direct integration", {
  ## The zeta recursion is four hand-differentiated lines. Nothing else in the
  ## package would catch a slip in the third or fourth of them.
  skip_on_cran()
  tn <- function(z) {
    d <- function(u) stats::dnorm(u - z) / stats::pnorm(z)
    m1 <- stats::integrate(function(u) u * d(u), 0, Inf, rel.tol = 1e-12)$value
    ck <- function(k) {
      stats::integrate(function(u) (u - m1)^k * d(u), 0, Inf, rel.tol = 1e-12)$value
    }
    m2 <- ck(2)
    c(ck(3) / m2^1.5, ck(4) / m2^2 - 3)
  }
  for (z in c(-6, -2, 0, 1.5, 4)) {
    got <- sfa:::.mr_tn_moments(z)
    expect_equal(c(got$g1, got$g2), tn(z), tolerance = 1e-6)
  }
})

test_that("the truncated normal tends to the exponential in the left tail", {
  m <- sfa:::.mr_tn_moments(-15)
  expect_lt(abs(m$g1 - 2), 0.03)
  expect_lt(abs(m$g2 - 6), 0.25)
})

test_that("the generalized exponential constants come out of the polygammas", {
  ## spec_test() carries these rounded to 1.61 and 4.08. The cumulant
  ## generating function log Gamma(3) + log Gamma(1-t) - log Gamma(3-t) gives
  ## them exactly, and the two must agree or one of the two files is wrong.
  expect_equal(unname(sfa:::.MR_GENEXP[["g1"]]), 1.61, tolerance = 1e-3)
  expect_equal(unname(sfa:::.MR_GENEXP[["g2"]]), 4.08, tolerance = 1e-3)
  expect_equal(unname(sfa:::.MR_GENEXP[["g1"]]),
    unname(sfa:::.SPEC_INEFF$genexponential[["g1"]]),
    tolerance = 1e-3
  )
  expect_equal(unname(sfa:::.MR_GENEXP[["g2"]]),
    unname(sfa:::.SPEC_INEFF$genexponential[["g2"]]),
    tolerance = 1e-3
  )
})

test_that("the generalized exponential constants match a direct sample", {
  skip_on_cran()
  set.seed(3)
  u <- -log(1 - sqrt(stats::runif(2e6)))
  uc <- u - mean(u)
  m2 <- mean(uc^2)
  expect_equal(mean(uc^3) / m2^1.5, unname(sfa:::.MR_GENEXP[["g1"]]), tolerance = 0.01)
  expect_equal(mean(uc^4) / m2^2 - 3, unname(sfa:::.MR_GENEXP[["g2"]]), tolerance = 0.05)
})

test_that("the kurtosis range agrees with brute force over R", {
  ## .mr_g2_range() solves for the stationary point analytically. A grid does
  ## the same job slowly, and disagreeing with it means the case analysis for
  ## the sign of gamma2(v) + gamma2(u) is wrong.
  R <- seq(0, 1, length.out = 20001)
  for (g2v in c(-6 / 5, 0, 6 / 5, 3)) {
    for (g2u in c(0.8692, 4.08, 6, -0.2429)) {
      got <- sfa:::.mr_g2_range(g2v, g2u, g2u)
      q <- g2v * (1 - R)^2 + g2u * R^2
      expect_equal(got, c(min(q), max(q)), tolerance = 1e-6)
    }
  }
  ## And with gamma2(u) itself free, as it is for the truncated normal.
  got <- sfa:::.mr_g2_range(3, -0.2429, 6)
  g <- outer(R, seq(-0.2429, 6, length.out = 501),
    function(r, gu) 3 * (1 - r)^2 + gu * r^2
  )
  expect_equal(got, c(min(g), max(g)), tolerance = 1e-4)
})

test_that("the tabulated composed-error moments in PP2021 come back out", {
  ## The header rows of Tables 1, 2, 4 and 5: gamma1(v-u) and gamma2(v-u) at
  ## SNR = 0.1 ... 3 with s_v = 1. These pin equations (5) and (10) themselves.
  snr <- c(0.1, 0.5, 1, 1.5, 2, 2.5, 3)
  R <- snr^2 / (1 + snr^2)
  hn <- sfa:::.mr_ineff_set("halfnormal")
  ex <- sfa:::.mr_ineff_set("exponential")
  expect_equal(round(-hn$g1max * R^1.5, 3),
    c(-0.001, -0.089, -0.352, -0.573, -0.712, -0.797, -0.850))
  expect_equal(round(hn$g2max * R^2, 3),
    c(0.000, 0.035, 0.217, 0.417, 0.556, 0.646, 0.704))
  expect_equal(round(-ex$g1max * R^1.5, 3),
    c(-0.002, -0.179, -0.707, -1.152, -1.431, -1.601, -1.708))
  expect_equal(round(ex$g2max * R^2, 3),
    c(0.001, 0.240, 1.500, 2.876, 3.840, 4.459, 4.860))
})

test_that("the ranges are the ones PP2021 states in words", {
  r <- moment_range(noise = c("normal", "uniform", "laplace"),
    inefficiency = c("halfnormal", "exponential", "truncatednormal", "gamma"))
  pick <- function(v, u) r[r$noise == v & r$inefficiency == u, ]

  ## "for the Normal-Half-Normal distributional pair we have that
  ##  gamma1(v-u) in (-1, 0)" and "bounded below by zero and above by
  ##  gamma2(u)", which for the half-normal is 0.869.
  x <- pick("normal", "halfnormal")
  expect_equal(x$g1_lo, -0.9953, tolerance = 1e-3)
  expect_equal(x$g1_hi, 0)
  expect_equal(x$g2_lo, 0)
  expect_equal(x$g2_hi, 0.8692, tolerance = 1e-3)

  ## "For the Exponential distribution ... gamma1(v-u) in (-2, 0)" and
  ## excess kurtosis up to 6.
  x <- pick("normal", "exponential")
  expect_equal(c(x$g1_lo, x$g1_hi, x$g2_lo, x$g2_hi), c(-2, 0, 0, 6))

  ## "as in the Exponential setting, gamma1(v-u) in (-2, 0) when inefficiency
  ##  is assumed to be distributed Truncated Normal", and the truncated normal
  ##  is the one u whose composed error can be platykurtic under normal noise.
  x <- pick("normal", "truncatednormal")
  expect_equal(c(x$g1_lo, x$g1_hi), c(-2, 0))
  expect_lt(x$g2_lo, 0)
  expect_equal(x$g2_hi, 6)

  ## The gamma inefficiency term can never suffer a Type II failure.
  x <- pick("normal", "gamma")
  expect_true(is.infinite(x$g1_lo) && is.infinite(x$g2_hi))

  ## Platykurtic noise drags the lower bound down to gamma2(v) itself.
  expect_equal(pick("uniform", "halfnormal")$g2_lo, -6 / 5)

  ## Leptokurtic noise: PP2021 section 2.2.2, the composed error DIPS below
  ## min(gamma2(v), gamma2(u)) before rising, so the floor is neither endpoint.
  x <- pick("laplace", "halfnormal")
  expect_lt(x$g2_lo, min(3, 0.8692))
  expect_equal(x$g2_lo, 3 * 0.8692 / (3 + 0.8692), tolerance = 1e-3)
  expect_equal(x$g2_hi, 3)
})

test_that("Student t noise enters through 6/(df - 4)", {
  r <- moment_range(noise = "t", inefficiency = "halfnormal", df = 10)
  expect_equal(r$g2_hi, 1)          # gamma2(v) = 6/6 = 1 > gamma2(u), so it is the max
  expect_error(moment_range(noise = "t", df = 4), "greater than 4")
  expect_error(moment_range(noise = "t", df = 3), "greater than 4")
  ## Large df is the normal case.
  a <- moment_range(noise = "t", inefficiency = "exponential", df = 1e6)
  b <- moment_range(noise = "normal", inefficiency = "exponential")
  expect_equal(a$g2_lo, b$g2_lo, tolerance = 1e-4)
  expect_equal(a$g2_hi, b$g2_hi, tolerance = 1e-4)
})

test_that("orientation flips the skewness interval and nothing else", {
  p <- moment_range(noise = "normal", inefficiency = "exponential")
  c_ <- moment_range(noise = "normal", inefficiency = "exponential",
    inefdec = "cost function")
  expect_equal(c(p$g1_lo, p$g1_hi), c(-2, 0))
  expect_equal(c(c_$g1_lo, c_$g1_hi), c(0, 2))
  expect_equal(c(p$g2_lo, p$g2_hi), c(c_$g2_lo, c_$g2_hi))
})

test_that("correctly specified data is inside the range, wrong u is outside", {
  skip_on_cran()
  set.seed(21)
  n <- 4000
  x <- stats::rnorm(n)

  ## Half-normal truth at a moderate SNR: every pair whose u can reach this
  ## much skewness should pass, the half-normal included.
  e <- stats::rnorm(n) - abs(stats::rnorm(n, 0, 1.2))
  r <- moment_range(e - mean(e), noise = "normal")
  expect_equal(r$verdict[r$inefficiency == "halfnormal"], "ok")

  ## Exponential truth at a high SNR: skewness beyond -1 refutes the
  ## half-normal outright and leaves the exponential standing.
  e <- stats::rnorm(n) - stats::rexp(n, 1 / 2.5)
  r <- moment_range(e - mean(e), noise = "normal")
  expect_true(grepl("type II", r$verdict[r$inefficiency == "halfnormal"]))
  expect_equal(r$verdict[r$inefficiency == "exponential"], "ok")
  ## gamma is unbounded, so it survives whatever the data does.
  expect_equal(r$verdict[r$inefficiency == "gamma"], "ok")
})

test_that("positively skewed residuals are a Type I failure for every pair", {
  set.seed(5)
  e <- stats::rnorm(500) + stats::rexp(500, 1 / 2)
  r <- moment_range(e, noise = "normal")
  expect_true(all(r$verdict == "type I (wrong skew)"))
  ## And under a cost reading the same residuals are fine.
  r2 <- moment_range(e, noise = "normal", inefdec = "cost function")
  expect_false(any(r2$verdict == "type I (wrong skew)"))
})

test_that("with no data it is a reference table", {
  r <- moment_range()
  expect_s3_class(r, "sfa_moment_range")
  expect_equal(nrow(r), 25L)
  expect_true(all(is.na(r$verdict)))
  expect_equal(attr(r, "nobs"), 0L)
  expect_output(print(r), "Attainable skewness")
})

test_that("a fitted sfm() object carries its own residuals and orientation", {
  skip_on_cran()
  d <- data_gen_cs(N = 400, rand = 4, sig_u = 1, sig_v = 0.5, cons = 1,
    beta1 = 0.5, beta2 = 0.5, a = 1, mu = 0.5)
  f <- sfm(y_pcs ~ x1 + x2, data = d, model_name = "NHN")
  r <- moment_range(f)
  expect_s3_class(r, "sfa_moment_range")
  expect_equal(attr(r, "nobs"), 400L)
  expect_true(attr(r, "production"))
  expect_true(is.finite(attr(r, "skewness")))
  expect_output(print(r), "sample skewness")
  ## The data really is normal-half-normal, so that pair must survive.
  expect_equal(r$verdict[r$noise == "normal" & r$inefficiency == "halfnormal"], "ok")
})

test_that("it refuses input it cannot use", {
  expect_error(moment_range(stats::rnorm(5)), "too few")
  expect_error(moment_range(list(a = 1)), "sfareg")
  expect_error(moment_range(noise = "cauchy"), "should be one of")
})
