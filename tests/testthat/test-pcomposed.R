## pcomposed()/dcomposed(): the CDF of the composed error (L5).
##
## This file is unusually well supplied with external standards, so almost
## nothing here is a self-comparison:
##
##   - lambda = 1 has an exact closed form, P(Q) = Phi(Q/(sqrt(2) sigma_u))^2,
##     derived by Amsler, Schmidt and Tsay (2019) for exactly this purpose.
##   - lambda = 0 collapses to the normal.
##   - the CDF must be the integral of dcomposed(), which is closed form.
##   - the paper prints a table of values, so the implementation can be checked
##     against numbers someone else published.

test_that("dcomposed is a density that integrates to one", {
  for (p in list(c(1, 0.4), c(0.3, 1.2), c(2, 2))) {
    su <- p[1]; sv <- p[2]
    for (inef in c(TRUE, FALSE)) {
      tot <- stats::integrate(
        function(x) dcomposed(x, su, sv, inefdec = inef),
        -Inf, Inf, rel.tol = 1e-12
      )$value
      expect_equal(tot, 1, tolerance = 1e-8)
    }
  }
  ## Production and cost are mirror images, which is the only difference
  ## between them.
  x <- seq(-4, 4, by = 0.5)
  expect_equal(
    dcomposed(x, 1, 0.5, inefdec = TRUE),
    dcomposed(-x, 1, 0.5, inefdec = FALSE)
  )
  expect_equal(
    dcomposed(x, 1, 0.5, log = TRUE),
    log(dcomposed(x, 1, 0.5))
  )
})

test_that("the quadrature reproduces the exact lambda = 1 closed form", {
  ## THE test in this file. Amsler et al. derive this case precisely so that
  ## numerical methods have a standard to be judged against, and it is the only
  ## non-trivial one that exists, so this is a comparison against an external
  ## truth rather than against another part of the package.
  su <- sv <- sqrt(0.5)              # lambda = 1, sigma = 1
  for (Q in c(-16, -12, -10, -6, -2, 0, 1)) {
    exact <- stats::pnorm(Q / (sqrt(2) * su))^2
    quad <- exp(pcomposed(Q, su, sv, inefdec = FALSE, log.p = TRUE))
    expect_equal(quad, exact, tolerance = 1e-9)
  }
  ## The log and non-log paths agree. (There is deliberately no lambda = 1
  ## fast path in the implementation -- see R/pcomposed.R for why. The closed
  ## form lives here, as the oracle, and nowhere else.)
  expect_equal(
    pcomposed(-3, su, sv, inefdec = FALSE),
    exp(pcomposed(-3, su, sv, inefdec = FALSE, log.p = TRUE))
  )
})

test_that("published table values are reproduced", {
  ## Amsler, Schmidt and Tsay (2019), Table 1: F(Q) with sigma = 1, cost
  ## frontier. Their figures are simulated at R = 10,000,000 and so carry
  ## Monte Carlo error of their own -- which is why the tolerance is 1% and
  ## not machine precision. Where an exact value exists (lambda = 1) this
  ## implementation is nearer to it than the published figure is, so a
  ## tighter tolerance here would be asserting that their simulation was
  ## exact.
  sc <- function(lam) {
    sv <- sqrt(1 / (1 + lam^2))
    c(su = lam * sv, sv = sv)
  }
  tab <- list(
    c(-16, 0, 6.38875e-58), c(-16, 0.25, 3.79523e-62), c(-16, 0.5, 6.34097e-73),
    c(-12, 0, 1.77648e-33), c(-12, 0.25, 4.48929e-36), c(-12, 0.5, 2.79129e-42),
    c(-12, 2, 9.8920e-161),
    c(-10, 0, 7.61985e-24), c(-10, 0.5, 3.47421e-30), c(-10, 2, 8.4227e-113)
  )
  for (e in tab) {
    p <- sc(e[2])
    got <- exp(pcomposed(e[1], p[["su"]], p[["sv"]],
      inefdec = FALSE, log.p = TRUE
    ))
    expect_equal(got, e[3], tolerance = 0.01)
  }
})

test_that("the CDF is the integral of the density", {
  su <- 1.2; sv <- 0.4
  q <- c(-5, -2, -0.5, 0, 0.5, 2, 5)
  for (inef in c(TRUE, FALSE)) {
    num <- vapply(q, function(z) {
      stats::integrate(function(t) dcomposed(t, su, sv, inefdec = inef),
        -Inf, z, rel.tol = 1e-12
      )$value
    }, numeric(1))
    expect_equal(pcomposed(q, su, sv, inefdec = inef), num, tolerance = 1e-8)
  }
})

test_that("it behaves like a distribution function", {
  su <- 1.2; sv <- 0.4
  q <- c(-6, -3, -1, 0, 1, 3, 6)
  lo <- pcomposed(q, su, sv)
  hi <- pcomposed(q, su, sv, lower.tail = FALSE)
  expect_true(all(lo >= 0 & lo <= 1))
  expect_true(all(diff(lo) > 0))
  expect_equal(lo + hi, rep(1, length(q)))
  expect_equal(log(lo), pcomposed(q, su, sv, log.p = TRUE))
  expect_length(pcomposed(q, su, sv), length(q))
})

test_that("the upper tail does not saturate, which is why L5 exists", {
  ## A Gaussian copula on the composed error needs qnorm(F(eps)). Forming the
  ## far tail as 1 - F returns exactly 1 in double precision and qnorm(1) is
  ## +Inf, taking the copula density with it. Computing the tail directly is
  ## the fix, and this is the assertion that it worked.
  su <- 1; sv <- 1
  far <- 30
  direct <- pcomposed(far, su, sv, lower.tail = FALSE, log.p = TRUE)
  expect_true(is.finite(direct))
  expect_lt(direct, log(1e-100))
  ## The complement really would have failed.
  expect_equal(1 - pcomposed(far, su, sv), 0)
  expect_true(is.finite(stats::qnorm(direct, log.p = TRUE)))
})

test_that("the boundary and the normal limit are exact", {
  q <- c(-3, -1, 0, 2)
  ## sigma_u = 0 is the no-inefficiency boundary that sfm() actually reports.
  expect_equal(pcomposed(q, 0, 0.7), stats::pnorm(q / 0.7))
  expect_equal(
    pcomposed(q, 0, 0.7, lower.tail = FALSE),
    stats::pnorm(q / 0.7, lower.tail = FALSE)
  )
  ## And approaching it from above agrees.
  expect_equal(pcomposed(q, 1e-8, 0.7), stats::pnorm(q / 0.7), tolerance = 1e-6)
})

test_that("simulation and quadrature agree to simulation error", {
  skip_on_cran()
  su <- 1.2; sv <- 0.4
  q <- c(-2, -0.5, 0, 1)
  a <- pcomposed(q, su, sv)
  b <- pcomposed(q, su, sv, method = "simulate", R = 2e6, seed = 1)
  expect_equal(a, b, tolerance = 5e-3)
  ## Seeded draws are reproducible, and must not disturb the caller's stream.
  set.seed(99)
  before <- stats::runif(1)
  set.seed(99)
  invisible(pcomposed(0, su, sv, method = "simulate", R = 1e4, seed = 5))
  expect_equal(stats::runif(1), before)
})

test_that("more nodes do not move the answer", {
  ## If the default were short of convergence this would drift.
  su <- 1.3; sv <- 0.5
  q <- c(-4, -1, 0.5, 3)
  expect_equal(
    pcomposed(q, su, sv, n_nodes = 128),
    pcomposed(q, su, sv, n_nodes = 512),
    tolerance = 1e-10
  )
})

test_that("the refusals name what is wrong", {
  expect_error(pcomposed(0, -1, 1), "sigma_u")
  expect_error(pcomposed(0, 1, 0), "sigma_v")
  expect_error(pcomposed(0, 1, 1, n_nodes = 4), "n_nodes")
  expect_error(pcomposed(0, 1, 1, method = "simulate", R = 10), "R")
  expect_error(pcomposed("a", 1, 1), "numeric")
  expect_error(dcomposed(0, 1, -1), "sigma_v")
})
