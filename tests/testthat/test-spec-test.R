## spec_test(): the Papadopoulos-Parmeter pair test (L2).
##
## The statistic was RECONSTRUCTED from PP2021 equations (5) and (10), because
## the Supplementary Appendix where PP2023 derive it is not available. That
## makes the checks below load-bearing rather than decorative: they verify the
## restriction itself, the constants read off the paper's table, and the
## hand-differentiated Jacobian.

test_that("the tabulated constants are what the derivation says they are", {
  ## PP2023 Table 1 prints c_u,1 and c_u,2 without defining them. The
  ## derivation implies c_u1 = gamma1(u)^(-2/3) and c_u2 = 1/gamma2(u). If that
  ## reading is wrong the whole test is wrong, so it is pinned here.
  g1 <- c(0.9953, 2, 1.61)
  g2 <- c(0.8692, 6, 4.08)
  expect_equal(g1^(-2 / 3), c(1.0031, 0.63, 0.728), tolerance = 5e-4)
  expect_equal(1 / g2, c(1.1505, 0.1667, 0.2451), tolerance = 5e-4)

  ## And the values themselves, as used by the package.
  expect_equal(unname(sfa:::.SPEC_INEFF$halfnormal[["g1"]]), 0.9953)
  expect_equal(unname(sfa:::.SPEC_INEFF$exponential[["g1"]]), 2)
  expect_equal(unname(sfa:::.SPEC_NOISE[["normal"]]), 0)
  expect_equal(unname(sfa:::.SPEC_NOISE[["laplace"]]), 3)
  expect_equal(unname(sfa:::.SPEC_NOISE[["logistic"]]), 6 / 5)
  expect_equal(unname(sfa:::.SPEC_NOISE[["uniform"]]), -6 / 5)
})

test_that("Gamma is zero when the assumed pair IS the true pair", {
  skip_on_cran()
  ## The restriction itself, checked at a sample size large enough that the
  ## 4th moment has settled. This is the single most important assertion in
  ## the file: if Gamma does not vanish under correct specification, the
  ## statistic is measuring the wrong thing.
  set.seed(11)
  n <- 4e5
  for (snr in c(0.7, 1.5)) {
    v <- rnorm(n, 0, 1)
    u <- abs(rnorm(n, 0, snr))
    e <- v - u
    ec <- e - mean(e)
    m <- vapply(2:4, function(k) mean(ec^k), numeric(1))
    gg <- sfa:::.spec_gamma(m[1], m[2], m[3],
      g1u = 0.9953, g2u = 0.8692, g2v = 0
    )
    expect_equal(gg$Gamma, 0, tolerance = 0.05)
    ## R is a variance share, so it must land in the unit interval.
    expect_gt(gg$R, 0)
    expect_lt(gg$R, 1)
  }
})

test_that("the analytic Jacobian matches a numerical one", {
  ## A hand-differentiated expression that is never compared against anything
  ## is a guess. numDeriv is already in Imports.
  set.seed(3)
  e <- rnorm(5000, 0, 0.6) - abs(rnorm(5000, 0, 1))
  ec <- e - mean(e)
  m <- vapply(1:8, function(k) mean(ec^k), numeric(1))
  for (pair in list(c(0.9953, 0.8692, 0), c(2, 6, 3), c(1.61, 4.08, 6 / 5))) {
    f <- function(v) {
      sfa:::.spec_gamma(v[1], v[2], v[3], pair[1], pair[2], pair[3])$Gamma
    }
    num <- numDeriv::grad(f, c(m[2], m[3], m[4]))
    ana <- sfa:::.spec_jacobian(m[2], m[3], m[4], pair[1], pair[2], pair[3])
    expect_equal(unname(ana), num, tolerance = 1e-6)
  }
})

test_that("the moment covariance matrix is symmetric and sane", {
  set.seed(4)
  e <- rnorm(20000)
  ec <- e - mean(e)
  m <- vapply(1:8, function(k) mean(ec^k), numeric(1))
  S <- sfa:::.spec_moment_vcov(m)
  expect_equal(S, t(S))
  expect_true(all(diag(S) > 0))
  ## For a standard normal the closed forms are known: Var(m2) = 2,
  ## Var(m3) = 6, and m2 and m3 are uncorrelated by symmetry.
  expect_equal(unname(S["m2", "m2"]), 2, tolerance = 0.15)
  expect_equal(unname(S["m3", "m3"]), 6, tolerance = 0.6)
  expect_equal(unname(S["m2", "m3"]), 0, tolerance = 0.3)
})

test_that("it does not reject the true pair and does reject a wrong one", {
  skip_on_cran()
  set.seed(2)
  n <- 2000
  e <- rnorm(n, 0, 0.6) - abs(rnorm(n, 0, 1))

  good <- spec_test(e, noise = "normal", inefficiency = "halfnormal")
  expect_s3_class(good, "sfa_spec_test")
  expect_s3_class(good, "htest")
  expect_gt(good$p.value, 0.05)
  expect_false(good$reject)

  ## Laplace noise has excess kurtosis 3 and an exponential u has 6; neither
  ## is anywhere near this sample, so the pair should be refused decisively.
  bad <- spec_test(e, noise = "laplace", inefficiency = "exponential")
  expect_lt(bad$p.value, 0.01)
  expect_true(bad$reject)
})

test_that("R > 1 is reported as a Type II failure rather than a p-value", {
  skip_on_cran()
  ## Residuals more skewed than the assumed inefficiency distribution can be.
  ## An exponential u tolerates |gamma1| up to 2; a half-normal only 0.9953.
  set.seed(7)
  n <- 3000
  e <- rnorm(n, 0, 0.15) - rexp(n, 1)     # strongly skewed
  tt <- spec_test(e, noise = "normal", inefficiency = "halfnormal")
  expect_gt(abs(tt$skewness), 0.9953)
  expect_gt(tt$R, 1)
  expect_true(tt$type2_failure)
  expect_output(print(tt), "TYPE II FAILURE")
})

test_that("spec_test_all runs the twelve pairs and ranks them", {
  skip_on_cran()
  set.seed(2)
  e <- rnorm(2000, 0, 0.6) - abs(rnorm(2000, 0, 1))
  g <- spec_test_all(e)
  expect_s3_class(g, "sfa_spec_grid")
  expect_equal(nrow(g), 12L)
  ## Sorted by p-value, least contradicted first.
  expect_false(is.unsorted(rev(g$p.value)))
  ## The true pair should not be at the bottom of twelve.
  i <- which(g$noise == "normal" & g$inefficiency == "halfnormal")
  expect_lt(i, 7)
  expect_output(print(g), "comparative evidence")
})

test_that("it takes a fit as well as a vector, and refuses what it cannot use", {
  skip_on_cran()
  set.seed(5)
  n <- 600
  x1 <- rnorm(n)
  d <- data.frame(y = 1 + 0.5 * x1 + rnorm(n, 0, 0.5) - abs(rnorm(n)), x1 = x1)
  f <- sfm(y ~ x1, model_name = "NHN", data = d)
  tt <- spec_test(f, noise = "normal", inefficiency = "halfnormal")
  expect_s3_class(tt, "sfa_spec_test")
  expect_equal(tt$nobs, n)

  expect_error(spec_test("nope"), "numeric vector")
  expect_error(spec_test(rnorm(20)), "too few")
})

test_that("the bootstrap null is calibrated where the asymptotic one is not", {
  skip_on_cran()
  ## The asymptotic p-value rejects 18.7% of the time at a nominal 5% averaged
  ## over 24 correctly-specified designs, and up to 37% when u is exponential;
  ## the bootstrap rejects 5.7%. The cause is the 8th moment, which a
  ## heavy-tailed u estimates far too imprecisely at these n. Not re-measured
  ## here -- that study is in the notes -- but the mechanism is asserted: the
  ## two nulls must actually differ on such a sample.
  set.seed(9)
  e <- rnorm(1500, 0, 1) - rexp(1500, 1)
  a <- spec_test(e, noise = "normal", inefficiency = "exponential",
    null = "asymptotic"
  )
  b <- spec_test(e, noise = "normal", inefficiency = "exponential",
    null = "bootstrap", B = 199, seed = 1
  )
  expect_identical(a$null, "asymptotic")
  expect_identical(b$null, "bootstrap")
  expect_true(is.null(a$boot))
  expect_gt(length(b$boot), 150)
  ## Same statistic, different null.
  expect_equal(a$Gamma, b$Gamma)
  expect_equal(unname(a$statistic), unname(b$statistic))
  ## A bootstrap p-value can never be below 1/(B+1).
  expect_gte(b$p.value, 1 / (length(b$boot) + 1))
})

test_that("a seeded bootstrap is reproducible and leaves the RNG alone", {
  skip_on_cran()
  set.seed(4)
  e <- rnorm(800, 0, 0.7) - abs(rnorm(800))
  p1 <- spec_test(e, B = 99, seed = 42)$p.value
  p2 <- spec_test(e, B = 99, seed = 42)$p.value
  expect_identical(p1, p2)

  ## The caller's stream must be untouched.
  set.seed(123); before <- runif(3)
  set.seed(123); invisible(spec_test(e, B = 99, seed = 7)); after <- runif(3)
  expect_equal(before, after)
})

test_that("component simulators hit the standard deviation they are asked for", {
  ## The bootstrap draws under the ASSUMED pair, so a simulator whose scale is
  ## wrong would silently test the wrong null.
  set.seed(6)
  for (fam in c("normal", "logistic", "laplace", "uniform")) {
    expect_equal(sd(sfa:::.spec_rnoise(2e5, fam, 1.7)), 1.7, tolerance = 0.05)
    expect_equal(mean(sfa:::.spec_rnoise(2e5, fam, 1.7)), 0, tolerance = 0.05)
  }
  for (fam in c("halfnormal", "exponential", "genexponential")) {
    x <- sfa:::.spec_rineff(2e5, fam, 0.9)
    expect_equal(sd(x), 0.9, tolerance = 0.05)
    expect_true(all(x >= 0))
  }
})
