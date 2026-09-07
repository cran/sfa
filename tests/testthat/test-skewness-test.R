## Formal wrong-skewness tests.

test_that("D'Agostino matches the reference implementation exactly", {
  ## Golden values from moments::agostino.test(), which exists for this
  ## statistic and nothing else, recorded 2026-08-29 on R 4.5.2. Pinned rather
  ## than computed live so the package does not take a Suggests dependency for
  ## a single validation test -- R CMD check raises "unstated dependencies in
  ## tests" otherwise, and CI fails the build on a WARNING.
  ##
  ## To re-derive, with moments installed:
  ##   set.seed(3)
  ##   for (nn in c(30, 60, 200, 1000)) {
  ##     e <- rnorm(nn) - abs(rnorm(nn, 0, 1.2))
  ##     print(moments::agostino.test(e))
  ##   }
  golden <- list(
    `30`   = c(z =  0.019343441854767732, p = 0.98456712881362896),
    `60`   = c(z = -1.7750652500210233,   p = 0.075887113954654195),
    `200`  = c(z = -0.79578274891458356,  p = 0.42615832035494572),
    `1000` = c(z = -3.733951714321849,    p = 0.00018849864797365612)
  )
  set.seed(3)
  for (nn in c(30, 60, 200, 1000)) {
    e <- rnorm(nn) - abs(rnorm(nn, 0, 1.2))
    mine <- skewness_test(e, test = "agostino", alternative = "two.sided")
    g <- golden[[as.character(nn)]]
    expect_equal(unname(mine$statistic), unname(g[["z"]]),
      tolerance = 1e-10, info = paste("n =", nn))
    expect_equal(mine$p.value, unname(g[["p"]]), tolerance = 1e-10)
  }
})

test_that("Coelli's M3T is the statistic it claims to be", {
  set.seed(5)
  e <- rnorm(300) - abs(rnorm(300))
  ec <- e - mean(e)
  m2 <- mean(ec^2); m3 <- mean(ec^3)
  z <- skewness_test(e, test = "coelli", alternative = "two.sided")
  expect_equal(unname(z$statistic), m3 / sqrt(6 * m2^3 / length(e)),
    tolerance = 1e-12)
  expect_equal(unname(z$estimate), m3 / m2^1.5, tolerance = 1e-12)
})

test_that("both tests hold their nominal size on symmetric residuals", {
  skip_on_cran()
  ## The reason "agostino" is the default: Coelli's asymptotic form is
  ## conservative in small samples. Loose bounds -- this is 600 reps, not the
  ## 4000 the documented table used -- but they would catch a sign error or a
  ## badly wrong transformation.
  set.seed(17)
  for (nn in c(50, 400)) {
    pa <- pc <- numeric(600)
    for (i in seq_len(600)) {
      e <- rnorm(nn)
      pa[i] <- skewness_test(e, "agostino", alternative = "less")$p.value
      pc[i] <- skewness_test(e, "coelli", alternative = "less")$p.value
    }
    expect_lt(abs(mean(pa < 0.05) - 0.05), 0.03, label = paste("agostino n =", nn))
    ## Coelli is allowed to be conservative but not anti-conservative.
    expect_lt(mean(pc < 0.05), 0.09)
  }
})

test_that("it detects the direction of skew", {
  set.seed(9)
  ## A strongly skewed design on purpose. The first attempt used
  ## rnorm(400) - abs(rnorm(400, 0, 1.5)), which gives a median p of 0.0075 but
  ## ranged up to 0.25 across seeds -- an assertion of p < 0.01 on that would
  ## have been another knife-edge test of exactly the kind this suite has been
  ## losing CI time to. Normal noise minus an exponential gives median
  ## p ~ 7e-24 and a worst case of 8e-14 over 20 seeds.
  neg <- rnorm(400, 0, 0.3) - rexp(400, 1)
  pos <- rnorm(400, 0, 0.3) + rexp(400, 1)
  expect_lt(skewness_test(neg)$p.value, 1e-8)
  expect_gt(skewness_test(pos)$p.value, 1 - 1e-8)
  expect_false(skewness_test(neg)$wrong_skew)
  expect_true(skewness_test(pos)$wrong_skew)
})

test_that("it reads an sfareg fit, and agrees with the fit's own diagnosis", {
  skip_on_cran()
  for (r in c(1, 4)) {
    d <- data_gen_cs(N = 100, rand = r, sig_u = 0.3, sig_v = 0.4, cons = 0.5,
                     beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.1)
    f <- suppressWarnings(sfm(y_pcs_e ~ x1 + x2, model_name = "NE", data = d))
    tt <- skewness_test(f)
    expect_s3_class(tt, "htest")
    ## The test and the fit must not disagree about the sign of m3.
    expect_identical(tt$wrong_skew, f$wrong_skew, info = paste("seed", r))
    expect_equal(tt$m3, f$residual_m3, tolerance = 1e-10)
    expect_equal(tt$nobs, nobs(f))
  }
})

test_that("it refuses what it cannot compute", {
  expect_error(skewness_test(list(1, 2)), "must be an")
  expect_error(skewness_test(rnorm(4), test = "agostino"), "not computable")
  expect_error(skewness_test(rep(1, 50)), "zero variance")
})
