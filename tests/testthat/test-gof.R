## Wang, Amsler and Schmidt (2011) goodness-of-fit tests.

.gof_data <- function(N, seed = 5) {
  as.data.frame(data_gen_cs(N = N, rand = seed, sig_u = 1, sig_v = 0.5,
    cons = 0.5, beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.1
  ))
}

test_that("the KS statistic is the usual two-sided one", {
  set.seed(1)
  p <- runif(500)
  expect_equal(.gof_ks(p), unname(stats::ks.test(p, "punif")$statistic),
    tolerance = 1e-12
  )
  ## a badly non-uniform PIT gives a large statistic
  expect_gt(.gof_ks(runif(500)^3), 0.2)
})

test_that("the chi-square statistic uses equiprobable cells", {
  ## A perfectly uniform PIT puts exactly n/k in every cell, so the statistic
  ## is 0 -- which is what "equiprobable" has to mean for E_j = n/k to hold.
  k <- 10
  p <- (seq_len(1000) - 0.5) / 1000
  expect_equal(.gof_chisq(p, k), 0, tolerance = 1e-8)
  ## everything crammed into one cell is the maximum, n*(k-1)
  expect_equal(.gof_chisq(rep(0.05, 1000), k), 1000 * (k - 1))
  ## boundary values do not fall outside the bins
  expect_true(is.finite(.gof_chisq(c(0, 1, 0.5), 5)))
})

test_that("gof_test() validates its inputs", {
  skip_on_cran()
  d <- .gof_data(200)
  f <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = d)
  expect_error(gof_test(f, data = d, cells = 2), "cells")
  expect_error(gof_test(f, data = d, cells = 500), "cells")
  expect_error(gof_test(structure(list(), class = "lm")), "sfareg")
  ## too many parameters for too few cells leaves no degrees of freedom
  expect_error(gof_test(f, data = d, test = "chisq", null = "asymptotic", cells = 3),
    "degrees of freedom"
  )
})

test_that("the asymptotic branch refuses to invent a KS p-value", {
  skip_on_cran()
  d <- .gof_data(300)
  f <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = d)
  expect_warning(g <- gof_test(f, data = d, null = "asymptotic"), "no valid asymptotic")
  expect_true(is.na(g$p.value[g$test == "ks"]))
  expect_true(is.finite(g$p.value[g$test == "chisq"]))
  expect_match(attr(g, "note"), "CONSERVATIVE")
})

test_that("the bootstrap runs, refits, and returns usable p-values", {
  skip_on_cran()
  d <- .gof_data(250)
  f <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = d)
  g <- gof_test(f, data = d, B = 39, seed = 1)
  expect_equal(nrow(g), 2L)
  expect_true(all(g$B > 30)) ## replications actually usable
  expect_true(all(g$p.value > 0 & g$p.value <= 1))
  ## the bootstrap p-value can never be exactly 0 with finite B
  expect_gte(min(g$p.value), 1 / 40)
  expect_true(is.matrix(attr(g, "boot")))
  ## reproducible
  g2 <- gof_test(f, data = d, B = 39, seed = 1)
  expect_equal(g$p.value, g2$p.value)
})

test_that("a correctly specified model is not rejected", {
  skip_on_cran()
  d <- .gof_data(300)
  f <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = d)
  g <- gof_test(f, data = d, B = 99, seed = 1)
  expect_true(all(g$p.value > 0.05))
})

test_that("misspecification is detected once n is large enough", {
  skip_on_cran()
  ## Power against the half-normal/exponential pair is about 0.4 at n = 200
  ## and essentially 1 by n = 800 (200 replications; see ?gof_test). A single
  ## draw at n = 300 therefore says nothing either way, so this asserts only
  ## the large-n half, where the answer is not a coin flip.
  d20 <- .gof_data(2000)
  p20 <- composed_cdf(sfm(y_pcs_e ~ x1 + x2, model_name = "NHN", data = d20), data = d20)
  expect_lt(suppressWarnings(stats::ks.test(p20, "punif")$p.value), 0.05)
  ## and the correctly specified fit on the SAME data is not flagged
  p20ok <- composed_cdf(sfm(y_pcs_e ~ x1 + x2, model_name = "NE", data = d20), data = d20)
  expect_gt(suppressWarnings(stats::ks.test(p20ok, "punif")$p.value), 0.10)
})
