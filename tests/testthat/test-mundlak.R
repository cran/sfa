## Mundlak adjustment terms for the true-effects panel models (L10).
##
## The helpers are pure functions of a formula and a data frame, so most of
## this is exact rather than statistical: which columns get built, what they
## contain, where they land in a piped formula, and what is refused. The one
## behavioural claim -- that the adjustment reduces heterogeneity bias -- is a
## Monte Carlo result recorded in the notes, not asserted here on one seed.

mund_panel <- function(seed = 1, N = 30, Tt = 5, rho = 1.5) {
  set.seed(seed)
  mu <- rnorm(N)
  id <- rep(seq_len(N), each = Tt)
  x <- mu[id] + rnorm(N * Tt, 0, 0.7)
  z <- rep(rnorm(N), each = Tt)      # time-invariant environmental factor
  a <- rho * mu + rnorm(N, 0, 0.3)
  data.frame(
    y = 1 + 0.5 * x + a[id] + rnorm(N * Tt, 0, 0.3) - abs(rnorm(N * Tt, 0, 0.6)),
    x = x, z = z, id = id, t = rep(seq_len(Tt), N)
  )
}

test_that("group means are the within-firm means, and named for their source", {
  d <- mund_panel()
  aug <- sfa:::.mundlak_augment(d, d$id, "x")
  expect_identical(aug$added, "mbar_x")
  expect_true("mbar_x" %in% names(aug$data))
  ## Exactly the firm mean, firm by firm.
  expect_equal(
    aug$data$mbar_x,
    as.numeric(ave(d$x, d$id, FUN = mean))
  )
  ## Constant within firm, by construction.
  expect_equal(
    length(unique(round(tapply(aug$data$mbar_x, d$id, sd), 10))), 1L
  )
  expect_equal(unname(tapply(aug$data$mbar_x, d$id, sd))[1], 0)
})

test_that("a variable already constant within firm gets NO mean column", {
  ## Its own mean IS itself, so including both would be exactly collinear.
  ## Karagiannis and Kellermann's time-invariant environmental factors are
  ## this case: they enter the auxiliary equation directly, not via a mean.
  d <- mund_panel()
  aug <- sfa:::.mundlak_augment(d, d$id, c("x", "z"))
  expect_identical(aug$added, "mbar_x")
  expect_identical(aug$dropped, "z")
  expect_false("mbar_z" %in% names(aug$data))
})

test_that("terms are inserted into the FRONTIER segment, never a variance one", {
  ## psfm()'s later pipe segments parameterize variances. A group mean landing
  ## there would silently change the model rather than the design matrix.
  f1 <- sfa:::.add_to_frontier(y ~ x1 + x2, "mbar_x1")
  expect_equal(deparse(f1[[3]]), "x1 + x2 + mbar_x1")

  f2 <- sfa:::.add_to_frontier(y ~ x1 | z1, c("mbar_x1"))
  s2 <- paste(deparse(f2), collapse = " ")
  expect_match(s2, "x1 \\+ mbar_x1")
  ## the variance segment is untouched
  expect_match(s2, "\\|\\s*z1\\s*$")

  f3 <- sfa:::.add_to_frontier(y ~ x1 | z1 | z2, "mbar_x1")
  s3 <- paste(deparse(f3), collapse = " ")
  expect_match(s3, "x1 \\+ mbar_x1")
  expect_match(s3, "\\|\\s*z1\\s*\\|\\s*z2\\s*$")
})

test_that("`~ .` means every frontier regressor, and only those", {
  d <- mund_panel()
  ## z1 is a variance determinant here, so it must not be averaged.
  out <- sfa:::.apply_mundlak(y ~ x + z | z, d, id = d$id, mundlak = ~.)
  expect_true("mbar_x" %in% out$added)
  expect_false("mbar_z" %in% out$added)   # z is constant within firm
  expect_identical(sfa:::.frontier_vars(y ~ x + z | z), c("x", "z"))
  expect_identical(sfa:::.frontier_vars(y ~ x1 | z1 | z2), "x1")
})

test_that("the refusals name the reason", {
  d <- mund_panel()
  expect_error(sfa:::.apply_mundlak(y ~ x, d, d$id, mundlak = "x"), "must be a formula")
  expect_error(sfa:::.apply_mundlak(y ~ x, d, d$id, mundlak = ~nope), "not in `data`")
  ## Nothing varies within firm -> the means would all be collinear.
  expect_error(sfa:::.apply_mundlak(y ~ z, d, d$id, mundlak = ~z),
    "collinear"
  )
  ## Naming a non-frontier variable still works but warns, because Mundlak's
  ## device is about the regressors the effect correlates with.
  expect_warning(sfa:::.apply_mundlak(y ~ x | z, d, d$id, mundlak = ~ x + t),
    "not frontier regressors"
  )
})

test_that("NULL is a clean no-op", {
  d <- mund_panel()
  out <- sfa:::.apply_mundlak(y ~ x, d, d$id, mundlak = NULL)
  expect_identical(out$formula, y ~ x)
  expect_identical(out$data, d)
  expect_length(out$added, 0L)
})

test_that("the between-groups collinearity warning is expected, and harmless", {
  skip_on_cran()
  ## EVERY Mundlak fit triggers this, and it is not a defect. mbar_x IS the
  ## firm mean of x, so in the BETWEEN-individual regression that supplies
  ## starting values the two are exactly collinear -- by construction, not by
  ## accident. The package drops it from the starting-value step only, falls
  ## back to pooled OLS starts, and estimates the likelihood with the full
  ## requested formula. Asserted here so that a user who reports the warning
  ## can be told it is expected, and so that a future change which makes it
  ## FATAL is caught.
  d <- mund_panel(N = 40, Tt = 5)
  expect_warning(
    f <- psfm(y ~ x, model_name = "TRE", data = d, individual = "id",
      halton_num = 50, rand.gtre = 7, mundlak = ~x
    ),
    "starting-value"
  )
  expect_true("mbar_x" %in% rownames(f$out))
  expect_true(is.finite(as.numeric(logLik(f))))
})

test_that("psfm() widens the frontier and reports the extra coefficient", {
  skip_on_cran()
  d <- mund_panel(N = 40, Tt = 5)
  f <- suppressWarnings(psfm(y ~ x, model_name = "TRE", data = d,
    individual = "id", halton_num = 50, rand.gtre = 7, mundlak = ~x
  ))
  expect_true("mbar_x" %in% rownames(f$out))
  expect_true("x" %in% rownames(f$out))
  ## The likelihood is unchanged -- this is a design-matrix change -- so the
  ## fit is still an ordinary sfareg with the usual machinery.
  expect_s3_class(f, "sfareg")
  expect_true(is.finite(as.numeric(logLik(f))))
})

test_that("mundlak is refused for the models that difference the effect out", {
  skip_on_cran()
  d <- mund_panel(N = 20, Tt = 4)
  for (m in c("TFE", "FD", "SSFE")) {
    expect_error(
      suppressWarnings(
        psfm(y ~ x, model_name = m, data = d, individual = "id", mundlak = ~x)
      ),
      "RANDOM"
    )
  }
})
