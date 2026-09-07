## Regression tests for the between-collinearity guard and the GTRE two-step
## standard errors. Both were reported from real use: psfm(model_name = "GTRE")
## with factor(year) on an unbalanced panel died inside plm::ercomp() with
## "system is exactly singular", because the guard could only drop whole
## formula TERMS and the offending columns were part of one factor term.

.mk_cohort_panel <- function(seed = 4, N = 60, y1 = 2000, y2 = 2015, ncoh = 4) {
  set.seed(seed)
  rows <- do.call(rbind, lapply(seq_len(N), function(i) {
    coh <- (i %% ncoh) + 1
    data.frame(id = i, year = (y1 + (coh - 1) * 3):y2)
  }))
  n <- nrow(rows)
  rows$x1 <- rnorm(n)
  rows$x2 <- rnorm(n)
  ai <- rnorm(N, 0, 0.4)[rows$id]
  hi <- abs(rnorm(N, 0, 0.5))[rows$id]
  rows$y <- 1 + 0.5 * rows$x1 + 0.3 * rows$x2 + ai - hi +
    rnorm(n, 0, 0.3) - abs(rnorm(n, 0, 0.6))
  rows
}

test_that("the cohort panel really is between-rank-deficient in only PART of a factor", {
  d <- .mk_cohort_panel()
  chk <- .check_collinearity(y ~ x1 + x2 + factor(year), d, "id")
  expect_true(chk$pooled_rank == chk$pooled_cols) ## full rank pooled
  expect_true(chk$between_rank < chk$between_cols) ## deficient between
  expect_true(length(chk$between_drop) > 0)
  ## every offending column belongs to factor(year), but NOT all of that
  ## factor's columns offend -- so term-level dropping is a no-op here.
  expect_true(all(grepl("^factor\\(year\\)", chk$between_drop)))
  expect_identical(
    .terms_for_columns(y ~ x1 + x2 + factor(year), d, chk$between_drop),
    character(0)
  )
})

test_that(".re_start_design() reduces at COLUMN granularity", {
  d <- .mk_cohort_panel()
  pd <- plm::pdata.frame(d, index = c("id", "year"))
  chk <- .check_collinearity(y ~ x1 + x2 + factor(year), d, "id")
  rd <- .re_start_design(y ~ x1 + x2 + factor(year), pd, chk$between_drop)
  expect_false(is.null(rd))
  ## the map points back at original model-matrix names, minus the dropped ones
  expect_true(all(c("x1", "x2") %in% unname(rd$map)))
  expect_false(any(chk$between_drop %in% unname(rd$map)))
  ## and the retained year dummies survive
  expect_true(any(grepl("^factor\\(year\\)", unname(rd$map))))
  ## the reduced between-design is now full rank, which is the whole point
  X <- stats::model.matrix(rd$formula, data = as.data.frame(rd$data))
  idv <- as.character(as.data.frame(plm::index(rd$data))[[1]])
  Z <- rowsum(as.data.frame(X), group = idv, reorder = FALSE)
  Z <- as.matrix(Z / as.numeric(table(idv)[rownames(Z)]))
  expect_equal(qr(Z)$rank, ncol(Z))
})

test_that(".remap_start_names() and .expand_start_se() align by name", {
  map <- c(.v1 = "x1", .v2 = "factor(year)2003")
  v <- c("(Intercept)" = 1, .v1 = 2, .v2 = 3)
  expect_identical(names(.remap_start_names(v, map)),
                   c("(Intercept)", "x1", "factor(year)2003"))
  se <- .expand_start_se(c(x1 = 0.1, x2 = 0.2), c("x1", "zzz", "x2"))
  expect_equal(unname(se), c(0.1, NA, 0.2))
  ## a column the starting-value fit could not identify gets NA, not a
  ## silently misaligned neighbour's SE
  expect_true(is.na(.expand_start_se(c(x1 = 0.1), c("zzz"))[["zzz"]]))
})

test_that("psfm() fits with a partially collinear factor instead of erroring", {
  skip_on_cran()
  d <- .mk_cohort_panel()
  for (mn in c("GTRE_SEQ1", "GTRE_SEQ2")) {
    fit <- suppressWarnings(
      psfm(y ~ x1 + x2 + factor(year), model_name = mn, data = d, individual = "id")
    )
    expect_s3_class(fit, "sfareg")
    expect_true(all(is.finite(fit$out[, "par"])))
    ## the year dummies are still in the reported model
    expect_true(sum(grepl("^factor\\(year\\)", rownames(fit$out))) > 0)
    ## slope estimates stay near truth despite the reduced starting values
    expect_equal(unname(fit$out["x1", "par"]), 0.5, tolerance = 0.15)
    expect_equal(unname(fit$out["x2", "par"]), 0.3, tolerance = 0.15)
  }
})

test_that("collinear_action is honoured for a partially collinear factor", {
  skip_on_cran()
  d <- .mk_cohort_panel()
  expect_error(
    psfm(y ~ x1 + x2 + factor(year), model_name = "GTRE_SEQ1", data = d,
         individual = "id", collinear_action = "error"),
    "Collinearity detected"
  )
  ## "warn_drop" cannot remove part of a factor term; it must say so rather
  ## than silently doing nothing.
  expect_warning(
    psfm(y ~ x1 + x2 + factor(year), model_name = "GTRE_SEQ1", data = d,
         individual = "id", collinear_action = "warn_drop"),
    "could not remove these"
  )
})

test_that(".gtre_two_step_se() matches the sampling sd of the moment inversion", {
  skip_on_cran()
  su <- 1; sv <- 0.5; n <- 2000; R <- 300
  set.seed(99)
  est <- se <- matrix(NA_real_, R, 2)
  for (r in seq_len(R)) {
    z <- rnorm(n, 0, sv) - abs(rnorm(n, 0, su))
    z <- z - mean(z)
    ts <- .gtre_two_step(z, z, 0)
    s5 <- .gtre_two_step_se(z, z, n, n, NA_real_)
    est[r, ] <- c(ts$gamma_uv, ts$sigmaSq_uv)
    se[r, ] <- c(s5[["gamma_uv"]], s5[["sigmaSq_uv"]])
  }
  ok <- stats::complete.cases(est) & stats::complete.cases(se)
  expect_gt(sum(ok), 0.9 * R)
  for (j in 1:2) {
    ratio <- mean(se[ok, j]) / stats::sd(est[ok, j])
    ## The pre-fix expressions gave 16x for gamma and 1.37x for sigmaSq.
    expect_gt(ratio, 0.85)
    expect_lt(ratio, 1.15)
  }
})

test_that(".gtre_two_step_se() returns NA at the wrong-skew boundary", {
  set.seed(1)
  z <- abs(rnorm(500)) ## positively skewed: m3 > 0, so min(0, m3) == 0
  s5 <- .gtre_two_step_se(z, z, 500, 500, 0.1)
  expect_true(all(is.na(s5)))
  expect_named(s5, c("gamma_uv", "sigmaSq_uv", "gamma_hr", "sigmaSq_hr", "beta_0"))
})
