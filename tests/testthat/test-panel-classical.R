## Cornwell-Schmidt-Sickles (1990) and Lee-Schmidt (1993).

.cl_panel <- function(N = 90, Tt = 8, seed = 3) {
  set.seed(seed)
  d <- expand.grid(t = seq_len(Tt), name = seq_len(N))
  d$year <- d$t
  n <- nrow(d)
  d$x1 <- runif(n, 1, 5)
  d$x2 <- runif(n, 1, 5)
  d
}

test_that("CSS recovers the frontier under firm-specific quadratic effects", {
  skip_on_cran()
  d <- .cl_panel()
  N <- max(d$name)
  th <- cbind(rnorm(N, 0, .5), rnorm(N, 0, .08), rnorm(N, 0, .008))
  d$y <- 0.6 * d$x1 + 0.3 * d$x2 +
    th[d$name, 1] + th[d$name, 2] * d$t + th[d$name, 3] * d$t^2 +
    rnorm(nrow(d), 0, .2)

  f <- psfm(y ~ x1 + x2, model_name = "CSS", data = d,
    individual = "name", time = "year")

  expect_s3_class(f, "sfareg")
  expect_identical(names(coef(f)), c("x1", "x2"))
  expect_true(all(abs(coef(f) - c(0.6, 0.3)) < 3 * f$std.errors))
  expect_equal(f$sigma_v, 0.2, tolerance = 0.05)
  expect_equal(dim(f$theta), c(N, 3L))
  ## Time-varying, and normalized within each period.
  expect_length(f$u_hat, nrow(d))
  expect_true(all(f$u_hat >= 0))
  expect_true(all(tapply(f$u_hat, d$year, min) < 1e-8))
  ## Not maximum likelihood, so no information criteria.
  expect_true(is.na(suppressWarnings(as.numeric(logLik(f)))))
})

test_that("LS recovers the common temporal pattern and the frontier", {
  skip_on_cran()
  d <- .cl_panel(seed = 9)
  N <- max(d$name)
  delta_true <- c(1, 1.4, 1.9, 2.3, 2.5, 2.4, 2.0, 1.5)
  ai <- -abs(rnorm(N, 0, .6))
  d$y <- 0.6 * d$x1 + 0.3 * d$x2 + delta_true[d$t] * ai[d$name] +
    rnorm(nrow(d), 0, .2)

  g <- psfm(y ~ x1 + x2, model_name = "LS", data = d,
    individual = "name", time = "year")

  expect_true(all(abs(coef(g) - c(0.6, 0.3)) < 3 * g$std.errors))
  expect_true(g$ls_converged)
  ## delta_1 is the normalization, so it is 1 exactly.
  expect_equal(unname(g$delta[1]), 1)
  ## The shape is what is identified; correlate rather than compare levels.
  expect_gt(cor(as.numeric(g$delta), delta_true), 0.99)
  expect_gt(abs(cor(g$alpha_i, ai)), 0.95)
  expect_true(all(tapply(g$u_hat, d$year, min) < 1e-8))
})

test_that("CSS refuses what it cannot identify", {
  skip_on_cran()
  d <- .cl_panel(N = 30, Tt = 3)
  d$y <- 0.6 * d$x1 + rnorm(nrow(d), 0, .2)
  ## T = 3 saturates the firm quadratic exactly.
  expect_error(
    psfm(y ~ x1 + x2, "CSS", d, individual = "name", time = "year"),
    "fewer than 4 periods"
  )

  d2 <- .cl_panel(N = 40, Tt = 8)
  d2$trend <- d2$t
  d2$y <- 0.6 * d2$x1 + rnorm(nrow(d2), 0, .2)
  ## A pure time trend is spanned by every firm's own quadratic.
  expect_error(
    psfm(y ~ x1 + trend, "CSS", d2, individual = "name", time = "year"),
    "not separately identified"
  )
})

test_that("both estimators drop the intercept into the firm effect", {
  skip_on_cran()
  d <- .cl_panel(N = 40, Tt = 6)
  d$y <- 2 + 0.6 * d$x1 + 0.3 * d$x2 - abs(rnorm(nrow(d), 0, .4)) +
    rnorm(nrow(d), 0, .2)
  for (m in c("CSS", "LS")) {
    f <- psfm(y ~ x1 + x2, m, d, individual = "name", time = "year")
    expect_false("(Intercept)" %in% names(coef(f)))
    expect_equal(nobs(f), nrow(d))
  }
})

## Schmidt-Sickles random effects (SSRE) and correlated random effects (SSCRE).

.ss_panel <- function(N = 80, Tt = 8, corr = 0, seed = 11) {
  set.seed(seed)
  d <- expand.grid(t = seq_len(Tt), name = seq_len(N))
  d$year <- d$t
  ## The firm effect, and regressors correlated with it by `corr`. corr = 0 is
  ## the world RE assumes; corr > 0 is the world it is wrong in.
  a_i <- -abs(rnorm(N, 0, 0.8))
  d$x1 <- runif(nrow(d), 1, 5) + corr * a_i[d$name]
  d$x2 <- runif(nrow(d), 1, 5)
  d$y <- 0.6 * d$x1 + 0.3 * d$x2 + a_i[d$name] + rnorm(nrow(d), 0, 0.25)
  attr(d, "alpha") <- a_i
  d
}

test_that("SSCRE reproduces the within estimator exactly (Mundlak)", {
  skip_on_cran()
  d <- .ss_panel()
  fe <- psfm(y ~ x1 + x2, "SSFE", d, individual = "name", time = "year")
  cre <- psfm(y ~ x1 + x2, "SSCRE", d, individual = "name", time = "year")
  ## Adding the group means to a random-effects fit recovers the within
  ## slopes. This is Mundlak (1978) and it is exact, not approximate.
  expect_equal(coef(cre)[c("x1", "x2")], coef(fe)[c("x1", "x2")],
    tolerance = 1e-6
  )
  expect_true(all(c(".mean_x1", ".mean_x2") %in% names(coef(cre))))
})

test_that("RE is biased under correlated effects where CRE is not", {
  skip_on_cran()
  ## x1 correlated with the firm effect: exactly what RE assumes away.
  d <- .ss_panel(corr = 1.5, seed = 4)
  re <- psfm(y ~ x1 + x2, "SSRE", d, individual = "name", time = "year")
  cre <- psfm(y ~ x1 + x2, "SSCRE", d, individual = "name", time = "year")
  err_re <- abs(coef(re)[["x1"]] - 0.6)
  err_cre <- abs(coef(cre)[["x1"]] - 0.6)
  expect_gt(err_re, err_cre)
  expect_lt(err_cre, 0.05)
  ## And the Mundlak term is what detects it -- significant here, unlike in
  ## the uncorrelated design above.
  expect_gt(abs(cre$t.values[[".mean_x1"]]), 2)
})

test_that("SSRE and SSCRE return per-observation scores and no logLik", {
  skip_on_cran()
  d <- .ss_panel()
  for (m in c("SSRE", "SSCRE")) {
    f <- psfm(y ~ x1 + x2, m, d, individual = "name", time = "year")
    expect_equal(nobs(f), nrow(d))
    expect_length(f$exp_u_hat, nrow(d))
    expect_true(all(f$u_hat >= 0))
    ## One firm is the benchmark, so its inefficiency is exactly zero.
    expect_equal(min(f$u_hat), 0)
    expect_length(f$alpha_hat, max(d$name))
    expect_true(is.na(suppressWarnings(as.numeric(logLik(f)))))
  }
})

test_that("SSCRE refuses a design with nothing to correct", {
  skip_on_cran()
  d <- .ss_panel(N = 40, Tt = 5)
  ## A regressor constant within firm is what RE exists to identify; with only
  ## such regressors there are no Mundlak means to form.
  d$z <- rep(rnorm(40), each = 5)
  expect_error(
    psfm(y ~ z, "SSCRE", d, individual = "name", time = "year"),
    "nothing for the Mundlak means to correct"
  )
})
