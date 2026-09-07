## start_from: seed a hard model from a simpler fit already in hand (G2/I3).

sf_data <- function(seed, n = 300) {
  set.seed(seed)
  x1 <- rnorm(n); x2 <- rnorm(n)
  y <- 2 + 0.5 * x1 + 0.5 * x2 + rnorm(n, 0, 0.3) - rgamma(n, shape = 2, scale = 0.4)
  data.frame(y = y, x1 = x1, x2 = x2)
}

test_that("matching is by NAME, so shared scale names carry and others do not", {
  skip_on_cran()
  d <- sf_data(1)
  ne <- sfm(y ~ x1 + x2, model_name = "NE", data = d)
  nhn <- sfm(y ~ x1 + x2, model_name = "NHN", data = d)

  ## NE and NG both call their scales sigv/sigu, so five of six carry.
  expect_message(
    sfm(y ~ x1 + x2, model_name = "NG", data = d, start_from = ne),
    "carried 5 of 6"
  )
  ## NHN calls them lambda/sigma. Those must NOT be poured into NG's sigv/sigu
  ## -- a positional copy would slide a variance ratio into a standard
  ## deviation. Only the three frontier coefficients are shared.
  expect_message(
    sfm(y ~ x1 + x2, model_name = "NG", data = d, start_from = nhn),
    "carried 3 of 6"
  )
  ## NHN -> NTN shares lambda and sigma, which is the documented idiom.
  expect_message(
    sfm(y ~ x1 + x2, model_name = "NTN", data = d, start_from = nhn),
    "carried 5 of 6"
  )
})

test_that("seeding never makes the fit worse", {
  skip_on_cran()
  ## "Never worse" is the only property of start_from that holds everywhere, so
  ## it is the only one asserted. The SIZE of any rescue does not travel: over
  ## seeds 1-40 on this machine, thirty-nine give an indistinguishable optimum
  ## either way and seed 4 gains 176.70 log-likelihood points. On the CI
  ## platforms seed 4 converges to the same optimum from both starts. Asserting
  ## the gain would pin a BLAS/optimizer path rather than the feature.
  ##
  ## Re-derived 2026-09-04. An earlier version of this comment said seeds 7 and
  ## 10 made the UNSEEDED fit fail outright in optim. That is no longer true:
  ## the I12 guard in opts.R catches a non-finite third-stage objective and
  ## falls back, and all forty seeds now converge from both starts. What
  ## remains is the quiet one -- seed 4 converging to a much worse optimum with
  ## nothing in the returned object to flag it.
  ##
  ## The error-tolerating branch below is KEPT even though nothing currently
  ## triggers it: it exists so that NG's start-value fragility cannot be misread
  ## as a start_from defect on some other platform, and that reason has not
  ## expired just because this machine no longer reproduces the failure.
  ll <- function(e) if (inherits(e, "error")) NA_real_ else -e$opt$value
  for (s in c(1, 2, 4)) {
    d <- sf_data(s)
    ne <- sfm(y ~ x1 + x2, model_name = "NE", data = d)
    a <- tryCatch(sfm(y ~ x1 + x2, model_name = "NG", data = d),
      error = function(e) e
    )
    b <- tryCatch(
      suppressMessages(
        sfm(y ~ x1 + x2, model_name = "NG", data = d, start_from = ne)
      ),
      error = function(e) e
    )
    expect_false(inherits(b, "error"))
    if (!is.na(ll(a))) expect_gte(ll(b), ll(a) - 1e-6)
  }
})

test_that("start_from rejects what is not a fit", {
  d <- sf_data(2)
  expect_error(
    sfm(y ~ x1 + x2, model_name = "NG", data = d, start_from = list(a = 1)),
    "must be a fitted"
  )
  expect_error(
    sfm(y ~ x1 + x2, model_name = "NG", data = d, start_from = 1:5),
    "must be a fitted"
  )
})
