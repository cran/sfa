## Observation weights (H6 / I5).

wt_data <- function(seed = 3, n = 600, het = FALSE) {
  set.seed(seed)
  x1 <- rnorm(n); x2 <- rnorm(n); z <- rnorm(n)
  su <- if (het) exp(-0.2 + 0.3 * z) else 0.8
  y <- 1 + 0.5 * x1 + 0.5 * x2 + rnorm(n, 0, 0.3) - abs(rnorm(n, 0, su))
  data.frame(y = y, x1 = x1, x2 = x2, z = z)
}

test_that("unit weights reproduce the unweighted fit", {
  skip_on_cran()
  d <- wt_data()
  a <- sfm(y ~ x1 + x2, model_name = "NHN", data = d)
  b <- sfm(y ~ x1 + x2, model_name = "NHN", data = d, weights = rep(1, nrow(d)))
  expect_equal(a$out[, "par"], b$out[, "par"], tolerance = 1e-6)
})

test_that("a weight of 2 is the same as including the row twice", {
  skip_on_cran()
  ## The defining property of FREQUENCY weights, and the strongest available
  ## check that they enter the likelihood correctly rather than merely being
  ## accepted. wscale = FALSE because duplicating rows changes n and rescaling
  ## would put the two objectives on different scales.
  d <- wt_data()
  n <- nrow(d)
  dup <- d[c(seq_len(n), seq_len(100)), ]
  w <- c(rep(2, 100), rep(1, n - 100))

  a <- sfm(y ~ x1 + x2, model_name = "NHN", data = dup)
  b <- sfm(y ~ x1 + x2, model_name = "NHN", data = d, weights = w, wscale = FALSE)
  expect_equal(unname(a$out[, "par"]), unname(b$out[, "par"]), tolerance = 1e-3)
})

test_that("the same holds on the heteroskedastic path", {
  skip_on_cran()
  d <- wt_data(het = TRUE)
  n <- nrow(d)
  a <- sfm(y ~ x1 + x2, model_name = "NHN", data = d, uhet = ~z)
  b <- sfm(y ~ x1 + x2, model_name = "NHN", data = d, uhet = ~z,
    weights = rep(1, n)
  )
  expect_equal(a$out[, "par"], b$out[, "par"], tolerance = 1e-6)

  dup <- d[c(seq_len(n), seq_len(100)), ]
  w <- c(rep(2, 100), rep(1, n - 100))
  c1 <- sfm(y ~ x1 + x2, model_name = "NHN", data = dup, uhet = ~z)
  c2 <- sfm(y ~ x1 + x2, model_name = "NHN", data = d, uhet = ~z,
    weights = w, wscale = FALSE
  )
  expect_equal(unname(c1$out[, "par"]), unname(c2$out[, "par"]), tolerance = 1e-3)
})

test_that("wscale rescales the objective without moving the estimates", {
  skip_on_cran()
  d <- wt_data()
  w <- runif(nrow(d), 0.5, 2)
  a <- sfm(y ~ x1 + x2, model_name = "NHN", data = d, weights = w, wscale = TRUE)
  b <- sfm(y ~ x1 + x2, model_name = "NHN", data = d, weights = w, wscale = FALSE)
  ## A common factor on the log-likelihood cannot move the argmax.
  expect_equal(unname(a$out[, "par"]), unname(b$out[, "par"]), tolerance = 1e-4)
  ## But it does move the log-likelihood itself, which is the point of it.
  expect_false(isTRUE(all.equal(a$opt$value, b$opt$value)))
})

test_that("bad weights are refused rather than recycled", {
  d <- wt_data(n = 300)
  n <- nrow(d)
  expect_error(
    sfm(y ~ x1 + x2, model_name = "NHN", data = d, weights = rep(1, 5)),
    "has length 5 but 300 observations"
  )
  expect_error(
    sfm(y ~ x1 + x2, model_name = "NHN", data = d, weights = c(-1, rep(1, n - 1))),
    "must be non-negative"
  )
  expect_error(
    sfm(y ~ x1 + x2, model_name = "NHN", data = d, weights = rep(0, n)),
    "cannot all be zero"
  )
  expect_error(
    sfm(y ~ x1 + x2, model_name = "NHN", data = d, weights = c(NA, rep(1, n - 1))),
    "must all be finite"
  )
  ## The robust divergences reweight observations themselves; a second set of
  ## weights on top makes neither interpretable.
  expect_error(
    sfm(y ~ x1 + x2, model_name = "NHN", data = d, weights = rep(1, n),
      robust = "mdpd"),
    "cannot be combined with a robust"
  )
})
