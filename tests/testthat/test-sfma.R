## sfma(): averaging over inefficiency distributions (L8).
##
## The optimiser is exact arithmetic and is tested as such -- against a
## brute-force grid, not against itself. The statistical claims are tested on
## data whose true inefficiency distribution is known, so "the right model gets
## the weight" is a checkable statement rather than a hope.

test_that("simplex projection lands on the simplex, exactly", {
  set.seed(1)
  for (i in 1:200) {
    v <- rnorm(sample(2:8, 1), 0, 2)
    p <- sfa:::.proj_simplex(v)
    expect_equal(sum(p), 1)
    expect_true(all(p >= 0))
    expect_length(p, length(v))
  }
  ## A point already on the simplex is left alone.
  w <- c(0.2, 0.5, 0.3)
  expect_equal(sfa:::.proj_simplex(w), w)
  ## Dominated coordinates are set to EXACTLY zero, not merely small. The
  ## weights are the output, so a zero weight has to be a real zero.
  p <- sfa:::.proj_simplex(c(5, -1, -2, -3))
  expect_identical(p[2:4], c(0, 0, 0))
  expect_equal(p[1], 1)
})

test_that("the QP finds the constrained optimum, checked against a grid", {
  ## The only honest check of an optimiser is a different method. S = 3 so the
  ## simplex can be enumerated finely.
  set.seed(2)
  n <- 40
  A <- matrix(rnorm(n * 3), n, 3)
  r <- rnorm(n)
  cv <- c(0.3, 0.1, 0.9)
  w <- sfa:::.simplex_qp(A, r, cv)

  expect_equal(sum(w), 1)
  expect_true(all(w >= 0))

  obj <- function(z) sum((A %*% z - r)^2) + sum(cv * z)
  g <- expand.grid(a = seq(0, 1, by = 0.005), b = seq(0, 1, by = 0.005))
  g <- g[g$a + g$b <= 1, ]
  g$c <- 1 - g$a - g$b
  best <- min(apply(as.matrix(g), 1, obj))
  ## The QP must be at least as good as the best grid point.
  expect_lte(obj(w), best + 1e-8)
})

test_that("a heavy penalty drives every weight onto the cheapest model", {
  ## The penalty term n^{1/2} log(n) k'w exists so that weight does not
  ## default to the largest model. Pushed to an extreme it must dominate.
  set.seed(3)
  A <- matrix(rnorm(50 * 3), 50, 3)
  r <- rnorm(50)
  w <- sfa:::.simplex_qp(A, r, c(0, 1e6, 1e6))
  expect_equal(w[1], 1, tolerance = 1e-6)
  expect_equal(sum(w[2:3]), 0, tolerance = 1e-6)
})

sfma_data <- function(seed = 3, n = 400, u = c("exp", "hn")) {
  u <- match.arg(u)
  set.seed(seed)
  x1 <- rnorm(n); x2 <- rnorm(n)
  ui <- if (u == "exp") rexp(n, 1) else abs(rnorm(n, 0, 1))
  data.frame(
    y = 1 + 0.5 * x1 + 0.5 * x2 + rnorm(n, 0, 0.4) - ui,
    x1 = x1, x2 = x2
  )
}

test_that("weights are a probability vector for every scheme", {
  skip_on_cran()
  d <- sfma_data()
  ## NTN is included so that "sfma" has the unique most general candidate its
  ## criterion requires; the other schemes do not care either way.
  for (wt in c("sfma", "tic", "aic", "bic", "equal")) {
    a <- sfma(y ~ x1 + x2, data = d, models = c("NHN", "NE", "NTN"),
      weights = wt, quiet = TRUE
    )
    expect_s3_class(a, "sfma")
    expect_equal(sum(a$weights), 1, tolerance = 1e-8)
    expect_true(all(a$weights >= 0))
    expect_length(a$weights, 3L)
    ## Averaged efficiency is a convex combination, so it cannot escape the
    ## range of the candidates it averages.
    E <- exp(-a$u_by_model)
    expect_true(all(a$efficiency <= apply(E, 1, max) + 1e-8))
    expect_true(all(a$efficiency >= apply(E, 1, min) - 1e-8))
  }
})

test_that("equal weighting really is equal", {
  skip_on_cran()
  d <- sfma_data()
  a <- sfma(y ~ x1 + x2, data = d, models = c("NHN", "NE", "NR"),
    weights = "equal", quiet = TRUE
  )
  expect_equal(unname(a$weights), rep(1 / 3, 3))
  expect_equal(a$u_hat, unname(rowMeans(a$u_by_model)))
})

test_that("the true distribution earns the weight", {
  skip_on_cran()
  ## The substantive claim, on data whose true inefficiency distribution is
  ## known. NHN/NE/NR all have k = 5, so "sfma" is refused here by design (see
  ## the degeneracy test below) and the information-criterion schemes are the
  ## ones on trial.
  d_exp <- sfma_data(u = "exp")
  for (wt in c("tic", "aic")) {
    a <- sfma(y ~ x1 + x2, data = d_exp, models = c("NHN", "NE", "NR"),
      weights = wt, quiet = TRUE
    )
    expect_identical(names(which.max(a$weights)), "NE")
  }

  ## And the other way, so a scheme that always answered "NE" could not pass.
  ## Compared against NE rather than NR: Rayleigh and half-normal are close in
  ## shape and at n = 400 the data does not reliably separate them (on this
  ## seed NR attains the higher log-likelihood, -392.30 against -392.48, so
  ## preferring it is correct behaviour, not a defect). NE is the candidate
  ## that is genuinely wrong here.
  d_hn <- sfma_data(u = "hn")
  for (wt in c("tic", "aic")) {
    a <- sfma(y ~ x1 + x2, data = d_hn, models = c("NHN", "NE", "NR"),
      weights = wt, quiet = TRUE
    )
    expect_gt(a$weights[["NHN"]], a$weights[["NE"]])
  }
})

test_that("the sfma criterion is REFUSED when the candidates tie on k", {
  skip_on_cran()
  ## The degeneracy that would otherwise produce a confident wrong answer:
  ## with equal dimensions the penalty is constant, the criterion collapses to
  ## "be closest to rho_full", and all the weight goes to whichever model was
  ## nominated as full. On exponential data that gave NHN 0.879 while TIC gave
  ## NE 0.934.
  d <- sfma_data(u = "exp")
  expect_error(
    sfma(y ~ x1 + x2, data = d, models = c("NHN", "NE", "NR"),
      weights = "sfma", quiet = TRUE),
    "tie at 5 parameters|degenerates"
  )
  ## With a genuinely more general candidate present it works.
  a <- sfma(y ~ x1 + x2, data = d, models = c("NHN", "NE", "NTN"),
    weights = "sfma", quiet = TRUE
  )
  expect_identical(a$criterion$full_model, "NTN")
})

test_that("the SFMA criterion starves the over-fitted model", {
  skip_on_cran()
  ## Their Theorem 1: weights on models that nest the truth go to zero. NTN
  ## nests NHN and carries an extra parameter, so on half-normal data it is
  ## the over-fitted candidate and should be squeezed out -- which is exactly
  ## what the penalty term is there to do.
  d <- sfma_data(u = "hn", n = 500)
  a <- sfma(y ~ x1 + x2, data = d, models = c("NHN", "NE", "NTN"),
    weights = "sfma", quiet = TRUE
  )
  expect_false(is.null(a$criterion))
  expect_identical(a$criterion$full_model, "NTN")
  expect_lt(a$weights[["NTN"]], 0.5)
})

test_that("it reports the disagreement that motivates averaging", {
  skip_on_cran()
  d <- sfma_data()
  a <- sfma(y ~ x1 + x2, data = d, models = c("NHN", "NE", "NR"),
    weights = "tic", quiet = TRUE
  )
  expect_equal(dim(a$u_by_model), c(nrow(d), 3L))
  expect_equal(a$nobs, nrow(d))
  expect_output(print(a), "spread across candidates")
  expect_output(print(a), "weight")
})

test_that("the refusals name the reason", {
  skip_on_cran()
  d <- sfma_data()
  expect_error(sfma(y ~ x1 + x2, data = d, models = "NHN"), "at least two")
  expect_error(sfma(y ~ x1 + x2, data = d, models = 1:3), "at least two")
  ## Two names that cannot both fit leaves nothing to average.
  expect_error(
    suppressMessages(sfma(y ~ x1 + x2, data = d,
      models = c("NOT_A_MODEL", "ALSO_NOT"), quiet = TRUE)),
    "fewer than two"
  )
})
