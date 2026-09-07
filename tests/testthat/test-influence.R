## influence_sfa(): the influence function as a specification diagnostic (L7).
##
## The claims tested are the ones a user acts on: that d is invariant to
## reparameterisation, that the raw sup-norm is NOT, and that the
## self-standardised sensitivity separates a specification with bounded
## influence from one without.

inf_data <- function(seed = 3, n = 200) {
  set.seed(seed)
  x <- runif(n, 1, 10)
  data.frame(y = 1 + 0.5 * log(x) + rnorm(n, 0, 0.3) - abs(rnorm(n, 0, 0.6)),
    x = log(x))
}

test_that("it returns a coherent object and refuses a fit it cannot use", {
  skip_on_cran()
  d <- inf_data()
  fit <- sfm(y ~ x, data = d, model_name = "NHN", keep_objective = TRUE)
  i <- influence_sfa(fit)

  expect_s3_class(i, "sfa_influence")
  expect_equal(dim(i$influence), c(nrow(d), length(fit$coefficients)))
  expect_identical(colnames(i$influence), names(fit$coefficients))
  expect_length(i$d, nrow(d))
  expect_true(all(i$d >= 0))
  expect_true(all(is.finite(i$influence)))

  ## The scores sum to zero at the optimum, so the influence function must
  ## average to zero. This is the check that the score matrix is the fit's own.
  expect_equal(unname(colMeans(i$influence)), rep(0, ncol(i$influence)),
    tolerance = 1e-4)

  ## Without keep_objective there is no likelihood to differentiate.
  plain <- sfm(y ~ x, data = d, model_name = "NHN")
  expect_error(influence_sfa(plain), "keep_objective|likelihood")
  expect_error(influence_sfa(lm(y ~ x, data = d)), "sfareg")
})

test_that("scale = TRUE is exactly the factor n", {
  skip_on_cran()
  d <- inf_data()
  fit <- sfm(y ~ x, data = d, model_name = "NHN", keep_objective = TRUE)
  a <- influence_sfa(fit, scale = TRUE)
  b <- influence_sfa(fit, scale = FALSE)
  expect_equal(a$influence, b$influence * nrow(d))
  ## d is defined with the n factor either way, so it must not move.
  expect_equal(a$d, b$d)
  expect_true(a$scaled)
  expect_false(b$scaled)
})

test_that("the case-influence quadratic form is invariant to reparameterisation", {
  ## An algebraic property of d = n * s\'Vs, tested as one: under theta -> A theta
  ## the scores become s A^{-1} and V becomes A V A\', so the quadratic form is
  ## unchanged while ||IF|| is not. Checked on synthetic matrices so it is exact
  ## arithmetic rather than a statement about one optimiser run.
  set.seed(11)
  n <- 50; p <- 4
  G <- matrix(rnorm(n * p), n, p)
  M <- matrix(rnorm(p * p), p, p); V <- crossprod(M)
  A <- diag(c(1, 10, 100, 0.01))

  d1 <- rowSums((G %*% V) * G) * n
  d2 <- rowSums((G %*% solve(A) %*% (A %*% V %*% t(A))) * (G %*% solve(A))) * n
  expect_equal(d1, d2)
  ## The raw sup-norm moves under the same change.
  expect_false(isTRUE(all.equal(
    max(sqrt(rowSums((G %*% V)^2))),
    max(sqrt(rowSums((G %*% solve(A) %*% (A %*% V %*% t(A)))^2))))))
})

test_that("rescaling a regressor leaves d essentially alone", {
  skip_on_cran()
  ## The same property on a real fit. Loosely: the two fits agree to six digits
  ## in the parameters and the log-likelihood, but the Hessian and scores are
  ## differentiated numerically, and a 100-fold change in one parameter\'s scale
  ## moves the finite-difference error by a few percent. Exact invariance is
  ## asserted above, where it is exact.
  d <- inf_data()
  d2 <- d; d2$x <- d2$x * 100

  a <- influence_sfa(sfm(y ~ x, data = d,  model_name = "NHN",
    keep_objective = TRUE))
  b <- influence_sfa(sfm(y ~ x, data = d2, model_name = "NHN",
    keep_objective = TRUE))

  expect_equal(a$sensitivity_std, b$sensitivity_std, tolerance = 0.1)
  expect_equal(sqrt(max(a$d)), a$sensitivity_std)
  expect_gt(cor(a$d, b$d), 0.99)
})

test_that("the raw sup-norm is not comparable across models, and the standardised one is", {
  skip_on_cran()
  ## The substantive claim of Stead, Wheat and Greene: a Student's t noise term
  ## has bounded influence where the normal does not, so contamination moves it
  ## far less. Measured as a RESPONSE to contamination, not as a level, because
  ## no level here has an absolute meaning.
  d <- inf_data()
  d2 <- d; d2$y[1] <- d2$y[1] - 5

  s <- function(dd, m) {
    influence_sfa(suppressWarnings(
      sfm(y ~ x, data = dd, model_name = m, keep_objective = TRUE)))
  }
  nhn_c <- s(d, "NHN");  nhn_x <- s(d2, "NHN")
  thn_c <- s(d, "tHN");  thn_x <- s(d2, "tHN")

  ## Contamination inflates the half-normal fit's sensitivity much more.
  r_nhn <- nhn_x$sensitivity_std / nhn_c$sensitivity_std
  r_thn <- thn_x$sensitivity_std / thn_c$sensitivity_std
  expect_gt(r_nhn, 2)
  expect_lt(r_thn, r_nhn)

  ## On clean data the two specifications are close on the standardised scale
  ## and wildly apart on the raw one -- the comparison the print method warns
  ## against making.
  expect_lt(abs(log(thn_c$sensitivity_std / nhn_c$sensitivity_std)), 0.5)
  expect_gt(thn_c$sensitivity / nhn_c$sensitivity, 5)

  ## The contaminated observation is the one the diagnostic points at.
  expect_identical(which.max(nhn_x$d), 1L)
})

test_that("print reports both scales and does not promise a smaller number", {
  skip_on_cran()
  d <- inf_data()
  i <- influence_sfa(sfm(y ~ x, data = d, model_name = "NHN",
    keep_objective = TRUE))
  expect_output(print(i), "self-standardised sensitivity")
  expect_output(print(i), "raw sup")
  expect_output(print(i), "Most influential observations")
  expect_output(print(i), "not comparable across models")
})
