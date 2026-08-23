## Cross-sectional entry point.
##
## The structural tests run everywhere; parameter recovery is behind
## skip_on_cran() because it needs a sample size that costs real time.

test_that("every sfm() model_name fits and returns a well-formed sfareg", {
  skip_on_cran()
  d <- cs_small(N = 250)
  specs <- list(
    NHN   = y_pcs     ~ x1 + x2,
    NR    = y_pcs     ~ x1 + x2,
    NE    = y_pcs_e   ~ x1 + x2,
    NTN   = y_pcs_tn  ~ x1 + x2,
    THT   = y_pcs_t   ~ x1 + x2,
    NU    = y_pcs     ~ x1 + x2,
    NGE   = y_pcs_e   ~ x1 + x2,
    NHN_Z = y_pcs_z   ~ x1 + x2 | z,
    NE_Z  = y_pcs_ez  ~ x1 + x2 | z
  )
  for (mn in names(specs)) {
    fit <- sfm(specs[[mn]], model_name = mn, data = d)
    expect_s3_class(fit, "sfareg")
    expect_identical(fit$model_name, mn)
    expect_true(all(is.finite(fit$coefficients)), info = mn)
    expect_true(is.finite(as.numeric(logLik(fit))), info = mn)
    ## out is the source of truth; the three vectors are its rows.
    expect_equal(unname(fit$coefficients), unname(fit$out[, "par"]), info = mn)
    expect_equal(nrow(fit$out), length(fit$coefficients), info = mn)
  }
})

test_that("NHN recovers its true parameters", {
  skip_on_cran()
  ## data_gen_cs(sig_u = 1, sig_v = 0.3) => lambda = 3.333, sigma = 1.044,
  ## beta = (0.5, 0.5, 0.5). See DATA_GENERATION_REFERENCE.md.
  fit <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = cs_small(N = 2000))
  cf  <- fit$coefficients
  ## NHN reports the (lambda, sigma) reparameterization, not raw sigmas.
  expect_equal(unname(cf["lambda"]), 1/0.3,          tolerance = 0.25)
  expect_equal(unname(cf["sigma"]),  sqrt(1 + 0.09), tolerance = 0.10)
  expect_equal(unname(cf["x1"]),     0.5,            tolerance = 0.05)
  expect_equal(unname(cf["x2"]),     0.5,            tolerance = 0.05)
})

test_that("NE recovers its true parameters", {
  skip_on_cran()
  ## NE reports raw sigv/sigu rather than the lambda/sigma reparameterization.
  fit <- sfm(y_pcs_e ~ x1 + x2, model_name = "NE", data = cs_small(N = 2000))
  cf  <- fit$coefficients
  expect_equal(unname(cf["sigv"]), 0.3, tolerance = 0.15)
  expect_equal(unname(cf["x1"]),   0.5, tolerance = 0.10)
  expect_equal(unname(cf["x2"]),   0.5, tolerance = 0.10)
})

test_that("the likelihood sign convention is right (a fit beats its own start)", {
  ## Every likelihood closure in this package must return the NEGATIVE summed
  ## log-likelihood, because bobyqa/psoptim/optim all minimize and none of the
  ## opts.R wrappers flip the sign. A flipped sign still "converges" and still
  ## returns believable-looking numbers -- it just maximizes nothing. The
  ## observable consequence is that the fitted parameters score WORSE than the
  ## starting values, which is what this checks.
  skip_on_cran()
  d   <- cs_small(N = 400)
  fit <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = d)
  X   <- cbind(1, d$x1, d$x2)
  nll <- function(p) {
    lam <- p[1]; sig <- p[2]; b <- p[-c(1,2)]
    e   <- as.numeric(d$y_pcs - X %*% b)
    -sum(log(2) - log(sig) + dnorm(e/sig, log = TRUE) +
           pnorm(-lam*e/sig, log.p = TRUE))
  }
  at_fit   <- nll(fit$coefficients)
  at_start <- nll(fit$start_v)
  expect_lte(at_fit, at_start + 1e-6)
  ## And the reported logLik really is the maximized log-likelihood.
  expect_equal(as.numeric(logLik(fit)), -at_fit, tolerance = 1e-4)
})

test_that("efficiency predictions are probabilities and track true inefficiency", {
  skip_on_cran()
  d   <- cs_small(N = 800)
  fit <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = d)
  expect_length(fit$exp_u_hat, nrow(d))
  expect_true(all(fit$exp_u_hat >= 0 & fit$exp_u_hat <= 1))
  ## Predicted efficiency should be strongly rank-correlated with exp(-u_true).
  expect_gt(cor(fit$exp_u_hat, exp(-d$u), method = "spearman"), 0.8)
})

test_that("inefdec = FALSE fits the cost orientation", {
  skip_on_cran()
  ## data_gen_cs() has no cost-oriented column (unlike data_gen_p(), which
  ## supplies c_gtre/c_tre/c_tfe/c_pcs), so build one from the returned error
  ## draws: a cost frontier composes as v + u rather than v - u.
  d      <- cs_small(N = 800)
  d$cost <- d$cons + 0.5*d$x1 + 0.5*d$x2 + d$v + d$u
  fit    <- sfm(cost ~ x1 + x2, model_name = "NHN", data = d, inefdec = FALSE)
  expect_true(all(is.finite(fit$coefficients)))
  expect_equal(unname(fit$coefficients["x1"]),     0.5,            tolerance = 0.05)
  expect_equal(unname(fit$coefficients["x2"]),     0.5,            tolerance = 0.05)
  expect_equal(unname(fit$coefficients["lambda"]), 1/0.3,          tolerance = 0.30)
  expect_equal(unname(fit$coefficients["sigma"]),  sqrt(1 + 0.09), tolerance = 0.10)
})

test_that("sfm() rejects an unknown model_name and a bad robust request", {
  d <- cs_small(N = 60)
  expect_error(sfm(y_pcs ~ x1 + x2, model_name = "NOT_A_MODEL", data = d))
  ## robust is deliberately NHN-only for now; it must say so rather than
  ## silently ignoring the argument.
  expect_error(sfm(y_pcs_e ~ x1 + x2, model_name = "NE", data = d, robust = "mdpd"),
               "only implemented for model_name")
})

test_that("robust divergence estimators run and stay close to MLE on clean data", {
  skip_on_cran()
  d <- cs_small(N = 500)
  mle <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = d)
  for (r in c("mlqe", "psi", "mdpd")) {
    fit <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = d, robust = r)
    expect_s3_class(fit, "sfareg")
    expect_true(all(is.finite(fit$coefficients)), info = r)
    ## No contamination present, so a robust fit should not wander far.
    expect_equal(unname(fit$coefficients["x1"]),
                 unname(mle$coefficients["x1"]), tolerance = 0.15, info = r)
  }
})
