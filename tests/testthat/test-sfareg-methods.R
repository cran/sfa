## The standard modelling generics on class "sfareg".

test_that("coef/vcov/logLik/nobs are mutually consistent", {
  skip_on_cran()
  d   <- cs_small(N = 300)
  fit <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = d)

  expect_equal(coef(fit), fit$out[, "par"])
  expect_equal(nobs(fit), nrow(d))
  expect_equal(dim(vcov(fit)), rep(length(coef(fit)), 2))
  expect_equal(rownames(vcov(fit)), names(coef(fit)))
  ## Standard errors are the square roots of the diagonal of vcov.
  expect_equal(unname(sqrt(diag(vcov(fit)))), unname(fit$std.errors), tolerance = 1e-8)
  ## t-values are par/se.
  expect_equal(unname(fit$t.values), unname(fit$coefficients/fit$std.errors), tolerance = 1e-8)

  ll <- logLik(fit)
  expect_s3_class(ll, "logLik")
  expect_equal(attr(ll, "df"), length(coef(fit)))
  expect_equal(attr(ll, "nobs"), nrow(d))
  ## AIC/BIC follow from logLik, so check they agree with the definition.
  expect_equal(AIC(fit), -2*as.numeric(ll) + 2*length(coef(fit)), tolerance = 1e-8)
  expect_equal(BIC(fit), -2*as.numeric(ll) + log(nrow(d))*length(coef(fit)), tolerance = 1e-8)
})

test_that("fitted + residuals reconstructs the response exactly", {
  skip_on_cran()
  d   <- cs_small(N = 300)
  fit <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = d)
  expect_length(fitted(fit), nrow(d))
  expect_length(residuals(fit), nrow(d))
  expect_equal(fitted(fit) + residuals(fit), d$y_pcs, tolerance = 1e-10, ignore_attr = TRUE)
})

test_that("predict() on the training data equals fitted()", {
  skip_on_cran()
  d   <- cs_small(N = 300)
  fit <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = d)
  expect_equal(unname(predict(fit)), unname(fitted(fit)), tolerance = 1e-10)
})

test_that("predict() honours newdata", {
  skip_on_cran()
  d   <- cs_small(N = 300)
  fit <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = d)
  nd  <- data.frame(x1 = c(0, 1, 2), x2 = c(1, 1, 1))
  p   <- predict(fit, newdata = nd)
  expect_length(p, 3)
  expect_true(all(is.finite(p)))
  ## The frontier is linear in x1, so equal steps in x1 give equal steps in p.
  expect_equal(diff(p)[1], diff(p)[2], tolerance = 1e-8)
  ## And the slope is the fitted x1 coefficient.
  expect_equal(unname(diff(p)[1]), unname(coef(fit)["x1"]), tolerance = 1e-8)
})

test_that("the generics work on panel fits too", {
  skip_on_cran()
  d   <- panel_small(t = 6, N = 40)
  fit <- fit_tfe_quietly(psfm(y_tfe ~ x1_w + x2_w, model_name = "TFE",
                              data = d, individual = "name"))
  expect_equal(nobs(fit), nrow(d))
  expect_true(is.finite(as.numeric(logLik(fit))))
  ## as.numeric() strips plm's pseries class off the stored response.
  expect_equal(as.numeric(fitted(fit) + residuals(fit)), as.numeric(d$y_tfe),
               tolerance = 1e-10)
})

test_that("nobs() works when the fit happens inside another function", {
  ## sfm() does not store the data on the fitted object, so nobs() has to
  ## re-evaluate the call's `data` argument. Resolving it against
  ## parent.frame() finds a global but not a local, which made BIC() -- which
  ## needs the "nobs" attribute -- return NA from inside any function.
  d <- cs_small(N = 200)
  fit_inside <- function(dat) sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = dat)
  local({
    dat <- d
    fit <- fit_inside(dat)
    expect_equal(nobs(fit), nrow(dat))
    expect_equal(attr(logLik(fit), "nobs"), nrow(dat))
    expect_true(is.finite(BIC(fit)))
  })
})

test_that("logLik() declines to answer for the non-likelihood estimators", {
  d <- panel_small()
  for (mn in c("SSFE", "GTRE_SEQ1", "GTRE_SEQ2")) {
    fit <- try(suppressWarnings(psfm(y_gtre ~ x1 + x2, model_name = mn,
                                     data = d, individual = "name")), silent = TRUE)
    if (inherits(fit, "try-error")) next
    expect_warning(ll <- logLik(fit))
    expect_true(is.na(ll), info = mn)
    expect_true(is.na(suppressWarnings(AIC(fit))), info = mn)
  }
})

test_that("print() and summary() run without error and return invisibly", {
  skip_on_cran()
  fit <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = cs_small(N = 200))
  expect_output(print(fit))
  expect_output(summary(fit))
  expect_invisible(print(fit))
})
