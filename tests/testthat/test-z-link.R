## The scale the variance-determinant linear predictor lives on.
##
## sfm() and ttsfm() historically put z'delta on the STANDARD DEVIATION,
## psfm() and every competitor put it on the VARIANCE. That is a factor of two
## in every delta and in every marginal effect, and it is a real trap when
## comparing a cross-sectional fit with a panel one. `z_link` lets a caller put
## them on the same footing; the default preserves sfm()'s own convention.

test_that("the default is unchanged and is the SD link", {
  skip_on_cran()
  d <- cs_small(N = 400)
  d$z <- runif(nrow(d), -1, 1)
  base <- suppressWarnings(sfm(y_pcs ~ x1 + x2 | z, model_name = "NHN_Z", data = d))
  sd_l <- suppressWarnings(sfm(y_pcs ~ x1 + x2 | z, model_name = "NHN_Z", data = d, z_link = "sd"))
  expect_equal(coef(base), coef(sd_l))
  expect_identical(base$z_spec$link, "sd")
})

test_that("the two links are reparameterizations, not different models", {
  skip_on_cran()
  d <- cs_small(N = 800)
  d$z <- runif(nrow(d), -1, 1)
  for (m in c("NHN_Z", "NE_Z")) {
    yv <- if (m == "NHN_Z") "y_pcs" else "y_pcs_e"
    f <- as.formula(paste(yv, "~ x1 + x2 | z"))
    a <- suppressWarnings(sfm(f, model_name = m, data = d, z_link = "sd"))
    b <- suppressWarnings(sfm(f, model_name = m, data = d, z_link = "var"))
    ## Same maximised likelihood: the model is the same, only delta's scale moves.
    expect_equal(as.numeric(logLik(a)), as.numeric(logLik(b)), tolerance = 1e-4, info = m)
    ## sigma_u = exp(eta) against sqrt(exp(eta)) = exp(eta/2): for the same
    ## fitted sigma_u the SD-link delta is HALF the variance-link delta.
    expect_equal(unname(coef(a)[["z"]]) / unname(coef(b)[["z"]]), 0.5, tolerance = 0.02, info = m)
    ## The frontier is untouched by the reparameterization.
    expect_equal(unname(coef(a)[["x1"]]), unname(coef(b)[["x1"]]), tolerance = 1e-3, info = m)
    expect_identical(b$z_spec$link, "var")
  }
})

test_that("marginal effects follow the link that was actually used", {
  skip_on_cran()
  d <- cs_small(N = 800)
  d$z <- runif(nrow(d), -1, 1)
  a <- suppressWarnings(sfm(y_pcs ~ x1 + x2 | z, model_name = "NHN_Z", data = d, z_link = "sd"))
  b <- suppressWarnings(sfm(y_pcs ~ x1 + x2 | z, model_name = "NHN_Z", data = d, z_link = "var"))
  ma <- marginal_effects(a, average = TRUE)
  mb <- marginal_effects(b, average = TRUE)
  ## delta halves and the derivative of sigma_u wrt z_k carries a compensating
  ## factor, so the MARGINAL EFFECT -- the thing with economic meaning -- is the
  ## same either way. If this ever fails, one of the two links is being
  ## differentiated with the other's formula.
  expect_equal(unname(ma[["dE_u.dz"]]), unname(mb[["dE_u.dz"]]), tolerance = 0.02)
})

test_that("an invalid link is refused", {
  d <- cs_small(N = 120)
  d$z <- runif(nrow(d), -1, 1)
  expect_error(sfm(y_pcs ~ x1 + x2 | z, model_name = "NHN_Z", data = d, z_link = "logit"),
               "should be one of")
})
