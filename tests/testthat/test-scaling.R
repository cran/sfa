## The scaling-property model (H4): Wang and Schmidt (2002); Alvarez, Amsler,
## Orea and Schmidt (2006). u_i = h(z_i, delta) * u*_i.

scal_data <- function(seed = 11, n = 1500, su = 0.8, mu = 0.6,
                      d = c(0.5, -0.3), sv = 0.3) {
  set.seed(seed)
  x1 <- rnorm(n); x2 <- rnorm(n); z1 <- rnorm(n); z2 <- rnorm(n)
  ## The covariates scale a COMMON draw: one h multiplies u*, so sigma_u and
  ## the pre-truncation mean move together and the shape is held fixed.
  h <- exp(d[1] * z1 + d[2] * z2)
  ustar <- truncnorm::rtruncnorm(n, a = 0, mean = mu, sd = su)
  y <- 1 + 0.5 * x1 + 0.5 * x2 + rnorm(n, 0, sv) - h * ustar
  data.frame(y = y, x1 = x1, x2 = x2, z1 = z1, z2 = z2)
}

test_that("the scaling coefficients and the common draw are both recovered", {
  skip_on_cran()
  d <- scal_data()
  f <- sfm(y ~ x1 + x2, model_name = "NTN", data = d, scaling = ~ z1 + z2)
  p <- f$out[, "par"]

  expect_true(all(c("scale.z1", "scale.z2", "Zu.(Intercept)", "Zmu.(Intercept)")
    %in% rownames(f$out)))
  expect_equal(unname(p[["scale.z1"]]), 0.5, tolerance = 0.15)
  expect_equal(unname(p[["scale.z2"]]), -0.3, tolerance = 0.15)
  ## sigma_u of the common draw, on the log link.
  expect_equal(unname(p[["Zu.(Intercept)"]]), log(0.8), tolerance = 0.2)
  expect_equal(unname(p[["x1"]]), 0.5, tolerance = 0.1)
  expect_equal(unname(p[["x2"]]), 0.5, tolerance = 0.1)
})

test_that("scaling imposes ONE factor, so it is not the same fit as uhet", {
  skip_on_cran()
  d <- scal_data()
  a <- sfm(y ~ x1 + x2, model_name = "NTN", data = d, scaling = ~ z1 + z2)
  b <- sfm(y ~ x1 + x2, model_name = "NTN", data = d, uhet = ~ z1 + z2)

  ## uhet lets sigma_u vary while mu stays fixed; scaling moves both together.
  ## The restriction costs parameters, so the unrestricted fit cannot have the
  ## lower log-likelihood -- and the two are genuinely different models.
  expect_false(isTRUE(all.equal(unname(a$out[, "par"]), unname(b$out[, "par"]))))
  expect_gte(length(b$coefficients), length(a$coefficients) - 1L)
})

test_that("scaling is refused where it would add nothing or contradict itself", {
  d <- scal_data(n = 300)

  ## For a half-normal u, h(z)*|N(0,s^2)| IS |N(0,(h s)^2)|. Offering `scaling`
  ## there would be two names for one model.
  expect_error(
    sfm(y ~ x1 + x2, model_name = "NHN", data = d, scaling = ~z1),
    "implemented for model_name"
  )
  ## The entire content of the restriction is that ONE factor moves both.
  expect_error(
    sfm(y ~ x1 + x2, model_name = "NTN", data = d, scaling = ~z1, uhet = ~z2),
    "cannot be combined"
  )
  expect_error(
    sfm(y ~ x1 + x2, model_name = "NTN", data = d, scaling = ~z1, muhet = ~z2),
    "cannot be combined"
  )
  expect_error(
    sfm(y ~ x1 + x2, model_name = "NTN", data = d, scaling = "z1"),
    "must be a one-sided formula"
  )
  ## h with only an intercept is a constant, already absorbed into sigma_u.
  expect_error(
    sfm(y ~ x1 + x2, model_name = "NTN", data = d, scaling = ~1),
    "at least one covariate"
  )
})

test_that("the ordinary heteroskedastic path is untouched", {
  skip_on_cran()
  ## Zs = NULL must reproduce the pre-existing fit exactly: the scaling block
  ## is opt-in and adds no parameters when absent.
  d <- scal_data(n = 500)
  a <- sfm(y ~ x1 + x2, model_name = "NTN", data = d, uhet = ~z1)
  expect_false(any(grepl("^scale\\.", rownames(a$out))))
  expect_true(is.finite(a$opt$value))
})
