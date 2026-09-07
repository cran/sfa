## simulation_se(): how much of a standard error is simulation noise (J4).

ss_data <- function() {
  data_gen_cs(N = 400, rand = 1, sig_u = 1, sig_v = 0.3, cons = 0.5,
    beta1 = 0.5, beta2 = 0.5, a = 1, mu = 0.5
  )
}

test_that("simulation noise is small for the slopes and larger for the scales", {
  skip_on_cran()
  d <- ss_data()
  f <- suppressWarnings(sfm(y_pcs_wb ~ x1 + x2, model_name = "NW", data = d, Nsim = 200))
  r <- suppressMessages(simulation_se(f, K = 6))

  expect_s3_class(r, "data.frame")
  expect_identical(r$parameter, names(f$coefficients))
  expect_true(all(r$simulation_se >= 0))
  expect_equal(attr(r, "seed_arg"), "sim_seed")
  expect_gte(attr(r, "n_refits"), 2L)

  ## Independently reproduces a documented property of these estimators: too
  ## few draws move everything EXCEPT the frontier slopes. So the slopes should
  ## carry markedly less simulation noise, relative to their sampling error,
  ## than the variance and shape parameters do.
  slope <- r$ratio[r$parameter %in% c("x1", "x2")]
  scale <- r$ratio[r$parameter %in% c("sigv", "sigu", "k")]
  expect_lt(max(slope), min(scale))
})

test_that("more draws mean less simulation noise", {
  skip_on_cran()
  ## The property that makes the number meaningful. If it did not fall with
  ## Nsim it would be measuring something else.
  d <- ss_data()
  lo <- suppressWarnings(sfm(y_pcs_wb ~ x1 + x2, model_name = "NW", data = d, Nsim = 200))
  hi <- suppressWarnings(sfm(y_pcs_wb ~ x1 + x2, model_name = "NW", data = d, Nsim = 800))
  a <- suppressMessages(simulation_se(lo, K = 6))
  b <- suppressMessages(simulation_se(hi, K = 6))
  expect_lt(median(b$simulation_se / a$simulation_se), 0.6)
})

test_that("it refuses fits that do no simulation", {
  d <- ss_data()
  g <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = d)
  ## A closed-form fit has no simulation variance. Reporting a column of zeros
  ## would look like a measurement; saying so is the honest answer.
  expect_error(simulation_se(g), "not estimated by simulation")

  f <- suppressWarnings(sfm(y_pcs_wb ~ x1 + x2, model_name = "NW", data = d, Nsim = 100))
  expect_error(simulation_se(f, K = 1), "at least two refits")
  expect_error(simulation_se(list()), "must be an \"sfareg\" fit")
})
