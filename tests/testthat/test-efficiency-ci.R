## Horrace and Schmidt (1996) intervals for individual inefficiency.
##
## The arithmetic is checked against an independent construction rather than
## against itself: the truncated-normal quantile is inverted directly, and in
## the far tail -- where the log-domain code exists -- against the exponential
## limit the truncated normal converges to. truncnorm::qtruncnorm() is
## deliberately NOT used as the reference: it underflows to a constant
## 3.1e-61 for both bounds once mu_star/sigma_star drops below about -6, which
## is exactly the regime worth testing.

test_that(".horrace_schmidt_ci inverts the truncated normal", {
  hs <- sfa:::.horrace_schmidt_ci
  for (ms in c(-1, -0.2, 0, 0.5, 2)) {
    for (ss in c(0.2, 1)) {
      r <- hs(ms, ss, level = 0.90)
      ## Independent inversion: F(u) = [Phi((u-m)/s) - Phi(-m/s)] / Phi(m/s).
      cdf <- function(u) {
        (pnorm((u - ms) / ss) - pnorm(-ms / ss)) / pnorm(ms / ss)
      }
      expect_equal(cdf(r$lower), 0.05, tolerance = 1e-6)
      expect_equal(cdf(r$upper), 0.95, tolerance = 1e-6)
    }
  }
})

test_that("the bounds stay finite where a naive implementation loses them", {
  hs <- sfa:::.horrace_schmidt_ci
  ## A very efficient unit drives mu_star/sigma_star strongly negative, so
  ## Phi(mu_star/sigma_star) underflows and 1 - c*A rounds to exactly 1.
  r <- hs(c(-5, -10, -20, -40), 0.5, level = 0.95)
  expect_true(all(is.finite(r$lower)))
  expect_true(all(is.finite(r$upper)))
  expect_true(all(r$lower >= 0))
  expect_true(all(r$upper > r$lower))
  ## As mu_star/sigma_star -> -Inf the posterior tends to an Exponential with
  ## rate |mu_star|/sigma_star^2, whose quantiles are known exactly. The
  ## approximation tightens as the ratio grows, so check the most extreme one.
  rate <- 40 / 0.5^2
  expect_equal(r$lower[4], -log1p(-0.025) / rate, tolerance = 1e-3)
  expect_equal(r$upper[4], -log(0.025) / rate, tolerance = 1e-3)
})

test_that("a wider level gives a wider interval, and level is validated", {
  hs <- sfa:::.horrace_schmidt_ci
  narrow <- hs(0.3, 0.8, level = 0.80)
  wide <- hs(0.3, 0.8, level = 0.99)
  expect_lt(narrow$upper - narrow$lower, wide$upper - wide$lower)
  expect_error(hs(0, 1, level = 1), "strictly between 0 and 1")
  expect_error(hs(0, 1, level = 0), "strictly between 0 and 1")
  expect_error(hs(0, 1, level = NA_real_), "strictly between 0 and 1")
})

test_that("efficiency_ci brackets the point predictor for every supported model", {
  skip_on_cran()
  d <- cs_small(N = 200)
  fits <- list(
    NHN = sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = d),
    NE = sfm(y_pcs ~ x1 + x2, model_name = "NE", data = d),
    NTN = sfm(y_pcs_tn ~ x1 + x2, model_name = "NTN", data = d)
  )
  for (nm in names(fits)) {
    ci <- efficiency_ci(fits[[nm]], level = 0.90)
    expect_equal(nrow(ci), nrow(d), info = nm)
    expect_named(ci, c("u_lower", "u_hat", "u_upper",
                       "te_lower", "te_hat", "te_upper"))
    expect_equal(attr(ci, "level"), 0.90)
    ## The JLMS predictor is the posterior MEAN, so it lies inside the
    ## posterior interval -- this is what ties the new output to the old one.
    expect_true(all(ci$u_lower <= ci$u_hat + 1e-8), info = nm)
    expect_true(all(ci$u_hat <= ci$u_upper + 1e-8), info = nm)
    ## Efficiency bounds are the inefficiency bounds mapped through exp(-u),
    ## which is decreasing, so the endpoints swap.
    expect_equal(ci$te_lower, exp(-ci$u_upper), tolerance = 1e-12)
    expect_equal(ci$te_upper, exp(-ci$u_lower), tolerance = 1e-12)
    expect_true(all(ci$te_lower >= 0 & ci$te_upper <= 1), info = nm)
  }
})

test_that("NHN_Z gets a per-observation posterior scale", {
  skip_on_cran()
  d <- cs_small(N = 200)
  fit <- sfm(y_pcs_z ~ x1 + x2 | z, model_name = "NHN_Z", data = d)
  ci <- efficiency_ci(fit)
  expect_equal(nrow(ci), nrow(d))
  expect_length(fit$u_posterior$mu_star, nrow(d))
  expect_length(fit$u_posterior$sigma_star, nrow(d))

  ## sigma_star must be built from the PER-OBSERVATION sigma_u rather than
  ## recycled from a single number. Test that by recomputing it from the fit's
  ## own z_spec, which is exact and does not depend on where the optimizer
  ## landed.
  ##
  ## The previous version asserted instead that sigma_star takes more than one
  ## distinct value at 10 decimal places. That is a much weaker claim than it
  ## looks: sigma_star = sigma_u sigma_v / sqrt(sigma_u^2 + sigma_v^2) is
  ## dominated by whichever scale is smaller, so it barely moves even when
  ## sigma_u moves a lot -- sd(sigma_star) is 2.1e-05 here. On CI's Windows the
  ## fit landed with just enough less variation to collapse it to one value at
  ## that rounding, and the test failed while nothing was wrong with the code.
  zs <- fit$z_spec
  eta <- as.numeric(zs$Z %*% zs$delta)
  sig_u <- if (identical(zs$link, "sd")) exp(eta) else sqrt(exp(eta))
  sig_v <- unname(coef(fit)[["sigma_v"]])
  expect_equal(fit$u_posterior$sigma_star,
               sig_u * sig_v / sqrt(sig_u^2 + sig_v^2),
               tolerance = 1e-8)
  ## A recycled scalar would only survive that if sigma_u were itself
  ## constant, so pin that it is not.
  expect_gt(stats::sd(sig_u), 0)
})

test_that("type selects the columns returned", {
  skip_on_cran()
  fit <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = cs_small(N = 150))
  expect_named(efficiency_ci(fit, type = "u"), c("u_lower", "u_hat", "u_upper"))
  expect_named(efficiency_ci(fit, type = "te"), c("te_lower", "te_hat", "te_upper"))
})

test_that("models without a truncated-normal posterior are refused by name", {
  skip_on_cran()
  ## NR: u ~ Rayleigh, whose posterior is not a truncated normal. (NGE would
  ## do as well statistically, but it fails to converge at this sample size.)
  fit <- sfm(y_pcs_r ~ x1 + x2, model_name = "NR", data = cs_small(N = 150))
  expect_null(fit$u_posterior)
  ## The refusal has to say which model it is and why, not just fail.
  expect_error(efficiency_ci(fit), "NR")
  expect_error(efficiency_ci(fit), "truncated normal")
  expect_error(efficiency_ci(structure(list(), class = "lm")), "sfareg")
})

test_that("the point predictors come from the fit, not a recomputation", {
  skip_on_cran()
  fit <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = cs_small(N = 150))
  ci <- efficiency_ci(fit)
  ## te_hat must be the Battese-Coelli score the fit already reported, so the
  ## interval and print()/summary() cannot drift apart.
  expect_equal(ci$te_hat, as.numeric(fit$exp_u_hat), tolerance = 1e-12)
})
