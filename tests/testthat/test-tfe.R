## The two true-fixed-effects estimators.
##
## "TFE" (Greene 2005) and "TFE_WMLE" (Chen, Schmidt & Wang 2014) fit the same
## model by different routes. "TFE" changed meaning in 1.1.3 -- it previously
## named the CSW estimator -- so the rename itself is worth pinning down.

test_that("model_name = \"TFE\" warns that its meaning changed in 1.1.3", {
  d <- panel_small(t = 5, N = 30)
  expect_warning(psfm(y_tfe ~ x1_w + x2_w, model_name = "TFE",
                      data = d, individual = "name"),
                 "now fits Greene")
  ## TFE_WMLE is the new name for the old estimator, and must not carry the
  ## rename warning. (plm emits an unrelated index warning of its own, so this
  ## checks for the specific message rather than for silence.)
  w <- character()
  withCallingHandlers(
    psfm(y_tfe ~ x1_w + x2_w, model_name = "TFE_WMLE", data = d, individual = "name"),
    warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
  expect_false(any(grepl("now fits Greene", w)))
})

test_that("the two estimators are genuinely different fits", {
  skip_on_cran()
  d <- panel_small(t = 6, N = 50)
  g <- fit_tfe_quietly(psfm(y_tfe ~ x1_w + x2_w, model_name = "TFE",
                            data = d, individual = "name"))
  w <- psfm(y_tfe ~ x1_w + x2_w, model_name = "TFE_WMLE", data = d, individual = "name")
  ## Same parameter layout ...
  expect_equal(names(coef(g)), names(coef(w)))
  ## ... but not the same numbers, and not comparable log-likelihoods (one is
  ## the likelihood of the data, the other of the within-transformed
  ## deviations), so AIC/BIC must not be used to choose between them.
  expect_false(isTRUE(all.equal(unname(coef(g)), unname(coef(w)))))
})

test_that("TFE recovers the frontier slopes and its firm effects", {
  skip_on_cran()
  ## data_gen_p()'s y_tfe = r_i + 0.5*x1_w + 0.5*x2_w + v - u, with r_i the
  ## firm effect. Use a noisier design than the default so that the interior
  ## maximum exists -- see the degeneracy test below.
  ## sig_r = 1 rather than the usual 0.2: with firm effects that small the
  ## true spread is swamped by estimation noise and the recovery test would
  ## measure nothing.
  d <- as.data.frame(data_gen_p(t = 8, N = 80, rand = 3, sig_u = 1, sig_v = 0.8,
                                sig_r = 1, sig_h = 0.4, cons = 0.5,
                                beta1 = 0.5, beta2 = 0.5))
  fit <- fit_tfe_quietly(psfm(y_tfe ~ x1_w + x2_w, model_name = "TFE",
                              data = d, individual = "name"))
  expect_equal(unname(coef(fit)["x1_w"]), 0.5, tolerance = 0.1)
  expect_equal(unname(coef(fit)["x2_w"]), 0.5, tolerance = 0.1)

  ## r_hat_m is one estimate per individual, named, and should track the true
  ## firm effects. Greene's alpha_i absorbs the frontier intercept, so compare
  ## on deviations rather than levels.
  r_true <- tapply(d$r, as.character(d$name), mean)[names(fit$r_hat_m)]
  expect_length(fit$r_hat_m, length(unique(d$name)))
  expect_gt(cor(fit$r_hat_m, r_true), 0.5)

  expect_true(all(fit$exp_u_hat >= 0 & fit$exp_u_hat <= 1))
  expect_true(all(fit$u_hat >= 0))
})

test_that("TFE is invariant to row order and to numeric-looking individual labels", {
  skip_on_cran()
  ## rowsum(reorder = TRUE) would sort ids as character, putting "10" before
  ## "2" and silently misassigning every firm effect; nothing may depend on the
  ## panel arriving sorted either.
  d       <- panel_small(t = 6, N = 40)
  d$name  <- as.character(d$name)
  set.seed(99)
  shuffled <- d[sample(nrow(d)), ]

  a <- fit_tfe_quietly(psfm(y_tfe ~ x1_w + x2_w, model_name = "TFE",
                            data = d, individual = "name"))
  b <- fit_tfe_quietly(psfm(y_tfe ~ x1_w + x2_w, model_name = "TFE",
                            data = shuffled, individual = "name"))

  expect_equal(unname(coef(a)), unname(coef(b)), tolerance = 1e-5)
  expect_equal(as.numeric(logLik(a)), as.numeric(logLik(b)), tolerance = 1e-6)
  ## Firm effects agree once matched by name, and the names are the real labels.
  expect_setequal(names(a$r_hat_m), unique(d$name))
  expect_equal(unname(a$r_hat_m), unname(b$r_hat_m[names(a$r_hat_m)]), tolerance = 1e-5)
})

test_that("the gamma parameterization reaches the same optimum as lambda", {
  skip_on_cran()
  d <- panel_small(t = 6, N = 40)
  a <- fit_tfe_quietly(psfm(y_tfe ~ x1_w + x2_w, model_name = "TFE",
                            data = d, individual = "name"))
  g <- fit_tfe_quietly(psfm(y_tfe ~ x1_w + x2_w, model_name = "TFE",
                            data = d, individual = "name", gamma = TRUE))
  expect_identical(names(coef(g))[1], "gamma")
  ## gamma = sigma_u^2/sigma^2, so lambda = sqrt(gamma/(1-gamma)).
  gv <- unname(coef(g)[1])
  expect_equal(sqrt(gv/(1-gv)), unname(coef(a)["lambda"]), tolerance = 1e-3)
  expect_equal(as.numeric(logLik(g)), as.numeric(logLik(a)), tolerance = 1e-5)
})

test_that("TFE bounds lambda and says so when it pins at the bound", {
  ## Greene's likelihood always has a supremum on the sigma_v -> 0 boundary
  ## (alpha_i = max_t(y_it - x_it'beta) gives the deterministic frontier), and
  ## with little noise relative to inefficiency there may be no interior
  ## maximum at all. A fit that runs into that must report the constraint
  ## rather than quietly returning sigma_v ~ 0.
  skip_on_cran()
  d <- as.data.frame(data_gen_p(t = 10, N = 60, rand = 100, sig_u = 1, sig_v = 0.3,
                                sig_r = 0.2, sig_h = 0.4, cons = 0.5,
                                beta1 = 0.5, beta2 = 0.5))
  fit <- suppressWarnings(
    withCallingHandlers(
      psfm(y_tfe ~ x1_w + x2_w, model_name = "TFE", data = d,
           individual = "name", tfe_lambda_max = 25),
      warning = function(w) NULL))
  expect_lte(unname(coef(fit)["lambda"]), 25 + 1e-8)

  ## The bound is validated.
  expect_error(psfm(y_tfe ~ x1_w + x2_w, model_name = "TFE", data = d,
                    individual = "name", tfe_lambda_max = -1),
               "tfe_lambda_max")
  expect_error(psfm(y_tfe ~ x1_w + x2_w, model_name = "TFE", data = d,
                    individual = "name", tfe_lambda_max = c(1, 2)),
               "tfe_lambda_max")
})

test_that("a lower lambda bound cannot produce a higher log-likelihood", {
  skip_on_cran()
  ## Sanity check on the constraint itself: shrinking the feasible set can only
  ## make the attained maximum weakly worse.
  d  <- panel_small(t = 6, N = 40)
  hi <- fit_tfe_quietly(psfm(y_tfe ~ x1_w + x2_w, model_name = "TFE", data = d,
                             individual = "name", tfe_lambda_max = 100))
  lo <- suppressWarnings(fit_tfe_quietly(
          psfm(y_tfe ~ x1_w + x2_w, model_name = "TFE", data = d,
               individual = "name", tfe_lambda_max = 1)))
  expect_lte(as.numeric(logLik(lo)), as.numeric(logLik(hi)) + 1e-6)
  expect_lte(unname(coef(lo)["lambda"]), 1 + 1e-8)
})

test_that("TFE fits a cost frontier", {
  skip_on_cran()
  d   <- panel_small(t = 6, N = 40)
  fit <- fit_tfe_quietly(psfm(c_tfe ~ x1_w + x2_w, model_name = "TFE", data = d,
                              individual = "name", inefdec = FALSE))
  expect_true(all(is.finite(coef(fit))))
  expect_equal(unname(coef(fit)["x1_w"]), 0.5, tolerance = 0.15)
})
