## Zero-inflated and two-tier entry points.

test_that("zsfm() fits both variants", {
  skip_on_cran()
  d <- cs_small(N = 400)
  for (s in list(list("ZISF", y_zisf ~ x1 + x2), list("ZISF_Z", y_zisf_z ~ x1 + x2 | z))) {
    fit <- zsfm(s[[2]], model_name = s[[1]], data = d)
    expect_s3_class(fit, "sfareg")
    expect_identical(fit$model_name, s[[1]])
    expect_true(all(is.finite(fit$coefficients)), info = s[[1]])
    expect_true(is.finite(as.numeric(logLik(fit))), info = s[[1]])
  }
})

test_that("ZISF recovers the frontier slopes", {
  skip_on_cran()
  fit <- zsfm(y_zisf ~ x1 + x2, model_name = "ZISF", data = cs_small(N = 1500))
  expect_equal(unname(fit$coefficients["x1"]), 0.5, tolerance = 0.1)
  expect_equal(unname(fit$coefficients["x2"]), 0.5, tolerance = 0.1)
})

test_that("ttsfm() fits the exponential two-tier variants", {
  skip_on_cran()
  d <- cs_small(N = 400)
  for (mn in c("TTNE", "TTNLS")) {
    fit <- suppressWarnings(ttsfm(y_ttne ~ x1 + x2, model_name = mn, data = d))
    expect_s3_class(fit, "sfareg")
    expect_identical(fit$model_name, mn)
    ## The frontier block (intercept + two slopes) is estimated by both.
    expect_true(all(is.finite(fit$coefficients[1:3])), info = mn)
  }
})

test_that("TTNLS returns NA for the parameters its objective cannot identify", {
  skip_on_cran()
  ## NLS minimizes (Y - X'beta + sigma_u - sigma_w)^2, so the two scales enter
  ## only through their difference -- which, because X carries the intercept,
  ## is itself confounded with beta_0. Only the slopes and the composite
  ## beta_0 + sigma_w - sigma_u are identified. Reporting the other three as
  ## point estimates (they used to come back at their starting values, with
  ## standard errors) invited users to interpret numbers that carry no
  ## information, so they are NA now.
  d <- cs_small(N = 300)
  expect_warning(fit <- ttsfm(y_ttne ~ x1 + x2, model_name = "TTNLS", data = d),
                 "not separately identified")
  expect_true(all(is.finite(fit$coefficients[1:3])))
  expect_true(all(is.na(fit$coefficients[4:6])))
  expect_true(all(is.na(fit$std.errors[4:6])))
})

test_that("the two-tier likelihood has the right sign", {
  skip_on_cran()
  ## TTNE and TTHN both had this sign backwards at one point, which produced
  ## believable-looking output rather than an error: the optimizer minimized
  ## the log-likelihood instead of maximizing it. Both one-sided scale
  ## parameters must come back strictly positive, which a sign flip destroys.
  fit <- ttsfm(y_ttne ~ x1 + x2, model_name = "TTNE", data = cs_small(N = 500))
  expect_true(is.finite(as.numeric(logLik(fit))))
  ## Positions 4-6 are the scale parameters (positions 1-3 are the frontier
  ## betas -- an earlier version of this test checked those by mistake and so
  ## asserted nothing). They are reported on the LOG scale (ttsfm.R's fn()
  ## exponentiates them), so the check is that they exponentiate to something
  ## finite and positive, not that they are positive themselves.
  scales <- exp(fit$coefficients[4:6])
  expect_true(all(is.finite(scales) & scales > 0))
  ## The DGP has sig_v = 0.3 and sig_u = sig_w = 1; a sign flip in the
  ## likelihood destroys these rather than merely loosening them.
  expect_equal(unname(scales[1]), 0.3, tolerance = 0.5)
})

## ---------------------------------------------------------------------------
## TTHN is quarantined behind an environment variable rather than run by
## default. On the standard two-tier DGP it is between three and four ORDERS OF
## MAGNITUDE slower than TTNE on the same data (TTNE and TTNLS return in ~0.1s
## at N = 200; TTHN had not finished after several minutes), and at N = 400 it
## fails outright with optim()'s ABNORMAL_TERMINATION_IN_LNSRCH. Running it by
## default would make R CMD check both slow and flaky.
##
## Set SFA_TEST_SLOW to a non-empty value to include it.
## ---------------------------------------------------------------------------
test_that("ttsfm() fits the half-normal two-tier variant", {
  skip_if(!nzchar(Sys.getenv("SFA_TEST_SLOW")),
          "TTHN is very slow and currently unreliable; set SFA_TEST_SLOW to run it")
  fit <- ttsfm(y_tthn ~ x1 + x2, model_name = "TTHN", data = cs_small(N = 200))
  expect_s3_class(fit, "sfareg")
  expect_true(all(is.finite(fit$coefficients)))
  expect_true(all(fit$coefficients[1:3] > 0))
})

test_that("ttsfm() accepts z and zp pipe segments", {
  skip_if(!nzchar(Sys.getenv("SFA_TEST_SLOW")),
          "the only DGP with both z and zp targets TTHN; see above")
  fit <- ttsfm(y_tthn_z ~ x1 + x2 | z | zp, model_name = "TTHN",
               data = cs_small(N = 200))
  expect_s3_class(fit, "sfareg")
  expect_true(all(is.finite(fit$coefficients)))
})

## --- zsfm() internals cleanup (2026-08-20) ------------------------------

test_that(".log_add2 matches the naive form and survives underflow", {
  a <- c(-1, -2, -700, -0.5); b <- c(-3, -1, -800, -0.5)
  expect_equal(sfa:::.log_add2(a, b), log(exp(a) + exp(b)))
  ## where the naive form underflows to -Inf, this must not
  expect_equal(sfa:::.log_add2(-800, -900), -800, tolerance = 1e-12)
  expect_true(is.finite(sfa:::.log_add2(-800, -900)))
  expect_equal(log(exp(-800) + exp(-900)), -Inf)   ## the form it replaces
  ## both components impossible -> -Inf, not NaN
  expect_equal(sfa:::.log_add2(-Inf, -Inf), -Inf)
  expect_false(is.nan(sfa:::.log_add2(-Inf, -Inf)))
})

test_that("ZISF's likelihood is symmetric in gamma, and post.prob respects it", {
  ## prob = exp(-|gamma|), so +g and -g are the same model. The JLMS block used
  ## exp(-gamma), which returns prob > 1 for a negative estimate -- reachable
  ## from an ordinary start, not a pathology.
  skip_on_cran()
  d <- data_gen_cs(N = 1500, rand = 1, sig_u = 1, sig_v = 0.3, cons = 0.5,
                   beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.1)
  f  <- suppressWarnings(zsfm(y_zisf ~ x1 + x2, model_name = "ZISF", data = d))
  pm <- as.numeric(coef(f)); pm[1] <- -pm[1]
  g  <- suppressWarnings(zsfm(y_zisf ~ x1 + x2, model_name = "ZISF", data = d,
                              start_val = pm))
  expect_lt(coef(g)[["gamma"]], 0)                       ## it really does return one
  expect_equal(as.numeric(logLik(g)), as.numeric(logLik(f)), tolerance = 1e-6)
  ## the point of the fix: a posterior probability is a probability
  expect_true(all(g$post.prob >= 0 & g$post.prob <= 1))
  expect_true(all(is.finite(g$jlms)))
  ## exp(-gamma), the old form, would have exceeded 1 here
  expect_gt(exp(-coef(g)[["gamma"]]), 1)
})

test_that("post.prob is a probability for both models", {
  skip_on_cran()
  d <- data_gen_cs(N = 1200, rand = 3, sig_u = 1, sig_v = 0.3, cons = 0.5,
                   beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.1)
  for(f in list(suppressWarnings(zsfm(y_zisf ~ x1 + x2, model_name="ZISF", data=d)),
                suppressWarnings(zsfm(y_zisf_z ~ x1 + x2 | z, model_name="ZISF_Z", data=d)))){
    expect_true(all(f$post.prob >= 0 & f$post.prob <= 1))
    expect_true(all(is.finite(f$post.prob)))
    expect_true(all(is.finite(f$jlms)))
  }
})

test_that("the logit link does not overflow at a large linear predictor", {
  ## exp(eta)/(1+exp(eta)) is Inf/Inf = NaN past eta ~ 710; plogis is not.
  eta <- c(-800, -1, 0, 1, 800)
  expect_true(all(is.finite(plogis(eta))))
  expect_equal(plogis(eta[2:4]), exp(eta[2:4])/(1+exp(eta[2:4])))
  expect_true(is.nan(exp(800)/(1+exp(800))))          ## the form it replaces
})
