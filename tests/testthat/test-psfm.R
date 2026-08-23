## Panel entry point.

test_that("every psfm() model_name fits and returns a well-formed sfareg", {
  skip_on_cran()
  d <- panel_small(t = 6, N = 50)
  specs <- list(
    list("TFE",           y_tfe  ~ x1_w + x2_w),
    list("TFE_WMLE",      y_tfe  ~ x1_w + x2_w),
    list("SSFE",          y_ssfe ~ x1   + x2),
    list("FD",            y_fd   ~ x_fd | z_fd),
    list("TRE",           y_tre  ~ x1   + x2),
    list("GTRE",          y_gtre ~ x1   + x2),
    list("PL80",          y_ssfe ~ x1   + x2),
    list("BC92",          y_bc92 ~ x1   + x2),
    list("K1990",         y_bc92 ~ x1   + x2),
    list("K1990modified", y_bc92 ~ x1   + x2)
  )
  for (s in specs) {
    ## GTRE routes through `estimator` since 1.1.4 and defaults to "fiml";
    ## pin it so this loop keeps exercising the simulated-ML path (GTRE_FML
    ## has its own coverage below).
    extra <- if (identical(s[[1]], "GTRE")) list(estimator = "sml") else list()
    fit <- fit_tfe_quietly(do.call(psfm, c(
      list(s[[2]], model_name = s[[1]], data = d, individual = "name"), extra
    )))
    expect_s3_class(fit, "sfareg")
    expect_identical(fit$model_name, s[[1]])
    expect_true(all(is.finite(fit$coefficients)), info = s[[1]])
    expect_equal(unname(fit$coefficients), unname(fit$out[, "par"]), info = s[[1]])
  }
})

test_that("psfm() accepts a plain data.frame as well as a pdata.frame", {
  skip_on_cran()
  skip_if_not_installed("plm")
  d  <- panel_small()
  pd <- plm::pdata.frame(d, index = c("name", "year"))
  a  <- psfm(y_ssfe ~ x1 + x2, model_name = "PL80", data = d,  individual = "name")
  b  <- psfm(y_ssfe ~ x1 + x2, model_name = "PL80", data = pd, individual = "name")
  expect_equal(unname(a$coefficients), unname(b$coefficients), tolerance = 1e-6)
})

test_that("PL80 and BC92 recover their true parameters", {
  skip_on_cran()
  ## y_ssfe uses time-invariant inefficiency (sig_u = 1), which is exactly
  ## PL80's assumption; y_bc92 adds eta = 0.1 decay.
  d <- panel_small(t = 8, N = 150)
  p <- psfm(y_ssfe ~ x1 + x2, model_name = "PL80", data = d, individual = "name")
  expect_equal(unname(p$coefficients["x1"]), 0.5, tolerance = 0.05)
  expect_equal(unname(p$coefficients["x2"]), 0.5, tolerance = 0.05)

  b <- psfm(y_bc92 ~ x1 + x2, model_name = "BC92", data = d, individual = "name")
  expect_equal(unname(b$coefficients["x1"]), 0.5, tolerance = 0.05)
  ## The decay parameter is reported as "time" (frontier's naming); the
  ## generator uses eta = 0.1.
  expect_equal(unname(b$coefficients["time"]), 0.1, tolerance = 0.06)
})

test_that("PL80 is the eta = 0 special case of BC92", {
  skip_on_cran()
  ## Both share one likelihood, differing only in B_it. On a time-invariant
  ## DGP, BC92's eta should be near zero and its fit close to PL80's.
  d <- panel_small(t = 8, N = 120)
  p <- psfm(y_ssfe ~ x1 + x2, model_name = "PL80", data = d, individual = "name")
  b <- psfm(y_ssfe ~ x1 + x2, model_name = "BC92", data = d, individual = "name")
  expect_lt(abs(unname(b$coefficients["time"])), 0.05)
  expect_equal(unname(b$coefficients["x1"]), unname(p$coefficients["x1"]), tolerance = 0.02)
  ## BC92 nests PL80, so it cannot fit worse.
  expect_gte(as.numeric(logLik(b)), as.numeric(logLik(p)) - 1e-4)
})

test_that("K1990 and K1990modified nest PL80's constant path", {
  skip_on_cran()
  d  <- panel_small(t = 8, N = 100)
  p  <- psfm(y_ssfe ~ x1 + x2, model_name = "PL80",  data = d, individual = "name")
  for (mn in c("K1990", "K1990modified")) {
    k <- psfm(y_ssfe ~ x1 + x2, model_name = mn, data = d, individual = "name")
    expect_true(all(is.finite(k$coefficients)), info = mn)
    expect_gte(as.numeric(logLik(k)), as.numeric(logLik(p)) - 1e-4)
  }
})

test_that("SSFE is deterministic and reports no log-likelihood", {
  d   <- panel_small()
  fit <- psfm(y_ssfe ~ x1 + x2, model_name = "SSFE", data = d, individual = "name")
  expect_null(fit$opt)
  ## Not an MLE, so logLik()/AIC()/BIC() must decline rather than invent a value.
  expect_warning(ll <- logLik(fit))
  expect_true(is.na(ll))
  expect_length(fit$alpha_hat, length(unique(d$name)))
  expect_true(all(fit$exp_u_hat > 0 & fit$exp_u_hat <= 1))
  ## Deterministic: refitting gives bit-identical numbers.
  again <- psfm(y_ssfe ~ x1 + x2, model_name = "SSFE", data = d, individual = "name")
  expect_identical(fit$coefficients, again$coefficients)
})

test_that("psfm() reports collinear between-individual designs rather than failing opaquely", {
  skip_on_cran()
  ## Time dummies are estimable in the pooled design but collapse onto the
  ## intercept once averaged within each individual, which is what breaks the
  ## random-effects starting-value regression.
  d          <- panel_small(t = 5, N = 40)
  d$year_fac <- factor(d$year)
  f          <- y_tre ~ x1 + x2 + year_fac

  expect_error(psfm(f, model_name = "TRE", data = d, individual = "name",
                    collinear_action = "error"), "collinear|rank")
  expect_warning(psfm(f, model_name = "TRE", data = d, individual = "name",
                      collinear_action = "start_only"))
})

test_that("psfm() rejects an unknown model_name", {
  d <- panel_small()
  expect_error(psfm(y_tre ~ x1 + x2, model_name = "NOPE", data = d, individual = "name"))
})

test_that("pipe-count validation rejects a formula the model cannot use", {
  d <- panel_small()
  ## TFE takes no pipe segments; TRE requires the explicit _Z name to use one.
  expect_error(psfm(y_tfe ~ x1_w + x2_w | z_fd, model_name = "TFE",
                    data = d, individual = "name"))
  expect_error(psfm(y_tre ~ x1 + x2 | z_gtre, model_name = "TRE",
                    data = d, individual = "name"))
})

## --- case-insensitive model names, and GTRE's estimator argument ------------

test_that("model_name is matched without regard to case", {
  skip_on_cran()
  d <- data_gen_cs(N = 150, rand = 3, sig_u = 1, sig_v = 0.3, cons = 0.5,
                   beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.5)
  vals <- vapply(c("NHN", "nhn", "Nhn", "nHn"), function(m)
    suppressWarnings(sfm(y_pcs ~ x1 + x2, model_name = m, data = d))$opt$value,
    numeric(1))
  expect_length(unique(signif(vals, 10)), 1L)
})

test_that("the THT / tHN pair is not confused by case folding", {
  expect_equal(sfa:::.match_model_name("tht", c("THT", "tHN")), "THT")
  expect_equal(sfa:::.match_model_name("THN", c("THT", "tHN")), "tHN")
  expect_equal(sfa:::.match_model_name("tHn", c("THT", "tHN")), "tHN")
})

test_that("an exact match beats a partial one", {
  ch <- c("TRE_Z", "GTRE_Z", "TRE", "GTRE", "GTRE_FML", "GTRE_SEQ1", "GTRE_SEQ2")
  ## "GTRE" prefixes four other names; it must resolve to itself
  expect_equal(sfa:::.match_model_name("gtre", ch), "GTRE")
  ## partial matching survives, where the partial is unambiguous ...
  expect_equal(sfa:::.match_model_name("GTRE_F", ch), "GTRE_FML")
  ## ... and an ambiguous one is still rejected rather than guessed at
  expect_error(sfa:::.match_model_name("GTRE_S", ch), "not a recognized choice")
})

test_that("an unrecognized name errors helpfully", {
  ch <- c("NHN", "NHN_Z", "NE")
  expect_error(sfa:::.match_model_name("nhnz", ch), "Did you mean")
  expect_error(sfa:::.match_model_name("nhnz", ch), "case does not matter")
  expect_error(sfa:::.match_model_name("xyzzy", ch), "not a recognized choice")
  expect_error(sfa:::.match_model_name(NA, ch), "must not be NA")
})

test_that("the unevaluated default vector still selects the first choice", {
  ch <- c("NHN", "NE", "NR")
  expect_equal(sfa:::.match_model_name(ch, ch), "NHN")
})

test_that("GTRE defaults to fiml and warns about the change", {
  skip_on_cran()
  pd <- panel_small(t = 6, N = 40)
  expect_warning(
    f <- psfm(y_gtre ~ x1 + x2, model_name = "GTRE", data = pd, individual = "name"),
    "defaults to estimator = \"fiml\""
  )
  expect_equal(f$model_name, "GTRE_FML")
})

test_that("estimator selects among the GTRE routes, and old names still work", {
  skip_on_cran()
  pd <- panel_small(t = 6, N = 40)
  f_sml <- psfm(y_gtre ~ x1 + x2, model_name = "GTRE", data = pd,
                individual = "name", estimator = "sml")
  expect_equal(f_sml$model_name, "GTRE")

  f_fiml <- psfm(y_gtre ~ x1 + x2, model_name = "GTRE", data = pd,
                 individual = "name", estimator = "fiml")
  f_old <- psfm(y_gtre ~ x1 + x2, model_name = "GTRE_FML", data = pd,
                individual = "name")
  ## the new route and the old name must be the same estimator
  expect_equal(f_fiml$coefficients, f_old$coefficients, tolerance = 1e-8)

  ## explicit estimator suppresses the default-change warning
  expect_silent(suppressMessages(
    psfm(y_gtre ~ x1 + x2, model_name = "GTRE", data = pd,
         individual = "name", estimator = "sml")))
})

test_that("estimator is refused where it has no meaning", {
  skip_on_cran()
  pd <- panel_small(t = 6, N = 40)
  expect_warning(
    psfm(y_tre ~ x1 + x2, model_name = "TRE", data = pd,
         individual = "name", estimator = "sml"),
    "only applies to model_name = \"GTRE\""
  )
})

test_that("an unbalanced panel falls back from fiml, but only when implicit", {
  skip_on_cran()
  pd <- panel_small(t = 6, N = 40)
  set.seed(1)
  ub <- pd[-sample(which(pd$name %in% unique(pd$name)[1:8]), 12), ]

  ## implicit default: warn and fall back rather than making GTRE unusable
  w <- capture_warnings(
    f <- psfm(y_gtre ~ x1 + x2, model_name = "GTRE", data = ub, individual = "name"))
  expect_true(any(grepl("unbalanced", w)))
  expect_equal(f$model_name, "GTRE")

  ## explicit choice: error rather than silently fitting something else
  expect_error(
    psfm(y_gtre ~ x1 + x2, model_name = "GTRE", data = ub,
         individual = "name", estimator = "fiml"),
    "requires a BALANCED panel"
  )
})
