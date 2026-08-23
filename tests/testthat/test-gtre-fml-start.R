## GTRE_FML starting values, and the two-step estimator they now draw on.
##
## Background: the FIML surface for the four-component GTRE carries a boundary
## optimum at sigma_h = 0, where the model collapses to TRE and the intercept
## absorbs the missing E[h] = sigma_h*sqrt(2/pi). Colombi (2010, sec. 3) and
## Colombi, Martini and Vittadini (2011, sec. 3.2) both recommend seeding the
## FIML search from the two-step moment estimates for exactly this reason.

test_that(".gtre_two_step reproduces the inline algebra GTRE_SEQ2 used to carry", {
  ## The moment inversion was lifted out of psfm()'s GTRE_SEQ2 branch verbatim.
  ## This pins the refactor: same numbers, to machine precision.
  set.seed(4)
  eps <- stats::rnorm(400) - abs(stats::rnorm(400))
  alp <- stats::rnorm(80)  - abs(stats::rnorm(80))
  b0  <- 1.25

  k      <- sqrt(pi/2)*(pi/(pi - 4))
  eps_2m <- mean(eps^2); eps_3m <- min(0, mean(eps^3))
  alp_2m <- mean(alp^2); alp_3m <- min(0, mean(alp^3))
  want <- list(
    gamma_uv   = min(1, 1/(eps_2m*(k*eps_3m)^(-2/3) + (2/pi))),
    sigmaSq_uv = eps_2m + (2/pi)*(k*eps_3m)^(2/3),
    gamma_hr   = min(1, 1/(alp_2m*(k*alp_3m)^(-2/3) + (2/pi))),
    sigmaSq_hr = alp_2m + (2/pi)*(k*alp_3m)^(2/3),
    beta_0     = b0 + sqrt(2/pi)*(k*eps_3m)^(1/3) + sqrt(2/pi)*(k*alp_3m)^(1/3))

  got <- sfa:::.gtre_two_step(eps, alp, b0)
  expect_equal(got[names(want)], want, tolerance = 1e-14)
})

test_that("the two-step split recovers the variance components it was given", {
  ## sigma_u = sqrt(gamma_uv * sigmaSq_uv) etc. On a large clean sample the
  ## moment estimator should land near the truth -- this is the mapping the
  ## GTRE_FML start relies on, so a sign or reciprocal slip would show here.
  ##
  ## The inputs are CENTRED first, which is the estimator's precondition: it
  ## reads raw second and third moments as central ones. That holds for the
  ## residuals and ranef() it is fed in psfm() but not for a raw draw of
  ## v - u, whose mean is -sigma_u*sqrt(2/pi) = -0.8 here. Skipping the
  ## centring saturates gamma_uv at 1 and drives sigma_v to 0.
  set.seed(11)
  n <- 2e5
  su <- 1; sv <- 0.3; sh <- 0.4; sr <- 0.2
  eps <- stats::rnorm(n, 0, sv) - abs(stats::rnorm(n, 0, su))
  alp <- stats::rnorm(n, 0, sr) - abs(stats::rnorm(n, 0, sh))
  ts  <- sfa:::.gtre_two_step(eps - mean(eps), alp - mean(alp), 0)
  expect_equal(sqrt(ts$sigmaSq_uv*ts$gamma_uv),       su, tolerance = 0.05)
  expect_equal(sqrt(ts$sigmaSq_uv*(1 - ts$gamma_uv)), sv, tolerance = 0.10)
  expect_equal(sqrt(ts$sigmaSq_hr*ts$gamma_hr),       sh, tolerance = 0.10)
  expect_equal(sqrt(ts$sigmaSq_hr*(1 - ts$gamma_hr)), sr, tolerance = 0.15)
})

test_that("wrong-skew samples are truncated rather than producing complex values", {
  ## m3 is truncated at 0 because a positively skewed sample carries no
  ## inefficiency signal, and (k*m3)^(1/3) with k < 0 would otherwise go
  ## complex. Feed it a deliberately wrong-skew sample.
  set.seed(5)
  z  <- abs(stats::rnorm(500))              ## positive skew
  ts <- sfa:::.gtre_two_step(z, z, 0)
  expect_true(all(vapply(ts, function(v) is.finite(v) && !is.complex(v), logical(1))))
  expect_equal(ts$beta_0, 0)                ## both shifts vanish
})

test_that("GTRE_FML records which start it chose, and never a worse one", {
  skip_on_cran()
  d <- as.data.frame(data_gen_p(t = 5, N = 60, rand = 9, sig_u = 1, sig_v = 0.3,
                                sig_r = 0.2, sig_h = 0.4, cons = 0.5,
                                beta1 = 0.5, beta2 = 0.5, eta = 0.1))
  f <- suppressWarnings(psfm(y_gtre ~ x1 + x2, model_name = "GTRE_FML",
                             data = d, individual = "name"))
  expect_equal(f$start_search$n_tried, 2L)
  expect_true(f$start_search$chosen %in% c("random-effects", "two-step"))
  ## The chosen candidate must be the better of the two on the likelihood.
  expect_equal(max(f$start_search$loglik),
               f$start_search$loglik[match(f$start_search$chosen,
                                           c("random-effects","two-step"))])
})

test_that("an explicit start_val still overrides both candidates", {
  skip_on_cran()
  d <- as.data.frame(data_gen_p(t = 5, N = 60, rand = 9, sig_u = 1, sig_v = 0.3,
                                sig_r = 0.2, sig_h = 0.4, cons = 0.5,
                                beta1 = 0.5, beta2 = 0.5, eta = 0.1))
  f <- suppressWarnings(psfm(y_gtre ~ x1 + x2, model_name = "GTRE_FML",
                             data = d, individual = "name",
                             start_val = c(0.5, 0.5, 0.5, 0.2, 0.3, 0.4, 1)))
  expect_null(f$start_search)
  expect_equal(length(coef(f)), 7L)
})

test_that("psfm() exposes maxit.nlminb and passes it through", {
  ## It used to be hard-coded (200 in the GTRE_FML branch, 500 elsewhere) with
  ## no way to raise it from the call.
  expect_true("maxit.nlminb" %in% names(formals(psfm)))
  expect_equal(eval(formals(psfm)$maxit.nlminb), 500)
  src <- deparse(body(psfm))
  expect_false(any(grepl("maxit\\.nlminb\\s*=\\s*[0-9]", src)))
})
