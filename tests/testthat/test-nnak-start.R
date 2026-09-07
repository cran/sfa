## NNAK's starting values.
##
## sigma_u -> 0 is a genuine attractor for the normal-Nakagami likelihood, and
## the old hard-coded start of sigma_u = sigma_v = 0.1 sat right next to it.
## On two of twelve samples at n = 3000 with a true sigma_u of 1 it converged
## to sigma_u = 0.0013 and 0.0000 -- inefficiency vanishing altogether -- for a
## log-likelihood 44 and 45 points worse than the moment start reaches. These
## tests pin the construction, one representative rescue, and the fallback.

test_that("the NNAK start is taken from the data, not from constants", {
  skip_on_cran()
  d <- cs_small(N = 400)
  ## Read the start from start_cs() itself: the fitted object's $start_v is
  ## overwritten by the three-stage scaffold and is not the initial value.
  sc <- sfa:::start_cs(y_pcs_nak ~ x1 + x2, d, c("(Intercept)", "x1", "x2"),
                       1, "NNAK", 3, FALSE, 0, NULL)

  sv <- sc$start_v
  expect_length(sv, 6)
  ## The two scales must be the moment estimates, not the old 0.1 constants.
  expect_false(isTRUE(all.equal(sv[1], 0.1)))
  expect_false(isTRUE(all.equal(sv[2], 0.1)))
  expect_true(all(sv[1:2] > 0))
  ## The shape still starts at 0.5, as it always did and as FronPy does.
  expect_equal(sv[3], 0.5)
})

test_that("the NNAK start matches the half-normal moment fit and shifts the intercept", {
  skip_on_cran()
  d <- cs_small(N = 400)
  sc <- sfa:::start_cs(y_pcs_nak ~ x1 + x2, d, c("(Intercept)", "x1", "x2"),
                       1, "NNAK", 3, FALSE, 0, NULL)

  ## Reconstruct the intended construction independently of start_cs().
  ols <- lm(y_pcs_nak ~ x1 + x2, data = d)
  e <- residuals(ols) - mean(residuals(ols))
  m2 <- mean(e^2)
  m3 <- mean(e^3)
  su <- (m3 / (sqrt(2 / pi) * (1 - 4 / pi)))^(1 / 3)
  sv_expect <- sqrt(m2 - su^2 * (1 - 2 / pi))
  skip_if(!is.finite(su) || su <= 0, "residuals wrong-skewed in this fixture")

  expect_equal(sc$start_v[1], sv_expect, tolerance = 1e-8)
  expect_equal(sc$start_v[2], su, tolerance = 1e-8)
  ## The intercept starts at the OLS intercept plus E[u], not at the raw one:
  ## OLS is biased downward by E[u] because the composed error has non-zero
  ## mean, and starting uncorrected leaves the frontier a whole E[u] too low.
  expect_equal(sc$start_v[4], unname(coef(ols)[1]) + su * sqrt(2 / pi),
               tolerance = 1e-8)
  expect_gt(sc$start_v[4], unname(coef(ols)[1]))
})

test_that("the moment start rescues a sample where the constant start collapsed", {
  skip_on_cran()
  ## rand = 701 is one of the two samples that drove sigma_u to ~0 under the
  ## old start. The truth is sigma_u = 1, m = 1, beta0 = 0.5.
  d <- data_gen_cs(N = 3000, rand = 701, sig_u = 1, sig_v = 0.3, cons = 0.5,
                   beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.5)
  fit <- sfm(y_pcs_nak ~ x1 + x2, model_name = "NNAK", data = d)
  p <- fit$out[, "par"]

  ## The failure mode being guarded against is sigma_u collapsing to zero.
  expect_gt(unname(p[2]), 0.5)
  expect_equal(unname(p[2]), 1, tolerance = 0.35)
  ## and the frontier intercept landing below zero instead of near 0.5.
  expect_gt(unname(p[4]), 0.25)
})

test_that("a wrong-skewed sample falls back instead of erroring", {
  skip_on_cran()
  ## Cost orientation against a production-frontier column reverses the sign of
  ## the third moment, so the half-normal moment inversion has no admissible
  ## solution. start_cs() must fall back to the constants rather than propagate
  ## a NaN start into the optimizer.
  d <- cs_small(N = 200)
  d$flipped <- -d$y_pcs_nak
  sc <- sfa:::start_cs(flipped ~ x1 + x2, d, c("(Intercept)", "x1", "x2"),
                       1, "NNAK", 3, FALSE, 0, NULL)
  expect_true(all(is.finite(sc$start_v)))
  ## No admissible moment solution here, so it must be the old constants.
  expect_equal(sc$start_v[1], 0.1)
  expect_equal(sc$start_v[2], 0.1)
  expect_equal(sc$start_v[3], 0.5)
})
