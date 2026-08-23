## tHN -- Student-t noise with a half-normal inefficiency term.
##
## Distinct from THT: there a single scale mixture is shared by both error
## components, so both are t with the same degrees of freedom and the composed
## error is a closed-form skew-t. Here v ~ sigma_v*t_nu and u ~ |N(0,sigma_u^2)|
## are independent with different tails, there is no closed form, and the
## density is a quadrature.

d_nhn <- function(e, sv, su){
  s <- sqrt(sv^2 + su^2)
  2/s*stats::dnorm(e/s)*stats::pnorm(-e*(su/sv)/s)
}

## Reference convolution with the panel boundaries placed AT the t peak, so the
## adaptive rule cannot miss a narrow spike. Plain integrate(0, Inf) does miss
## it once sigma_v is small, which is why this helper exists.
d_ref <- function(e, sv, su, nu) vapply(e, function(ei){
  f   <- function(u) stats::dt((ei + u)/sv, df = nu)/sv *
                     sqrt(2/pi)/su*exp(-(u/su)^2/2)
  brk <- sort(unique(pmax(0, c(0, -ei + sv*c(-40,-8,-2,0,2,8,40), su*c(.5,1,2,4,8)))))
  sum(vapply(seq_along(brk), function(k)
    stats::integrate(f, brk[k], if (k == length(brk)) Inf else brk[k + 1],
                     rel.tol = 1e-11, subdivisions = 2000)$value, numeric(1)))
}, numeric(1))

test_that("tHN density integrates to one", {
  ## Integrate over (-Inf, Inf), not a wide finite range. On a finite range the
  ## adaptive rule can miss the bulk entirely -- integrate(-200, 200) returns
  ## 0.983 at sigma_v = 0.05, sigma_u = 1 because the density lives in about
  ## [-5, 2] and the initial subdivision never looks there. That is a property
  ## of integrate(), not of the density.
  for (p in list(c(0.3, 1, 5), c(0.1, 0.3, 3), c(0.05, 1, 10), c(0.3, 1, 200),
                 c(1, 1, 2.05))) {
    m <- stats::integrate(function(e) exp(sfa:::.log_d_thn(e, p[1], p[2], p[3])),
                          -Inf, Inf)$value
    expect_equal(m, 1, tolerance = 1e-6)
  }
})

test_that("tHN quadrature stays accurate as lambda grows", {
  ## THE regression test for this model. A FIXED node count (the reference
  ## implementation used 96) is fine at lambda ~ 3 and badly wrong later: 4e-2
  ## relative error at lambda = 20 and 6e-1 at lambda = 62 -- and the model
  ## really does visit lambda in the 60s on real data. .thn_m_for() therefore
  ## scales the node count with lambda.
  ##
  ## Note that integrating the density to 1 does NOT catch this: the error
  ## redistributes across e and integrates away, so total mass still reads
  ## 1.000 while the density is 40% wrong pointwise. Only a pointwise check
  ## against an independent reference finds it.
  skip_on_cran()
  for (lam in c(0.5, 2, 5, 12, 20, 40, 70)) {
    su <- 1; sv <- su/lam
    for (nu in c(2.05, 5, 30)) {
      ee <- c(-4, -2, -1, -0.3, 0, 0.4, 1, 3)*max(sv, su)
      r  <- d_ref(ee, sv, su, nu)
      v  <- exp(sfa:::.log_d_thn(ee, sv, su, nu))
      ok <- is.finite(r) & r > 1e-290
      expect_lt(max(abs(v[ok]/r[ok] - 1)), 1e-4)
    }
  }
})

test_that("tHN node count is bounded and survives a degenerate lambda", {
  ## sigma_v -> 0 arrives as Inf/NaN. Before it was clamped, the node-count
  ## arithmetic overflowed the integer range and as.integer() returned NA,
  ## which surfaced far away as "missing value where TRUE/FALSE needed".
  ## Non-finite lambda is clamped to the safe minimum rather than the cap:
  ## the density is refused outright above .THN_LAMBDA_MAX, so there is no
  ## point paying for 4096 nodes on the way there.
  expect_equal(sfa:::.thn_m_for(NaN), 96L)
  expect_equal(sfa:::.thn_m_for(Inf), 96L)
  expect_equal(sfa:::.thn_m_for(0), 96L)
  expect_lte(sfa:::.thn_m_for(1e18), 4096L)
  expect_gt(sfa:::.thn_m_for(60), sfa:::.thn_m_for(3))
  ## and the degenerate region is refused rather than quadratured
  expect_true(all(sfa:::.log_d_thn(c(-1, 0, 1), 1e-4, 1, 5) == -1e6))
  expect_true(all(is.finite(sfa:::.log_d_thn(c(-1, 0, 1), 1e-300, 1, 5))))
})

test_that("tHN converges to normal-half-normal as nu grows", {
  ## The strongest available check: NHN is already in the package, so it is a
  ## free ground truth. Convergence is O(1/nu) and must be compared where the
  ## data are -- in the far tail the t is polynomial and the normal is
  ## exponential, so the relative gap there is large by construction, not by
  ## error.
  sv <- 0.3; su <- 1
  ee <- seq(-3.2, 0.35, length.out = 25)          # central ~98% of v - u
  prev <- Inf
  for (nu in c(30, 200, 1e4, 1e6)) {
    err <- max(abs(exp(sfa:::.log_d_thn(ee, sv, su, nu))/d_nhn(ee, sv, su) - 1))
    expect_lt(err, prev)                          # monotone convergence
    prev <- err
  }
  expect_lt(prev, 1e-5)
})

test_that("tHN efficiency converges to Battese-Coelli as nu grows", {
  bc <- function(e, sv, su){
    s2 <- sv^2 + su^2; mu <- -e*su^2/s2; ss <- sqrt(sv^2*su^2/s2)
    exp(-mu + ss^2/2)*stats::pnorm(mu/ss - ss)/stats::pnorm(mu/ss)
  }
  ee <- c(-3, -1.5, -0.5, 0, 0.5, 1.5)
  expect_equal(sfa:::.thn_eff(ee, 0.3, 1, 1e6)$exp_u_hat, bc(ee, 0.3, 1),
               tolerance = 1e-4)
  ## and it is a genuine efficiency: in [0,1], and u_hat non-negative
  eff <- sfa:::.thn_eff(ee, 0.3, 1, 5)
  expect_true(all(eff$exp_u_hat >= 0 & eff$exp_u_hat <= 1))
  expect_true(all(eff$u_hat >= 0))
})

test_that("data_gen_cs() provides a tHN column with independent t and half-normal", {
  d <- data_gen_cs(N = 300, rand = 1, sig_u = 1, sig_v = 0.3, cons = 0.5,
                   beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.1)
  expect_true(all(c("v_thn", "u_thn", "y_pcs_thn") %in% names(d)))
  expect_true(all(d$u_thn >= 0))
  expect_equal(d$y_pcs_thn, 0.5 + 0.5*d$x1 + 0.5*d$x2 + d$v_thn - d$u_thn)
  ## must NOT be the THT columns -- different model, different draws
  expect_false(isTRUE(all.equal(d$v_thn, d$v_t)))
  expect_false(isTRUE(all.equal(d$u_thn, d$u_t)))
})

test_that("tHN parameter labels match the likelihood's x[] positions", {
  ## THT's labels were once transposed relative to its likelihood, so every fit
  ## reported each scale under the other name until a convergence sweep caught
  ## it. tHN is (sigma_v, sigma_u, nu) -- conventional order, NOT THT's.
  ## Asserted by moving ONE parameter at a time and checking the log-likelihood
  ## responds the way that parameter should.
  skip_on_cran()
  set.seed(9)
  e  <- stats::rnorm(400, 0, 0.3) - abs(stats::rnorm(400, 0, 1))
  ll <- function(sv, su, nu) sum(sfa:::.log_d_thn(e, sv, su, nu))
  ## truth is sigma_v = 0.3, sigma_u = 1: the likelihood must prefer those to
  ## the swapped assignment, which is what a transposition would produce.
  expect_gt(ll(0.3, 1.0, 5), ll(1.0, 0.3, 5))
  ## and the fitted object must label them in that order
  d   <- data_gen_cs(N = 150, rand = 5, sig_u = 1, sig_v = 0.3, cons = 0.5,
                     beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.1)
  fit <- suppressWarnings(sfm(y_pcs_thn ~ x1 + x2, data = d, model_name = "tHN"))
  ## $out is stored transposed (3 columns par/st_err/t-val), so the parameter
  ## names are its ROWNAMES.
  expect_identical(colnames(fit$out), c("par", "st_err", "t-val"))
  expect_identical(rownames(fit$out)[1:3], c("sigv", "sigu", "nu"))
  expect_identical(names(coef(fit))[1:3], c("sigv", "sigu", "nu"))
})

test_that("tHN fit returns efficiency, diagnostics and working S3 methods", {
  skip_on_cran()
  d   <- data_gen_cs(N = 150, rand = 5, sig_u = 1, sig_v = 0.3, cons = 0.5,
                     beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.1)
  fit <- suppressWarnings(sfm(y_pcs_thn ~ x1 + x2, data = d, model_name = "tHN"))
  expect_s3_class(fit, "sfareg")
  ## the no-efficiency-predictor fallback must not have overwritten the object
  expect_false(is.null(fit$exp_u_hat))
  expect_length(fit$exp_u_hat, nrow(d))
  expect_true(all(fit$exp_u_hat >= 0 & fit$exp_u_hat <= 1, na.rm = TRUE))
  expect_true(all(fit$u_hat >= 0, na.rm = TRUE))
  ## multi-start diagnostic and boundary flag are both reported
  expect_false(is.null(fit$thn_starts))
  expect_true(is.numeric(fit$thn_starts$n_distinct))
  expect_true(is.logical(fit$thn_sigma_u_at_bound))
  ## standard generics
  expect_length(coef(fit), nrow(fit$out))
  expect_true(is.finite(as.numeric(logLik(fit))))
  expect_equal(nobs(fit), nrow(d))
  expect_true(is.finite(AIC(fit)))
  expect_length(residuals(fit), nrow(d))
  expect_length(fitted(fit), nrow(d))
  expect_output(print(fit))
})
