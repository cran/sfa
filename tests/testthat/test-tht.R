## THT -- the skew-t stochastic frontier of Tancredi (2002).
##
## Two independent defects lived here and neither announced itself:
##
##  1. The likelihood evaluated dt() at the RAW residual with no 1/omega
##     Jacobian, pinning the skew-t scale at 1. By Azzalini's lemma
##     2*f(e)*G(w(e)) is a valid density for any symmetric f and odd w, so the
##     wrong version still integrated to 1 and still produced believable fits.
##     sigma_u and sigma_v were left identified only through the skewing term,
##     and the degrees of freedom absorbed the scale mismatch.
##
##  2. data_gen_cs()'s y_pcs_t column drew the two error components as two
##     INDEPENDENT rt() variates. That shares the degrees of freedom but not
##     the mixing variable; Tancredi eq (5) needs a single common
##     Gamma(a/2, a/2). The composed error was therefore not skew-t at all.
##
## The tests below pin the density against an independent reference, and pin
## the efficiency predictor against its two known limits.

skew_t_logpdf <- function(e, sig_u, sig_v, a) {
  om <- sqrt(sig_v^2 + sig_u^2); z <- e/om
  log(2) - log(om) + stats::dt(z, df = a, log = TRUE) +
    stats::pt(-z*(sig_u/sig_v)*sqrt((a + 1)/(z^2 + a)), df = a + 1, log.p = TRUE)
}

test_that("THT log-density matches Tancredi (2002) eq (4) via sn::dst", {
  skip_if_not_installed("sn")
  e <- c(-3, -2, -1, -0.4, 0, 0.3, 1.1, 2.5)
  for (p in list(c(1, 0.3, 5), c(0.8, 0.4, 6), c(0.5, 0.5, 12))) {
    expect_equal(
      skew_t_logpdf(e, p[1], p[2], p[3]),
      sn::dst(e, xi = 0, omega = sqrt(p[2]^2 + p[1]^2),
              alpha = -p[1]/p[2], nu = p[3], log = TRUE),
      tolerance = 1e-10)
  }
})

test_that("THT density integrates to one and is NOT scale-1", {
  ## Deliberately sig_u = 2, sig_v = 0.8, so omega = 2.154 is well away from 1.
  ## At the convergence sweep's own settings (sig_u = 1, sig_v = 0.3) omega is
  ## 1.044, so dropping the scale is nearly a no-op there -- which is exactly
  ## why the defect survived so long and why this test must not use those
  ## values.
  su <- 2; sv <- 0.8; a <- 5; om <- sqrt(sv^2 + su^2)
  f <- function(e) exp(skew_t_logpdf(e, su, sv, a))
  expect_equal(stats::integrate(f, -Inf, Inf)$value, 1, tolerance = 1e-6)
  ## The old bug: dt() at the raw residual. Also integrates to 1 -- which is
  ## why it was invisible -- but is a different density. If these ever agree,
  ## the scale has been dropped again.
  wrong <- function(e)
    2*stats::dt(e, df = a) *
      stats::pt(-(e/om)*(su/sv)*sqrt((a + 1)/((e/om)^2 + a)), df = a + 1)
  expect_equal(stats::integrate(wrong, -Inf, Inf)$value, 1, tolerance = 1e-6)
  expect_gt(max(abs(log(f(c(-2, 0, 2))) - log(wrong(c(-2, 0, 2))))), 0.5)
})

test_that("data_gen_cs()'s y_pcs_st composed error is skew-t, y_pcs_t is not", {
  skip_on_cran()
  skip_if_not_installed("sn")
  d  <- data_gen_cs(N = 3e5, rand = 7, sig_u = 1, sig_v = 0.3, cons = 0.5,
                    beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.1)
  pr <- c(0.05, 0.25, 0.5, 0.75, 0.95, 0.99)
  q  <- sn::qst(pr, xi = 0, omega = sqrt(0.3^2 + 1^2), alpha = -1/0.3, nu = 5)
  err_st  <- max(abs(unname(stats::quantile(d$v_st - d$u_st, pr)) - q))
  err_old <- max(abs(unname(stats::quantile(d$v_t  - d$u_t,  pr)) - q))
  expect_lt(err_st, 0.02)      # shared mixing: Monte Carlo error only
  expect_gt(err_old, 0.05)     # independent draws: genuinely a different law
})

test_that("adding the skew-t columns left every earlier column untouched", {
  skip_on_cran()
  ## The new draws must stay at the END of data_gen_cs(). Inserting them
  ## earlier advances the RNG stream and silently changes every column after.
  d <- data_gen_cs(N = 200, rand = 1, sig_u = 1, sig_v = 0.3, cons = 0.5,
                   beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.1)
  expect_true(all(c("lam_st", "u_st", "v_st", "y_pcs_st") %in% names(d)))
  expect_true(all(d$u_st >= 0))
  expect_true(all(d$lam_st > 0))
  expect_equal(d$y_pcs_st, 0.5 + 0.5*d$x1 + 0.5*d$x2 + d$v_st - d$u_st)
})

test_that("THT conditional inefficiency reduces to Jondrow as a -> Inf", {
  ## Tancredi eq (7) is a Student-t truncated to [0, Inf) with df = a+1,
  ## location -e*su^2/om^2 and scale^2 (a + e^2/om^2)*sv^2*su^2/(om^2*(a+1)).
  ## As a -> Inf that must become Jondrow et al.'s N+(mu*, sigma*^2).
  su <- 1; sv <- 0.3; om2 <- sv^2 + su^2
  trt <- function(z, e, a) {
    mu <- -e*su^2/om2
    s  <- sqrt((a + e^2/om2)*sv^2*su^2/(om2*(a + 1)))
    stats::dt((z - mu)/s, df = a + 1)/(s*stats::pt(mu/s, df = a + 1))
  }
  jondrow <- function(z, e) {
    mu <- -e*su^2/om2; s <- sqrt(sv^2*su^2/om2)
    stats::dnorm((z - mu)/s)/(s*stats::pnorm(mu/s))
  }
  zz <- c(0.05, 0.5, 1.5, 3)
  for (e in c(-2, -0.5, 1))
    expect_equal(trt(zz, e, 1e7), jondrow(zz, e), tolerance = 1e-5)
  ## and it is a proper density for finite a
  for (e in c(-3, -1, 0.5))
    expect_equal(stats::integrate(function(z) trt(z, e, 5), 0, Inf)$value,
                 1, tolerance = 1e-6)
})

test_that("the THT efficiency quadrature reproduces integrate()", {
  ## sfm() computes E[exp(-u)|e] by a 256-node Gauss-Legendre rule on the
  ## probability scale, and E[u|e] by the closed-form truncated-t mean. Both
  ## must match direct numerical integration.
  su <- 1; sv <- 0.3; a <- 5; om2 <- sv^2 + su^2; nu <- a + 1
  e  <- c(-4, -2, -1, -0.3, 0, 0.5, 1.5, 3, 5)
  mu <- -e*su^2/om2
  s  <- sqrt((a + e^2/om2)*sv^2*su^2/(om2*(a + 1)))
  cc <- -mu/s

  gl <- sfa:::.gauss_legendre_01(256L)
  expect_equal(sum(gl$weights), 1, tolerance = 1e-12)
  expect_true(all(gl$nodes > 0 & gl$nodes < 1))
  ## exact on polynomials up to degree 2n-1
  for (k in c(1, 7, 31)) expect_equal(sum(gl$weights*gl$nodes^k), 1/(k + 1), tolerance = 1e-12)

  P  <- outer(stats::pt(cc, df = nu), gl$nodes, function(f, p) f + p*(1 - f))
  Zq <- pmax(mu + s*stats::qt(P, df = nu), 0)
  W  <- matrix(gl$weights, nrow = length(mu), ncol = 256L, byrow = TRUE)
  te_quad <- rowSums(W*exp(-Zq))
  u_closed <- mu + s*((nu + cc^2)/(nu - 1))*stats::dt(cc, df = nu)/
    stats::pt(cc, df = nu, lower.tail = FALSE)

  ref <- function(g) vapply(seq_along(mu), function(i)
    stats::integrate(function(z) g(z)*stats::dt((z - mu[i])/s[i], df = nu)/s[i],
                     0, Inf, rel.tol = 1e-12, subdivisions = 2000)$value /
      stats::pt(mu[i]/s[i], df = nu), numeric(1))

  expect_equal(te_quad,  ref(function(z) exp(-z)), tolerance = 1e-6)
  expect_equal(u_closed, ref(function(z) z),       tolerance = 1e-10)
})

test_that("THT reports efficiency, and outliers do NOT drive it to one", {
  skip_on_cran()
  d   <- cs_small(N = 400)
  fit <- sfm(y_pcs_st ~ x1 + x2, data = d, model_name = "THT")
  expect_s3_class(fit, "sfareg")
  expect_false(is.null(fit$exp_u_hat))
  expect_length(fit$exp_u_hat, nrow(d))
  expect_true(all(fit$exp_u_hat >= 0 & fit$exp_u_hat <= 1, na.rm = TRUE))
  expect_true(all(fit$u_hat >= 0, na.rm = TRUE))
  expect_false(is.null(fit$sd_exp_u_hat))

  ## Tancredi section 2.2: under the half-normal model a large POSITIVE
  ## residual sends E[exp(-u)|e] to 1; under the skew-t it does not, because
  ## such a point is read as an outlier rather than as evidence of full
  ## efficiency. Check the model's own predictor turns over.
  su <- 1; sv <- 0.3; a <- 5; om2 <- sv^2 + su^2
  te <- function(e) {
    mu <- -e*su^2/om2
    s  <- sqrt((a + e^2/om2)*sv^2*su^2/(om2*(a + 1)))
    stats::integrate(function(z) exp(-z)*stats::dt((z - mu)/s, df = a + 1)/s,
                     0, Inf)$value / stats::pt(mu/s, df = a + 1)
  }
  expect_lt(te(3), te(1))   # turns over instead of climbing to 1
})
