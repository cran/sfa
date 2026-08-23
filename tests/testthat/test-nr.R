## The normal-Rayleigh model, sfm(model_name = "NR").
##
## The first block is the important one. NR was long described in this package
## as "an alternative closed-form derivation of the same normal/half-normal
## composed error" as NHN, and was tested against NHN's DGP column on that
## basis. It is a different family, and these tests pin that down so the claim
## cannot come back.

## --- the density -------------------------------------------------------

## sfm.R's NR likelihood, transcribed. If the implementation changes, this
## copy must be updated to match, and the tests below then re-check it against
## an independent numerical convolution rather than against itself.
.nr_logden <- function(eps, sigv, sigu){
  sigma <- sqrt(2*sigv^2 + sigu^2)
  z     <- (eps*sigu/sigv)/sigma
  log(sigv) - 2*log(sigma) - 0.5*(eps/sigv)^2 + 0.5*z^2 +
    log(pmax(sqrt(2/pi)*exp(-0.5*z^2) - z*(1 - gsl::erf(z/sqrt(2))),
             .Machine$double.xmin))
}

test_that("the NR density is the normal-Rayleigh convolution, not the half-normal one", {
  sv <- 0.4; su <- 1.1
  ## u ~ Rayleigh(scale th), eps = v - u, integrated numerically.
  conv_ray <- function(e, th)
    stats::integrate(function(u) stats::dnorm((e+u)/sv)/sv * (u/th^2)*exp(-u^2/(2*th^2)),
                     0, Inf, rel.tol = 1e-12)$value
  conv_hn <- function(e)
    stats::integrate(function(u) stats::dnorm((e+u)/sv)/sv * sqrt(2/pi)/su*exp(-u^2/(2*su^2)),
                     0, Inf, rel.tol = 1e-12)$value

  es <- c(-3, -1.5, -0.4, 0, 0.3, 1, 2.5)
  coded <- exp(.nr_logden(es, sv, su))

  ## The scale convention: sigma_u is the SECOND RAW MOMENT of u, so the
  ## Rayleigh scale is sigma_u/sqrt(2). This is the assertion that fixes the
  ## meaning of sigma_u for the whole model.
  expect_equal(coded, vapply(es, conv_ray, 0, th = su/sqrt(2)), tolerance = 1e-6)
  ## and it is emphatically NOT the half-normal convolution
  expect_gt(max(abs(coded/vapply(es, conv_hn, 0) - 1)), 0.5)
})

test_that("the NR density integrates to one", {
  ## Finite range: the upper tail is floored at double.xmin by the pmax, which
  ## is deliberate -- the bracket is an O(exp(-z^2/4)) difference of two
  ## O(exp(-z^2/2)) terms and cancels completely past z ~ 35. Well beyond any
  ## region carrying mass.
  expect_equal(stats::integrate(function(e) exp(.nr_logden(e, 0.4, 1.1)),
                                -40, 15, rel.tol = 1e-10)$value,
               1, tolerance = 1e-7)
})

test_that("NR and NHN are different families, so neither can call the other", {
  ## Both are two-scale location families, so the shape of the composed error
  ## is pinned by ONE number: the standardized skewness the inefficiency
  ## contributes. Those constants differ, and no reparameterization moves a
  ## standardized moment -- which is why NR cannot be implemented as NHN under
  ## a parameter transformation.
  set.seed(11); n <- 2e6
  sk <- function(x) mean((x - mean(x))^3)/stats::sd(x)^3
  expect_equal(sk(abs(stats::rnorm(n))), 0.9953, tolerance = 1e-2)
  expect_equal(sk(sqrt(stats::rnorm(n)^2 + stats::rnorm(n)^2)),
               2*sqrt(pi)*(pi-3)/(4-pi)^1.5, tolerance = 1e-2)   ## 0.6311
})

## --- the DGP column ----------------------------------------------------

test_that("data_gen_cs() gains a Rayleigh column on the documented convention", {
  skip_on_cran()
  d <- data_gen_cs(N = 2e5, rand = 1, sig_u = 1, sig_v = 0.3, cons = 0.5,
                   beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.1)
  expect_true(all(c("u_r","y_pcs_r") %in% names(d)))
  expect_true(all(d$u_r > 0))
  expect_equal(mean(d$u_r^2), 1,            tolerance = 0.02)  ## E[u^2] = sig_u^2
  expect_equal(mean(d$u_r),   sqrt(pi)/2,   tolerance = 0.02)  ## E[u]
  expect_equal(stats::var(d$u_r), 1 - pi/4, tolerance = 0.02)  ## Var(u)
})

test_that("adding the Rayleigh draws did not disturb the existing RNG stream", {
  ## The columns are appended at the END of data_gen_cs() for exactly this
  ## reason: any draw inserted earlier shifts every column generated after it
  ## and silently changes every other model's test data. These are the values
  ## y_pcs had before u_r existed.
  d <- data_gen_cs(N = 5, rand = 42, sig_u = 1, sig_v = 0.3, cons = 0.5,
                   beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.1)
  expect_equal(d$y_pcs,
               c(3.5565605040, 2.6137019850, 1.8333047135,
                 1.2377404678, 3.1733961231),
               tolerance = 1e-8)
})

## --- the starting values -----------------------------------------------

test_that(".nr_start inverts the Rayleigh moment equations", {
  set.seed(4); n <- 2e5
  sv <- 0.3; su <- 1
  e  <- stats::rnorm(n, 0, sv) - (su/sqrt(2))*sqrt(stats::rnorm(n)^2 + stats::rnorm(n)^2)
  st <- sfa:::.nr_start(e - mean(e), 0, c(0.5, 0.5))
  expect_equal(st[1], sv, tolerance = 0.05)   ## sigma_v
  expect_equal(st[2], su, tolerance = 0.05)   ## sigma_u
  ## the intercept is lifted by E[u] = sigma_u*sqrt(pi)/2
  expect_equal(st[3], st[2]*sqrt(pi)/2, tolerance = 1e-10)
  expect_equal(st[4:5], c(0.5, 0.5))
})

test_that(".nr_start declines wrongly skewed residuals rather than guessing", {
  ## Positive skew: the moment equations have no admissible solution, and
  ## start_cs() must fall back to the old flat start instead of taking a
  ## complex root or a negative variance.
  set.seed(5)
  expect_null(sfa:::.nr_start(abs(stats::rnorm(500)), 0, c(0.5, 0.5)))
  expect_null(sfa:::.nr_start(numeric(0), 0, 0.5))
  expect_null(sfa:::.nr_start(c(1, NA, 2), 0, 0.5))
})

test_that("start_cs() hands NR the moment start, not the flat 0.1", {
  skip_on_cran()
  d <- data_gen_cs(N = 2000, rand = 3, sig_u = 1, sig_v = 0.3, cons = 0.5,
                   beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.1)
  s <- start_cs(y_pcs_r ~ x1 + x2, data_orig = d,
                x_vars_vec = c("(Intercept)","x1","x2"), intercept = 1,
                model_name = "NR", n_x_vars = 3, start_val = NA,
                n_z_vars = 0, z_vars = NA)$start_v
  expect_false(isTRUE(all.equal(s[1:2], c(0.1, 0.1))))   ## the old flat start
  expect_equal(s[1], 0.3, tolerance = 0.15)
  expect_equal(s[2], 1.0, tolerance = 0.25)
})

## --- end to end --------------------------------------------------------

test_that("NR recovers its own DGP and beats the truth's likelihood", {
  skip_on_cran()
  ll <- function(p, y, X) sum(.nr_logden(as.numeric(y - X %*% p[-(1:2)]), p[1], p[2]))
  tru <- c(0.3, 1, 0.5, 0.5, 0.5)
  beat <- 0L
  for(r in 1:5){
    d <- data_gen_cs(N = 4000, rand = 700 + r, sig_u = 1, sig_v = 0.3, cons = 0.5,
                     beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.1)
    f <- sfm(y_pcs_r ~ x1 + x2, model_name = "NR", data = d)
    p <- as.numeric(coef(f))
    expect_equal(p[1], 0.3, tolerance = 0.15)
    expect_equal(p[2], 1.0, tolerance = 0.15)
    expect_equal(p[3], 0.5, tolerance = 0.25)
    ## A fit BELOW the truth is an optimizer failure -- the symptom that made
    ## NR fail every convergence test from the flat start.
    X <- cbind(1, d$x1, d$x2)
    if(ll(p, d$y_pcs_r, X) >= ll(tru, d$y_pcs_r, X)) beat <- beat + 1L
  }
  expect_equal(beat, 5L)
})

test_that("NR still reports efficiency and standard errors", {
  skip_on_cran()
  d <- data_gen_cs(N = 1500, rand = 8, sig_u = 1, sig_v = 0.3, cons = 0.5,
                   beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.1)
  f <- sfm(y_pcs_r ~ x1 + x2, model_name = "NR", data = d)
  expect_true(all(is.finite(f$exp_u_hat)))
  expect_true(all(f$exp_u_hat >= 0))
  expect_true(all(is.finite(f$std.errors)))
})
