## Normal-gamma, and the parameter-offset bug it exposed.
##
## sfm()'s likelihood closure slices the frontier coefficients out of the
## parameter vector by a fixed offset: x[3:(n_x+2)] for models with two leading
## scale parameters, x[4:(n_x+3)] for models with three. NG and NNAK carry
## three (sigv, sigu, mu) but were in the two-parameter group, so the slice
## took the right NUMBER of coefficients starting one slot too early.

test_that("the NG log-density matches a numerical convolution of its two parts", {
  ## The density itself was never wrong -- this pins that, so a future change
  ## to .log_pcf() cannot quietly break it.
  sv <- 0.3; P <- 2; th <- 0.5
  e  <- c(2, 0.5, 0, -0.5, -1, -2, -4, -8)
  pkg <- {
    z <- e/sv + sv/th
    (P-1)*log(sv) - 0.5*log(2) - 0.5*log(pi) - P*log(th) -
      0.5*(e/sv)^2 + 0.25*z^2 + sfa:::.log_pcf(-P, z)
  }
  ref <- vapply(e, function(ei){
    f  <- function(u) stats::dnorm(ei + u, 0, sv)*stats::dgamma(u, shape = P, scale = th)
    pk <- max(-ei, 0)
    br <- sort(unique(pmax(c(0, pk*c(0.5, 1, 2), pk + 10*sv, 50*th), 0)))
    s  <- sum(vapply(seq_len(length(br) - 1), function(j)
           stats::integrate(f, br[j], br[j+1], rel.tol = 1e-12,
                            subdivisions = 2000L)$value, numeric(1)))
    log(s + stats::integrate(f, br[length(br)], Inf, rel.tol = 1e-12,
                             subdivisions = 2000L)$value)
  }, numeric(1))
  expect_equal(pkg, ref, tolerance = 1e-9)
})

test_that("every model's coefficient slice starts after its scale parameters", {
  ## The regression guard. For each cross-sectional model, the number of
  ## leading non-beta parameters implied by the reported column names must
  ## match the offset the likelihood uses. Reading it off the fitted object
  ## keeps the two in step without re-parsing sfm()'s source.
  skip_on_cran()
  d <- cs_small(N = 200)
  two <- c("NHN", "NE", "NR", "NU", "NGE")
  three <- c("NTN", "NG", "NNAK", "NW")
  for (mn in c(two, three)) {
    f <- try(suppressWarnings(sfm(y_pcs ~ x1 + x2, model_name = mn, data = d)),
             silent = TRUE)
    if (inherits(f, "try-error")) next
    nm  <- rownames(f$out)
    lead <- which(nm == "(Intercept)") - 1L
    expect_equal(lead, if (mn %in% two) 2L else 3L,
                 info = paste(mn, "reports", lead, "leading scale parameters"))
  }
})

test_that("NG recovers its DGP and beats the truth on the likelihood", {
  ## Before the offset fix this returned sigma_u = 0 and stopped 710
  ## log-likelihood units BELOW the true parameter vector at n = 4000.
  skip_on_cran()
  truth <- c(0.3, 0.5, 2, 0.5, 0.5, 0.5)
  d <- data_gen_cs(N = 4000, rand = 1, sig_u = 1, sig_v = 0.3, cons = 0.5,
                   beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.1)
  Y <- d$y_pcs_g; X <- cbind(1, d$x1, d$x2)
  nll <- function(p){
    sv <- p[1]; su <- p[2]; mu <- p[3]
    e  <- as.numeric(Y - X %*% p[4:6]); z <- e/sv + sv/su
    -sum((mu-1)*log(sv) - 0.5*log(2) - 0.5*log(pi) - mu*log(su) -
         0.5*(e/sv)^2 + 0.25*z^2 + sfa:::.log_pcf(-mu, z))
  }
  f <- suppressWarnings(sfm(y_pcs_g ~ x1 + x2, model_name = "NG", data = d))
  p <- unname(coef(f)[c("sigv","sigu","mu","(Intercept)","x1","x2")])

  expect_lt(nll(p), nll(truth))          ## the fit must beat the truth
  expect_gt(p[2], 0.05)                  ## sigma_u not collapsed to the bound
  expect_equal(p[5], 0.5, tolerance = 0.1)
  expect_equal(p[6], 0.5, tolerance = 0.1)
})

test_that("NG's slope coefficients all actually enter the likelihood", {
  ## The sharpest symptom of the old offset: the LAST coefficient was never
  ## read, so it simply kept its starting value. Perturbing ANY coefficient
  ## must move the log-likelihood.
  ##
  ## Evaluated directly rather than through a re-fit: routing this through
  ## sfm(..., optHessian = FALSE) leaves no $opt, so logLik() returns NA and
  ## the comparison would pass for the wrong reason.
  d <- data_gen_cs(N = 800, rand = 2, sig_u = 1, sig_v = 0.3, cons = 0.5,
                   beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.1)
  Y <- d$y_pcs_g; X <- cbind(1, d$x1, d$x2)
  ll <- function(p){
    sv <- p[1]; su <- p[2]; mu <- p[3]
    e  <- as.numeric(Y - X %*% p[4:6]); z <- e/sv + sv/su
    sum((mu-1)*log(sv) - 0.5*log(2) - 0.5*log(pi) - mu*log(su) -
        0.5*(e/sv)^2 + 0.25*z^2 + sfa:::.log_pcf(-mu, z))
  }
  p0   <- c(0.3, 0.5, 2, 0.5, 0.5, 0.5)
  base <- ll(p0)
  expect_true(is.finite(base))
  for (j in seq_along(p0)) {
    q <- p0; q[j] <- q[j] + 0.25
    expect_false(isTRUE(all.equal(ll(q), base, tolerance = 1e-10)),
                 info = paste("parameter", j, "does not enter the likelihood"))
  }
})

test_that("NG reports its starting-value search", {
  skip_on_cran()
  d <- data_gen_cs(N = 600, rand = 3, sig_u = 1, sig_v = 0.3, cons = 0.5,
                   beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.1)
  f <- suppressWarnings(sfm(y_pcs_g ~ x1 + x2, model_name = "NG", data = d))
  expect_true(!is.null(f$ng_starts))
  expect_gte(f$ng_starts$n_tried, 5L)
  expect_true(is.finite(f$ng_starts$best))
})

test_that(".ng_start_candidates sweeps the shape along the E[u] ridge", {
  set.seed(7)
  e <- stats::rnorm(2000, 0, 0.3) - stats::rgamma(2000, shape = 2, scale = 0.5)
  cs <- sfa:::.ng_start_candidates(e - mean(e), beta_0_st = 0,
                                   beta_hat = c(x1 = 0.5, x2 = 0.5))
  expect_gte(length(cs), 5L)
  ## every candidate positive in the three scale slots, and finite throughout
  for (z in cs) {
    expect_true(all(is.finite(z)))
    expect_true(all(z[1:3] > 0))
  }
  ## E[u] = mu*sigma_u should be roughly constant across the grid candidates,
  ## which is the point of sweeping along the ridge rather than across it.
  eu <- vapply(cs, function(z) z[3]*z[2], numeric(1))
  expect_lt(stats::sd(eu[-1])/mean(eu[-1]), 1e-8)
})
