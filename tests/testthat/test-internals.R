## Unit tests for the numerical helpers in matrix_utils.R.
##
## These are the pieces most worth testing directly: they are shared by several
## model branches, they are where the package's real numerical content lives,
## and each one below has an independent reference to check against, so a wrong
## answer is detectable rather than merely plausible.

test_that(".gauss_hermite_nodes integrates known moments exactly", {
  gh <- .gauss_hermite_nodes(40L)
  expect_named(gh, c("nodes", "weights"))
  ## The rule integrates f(x)exp(-x^2); with x = z/sqrt(2) it gives E[f(Z)]
  ## for Z ~ N(0,1). Check the first four moments (1, 0, 1, 3).
  z <- gh$nodes * sqrt(2)
  w <- gh$weights / sqrt(pi)
  expect_equal(sum(w), 1, tolerance = 1e-12)
  expect_equal(sum(w * z), 0, tolerance = 1e-10)
  expect_equal(sum(w * z^2), 1, tolerance = 1e-10)
  expect_equal(sum(w * z^4), 3, tolerance = 1e-10)
})

test_that(".log_mvn_cdf_rank1 matches Monte Carlo for I + c*11'", {
  skip_on_cran()
  gh <- .gauss_hermite_nodes(64L)
  set.seed(11)
  for (cc in c(0.25, 1, 4)) {
    upper <- c(0.3, -0.5, 1.2, 0.0)
    Tn    <- length(upper)
    ## Z = xi + sqrt(c)*eta reproduces the covariance I + c*11'
    draws <- matrix(rnorm(2e5 * Tn), ncol = Tn) + sqrt(cc) * rnorm(2e5)
    mc    <- mean(apply(sweep(draws, 2, upper, "<="), 1, all))
    expect_equal(exp(.log_mvn_cdf_rank1(upper, cc, gh)), mc, tolerance = 0.01)
  }
})

test_that(".log_within_mvn_density matches the explicit multivariate normal", {
  ## Sigma = sigma^2 (I - 11'/T) is singular, so the density is the one for the
  ## first T-1 deviations. Compare against a direct evaluation on that block.
  e_full <- c(0.4, -0.2, 0.9, -1.1)
  Tn     <- length(e_full)
  e_dev  <- e_full - mean(e_full)
  sigma2 <- 1.7
  Sig    <- sigma2 * (diag(Tn) - matrix(1/Tn, Tn, Tn))
  S11    <- Sig[1:(Tn-1), 1:(Tn-1)]
  x1     <- e_dev[1:(Tn-1)]
  direct <- -0.5 * (Tn-1) * log(2*pi) - 0.5 * determinant(S11, logarithm = TRUE)$modulus -
              0.5 * as.numeric(t(x1) %*% solve(S11) %*% x1)
  expect_equal(.log_within_mvn_density(x1, sigma2), as.numeric(direct), tolerance = 1e-10)
})

test_that(".log_pnorm_diff equals log(Phi(b)-Phi(a)) where that is computable", {
  a <- c(-2, -1, 0, 0.5, -4)
  b <- c(-1,  0, 1, 2.0, -3)
  expect_equal(.log_pnorm_diff(a, b), log(pnorm(b) - pnorm(a)), tolerance = 1e-10)
})

test_that(".log_pnorm_diff stays finite deep in a tail where the naive form underflows", {
  ## Both limits far out in the same tail: pnorm(b)-pnorm(a) rounds to 0 and the
  ## naive log is -Inf. This is the case the helper exists for.
  expect_true(pnorm(-40) - pnorm(-41) == 0)
  expect_true(is.finite(.log_pnorm_diff(-41, -40)))
  expect_lt(.log_pnorm_diff(-41, -40), 0)
})

test_that(".te_battese_coelli and .jlms_u are in range and move the right way", {
  mu_star    <- c(-1, 0, 1, 3)
  sigma_star <- 0.5
  te <- .te_battese_coelli(mu_star, sigma_star)
  u  <- .jlms_u(mu_star, sigma_star)
  expect_true(all(te >= 0 & te <= 1))
  expect_true(all(u >= 0))
  ## More expected inefficiency => lower efficiency score, and E[u|e] rising.
  expect_true(all(diff(te) < 0))
  expect_true(all(diff(u)  > 0))
  ## Against a direct numerical integral of the truncated normal posterior.
  f  <- function(x) dnorm(x, 1, sigma_star) / pnorm(1/sigma_star)
  expect_equal(.te_battese_coelli(1, sigma_star),
               integrate(function(x) exp(-x)*f(x), 0, Inf)$value, tolerance = 1e-6)
  expect_equal(.jlms_u(1, sigma_star),
               integrate(function(x) x*f(x), 0, Inf)$value, tolerance = 1e-6)
})

test_that(".grad_nhn is the analytic gradient of the NHN negative log-likelihood", {
  skip_if_not_installed("numDeriv")
  set.seed(3)
  n    <- 120
  X    <- cbind(1, rnorm(n), rnorm(n))
  beta <- c(0.5, 0.5, 0.5)
  eps  <- rnorm(n, 0, 0.3) - abs(rnorm(n, 0, 1))
  Y    <- as.numeric(X %*% beta + eps)

  ## The objective every minimizer in this package is handed: the NEGATIVE
  ## summed log-likelihood.
  nll <- function(p) {
    lam <- p[1]; sig <- p[2]; b <- p[-c(1,2)]
    e   <- as.numeric(Y - X %*% b)
    -sum(log(2) - log(sig) + dnorm(e/sig, log = TRUE) +
           pnorm(-lam*e/sig, log.p = TRUE))
  }
  for (p in list(c(3.3, 1.05, 0.5, 0.5, 0.5), c(1.0, 0.8, 0.2, 0.7, 0.3))) {
    expect_equal(.grad_nhn(p, Y, X, inefdec_n = 1),
                 numDeriv::grad(nll, p), tolerance = 1e-6)
  }
})

test_that(".grad_nhn handles the cost frontier orientation too", {
  skip_if_not_installed("numDeriv")
  set.seed(4)
  n    <- 100
  X    <- cbind(1, rnorm(n))
  Y    <- as.numeric(X %*% c(0.5, 0.5) - rnorm(n, 0, 0.3) + abs(rnorm(n, 0, 1)))
  nll <- function(p) {
    lam <- p[1]; sig <- p[2]; b <- p[-c(1,2)]
    e   <- -as.numeric(Y - X %*% b)
    -sum(log(2) - log(sig) + dnorm(e/sig, log = TRUE) +
           pnorm(-lam*e/sig, log.p = TRUE))
  }
  p <- c(2.5, 1.1, 0.4, 0.6)
  expect_equal(.grad_nhn(p, Y, X, inefdec_n = -1),
               numDeriv::grad(nll, p), tolerance = 1e-6)
})

test_that(".grad_nhn degrades to zeros rather than NaN outside the parameter space", {
  ## optim()'s numerical Hessian steps outside the optimizer's own bounds, so
  ## the gradient must survive a non-positive lambda or sigma.
  X <- cbind(1, rnorm(20)); Y <- rnorm(20)
  expect_equal(.grad_nhn(c(-1, 1, 0, 0), Y, X), rep(0, 4))
  expect_equal(.grad_nhn(c(1, 0, 0, 0),  Y, X), rep(0, 4))
  expect_equal(.grad_nhn(c(NA, 1, 0, 0), Y, X), rep(0, 4))
})

test_that(".gsum sums by group independently of row order and label ordering", {
  ## rowsum(reorder = TRUE) sorts group labels as CHARACTER, so "10" lands
  ## before "2". This helper must not depend on that, nor on the data arriving
  ## grouped.
  set.seed(5)
  gid <- sample(rep(1:12, each = 4))
  x   <- rnorm(length(gid))
  expect_equal(.gsum(x, gid, 12), as.numeric(tapply(x, gid, sum)))
  ## Same data, permuted rows: identical answer.
  p <- sample(length(gid))
  expect_equal(.gsum(x[p], gid[p], 12), .gsum(x, gid, 12))
})

test_that(".build_Bit gives the documented time paths", {
  t_i <- 1:5
  expect_equal(.build_Bit("PL80", t_i), rep(1, 5))
  expect_equal(.build_Bit("BC92", t_i, Tref = 5, decay_par = 0.1),
               exp(-0.1 * (t_i - 5)))
  expect_equal(.build_Bit("K1990", t_i, decay_par = c(0.1, 0.02)),
               1/(1 + exp(0.1*t_i + 0.02*t_i^2)))
  expect_equal(.build_Bit("K1990modified", t_i, decay_par = c(0.1, 0.02)),
               1 + 0.1*(t_i - 5) + 0.02*(t_i - 5)^2)
  expect_error(.build_Bit("nope", t_i), "Unknown decay model_name")
})

test_that(".tfe_alpha_profile solves the same problem as a nested optimize()", {
  ## The concentrated firm effects underpinning psfm(model_name = "TFE"). A
  ## nested 1-D optimize() per firm is the obvious-but-slower way to get the
  ## same quantity, so it makes a good independent reference.
  set.seed(6)
  N <- 25; Tt <- 6
  gid <- rep(seq_len(N), each = Tt)
  r   <- rnorm(N * Tt) + rep(rnorm(N), each = Tt)

  ll_obs <- function(e, s, lam)
    log(2) - log(s) + dnorm(e/s, log = TRUE) + pnorm(-lam*e/s, log.p = TRUE)

  for (par in list(c(lam = 1, sig = 1), c(lam = 3.3, sig = 1.05), c(lam = 0.4, sig = 2))) {
    got <- .tfe_alpha_profile(r, gid, N, par[["lam"]], par[["sig"]])
    ref <- vapply(split(r, gid), function(rg)
      optimize(function(a) -sum(ll_obs(rg - a, par[["sig"]], par[["lam"]])),
               interval = c(min(rg) - 20*par[["sig"]], max(rg) + 20*par[["sig"]]),
               tol = 1e-12)$minimum, numeric(1))
    expect_equal(unname(got), unname(ref), tolerance = 1e-5)
  }
})

test_that(".tfe_alpha_profile returns a true stationary point", {
  ## Independent of any reference implementation: the score must vanish.
  set.seed(7)
  N <- 15; gid <- rep(seq_len(N), each = 5)
  r <- rnorm(N * 5); lam <- 2.2; sig <- 0.9
  alpha <- .tfe_alpha_profile(r, gid, N, lam, sig)
  e <- r - alpha[gid]
  a <- -(lam/sig) * e
  m <- exp(dnorm(a, log = TRUE) - pnorm(a, log.p = TRUE))
  score <- .gsum(e/sig^2 + (lam/sig)*m, gid, N)
  expect_lt(max(abs(score)), 1e-7)
})

test_that(".tfe_alpha_profile rejects a parameter outside the parameter space", {
  gid <- rep(1:3, each = 4); r <- rnorm(12)
  expect_null(.tfe_alpha_profile(r, gid, 3, lambda = -1, sigma = 1))
  expect_null(.tfe_alpha_profile(r, gid, 3, lambda = 1,  sigma = 0))
  expect_null(.tfe_alpha_profile(r, gid, 3, lambda = NA, sigma = 1))
})

test_that(".SFA_CONSTANTS keeps exp() and pnorm() arguments in a safe range", {
  expect_true(is.finite(exp(.SFA_CONSTANTS$EXP_CLIP_UPPER)))
  expect_gt(.SFA_CONSTANTS$MIN_POSITIVE, 0)
  expect_true(is.finite(pnorm(.SFA_CONSTANTS$CLIP_Z1_UPPER)))
  expect_equal(pnorm(.SFA_CONSTANTS$CLIP_Z1_UPPER), 1)
})
