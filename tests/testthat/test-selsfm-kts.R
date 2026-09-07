## Kumbhakar, Tsionas and Sipilainen (2009): joint technology choice and
## technical efficiency, selsfm(model_name = "kts").
##
## The load-bearing checks here do not fit anything. A fit at this model's
## sample sizes costs minutes, and two cheaper properties pin the likelihood
## far more tightly than one fit would: that it is a proper density, and that
## the true parameters maximise it.

kts_par <- list(b0 = c(1.0, 0.5), b1 = c(1.4, 0.4), sv0 = 0.3, sv1 = 0.35,
  su0 = 0.8, su1 = 0.6, gam = c(-0.2, 0.6), del = 0.9)

kts_theta <- function(p = kts_par) {
  c(p$b0, p$b1, log(p$sv0), log(p$sv1), log(p$su0), log(p$su1), p$gam, p$del)
}

## The generative process the paper's equations (8) and (9) imply. The SAME u
## drives the technology choice and enters output -- that is the "dual role",
## and drawing two independent u's instead silently fits a different model.
kts_sim <- function(n, seed, p = kts_par, nq = 200L) {
  set.seed(seed)
  nd <- sfa:::.kts_nodes(nq)
  x <- stats::rnorm(n); z <- stats::rnorm(n)
  X <- cbind(1, x); Z <- cbind(1, z)
  a <- as.numeric(Z %*% p$gam)
  A <- function(su) {
    as.numeric(exp(stats::pnorm(outer(a, rep(1, length(nd$z))) +
      p$del * outer(rep(1, n), su * nd$z), log.p = TRUE)) %*% nd$w)
  }
  A0 <- A(p$su0); A1 <- A(p$su1)
  phi <- pmin(pmax(A0 / (1 + A0 - A1), 0), 1)
  u <- ifelse(stats::runif(n) < phi, abs(stats::rnorm(n, 0, p$su1)),
    abs(stats::rnorm(n, 0, p$su0)))
  I <- stats::rbinom(n, 1, stats::pnorm(a + p$del * u))
  mu <- ifelse(I == 1, as.numeric(X %*% p$b1), as.numeric(X %*% p$b0))
  data.frame(y = mu - u + stats::rnorm(n, 0, ifelse(I == 1, p$sv1, p$sv0)),
    x = x, z = z, I = I)
}

test_that("the joint density of (y, I) integrates to 1", {
  ## The check that a likelihood is a likelihood. Summed over both regimes,
  ## because p(y, I) is a joint density of a continuous and a discrete part.
  nd <- sfa:::.kts_nodes(96L)
  th <- kts_theta()
  X <- matrix(c(1, 1.2), 1, 2); Z <- matrix(c(1, 0.7), 1, 2)
  yg <- seq(-14, 10, length.out = 4001); dy <- diff(yg)[1]
  n <- length(yg)
  tot <- 0
  for (I in c(0, 1)) {
    ll <- sfa:::.kts_loglik(th, y = yg,
      X0 = X[rep(1, n), , drop = FALSE], X1 = X[rep(1, n), , drop = FALSE],
      Z = Z[rep(1, n), , drop = FALSE], I = rep(I, n), nd = nd)
    tot <- tot + sum(exp(ll)) * dy
  }
  expect_equal(tot, 1, tolerance = 1e-3)
})

test_that("the simulator and the likelihood describe the same model", {
  skip_on_cran()
  ## If these disagree, a recovery failure would be blamed on the estimator.
  ## Model-implied P(I = 1) is phi from equation (9); the simulator's is the
  ## realised share.
  p <- kts_par
  set.seed(4); n <- 20000
  nd <- sfa:::.kts_nodes(200L)
  z <- stats::rnorm(n); a <- as.numeric(cbind(1, z) %*% p$gam)
  A <- function(su) {
    as.numeric(exp(stats::pnorm(outer(a, rep(1, length(nd$z))) +
      p$del * outer(rep(1, n), su * nd$z), log.p = TRUE)) %*% nd$w)
  }
  A0 <- A(p$su0); A1 <- A(p$su1)
  phi <- pmin(pmax(A0 / (1 + A0 - A1), 0), 1)
  d <- kts_sim(n, 4, p)
  expect_equal(mean(d$I), mean(phi), tolerance = 0.02)
})

test_that("the true delta maximises the likelihood, and sharply", {
  skip_on_cran()
  ## The check that would have caught the shape bug this file was first written
  ## with. Under that bug the choice probability was evaluated at ONE
  ## quadrature node, so the integral had no delta-dependence at all and the
  ## surface was flat to 0.0007 per observation; the argmax then wandered
  ## between grid points and looked like weak identification. It is not weak:
  ## corrected, the spread is ~0.069 per observation, a hundredfold larger.
  d <- kts_sim(30000, 7)
  X <- cbind(1, d$x); Z <- cbind(1, d$z)
  nd <- sfa:::.kts_nodes(150L)
  mean_ll <- function(del) {
    th <- kts_theta(); th[length(th)] <- del
    mean(sfa:::.kts_loglik(th, d$y, X, X, Z, d$I, nd))
  }
  grid <- c(0.1, 0.5, 0.9, 1.3, 1.8)          # truth is the third
  v <- vapply(grid, mean_ll, numeric(1))
  expect_identical(which.max(v), 3L)
  ## Strictly decreasing away from the truth in both directions.
  expect_true(all(diff(v[1:3]) > 0))
  expect_true(all(diff(v[3:5]) < 0))
  ## And the surface has real curvature, which is what makes delta estimable.
  expect_gt(max(v) - min(v), 0.03)
})

test_that("the selection factor varies across the quadrature grid", {
  ## A direct guard on the bug itself, independent of any fit. If the choice
  ## probability is constant in u then delta cannot be identified through the
  ## integral, and the model reduces to one without a dual role.
  d <- kts_sim(200, 2)
  X <- cbind(1, d$x); Z <- cbind(1, d$z)
  nd <- sfa:::.kts_nodes(24L)
  th <- kts_theta()
  ## Two deltas that differ only in the coupling must give different densities.
  th_a <- th; th_a[length(th_a)] <- 0
  th_b <- th; th_b[length(th_b)] <- 2
  la <- sfa:::.kts_loglik(th_a, d$y, X, X, Z, d$I, nd)
  lb <- sfa:::.kts_loglik(th_b, d$y, X, X, Z, d$I, nd)
  expect_gt(mean(abs(la - lb)), 1e-3)
})

test_that("selsfm(model_name = 'kts') returns a coherent fit", {
  skip_on_cran()
  skip_if(!nzchar(Sys.getenv("SFA_TEST_SLOW")),
    "the KTS likelihood is a quadrature per observation; set SFA_TEST_SLOW")
  d <- kts_sim(400, 3)
  f <- suppressWarnings(selsfm(selection = I ~ z, frontier = y ~ x, data = d,
    model_name = "kts", n_nodes = 32))
  expect_s3_class(f, "sfareg")
  expect_identical(f$model_name, "KTS")
  expect_true(all(is.finite(f$coefficients)))
  expect_true(all(c("t0.(Intercept)", "t1.(Intercept)", "sigma_u0", "sigma_u1",
    "delta") %in% rownames(f$out)))
  ## Scales are reported on the natural scale, so they must be positive.
  expect_true(all(f$out[c("sigma_v0", "sigma_v1", "sigma_u0", "sigma_u1"),
    "par"] > 0))
  expect_length(f$efficiency, nrow(d))
  expect_true(all(f$efficiency > 0 & f$efficiency <= 1))
})

test_that("it refuses a sample that cannot support two frontiers", {
  skip_on_cran()
  d <- kts_sim(60, 5)
  d$I[d$I == 1][-(1:2)] <- 0            # leave 2 in one regime
  expect_error(
    selsfm(selection = I ~ z, frontier = y ~ x, data = d, model_name = "kts"),
    "SEPARATE frontier to each regime|more observations than it has"
  )
})
