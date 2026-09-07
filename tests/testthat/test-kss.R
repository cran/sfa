## Kneip, Sickles and Song (2012).

.kss_dgp <- function(N = 110, Tt = 10, L = 2, seed = 21, sd_v = 0.2) {
  set.seed(seed)
  d <- expand.grid(t = seq_len(Tt), name = seq_len(N))
  d$year <- d$t
  n <- nrow(d)
  d$x1 <- runif(n, 1, 5)
  d$x2 <- runif(n, 1, 5)
  tt <- seq_len(Tt) / Tt
  G <- cbind(sin(pi * tt), cos(2 * pi * tt), tt^2 - mean(tt^2))[, seq_len(L), drop = FALSE]
  TH <- matrix(rnorm(N * L, 0, c(.8, .4, .25)[seq_len(L)]), N, L, byrow = TRUE)
  a <- rowSums(TH[d$name, , drop = FALSE] * G[d$t, , drop = FALSE])
  d$y <- 0.6 * d$x1 + 0.3 * d$x2 + a + rnorm(n, 0, sd_v)
  attr(d, "alpha") <- a
  d
}

test_that("KSS recovers the rank of the factor space", {
  skip_on_cran()
  for (L in 1:3) {
    d <- .kss_dgp(L = L, seed = 20 + L)
    f <- psfm(y ~ x1 + x2, "KSS", d, individual = "name", time = "year")
    expect_identical(f$kss$L, L)
    expect_true(all(abs(coef(f) - c(0.6, 0.3)) < 3 * f$std.errors))
    expect_equal(f$sigma_v, 0.2, tolerance = 0.05)
    ## Against the PERIOD-CENTRED truth: KSS's normalization is
    ## sum_i u_i(t) = 0, so its alpha_hat carries no period mean.
    ac <- attr(d, "alpha") - stats::ave(attr(d, "alpha"), d$year)
    expect_gt(cor(f$alpha_hat, ac), 0.95)
    expect_true(f$kss$converged)
  }
})

test_that("kss_L overrides the criterion and kss_smooth the GCV", {
  skip_on_cran()
  d <- .kss_dgp(L = 2)
  f4 <- psfm(y ~ x1 + x2, "KSS", d, individual = "name", time = "year", kss_L = 4)
  expect_identical(f4$kss$L, 4L)
  expect_equal(dim(f4$kss$basis), c(10L, 4L))
  f0 <- psfm(y ~ x1 + x2, "KSS", d,
    individual = "name", time = "year", kss_smooth = 0
  )
  expect_equal(f0$kss$kappa, 0)
  ## An unsmoothed fit is still a valid factor model, just a rougher one.
  expect_true(all(abs(coef(f0) - c(0.6, 0.3)) < 3 * f0$std.errors))
})

test_that("the estimated basis spans the true factor space", {
  skip_on_cran()
  d <- .kss_dgp(L = 2, seed = 33)
  f <- psfm(y ~ x1 + x2, "KSS", d, individual = "name", time = "year")
  tt <- seq_len(10) / 10
  Gt <- qr.Q(qr(cbind(sin(pi * tt), cos(2 * pi * tt))))
  G <- f$kss$basis
  ## Mean squared canonical correlation between the two 2-dim subspaces; 1
  ## means identical span. The BASIS is not identified, only the space it
  ## spans, so this is the right thing to check rather than g_1 against g_1.
  align <- sum(diag(t(G) %*% Gt %*% t(Gt) %*% G)) / ncol(G)
  expect_gt(align, 0.95)
})

test_that("KSS nests CSS's and LS's structures", {
  skip_on_cran()
  ## Lee-Schmidt is the one-factor case: KSS should find L = 1 and track it.
  set.seed(9)
  N <- 110; Tt <- 10
  d <- expand.grid(t = seq_len(Tt), name = seq_len(N))
  d$year <- d$t
  d$x1 <- runif(nrow(d), 1, 5)
  d$x2 <- runif(nrow(d), 1, 5)
  delta <- c(1, 1.4, 1.9, 2.3, 2.5, 2.4, 2.0, 1.5, 1.2, 1.0)
  ai <- -abs(rnorm(N, 0, 0.6))
  a <- delta[d$t] * ai[d$name]
  d$y <- 0.6 * d$x1 + 0.3 * d$x2 + a + rnorm(nrow(d), 0, 0.2)

  fk <- psfm(y ~ x1 + x2, "KSS", d, individual = "name", time = "year")
  fl <- psfm(y ~ x1 + x2, "LS", d, individual = "name", time = "year")
  expect_identical(fk$kss$L, 1L)

  ## Compare on u_hat, not alpha_hat. KSS centres by period (its equation 74
  ## imposes sum_i u_i(t) = 0) and LS does not, so their alpha_hat differ by a
  ## per-period constant BY CONSTRUCTION -- against the raw truth KSS scores
  ## 0.934 and LS 0.995, and against the period-centred truth the two swap
  ## places exactly. u_hat is a within-period contrast and so is invariant to
  ## that shift; it is also what the model is actually for.
  u_true <- stats::ave(a, d$year, FUN = max) - a
  expect_gt(cor(fk$u_hat, u_true), 0.97)
  expect_gt(cor(fk$u_hat, fl$u_hat), 0.98)
  ## LS is correctly specified for this DGP and KSS has to estimate the basis,
  ## so LS should be the more accurate of the two -- but not by much.
  expect_lt(
    sqrt(mean((fk$u_hat - u_true)^2)) - sqrt(mean((fl$u_hat - u_true)^2)),
    0.05
  )
})

test_that("KSS refuses what it is not defined on", {
  skip_on_cran()
  d <- .kss_dgp(N = 40, Tt = 6)
  expect_error(
    psfm(y ~ x1 + x2, "KSS", d[-3, ], individual = "name", time = "year"),
    "BALANCED panel"
  )
  d2 <- .kss_dgp(N = 40, Tt = 2)
  expect_error(
    psfm(y ~ x1 + x2, "KSS", d2, individual = "name", time = "year"),
    "at least 3 periods"
  )
})

test_that("KSS returns per-observation scores and no logLik", {
  skip_on_cran()
  d <- .kss_dgp(N = 60, Tt = 8)
  f <- psfm(y ~ x1 + x2, "KSS", d, individual = "name", time = "year")
  expect_equal(nobs(f), nrow(d))
  expect_length(f$u_hat, nrow(d))
  expect_true(all(f$u_hat >= 0))
  expect_true(all(tapply(f$u_hat, d$year, min) < 1e-8))
  expect_true(all(f$exp_u_hat > 0 & f$exp_u_hat <= 1))
  expect_true(is.na(suppressWarnings(as.numeric(logLik(f)))))
})

test_that("KSS finds CSS's polynomial space without being told the basis", {
  skip_on_cran()
  ## A CSS design: every firm effect lies in span{1, t, t^2}, so the true
  ## factor space is three-dimensional AND polynomial. KSS is given neither
  ## fact.
  set.seed(5)
  N <- 120
  Tt <- 10
  d <- expand.grid(t = seq_len(Tt), name = seq_len(N))
  d$year <- d$t
  d$x1 <- runif(nrow(d), 1, 5)
  d$x2 <- runif(nrow(d), 1, 5)
  th <- cbind(rnorm(N, 0, .5), rnorm(N, 0, .10), rnorm(N, 0, .03))
  a <- th[d$name, 1] + th[d$name, 2] * d$t + th[d$name, 3] * d$t^2
  d$y <- 0.6 * d$x1 + 0.3 * d$x2 + a + rnorm(nrow(d), 0, 0.2)

  fk <- psfm(y ~ x1 + x2, "KSS", d, individual = "name", time = "year")
  fc <- psfm(y ~ x1 + x2, "CSS", d, individual = "name", time = "year")

  expect_identical(fk$kss$L, 3L)
  ## The estimated basis should lie inside the polynomial span, even though
  ## the individual functions are not identified (only the space is).
  tt <- seq_len(Tt)
  P <- qr.Q(qr(cbind(1, tt, tt^2)))
  G <- fk$kss$basis
  expect_gt(sum(diag(t(G) %*% P %*% t(P) %*% G)) / ncol(G), 0.99)

  ## And the two should agree on the scores, CSS being correctly specified here.
  u_true <- stats::ave(a, d$year, FUN = max) - a
  expect_gt(cor(fk$u_hat, fc$u_hat), 0.98)
  expect_gt(cor(fk$u_hat, u_true), 0.98)
})

test_that("the rank criterion tracks the signal, not the arithmetic", {
  skip_on_cran()
  ## A third factor at the noise floor should NOT be selected, and one clearly
  ## above it should be. This is the property that broke when the criterion was
  ## scored on the smoothed residuals instead of the raw ones -- that version
  ## returned L_max regardless.
  mk <- function(sd3) {
    set.seed(5)
    N <- 120; Tt <- 10
    d <- expand.grid(t = seq_len(Tt), name = seq_len(N))
    d$year <- d$t
    d$x1 <- runif(nrow(d), 1, 5)
    d$x2 <- runif(nrow(d), 1, 5)
    th <- cbind(rnorm(N, 0, .5), rnorm(N, 0, .10), rnorm(N, 0, sd3))
    a <- th[d$name, 1] + th[d$name, 2] * d$t + th[d$name, 3] * d$t^2
    d$y <- 0.6 * d$x1 + 0.3 * d$x2 + a + rnorm(nrow(d), 0, 0.2)
    d
  }
  weak <- psfm(y ~ x1 + x2, "KSS", mk(0.004),
    individual = "name", time = "year"
  )
  strong <- psfm(y ~ x1 + x2, "KSS", mk(0.080),
    individual = "name", time = "year"
  )
  expect_identical(weak$kss$L, 2L)
  expect_identical(strong$kss$L, 3L)
  expect_lt(weak$kss$L, strong$kss$L)
})

test_that("the dimension criterion is not let loose on a short panel", {
  skip_on_cran()
  ## Bai-Ng needs T large relative to L. Allowing the criterion up to T - 1
  ## lets the factors span nearly the whole time dimension, V(L) collapses,
  ## and it returns the cap every time -- measured at 10 of 10 designs for
  ## both T = 6 and T = 8 before the floor(T/2) cap was imposed.
  for (Tt in c(8, 10, 15)) {
    for (L in 1:2) {
      d <- .kss_dgp(N = 100, Tt = Tt, L = L, seed = 40 + L)
      f <- suppressWarnings(
        psfm(y ~ x1 + x2, "KSS", d, individual = "name", time = "year")
      )
      expect_identical(f$kss$L, L,
        info = paste("T =", Tt, "true L =", L)
      )
    }
  }
})

test_that("KSS warns when the criterion lands on its own cap", {
  skip_on_cran()
  ## T = 6 caps the automatic search at 3. A design whose true rank is 3
  ## therefore selects the cap, and that is not a free choice.
  d <- .kss_dgp(N = 100, Tt = 6, L = 3, seed = 1)
  expect_warning(
    psfm(y ~ x1 + x2, "KSS", d, individual = "name", time = "year"),
    "largest it is allowed to consider"
  )
  ## An explicit kss_L is the user's own call and is honoured above the
  ## automatic cap, without the warning.
  f <- psfm(y ~ x1 + x2, "KSS", d,
    individual = "name", time = "year", kss_L = 5
  )
  expect_identical(f$kss$L, 5L)
})
