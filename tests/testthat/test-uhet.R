## Kim and Schmidt (2008) two-step test.

.ks_data <- function(n, seed, rho = 0.9, delta = 0, dgp = c("NHN", "NE")) {
  dgp <- match.arg(dgp)
  set.seed(seed)
  x1 <- rnorm(n); x2 <- rnorm(n)
  z <- rho * x1 + sqrt(1 - rho^2) * rnorm(n)
  us <- if (dgp == "NHN") abs(rnorm(n, 0, 1)) else rexp(n, rate = 1)
  y <- 0.5 + 0.5 * x1 + 0.5 * x2 + rnorm(n, 0, 0.5) - exp(z * delta) * us
  data.frame(y, x1, x2, z)
}

test_that("the two tests agree when z is independent of x, and differ when it is not", {
  skip_on_cran()
  ## G = E[z grad_psi b] is zero when z is independent of x, and that is
  ## exactly when Kim and Schmidt say the naive test is valid.
  di <- .ks_data(1500, 4, rho = 0)
  fi <- sfm(y ~ x1 + x2, model_name = "NHN", data = di, keep_objective = TRUE)
  ti <- uhet_test(fi, ~z, data = di)
  expect_equal(ti$statistic[1], ti$statistic[2], tolerance = 0.35)

  dc <- .ks_data(1500, 4, rho = 0.9)
  fc <- sfm(y ~ x1 + x2, model_name = "NHN", data = dc, keep_objective = TRUE)
  tc <- uhet_test(fc, ~z, data = dc)
  ## with z correlated with x the naive standard error is far too large
  expect_gt(attr(tc, "se_naive")[["z"]], 1.2 * attr(tc, "se_corrected")[["z"]])
})

test_that("the output is well formed", {
  skip_on_cran()
  d <- .ks_data(600, 7)
  f <- sfm(y ~ x1 + x2, model_name = "NHN", data = d, keep_objective = TRUE)
  tt <- uhet_test(f, ~z, data = d)
  expect_equal(nrow(tt), 2L)
  expect_named(tt, c("test", "statistic", "df", "p.value"))
  expect_equal(tt$df, c(1L, 1L))
  expect_true(all(tt$statistic >= 0))
  expect_equal(tt$p.value, pchisq(tt$statistic, 1, lower.tail = FALSE))
  expect_named(attr(tt, "gamma"), "z")
  ## a matrix z works as well as a formula, and two columns give 2 df
  tt2 <- uhet_test(f, cbind(z = d$z, w = d$x2), data = d)
  expect_equal(tt2$df, c(2L, 2L))
})

test_that("it has power against a real dependence", {
  skip_on_cran()
  d <- .ks_data(1500, 9, rho = 0.5, delta = 0.4)
  f <- sfm(y ~ x1 + x2, model_name = "NHN", data = d, keep_objective = TRUE)
  tt <- uhet_test(f, ~z, data = d)
  expect_lt(tt$p.value[tt$test == "two-step (corrected)"], 0.05)
})

test_that("it works for the exponential model too", {
  skip_on_cran()
  d <- .ks_data(1200, 11, rho = 0.5, dgp = "NE")
  f <- sfm(y ~ x1 + x2, model_name = "NE", data = d, keep_objective = TRUE)
  expect_s3_class(uhet_test(f, ~z, data = d), "data.frame")
})

test_that("unsupported and degenerate cases are refused, not answered", {
  skip_on_cran()
  d <- .ks_data(500, 3)
  ## the correction is not distribution-free
  fr <- sfm(y ~ x1 + x2, model_name = "NR", data = d, keep_objective = TRUE)
  expect_error(uhet_test(fr, ~z, data = d), "not distribution-free")
  ## the score is needed, so keep_objective is not optional
  fn <- sfm(y ~ x1 + x2, model_name = "NHN", data = d)
  expect_error(uhet_test(fn, ~z, data = d), "keep_objective")
  ## z must contain something other than a constant
  f <- sfm(y ~ x1 + x2, model_name = "NHN", data = d, keep_objective = TRUE)
  expect_error(uhet_test(f, ~1, data = d), "at least one non-constant")
  expect_error(uhet_test(f, ~z, data = d[1:100, ]), "rows")
  expect_error(uhet_test(structure(list(), class = "lm"), ~z, data = d), "sfareg")
})

test_that("a wrong-skew fit is refused: there is no inefficiency to explain", {
  skip_on_cran()
  set.seed(15)
  n <- 400
  x1 <- rnorm(n); x2 <- rnorm(n); z <- rnorm(n)
  ## cost data read as production -> sigma_u collapses to the boundary
  y <- 0.5 + 0.5 * x1 + 0.5 * x2 + rnorm(n, 0, 0.5) + abs(rnorm(n, 0, 1))
  d <- data.frame(y, x1, x2, z)
  f <- suppressWarnings(sfm(y ~ x1 + x2, model_name = "NHN", data = d,
    keep_objective = TRUE
  ))
  skip_if_not(isTRUE(f$wrong_skew) || isTRUE(f$sigma_u_at_bound))
  expect_error(uhet_test(f, ~z, data = d), "wrong-skew|singular")
})

test_that("the first-step score matches an analytic half-normal score", {
  ## The correction leans entirely on the score, so this pins it. estfun()
  ## differences the likelihood numerically; the closed form is written out
  ## here independently.
  skip_on_cran()
  d <- .ks_data(800, 21, rho = 0.5)
  f <- sfm(y ~ x1 + x2, model_name = "NHN", data = d, keep_objective = TRUE)
  psi <- f$out[, "par"]
  lam <- psi[["lambda"]]; sg <- psi[["sigma"]]
  X <- model.matrix(y ~ x1 + x2, d)
  eps <- as.numeric(d$y - X %*% psi[colnames(X)])
  zz <- -eps * lam / sg
  m <- exp(dnorm(zz, log = TRUE) - pnorm(zz, log.p = TRUE))
  Sa <- cbind(
    lambda = m * (-eps / sg),
    sigma = -1 / sg + eps^2 / sg^3 + m * (eps * lam / sg^2),
    X * (eps / sg^2 + (lam / sg) * m)
  )
  colnames(Sa)[3:5] <- colnames(X)
  Sn <- estfun.sfareg(f)[, colnames(Sa)]
  expect_equal(unname(Sa), unname(Sn), tolerance = 1e-6)
  ## and the scores sum to zero at the MLE
  expect_lt(max(abs(colMeans(Sn))), 1e-5)
})
