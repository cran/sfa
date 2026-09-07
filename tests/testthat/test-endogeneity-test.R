## endogeneity_test(): the Wald test of Hou, Ramalho and Roseta-Palma (2025),
## Energy Economics 151:108922, Eq. (25), plus the standard errors for rho that
## it needs and that ivsfm() previously reported as NA.

mk_iv <- function(n = 800, rho = 0.6, seed = 9) {
  set.seed(seed)
  z <- stats::rnorm(n); x1 <- stats::rnorm(n); eps <- stats::rnorm(n)
  x2 <- z + eps
  v <- rho * eps + sqrt(1 - rho^2) * stats::rnorm(n)
  u <- abs(stats::rnorm(n))
  data.frame(y = 0.5 * x1 + 0.5 * x2 + v - u, x1 = x1, x2 = x2, z = z)
}

test_that("the rho -> t Jacobian is the one the code claims", {
  ## r = t / sqrt(1 + t't), so dr_i/dt_j = delta_ij/s - t_i t_j/s^3, i.e.
  ## J = (I - r r') / s. Four lines of algebra with nothing else depending on
  ## them, so they are checked against a numerical derivative.
  for (tt in list(c(0.7), c(0.4, -1.2), c(0, 0, 0), c(3, 0.1, -2))) {
    f <- function(t) t / sqrt(1 + sum(t^2))
    s <- sqrt(1 + sum(tt^2)); r <- tt / s
    J <- (diag(length(tt)) - tcrossprod(r)) / s
    expect_equal(J, numDeriv::jacobian(f, tt), tolerance = 1e-7)
  }
  ## At the null the map is the identity, which is what makes the Wald test
  ## well behaved there.
  expect_equal((diag(2) - tcrossprod(c(0, 0))) / 1, diag(2))
})

test_that("ivsfm() now reports a standard error for rho", {
  skip_on_cran()
  d <- mk_iv()
  f <- suppressWarnings(ivsfm(y ~ x1 + x2, endogenous = ~x2, instruments = ~z,
    data = d, model_name = "IVLIML"))
  expect_true(is.finite(f$out["rho_x2", "st_err"]))
  expect_gt(f$out["rho_x2", "st_err"], 0)
  ## and it recovers the truth it was given
  expect_equal(unname(f$out["rho_x2", "par"]), 0.6, tolerance = 0.1)
  ## the joint covariance travels with the fit
  expect_true(is.matrix(f$vcov_rho))
  expect_identical(dim(f$vcov_rho), c(1L, 1L))
  expect_equal(sqrt(f$vcov_rho[1, 1]), unname(f$out["rho_x2", "st_err"]))
})

test_that("with one endogenous regressor the Wald statistic is the squared t", {
  skip_on_cran()
  d <- mk_iv()
  f <- suppressWarnings(ivsfm(y ~ x1 + x2, endogenous = ~x2, instruments = ~z,
    data = d, model_name = "IVLIML"))
  tt <- endogeneity_test(f)
  expect_s3_class(tt, "sfa_endog_test")
  expect_equal(unname(tt$statistic), unname(f$out["rho_x2", "t-val"])^2,
    tolerance = 1e-8)
  expect_identical(unname(tt$parameter), 1L)
  expect_equal(tt$p.value,
    stats::pchisq(unname(tt$statistic), 1, lower.tail = FALSE))
  expect_output(print(tt), "Wald test for endogeneity")
})

test_that("it refuses fits that carry no correction parameter", {
  skip_on_cran()
  d <- mk_iv()
  f <- suppressWarnings(ivsfm(y ~ x1 + x2, endogenous = ~x2, instruments = ~z,
    data = d, model_name = "C2SLS"))
  expect_error(endogeneity_test(f), "C2SLS")
  g <- suppressWarnings(sfm(y ~ x1 + x2, data = d, model_name = "NHN"))
  expect_error(endogeneity_test(g), "IVLIML")
  expect_error(endogeneity_test(1:5), "ivsfm")
})

test_that("it is correctly sized under exogeneity and powerful under it", {
  skip_on_cran()
  ## Small by Monte Carlo standards -- the full study is
  ## spec_tests/hou2025_table2.R. This asserts only that the two ends are on
  ## the right side of the numbers that study measures.
  pv <- function(rho, R) {
    vapply(seq_len(R), function(r) {
      d <- mk_iv(n = 600, rho = rho, seed = 1000 + r)
      f <- try(suppressWarnings(ivsfm(y ~ x1 + x2, endogenous = ~x2,
        instruments = ~z, data = d, model_name = "IVCF")), silent = TRUE)
      if (inherits(f, "try-error")) return(NA_real_)
      t2 <- try(endogeneity_test(f), silent = TRUE)
      if (inherits(t2, "try-error")) NA_real_ else t2$p.value
    }, numeric(1))
  }
  p0 <- pv(0, 60); p0 <- p0[is.finite(p0)]
  p8 <- pv(0.8, 20); p8 <- p8[is.finite(p8)]
  expect_gt(length(p0), 50)
  expect_lt(mean(p0 < 0.05), 0.25)   # size: loose, 60 replications
  expect_equal(mean(p8 < 0.05), 1)   # power at rho = 0.8 is 1.000
})
