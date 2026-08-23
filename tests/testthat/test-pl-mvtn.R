## Pitt and Lee (1981) Model III -- multivariate truncated normal panel model.

mvtn_data <- function(N = 50, TT = 4, seed = 5) {
  as.data.frame(data_gen_p(t = TT, N = N, rand = seed, sig_u = 1, sig_v = 0.3,
                           sig_r = 0.2, sig_h = 0.4, cons = 0.5,
                           beta1 = 0.5, beta2 = 0.5))
}

test_that("the DGP column is a valid negative-orthant draw", {
  d <- mvtn_data()
  expect_true(all(c("u_mvtn", "y_pl_mvtn") %in% names(d)))
  expect_true(all(d$u_mvtn <= 1e-9))
  ## correlated within firm across periods (attenuated from rho by truncation)
  m <- matrix(d$u_mvtn, nrow = 4)
  expect_gt(mean(cor(t(m))[upper.tri(diag(4))]), 0.1)
})

test_that(".pl_mvtn_parts rejects inadmissible parameters", {
  expect_null(sfa:::.pl_mvtn_parts(-1, 0.5, 0.3, 4))     # sigma_u <= 0
  expect_null(sfa:::.pl_mvtn_parts(1, 0.5, -0.3, 4))     # sigma_v <= 0
  expect_null(sfa:::.pl_mvtn_parts(1, 1, 0.3, 4))        # rho = 1 is degenerate
  ## rho = 0.999 is near-degenerate but Sigma is still PD, so it is admissible
  expect_type(sfa:::.pl_mvtn_parts(1, 0.999, 0.3, 4), "list")
  expect_null(sfa:::.pl_mvtn_parts(1, -0.9, 0.3, 4))     # rho below -1/(T-1)
  p <- sfa:::.pl_mvtn_parts(1, 0.5, 0.3, 4)
  expect_type(p, "list")
  ## closed-form pieces must agree with brute force
  expect_equal(p$logdet_Sig, as.numeric(determinant(p$Sigma, logarithm = TRUE)$modulus),
               tolerance = 1e-10)
  Qinv <- solve(p$Sigma) + diag(1 / 0.3^2, 4)
  expect_equal(p$Q, solve(Qinv), tolerance = 1e-8)
  expect_equal(p$logdet_Q, as.numeric(determinant(p$Q, logarithm = TRUE)$modulus),
               tolerance = 1e-8)
})

test_that("the likelihood returns the penalty rather than erroring on bad input", {
  d <- mvtn_data(20, 4)
  Y <- split(d$y_pl_mvtn, d$name)
  X <- lapply(split(seq_len(nrow(d)), d$name), function(i) as.matrix(d[i, c("x1","x2")]))
  bad <- c(-1, 1, 0.5, 0.5, 0.5)          # sigma_v < 0
  expect_equal(sfa:::.pl_mvtn_nll(bad, Y, X, 4), sfa:::.SFA_CONSTANTS$MAX_VALUE)
  ok <- sfa:::.pl_mvtn_nll(c(0.3, 1, 0.5, 0.5, 0.5), Y, X, 4)
  expect_true(is.finite(ok))
})

test_that("the likelihood is maximized at the truth", {
  skip_on_cran()
  set.seed(20); N <- 100; TT <- 4
  SU <- 0.8; RHO <- 0.5; SV <- 0.3; B <- c(0.5, 0.5)
  Sig <- SU^2 * ((1 - RHO) * diag(TT) + RHO)
  U <- tmvtnorm::rtmvnorm(N, rep(0, TT), Sig, lower = rep(-Inf, TT),
                          upper = rep(0, TT), algorithm = "gibbs",
                          burn.in.samples = 100)
  Y <- X <- vector("list", N)
  for (i in seq_len(N)) {
    Xi <- cbind(runif(TT, 1, 3), runif(TT, 1, 3))
    Y[[i]] <- as.numeric(Xi %*% B) + U[i, ] + rnorm(TT, 0, SV)
    X[[i]] <- Xi
  }
  truth <- c(SV, SU, RHO, B)
  at_truth <- sfa:::.pl_mvtn_nll(truth, Y, X, TT)
  ## perturbing any single parameter must not improve the fit
  for (j in seq_along(truth)) {
    for (delta in c(-0.15, 0.15)) {
      p <- truth; p[j] <- p[j] + delta
      expect_gt(sfa:::.pl_mvtn_nll(p, Y, X, TT), at_truth)
    }
  }
})

test_that("psfm() fits it and reports the right parameters", {
  skip_on_cran()
  d <- mvtn_data(60, 4, seed = 5)
  f <- psfm(y_pl_mvtn ~ x1 + x2, model_name = "PL80_MVTN", data = d,
            individual = "name")
  expect_s3_class(f, "sfareg")
  expect_equal(f$model_name, "PL80_MVTN")
  expect_identical(names(f$coefficients),
                   c("sigv", "sigu", "rho", "(Intercept)", "x1", "x2"))
  expect_true(all(is.finite(f$coefficients)))
  expect_true(is.finite(as.numeric(logLik(f))))
  ## rho must stay inside the PD region
  expect_gt(f$coefficients[["rho"]], -1 / (4 - 1))
  expect_lt(f$coefficients[["rho"]], 1)
  ## efficiency predictions are on (0, 1]
  expect_true(all(f$exp_u_hat > 0 & f$exp_u_hat <= 1, na.rm = TRUE))
  expect_true(all(f$u_hat >= 0, na.rm = TRUE))
})

test_that("an unbalanced panel and T = 1 are refused with actionable messages", {
  skip_on_cran()
  d <- mvtn_data(40, 4)
  expect_error(
    psfm(y_pl_mvtn ~ x1 + x2, model_name = "PL80_MVTN", data = d[-(1:3), ],
         individual = "name"),
    "requires a BALANCED panel"
  )
  d1 <- mvtn_data(40, 4)
  d1 <- d1[d1$year == 1, ]
  expect_error(
    psfm(y_pl_mvtn ~ x1 + x2, model_name = "PL80_MVTN", data = d1,
         individual = "name"),
    "at least two periods"
  )
})
