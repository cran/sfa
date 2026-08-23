## The two data generators are the package's testing substrate -- every
## correctness claim about an estimator is a claim about recovering the
## parameters these put in. Their structure is therefore worth pinning.

test_that("data_gen_cs() is reproducible and returns the documented columns", {
  a <- data_gen_cs(N = 50, rand = 1, sig_u = 1, sig_v = 0.3, cons = 0.5,
                   beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.5)
  b <- data_gen_cs(N = 50, rand = 1, sig_u = 1, sig_v = 0.3, cons = 0.5,
                   beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.5)
  expect_identical(a, b)
  expect_equal(nrow(a), 50)
  ## Columns the model tests in this suite depend on. u_tn and the two-tier
  ## columns were added after each was found missing, so they are named here
  ## explicitly to stop them being dropped again.
  expect_true(all(c("y_pcs", "y_pcs_e", "y_pcs_tn", "y_pcs_t", "y_pcs_z",
                    "y_pcs_ez", "y_ttne", "y_tthn", "y_tthn_z", "y_zisf",
                    "y_zisf_z", "u", "u_tn", "z", "zp") %in% names(a)))
})

test_that("data_gen_p() is reproducible and returns the documented columns", {
  a <- data_gen_p(t = 4, N = 20, rand = 1, sig_u = 1, sig_v = 0.3, sig_r = 0.2,
                  sig_h = 0.4, cons = 0.5, beta1 = 0.5, beta2 = 0.5)
  b <- data_gen_p(t = 4, N = 20, rand = 1, sig_u = 1, sig_v = 0.3, sig_r = 0.2,
                  sig_h = 0.4, cons = 0.5, beta1 = 0.5, beta2 = 0.5)
  expect_identical(a, b)
  expect_equal(nrow(a), 80)
  expect_equal(length(unique(a$name)), 20)
  expect_true(all(c("y_gtre", "y_tre", "y_tfe", "y_fd", "y_ssfe", "y_bc92",
                    "y_gtre_zz", "y_tre_z", "h_z", "u_inv") %in% names(a)))
})

test_that("data_gen_p()'s time-invariant draws really are constant within a firm", {
  ## u_inv (and h) are drawn once per firm; SSFE and PL80 both assume it. A
  ## per-observation draw here would silently test the wrong model.
  d <- as.data.frame(data_gen_p(t = 5, N = 30, rand = 2, sig_u = 1, sig_v = 0.3,
                                sig_r = 0.2, sig_h = 0.4, cons = 0.5,
                                beta1 = 0.5, beta2 = 0.5))
  within_sd <- function(x) tapply(x, d$name, function(v) diff(range(v)))
  expect_true(all(within_sd(d$u_inv) == 0))
  expect_true(all(within_sd(d$r) == 0))
  expect_true(all(within_sd(d$h) == 0))
  ## u, by contrast, is meant to vary within a firm (GTRE/TRE/TFE/FD use it).
  expect_true(any(within_sd(d$u) > 0))
})

test_that("data_gen_p()'s BC92 column decays at the requested rate", {
  ## y_bc92 uses B_it = exp(-eta*(t - Tref)); with eta = 0 it must reduce to
  ## the time-invariant case.
  args <- list(t = 6, N = 40, rand = 4, sig_u = 1, sig_v = 0.3, sig_r = 0.2,
               sig_h = 0.4, cons = 0.5, beta1 = 0.5, beta2 = 0.5)
  d0 <- as.data.frame(do.call(data_gen_p, c(args, list(eta = 0))))
  d1 <- as.data.frame(do.call(data_gen_p, c(args, list(eta = 0.5))))
  expect_false(isTRUE(all.equal(d0$y_bc92, d1$y_bc92)))
  ## At eta = 0 the BC92 outcome coincides with the time-invariant one.
  expect_equal(d0$y_bc92, d0$y_ssfe, tolerance = 1e-10)
})

test_that("the generators respond to their scale parameters", {
  big <- data_gen_cs(N = 4000, rand = 7, sig_u = 3, sig_v = 0.3, cons = 0.5,
                     beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.5)
  small <- data_gen_cs(N = 4000, rand = 7, sig_u = 0.5, sig_v = 0.3, cons = 0.5,
                       beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.5)
  expect_gt(sd(big$u), sd(small$u))
  ## u is half-normal, so its mean is sigma_u*sqrt(2/pi).
  expect_equal(mean(big$u), 3*sqrt(2/pi), tolerance = 0.1)
})
