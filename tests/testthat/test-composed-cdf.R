## The general composed-error CDF (prerequisite for goodness-of-fit testing).

test_that("the general path reproduces pcomposed() for the half-normal", {
  su <- 1.3; sv <- 0.7
  par <- c(lambda = su / sv, sigma = sqrt(su^2 + sv^2))
  qs <- c(-8, -4, -2, -1, 0, 1, 2, 4)
  for (dec in c(TRUE, FALSE)) {
    a <- pcomposed(qs, sigma_u = su, sigma_v = sv, inefdec = dec)
    b <- pcomposed_model(qs, "NHN", par, inefdec = dec)
    expect_equal(b, a, tolerance = 1e-10)
  }
  ## and in the far tail, where forming 1 - F would have saturated
  la <- pcomposed(c(-12, -20, -30), su, sv, TRUE, log.p = TRUE)
  lb <- pcomposed_model(c(-12, -20, -30), "NHN", par, TRUE, log.p = TRUE)
  expect_equal(lb, la, tolerance = 1e-8)
  expect_lt(lb[3], -200)
})

test_that("it is a proper CDF for every supported model", {
  cases <- list(
    NHN = c(lambda = 2, sigma = 1.2), NTN = c(lambda = 2, sigma = 1.5, mu = -0.4),
    NE = c(sigv = 0.6, sigu = 1.1), NR = c(sigv = 0.5, sigu = 1.2),
    NU = c(sigv = 0.4, theta = 1.8), NGE = c(sigv = 0.5, sigu = 0.9),
    NLN = c(sigv = 0.4, sigu = 0.8, mu = -0.2), NW = c(sigv = 0.4, sigu = 1.3, k = 1.8),
    NG = c(sigv = 0.5, sigu = 0.5, mu = 2), NNAK = c(sigv = 0.5, sigu = 1, mu = 1),
    TSL = c(sigv = 0.5, sigu = 0.9, lambda = 2),
    THT = c(sigu = 1.2, sigv = 0.5, a = 9), tHN = c(sigv = 0.35, sigu = 0.95, nu = 6)
  )
  qs <- seq(-12, 6, by = 0.5)
  for (mn in names(cases)) {
    p <- pcomposed_model(qs, mn, cases[[mn]])
    expect_true(all(is.finite(p)), info = mn)
    expect_true(all(p >= 0 & p <= 1), info = mn)
    expect_false(is.unsorted(p), info = mn) ## monotone
    expect_lt(p[1], 0.02, label = mn)
    expect_gt(p[length(p)], 0.98, label = mn)
    ## upper tail is computed directly, not as 1 - F
    up <- pcomposed_model(qs, mn, cases[[mn]], lower.tail = FALSE)
    expect_equal(up + p, rep(1, length(qs)), tolerance = 1e-6, info = mn)
  }
})

test_that("THT is a scale mixture, NOT an independent t-plus-half-normal", {
  ## Tancredi's composed error divides v and u by the SAME sqrt(V/a). The
  ## independent-convolution reading is a different distribution, and this
  ## test pins the difference so the two models cannot be merged by accident.
  su <- 1.2; sv <- 0.5; a <- 9
  set.seed(9); R <- 3e5
  W <- rgamma(R, shape = a / 2, rate = a / 2)
  mix <- (rnorm(R, 0, sv) - abs(rnorm(R, 0, su))) / sqrt(W)   ## skew-t: correct
  ind <- sv * rt(R, df = a) - abs(rnorm(R, 0, su))            ## tHN: different
  qs <- as.numeric(quantile(mix, c(.05, .25, .5, .75, .95)))
  thy <- pcomposed_model(qs, "THT", c(sigu = su, sigv = sv, a = a))
  emp_mix <- vapply(qs, function(q) mean(mix <= q), numeric(1))
  emp_ind <- vapply(qs, function(q) mean(ind <= q), numeric(1))
  expect_equal(thy, emp_mix, tolerance = 0.01)
  expect_gt(max(abs(thy - emp_ind)), 0.02) ## the two really do differ
})

test_that("bad input is refused clearly", {
  expect_error(pcomposed_model(0, "NHN", c(1, 2)), "NAMED")
  expect_error(pcomposed_model(0, "NOPE", c(sigv = 1, sigu = 1)), "no composed-error CDF")
  expect_error(pcomposed_model(0, "NHN", c(lambda = 2, sigma = 1), n_nodes = 4), "n_nodes")
  expect_error(pcomposed_model("a", "NHN", c(lambda = 2, sigma = 1)), "numeric")
  ## a parameter the model does not carry
  expect_error(pcomposed_model(0, "NE", c(sigv = 1)), "not found")
})

test_that("composed_cdf() maps a fit's residuals onto a uniform", {
  skip_on_cran()
  d <- as.data.frame(data_gen_cs(N = 400, rand = 5, sig_u = 1, sig_v = 0.5,
    cons = 0.5, beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.1
  ))
  f <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = d)
  p <- composed_cdf(f, data = d)
  expect_length(p, nrow(d))
  expect_true(all(p > 0 & p < 1))
  ## correctly specified, so the PIT should be roughly uniform
  expect_gt(stats::ks.test(p, "punif")$p.value, 0.01)
  expect_error(composed_cdf(f), NA) ## recoverable from the call
})
