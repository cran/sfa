## bread() and estfun(), which are all `sandwich` needs to work on these fits.

.sw_fit <- function(N = 400) {
  d <- data_gen_cs(N = N, rand = 1, sig_u = 0.8, sig_v = 0.3, cons = 0.5,
                   beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.1)
  list(fit = sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = d,
                 keep_objective = TRUE), data = d)
}

test_that("estfun agrees with the ANALYTIC gradient", {
  skip_on_cran()
  ## The scores are taken by central differences, so they need checking against
  ## something independent. NHN is the one model here with a hand-derived
  ## gradient, which makes it the right place to do it.
  o <- .sw_fit()
  X <- stats::model.matrix(y_pcs ~ x1 + x2, o$data)
  Y <- as.matrix(o$data$y_pcs)
  ## .grad_nhn is the gradient of the NEGATIVE summed log-likelihood; estfun
  ## holds contributions to the POSITIVE one.
  g_analytic <- sfa:::.grad_nhn(o$fit$opt$par, Y, X, inefdec_n = 1)
  g_numeric <- -colSums(sfa:::estfun.sfareg(o$fit))
  expect_equal(unname(g_analytic), unname(g_numeric), tolerance = 1e-5)
})

test_that("the scores vanish at the optimum, which is what makes them scores", {
  skip_on_cran()
  o <- .sw_fit()
  S <- sfa:::estfun.sfareg(o$fit)
  expect_equal(dim(S), c(nobs(o$fit), length(coef(o$fit))))
  expect_identical(colnames(S), names(coef(o$fit)))
  expect_true(all(is.finite(S)))
  ## Judged on the scale of the scores themselves, not absolutely: a column
  ## whose entries are of order 1 and whose sum is of order 1e-4 has vanished.
  rel <- abs(colSums(S)) / (apply(S, 2, stats::sd) * sqrt(nrow(S)))
  expect_true(all(rel < 1e-3), info = paste(signif(rel, 3), collapse = ", "))
})

test_that("bread is n times the covariance", {
  skip_on_cran()
  o <- .sw_fit(N = 200)
  expect_equal(sfa:::bread.sfareg(o$fit), vcov(o$fit) * nobs(o$fit))
})

test_that("sandwich and vcovCL produce usable covariances", {
  skip_on_cran()
  skip_if_not_installed("sandwich")
  o <- .sw_fit()
  V_rob <- sandwich::sandwich(o$fit)
  V_cl <- sandwich::vcovCL(o$fit, cluster = rep(seq_len(40), each = 10))
  for (V in list(V_rob, V_cl)) {
    expect_equal(dim(V), rep(length(coef(o$fit)), 2))
    expect_true(all(is.finite(V)))
    expect_true(all(diag(V) > 0))
    ## symmetric, as any covariance must be
    expect_equal(V, t(V), tolerance = 1e-10)
  }
  ## Robust and classical should be in the same ballpark on well-specified
  ## data -- an order-of-magnitude gap would mean something is wrong.
  r <- sqrt(diag(V_rob)) / sqrt(diag(vcov(o$fit)))
  expect_true(all(r > 0.3 & r < 3), info = paste(signif(r, 3), collapse = ", "))
})

test_that("it refuses fits it cannot build scores for", {
  skip_on_cran()
  d <- data_gen_cs(N = 150, rand = 1, sig_u = 0.8, sig_v = 0.3, cons = 0.5,
                   beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.1)
  ## keep_objective defaults to FALSE, so there is no likelihood to difference.
  f0 <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = d)
  expect_error(sfa:::estfun.sfareg(f0), "keep_objective")
  ## The robust divergences do not maximise a log-likelihood at all.
  fr <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = d,
            robust = "mdpd", keep_objective = TRUE)
  expect_error(sfa:::estfun.sfareg(fr), "do not maximise a log-likelihood")
})

test_that("the sandwich methods are REGISTERED, not merely defined", {
  ## This is the test that would have caught a real shipped defect. Both
  ## methods existed and both worked under devtools::load_all(), which exports
  ## everything -- so every test passed while the INSTALLED package failed with
  ## "no applicable method for 'estfun' applied to an object of class sfareg",
  ## taking vcovCL() and coeftest() down with it.
  ##
  ## `sandwich` is in Suggests, so its generics do not exist when sfa loads and
  ## a plain S3method() directive cannot work. The `pkg::generic` form defers
  ## registration until sandwich is loaded. Asserting on NAMESPACE is crude but
  ## it is the only check that survives load_all(): a dispatch test here would
  ## pass whether or not the directive is present.
  ns <- readLines(system.file("NAMESPACE", package = "sfa"), warn = FALSE)
  if (!length(ns)) ns <- readLines(testthat::test_path("..", "..", "NAMESPACE"), warn = FALSE)
  ns <- paste(ns, collapse = "\n")

  expect_match(ns, "S3method\\(sandwich::bread, *sfareg\\)")
  expect_match(ns, "S3method\\(sandwich::estfun, *sfareg\\)")
})
