## sigma_u on the zero boundary under wrong skew.
##
## This is the CORRECT maximum likelihood estimate, not a numerical failure
## (Waldman 1982; Olson, Schmidt and Waldman 1980 call it the Type I failure),
## so the package must report it rather than bound sigma_u away from zero. A
## bound would corrupt precisely the samples where the boundary is the answer.
## Before this was added the MLE path said nothing at all -- the user saw only
## a run of uninformative "NaNs produced" warnings from the optimizer.

test_that("the boundary case is flagged, and only when skew is wrong", {
  skip_on_cran()
  ## lambda = 0.75 at N = 100 puts a useful fraction of samples in the wrong
  ## skew region. Seed 4 is one of them; seeds 1-3 are not.
  flagged <- list()
  for (r in 1:6) {
    d <- data_gen_cs(N = 100, rand = r, sig_u = 0.3, sig_v = 0.4, cons = 0.5,
                     beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.1)
    f <- suppressWarnings(sfm(y_pcs_e ~ x1 + x2, model_name = "NE", data = d))
    expect_false(is.null(f$wrong_skew), info = paste("seed", r))
    expect_false(is.null(f$sigma_u_at_bound), info = paste("seed", r))
    ## The two must agree with the residual moment the fit actually saw.
    expect_equal(f$wrong_skew, f$residual_m3 >= 0, info = paste("seed", r))
    ## sigma_u only ever collapses under wrong skew; a correctly skewed sample
    ## reaching the boundary would be a genuine numerical failure.
    if (isTRUE(f$sigma_u_at_bound)) {
      expect_true(f$wrong_skew, info = paste("seed", r))
      flagged[[length(flagged) + 1L]] <- r
    }
  }
  expect_gt(length(flagged), 0)
})

test_that("the warning fires on the boundary case and stays quiet otherwise", {
  skip_on_cran()
  grab <- function(r) {
    d <- data_gen_cs(N = 100, rand = r, sig_u = 0.3, sig_v = 0.4, cons = 0.5,
                     beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.1)
    w <- character(0)
    f <- withCallingHandlers(
      suppressMessages(sfm(y_pcs_e ~ x1 + x2, model_name = "NE", data = d)),
      warning = function(x) {
        m <- conditionMessage(x)
        if (!grepl("NaNs produced", m)) w <<- c(w, m)
        invokeRestart("muffleWarning")
      })
    list(fit = f, warns = w)
  }
  bad <- grab(4)
  expect_true(isTRUE(bad$fit$sigma_u_at_bound))
  expect_true(any(grepl("Olson, Schmidt and Waldman", bad$warns)))
  expect_true(any(grepl("boundary IS the maximum likelihood estimate", bad$warns)))

  good <- grab(1)
  expect_false(isTRUE(good$fit$sigma_u_at_bound))
  expect_false(any(grepl("Olson, Schmidt and Waldman", good$warns)))
})

test_that("the boundary threshold is not knife-edge", {
  skip_on_cran()
  ## The point of this test is the MARGIN, not the classification. With
  ## rel_tol = 1e-3 the collapsed sample cleared the threshold by a factor of
  ## 1.16 while the interior ones cleared it by ~470 -- and the optimizer's
  ## stopping point moves by more than 16% between BLAS implementations, so it
  ## passed locally and failed R CMD check on CI's macOS. Every sample must now
  ## sit at least 5x clear of the threshold on whichever side it belongs.
  ratios <- vapply(1:6, function(r) {
    d <- data_gen_cs(N = 100, rand = r, sig_u = 0.3, sig_v = 0.4, cons = 0.5,
                     beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.1)
    f <- suppressWarnings(sfm(y_pcs_e ~ x1 + x2, model_name = "NE", data = d))
    X <- stats::model.matrix(y_pcs_e ~ x1 + x2, d)
    unname(coef(f)[["sigu"]]) / stats::sd(stats::lm.fit(X, d$y_pcs_e)$residuals)
  }, numeric(1))
  tol <- eval(formals(sfa:::.wrong_skew_boundary)$rel_tol)
  margin <- ifelse(ratios < tol, tol / ratios, ratios / tol)
  expect_true(all(margin > 5), info = paste(round(margin, 1), collapse = ", "))
  ## And the design still contains both kinds of sample, or it tests nothing.
  expect_true(any(ratios < tol))
  expect_true(any(ratios > tol))
})

test_that("the helper refuses to guess when it has nothing to go on", {
  ## Degenerate inputs must return NA rather than a confident FALSE.
  z <- sfa:::.wrong_skew_boundary(c(1, 2), 0.5, 1, "NE")
  expect_true(is.na(z$wrong_skew))
  z2 <- sfa:::.wrong_skew_boundary(rnorm(50), NA_real_, 1, "NE")
  expect_true(is.na(z2$at_bound))
  z3 <- sfa:::.wrong_skew_boundary(rnorm(50), 0.5, 0, "NE")
  expect_true(is.na(z3$at_bound))
})
