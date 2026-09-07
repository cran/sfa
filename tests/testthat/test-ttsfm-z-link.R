## ttsfm()'s variance-determinant scale.
##
## Two-tier models put z'delta on the standard deviation, like sfm() and unlike
## psfm() and the competing packages -- gap C1. `z_link` lets a caller align
## them; the default preserves ttsfm()'s own convention.
##
## Every fit here sets rand.psoptim. That is defensive rather than load-bearing
## -- PSopt defaults to FALSE, so the particle swarm does not actually run --
## but a fit that DOES enable it draws from the session RNG and then returns a
## different answer every session (gap I11). Comparing two builds with unseeded
## fits produces pure noise; that mistake cost three wrong conclusions before it
## was spotted. Do not remove the seed.
##
## TTHN is quarantined behind SFA_TEST_SLOW, matching test-zsfm-ttsfm.R. The
## seed does NOT protect it: its instability is in optim()'s final stage, which
## is deterministic given the data but follows a different trajectory under a
## different BLAS. Run ungated, it passed locally and failed R CMD check on
## macOS and Windows with "ABNORMAL_TERMINATION_IN_LNSRCH" -- red CI from
## 2026-08-28 until this gate was added.

tt_data <- function(N = 800) {
  d <- data_gen_cs(N = N, rand = 2, sig_u = 0.6, sig_v = 0.3, cons = 0.5,
                   beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.1)
  set.seed(99)
  d$z <- runif(nrow(d), -1, 1)
  d$zp <- runif(nrow(d), -1, 1)
  d
}

test_that("the default is the SD link and is unchanged", {
  skip_on_cran()
  d <- tt_data()
  a <- suppressWarnings(ttsfm(y_ttne ~ x1 + x2 | z | zp, model_name = "TTNE",
                              data = d, rand.psoptim = 7L))
  b <- suppressWarnings(ttsfm(y_ttne ~ x1 + x2 | z | zp, model_name = "TTNE",
                              data = d, rand.psoptim = 7L, z_link = "sd"))
  expect_equal(coef(a), coef(b))
  expect_equal(as.numeric(logLik(a)), as.numeric(logLik(b)))
})

## Factored out so TTNE can run by default and TTHN only under SFA_TEST_SLOW.
.tt_link_check <- function(d, m) {
    sd_l <- suppressWarnings(ttsfm(y_ttne ~ x1 + x2 | z | zp, model_name = m,
                                   data = d, rand.psoptim = 7L, z_link = "sd"))
    vr_l <- suppressWarnings(ttsfm(y_ttne ~ x1 + x2 | z | zp, model_name = m,
                                   data = d, rand.psoptim = 7L, z_link = "var"))
    ## Same model, so the maximised likelihood must agree.
    expect_equal(as.numeric(logLik(sd_l)), as.numeric(logLik(vr_l)),
                 tolerance = 1e-4, info = m)
    ## sigma = exp(eta) against exp(eta/2): the SD-link delta is HALF.
    expect_equal(unname(coef(sd_l)[["z"]]) / unname(coef(vr_l)[["z"]]), 0.5,
                 tolerance = 0.02, info = m)
    ## The second tier moves the same way.
    expect_equal(unname(coef(sd_l)[["zp"]]) / unname(coef(vr_l)[["zp"]]), 0.5,
                 tolerance = 0.02, info = m)
    ## The frontier is untouched by the reparameterization.
    expect_equal(unname(coef(sd_l)[["x1"]]), unname(coef(vr_l)[["x1"]]),
                 tolerance = 1e-3, info = m)
    invisible(TRUE)
}

test_that("the two links are reparameterizations and delta halves exactly", {
  skip_on_cran()
  .tt_link_check(tt_data(), "TTNE")
})

test_that("the same holds for the half-normal two-tier model", {
  ## See the file header: TTHN is not stable enough to run by default.
  skip_if(!nzchar(Sys.getenv("SFA_TEST_SLOW")),
          "TTHN is slow and numerically fragile; set SFA_TEST_SLOW to run it")
  .tt_link_check(tt_data(N = 300), "TTHN")
})

test_that("an invalid link is refused", {
  d <- tt_data(N = 200)
  expect_error(ttsfm(y_ttne ~ x1 + x2 | z | zp, model_name = "TTNE", data = d,
                     z_link = "logit"), "should be one of")
})
