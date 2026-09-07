## A persistent scale on the zero boundary in a GTRE fit.

test_that("the helper judges a collapse relative to scale", {
  expect_true(sfa:::.panel_scale_at_bound(1e-7, 1.0))
  expect_false(sfa:::.panel_scale_at_bound(0.4, 1.0))
  ## Relative, not absolute: 0.005 is a collapse against a scale of 1 but not
  ## against a scale of 0.05.
  expect_true(sfa:::.panel_scale_at_bound(0.005, 1.0))
  expect_false(sfa:::.panel_scale_at_bound(0.005, 0.05))
  ## Degenerate inputs return NA rather than a confident answer.
  expect_true(is.na(sfa:::.panel_scale_at_bound(NA_real_, 1)))
  expect_true(is.na(sfa:::.panel_scale_at_bound(0.1, 0)))
  expect_true(is.na(sfa:::.panel_scale_at_bound(0.1, NA_real_)))
})

test_that("GTRE reports a collapsed persistent scale, in both directions", {
  skip_on_cran()
  ## Two seeds that collapse in OPPOSITE directions -- 1005 sends sigh to the
  ## boundary, 1001 sends sigr. Both are genuine maximum likelihood solutions:
  ## on seed 1005 the boundary fit beats the TRUE parameter vector by 3.66 log
  ## units. See notes/sigh_investigation_2026-08-29.md.
  grab <- function(seed) {
    d <- data_gen_p(t = 6, N = 200, rand = seed, sig_u = 1, sig_v = 0.3,
                    sig_r = 0.2, sig_h = 0.4, cons = 0.5, beta1 = 0.5,
                    beta2 = 0.5)
    pd <- plm::pdata.frame(d, index = c("name", "year"))
    w <- character(0)
    f <- withCallingHandlers(
      psfm(y_gtre ~ x1 + x2, "GTRE", pd, individual = "name",
           halton_num = 50, rand.gtre = 7L, estimator = "sml"),
      warning = function(x) {
        w <<- c(w, conditionMessage(x))
        invokeRestart("muffleWarning")
      })
    list(fit = f, warns = w)
  }

  a <- grab(1005)
  expect_true(isTRUE(a$fit$sigh_at_bound))
  expect_false(isTRUE(a$fit$sigr_at_bound))
  expect_true(any(grepl("sigh has collapsed", a$warns)))
  ## The warning must name the OTHER scale, because they have to be read
  ## together -- the variation is reclassified, not lost.
  expect_true(any(grepl("sigr = ", a$warns)))

  b <- grab(1001)
  expect_true(isTRUE(b$fit$sigr_at_bound))
  expect_false(isTRUE(b$fit$sigh_at_bound))
  expect_true(any(grepl("sigr has collapsed", b$warns)))
  expect_true(any(grepl("sigh = ", b$warns)))
})

test_that("the warning explains that the boundary can be the correct MLE", {
  skip_on_cran()
  d <- data_gen_p(t = 6, N = 200, rand = 1005, sig_u = 1, sig_v = 0.3,
                  sig_r = 0.2, sig_h = 0.4, cons = 0.5, beta1 = 0.5, beta2 = 0.5)
  pd <- plm::pdata.frame(d, index = c("name", "year"))
  expect_warning(
    psfm(y_gtre ~ x1 + x2, "GTRE", pd, individual = "name", halton_num = 50,
         rand.gtre = 7L, estimator = "sml"),
    "CORRECT maximum likelihood estimate"
  )
})

test_that("psfm_bootstrap warns when bootstrapping a boundary fit", {
  skip_on_cran()
  ## A parametric bootstrap resamples from the FITTED model, so a collapsed
  ## persistent scale means resampling from a DGP without that component --
  ## and the bootstrap is inconsistent on the boundary of the parameter space
  ## regardless. The interval would understate the uncertainty in the
  ## persistent split rather than represent it.
  d <- data_gen_p(t = 6, N = 200, rand = 1005, sig_u = 1, sig_v = 0.3,
                  sig_r = 0.2, sig_h = 0.4, cons = 0.5, beta1 = 0.5, beta2 = 0.5)
  pd <- plm::pdata.frame(d, index = c("name", "year"))
  f <- suppressWarnings(
    psfm(y_gtre ~ x1 + x2, "GTRE", pd, individual = "name", halton_num = 50,
         rand.gtre = 7L, estimator = "sml"))
  expect_true(isTRUE(f$sigh_at_bound))

  w <- character(0)
  invisible(withCallingHandlers(
    try(psfm_bootstrap(f, B = 2, numCores = 1, individual = "name"),
        silent = TRUE),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    }))
  expect_true(any(grepl("persistent scale on the zero boundary", w)))
})

## ---------------------------------------------------------------------------
## The same check on the DEFAULT estimator.
##
## psfm(model_name = "GTRE") returns model_name = "GTRE_FML" under its default
## estimator = "fiml", and that branch packs its results separately -- so the
## boundary block, which keyed on model_name == "GTRE", never ran for the path
## most users take. The FIML archive puts the collapse rate at 19.3% at
## N = 100, T = 10, so this was silent on roughly one fit in five.
## ---------------------------------------------------------------------------

fml_fit <- function(seed, sig_r = 0.2, sig_h = 0.4) {
  d <- data_gen_p(t = 10, N = 100, rand = seed, sig_u = 1, sig_v = 0.3,
                  sig_r = sig_r, sig_h = sig_h, cons = 0.5, beta1 = 0.5,
                  beta2 = 0.5)
  w <- character(0)
  f <- withCallingHandlers(
    psfm(y_gtre ~ x1 + x2, "GTRE", as.data.frame(d), individual = "name",
         estimator = "fiml"),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  list(fit = f, warns = w)
}

## ---------------------------------------------------------------------------
## WHY THIS TEST DOES NOT PIN A SEED THAT COLLAPSES (2026-08-30)
##
## It used to, and CI went red: seeds 5 and 6 from the FIML archive collapsed
## sigh and sigr respectively here and on ubuntu-devel and macOS, but seed 6
## did NOT collapse on ubuntu-release or windows-release. The old comment
## defended the assertion on the grounds that the collapse cleared the 1%
## threshold by three orders of magnitude -- but that argues about precision AT
## an optimum, not about which optimum a different BLAS walks to. A collapse
## sample is a flat, wrong-skew region of the likelihood; which basin the
## optimizer lands in is exactly what is not portable.
##
## So the two directions are now tested by the means each one admits:
##
##   sigh: forced by the DGP. Generating with sig_h = 0 -- no persistent
##         inefficiency at all, which is a case users really do hit -- puts
##         sigh on the boundary in 7 of 8 seeds, at 1e-3 or below. That is a
##         property of the data, not of the optimizer.
##
##   sigr: NOT forceable, and the asymmetry is the finding. Generating with
##         sig_r = 0 collapses sigr in only 1 of 8 seeds; the MLE returns
##         sigr = 0.06-0.22 against a true 0. A symmetric sigr absorbs part of
##         the between-firm variance at almost no likelihood cost, because the
##         r-versus-h split is the weakly identified thing. More firms shrink
##         it (0.18 -> 0.015 from N = 100 to 300) without driving it to the
##         bound. So the sigr direction is covered where it is deterministic:
##         .report_panel_boundary() itself, below.
## ---------------------------------------------------------------------------

test_that("the default FIML estimator reports a collapsed persistent scale", {
  skip_on_cran()
  ## sig_h = 0: the DGP has no persistent inefficiency, so the boundary is the
  ## right answer on every platform rather than one basin among several.
  a <- fml_fit(3, sig_r = 0.4, sig_h = 0)
  expect_identical(a$fit$model_name, "GTRE_FML")
  expect_true(isTRUE(a$fit$sigh_at_bound))
  expect_false(isTRUE(a$fit$sigr_at_bound))
  expect_true(any(grepl("sigh has collapsed", a$warns)))
  ## The warning must name the other scale, as in the SML case above.
  expect_true(any(grepl("sigr = ", a$warns)))
})

test_that("the boundary reporter flags and explains a collapsed sigr", {
  ## The sigr direction, deterministically. .report_panel_boundary() is what
  ## both estimator branches call, so this covers the reporting contract
  ## without asking an optimizer to land in a particular basin.
  w <- character(0)
  res <- withCallingHandlers(
    sfa:::.report_panel_boundary(list(), "GTRE_FML",
                                 sig_h = 0.52, sig_r = 1.1e-06, ref = 1.04),
    warning = function(x) {
      w <<- c(w, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  expect_true(isTRUE(res$sigr_at_bound))
  expect_false(isTRUE(res$sigh_at_bound))
  expect_true(any(grepl("sigr has collapsed", w)))
  expect_true(any(grepl("GTRE_FML", w)))
  ## Names the absorbing scale and its value, so the two are read together.
  expect_true(any(grepl("sigh = 0.52", w)))
  expect_true(any(grepl("CORRECT maximum likelihood estimate", w)))

  ## And the mirror image: only one of the two can be reported at a time.
  w2 <- character(0)
  res2 <- withCallingHandlers(
    sfa:::.report_panel_boundary(list(), "GTRE",
                                 sig_h = 9.1e-04, sig_r = 0.31, ref = 1.04),
    warning = function(x) {
      w2 <<- c(w2, conditionMessage(x))
      invokeRestart("muffleWarning")
    })
  expect_true(isTRUE(res2$sigh_at_bound))
  expect_false(isTRUE(res2$sigr_at_bound))
  expect_true(any(grepl("sigh has collapsed", w2)))
  expect_true(any(grepl("sigr = 0.31", w2)))

  ## An interior fit sets both flags to FALSE and warns about neither.
  expect_silent(res3 <- sfa:::.report_panel_boundary(
    list(), "GTRE_FML", sig_h = 0.40, sig_r = 0.20, ref = 1.04))
  expect_false(res3$sigh_at_bound)
  expect_false(res3$sigr_at_bound)
})

test_that("an interior FIML fit reports FALSE rather than nothing at all", {
  skip_on_cran()
  ## The distinction that matters: absent is not the same as FALSE. Before the
  ## fix these fields did not exist on a FIML fit, so isTRUE(NULL) was FALSE
  ## and the condition read as "no collapse" whether or not there was one.
  a <- fml_fit(1)
  expect_false(is.null(a$fit$sigh_at_bound))
  expect_false(is.null(a$fit$sigr_at_bound))
  expect_false(a$fit$sigh_at_bound)
  expect_false(a$fit$sigr_at_bound)
})

test_that("the boundary warning says the intercept moved too", {
  skip_on_cran()
  ## E[y] = beta0 - E[h] - E[u], so collapsing sigma_h removes
  ## E[h] = sigh*sqrt(2/pi) from the error and the intercept absorbs it. A user
  ## reading a level off coef() has to be told that, not only that the split is
  ## unidentified.
  a <- fml_fit(5)
  expect_true(any(grepl("frontier intercept has shifted", a$warns)))
})

test_that("psfm_bootstrap accepts a fit from the default estimator", {
  skip_on_cran()
  ## psfm(model_name = "GTRE") returns "GTRE_FML", which psfm_bootstrap()
  ## rejected outright -- so the natural two-line sequence (fit, then
  ## bootstrap) failed on the default path.
  d <- data_gen_p(t = 6, N = 40, rand = 11, sig_u = 1, sig_v = 0.3,
                  sig_r = 0.2, sig_h = 0.4, cons = 0.5, beta1 = 0.5, beta2 = 0.5)
  f <- suppressWarnings(
    psfm(y_gtre ~ x1 + x2, "GTRE", as.data.frame(d), individual = "name",
         estimator = "fiml"))
  expect_identical(f$model_name, "GTRE_FML")

  ## No `inefdec` here on purpose: psfm() defaults it to TRUE and never asked
  ## the user about it, so psfm_bootstrap() must recover it from the fit
  ## rather than demand it back.
  b <- suppressWarnings(
    psfm_bootstrap(f, BOOT = 2, numCores = 1, individual = "name"))
  ## The FIML layout is (beta..., sigr, sigv, sigh, sigu), so the bootstrap
  ## must carry one column per parameter, named from this model's own $out --
  ## reading the GTRE row order here would silently mislabel every column.
  expect_equal(nrow(b$boot_par), 2L)
  expect_identical(colnames(b$boot_par)[seq_len(nrow(f$out))], rownames(f$out))
  expect_true(all(c("sigr", "sigv", "sigh", "sigu") %in% colnames(b$boot_par)))
  ## $H is returned by this model, so the persistent-score draws come with it.
  expect_false(is.null(b$boot_eff_h))
  expect_equal(nrow(b$boot_eff_h), 2L)
})


test_that("psfm_bootstrap recovers inefdec from the fit rather than demanding it", {
  skip_on_cran()
  ## Omitting `inefdec` used to fail inside parallel::clusterExport(): the
  ## promise was forced during serialization, so the error arrived from
  ## postNode/sendData/serialize with nothing naming the user's call. And
  ## requiring it was the wrong contract -- a user who restates it can restate
  ## it WRONG, and bootstrapping a cost frontier against a production fit is a
  ## silent error, not a loud one.
  d <- as.data.frame(
    data_gen_p(t = 6, N = 40, rand = 11, sig_u = 1, sig_v = 0.3, sig_r = 0.2,
               sig_h = 0.4, cons = 0.5, beta1 = 0.5, beta2 = 0.5))

  ## Cost frontier: the fit NAMES inefdec = FALSE, so the bootstrap must pick
  ## that up. Getting it wrong would resample from a differently oriented
  ## frontier and report perfectly ordinary-looking intervals.
  g <- suppressWarnings(
    psfm(c_gtre ~ x1 + x2, "GTRE", d, individual = "name",
         estimator = "fiml", inefdec = FALSE))
  expect_true("inefdec" %in% names(as.list(g$call)))
  bg <- suppressWarnings(
    psfm_bootstrap(g, BOOT = 2, numCores = 1, individual = "name"))
  expect_equal(nrow(bg$boot_par), 2L)

  ## The other required arguments now say which one is missing, in place of an
  ## error from inside the cluster.
  f <- suppressWarnings(
    psfm(y_gtre ~ x1 + x2, "GTRE", d, individual = "name", estimator = "fiml"))
  expect_error(psfm_bootstrap(f, BOOT = 2, individual = "name"),
               "`numCores` has no default")
})

test_that("print() and summary() say when a scale is on the boundary", {
  ## Tested on the flags rather than through a fit: the reporting logic is what
  ## is under test, and a warning at fit time does not survive saveRDS(), a new
  ## session, or a loop that suppressed warnings -- which is exactly why the
  ## printed output has to carry it.
  fake <- structure(
    list(out = matrix(0, 1, 3), formula = y ~ x, total_time = "0s",
         sigh_at_bound = TRUE, sigr_at_bound = FALSE),
    class = "sfareg")
  expect_output(sfa:::.sfa_report_boundary(fake), "sigh is on the zero boundary")
  expect_output(sfa:::.sfa_report_boundary(fake), "sigr has absorbed its variance")
  expect_output(sfa:::.sfa_report_boundary(fake), "absorbed its mean")

  ## The other direction names the other scale, not the same one twice.
  fake$sigh_at_bound <- FALSE; fake$sigr_at_bound <- TRUE
  expect_output(sfa:::.sfa_report_boundary(fake), "sigr is on the zero boundary")
  expect_output(sfa:::.sfa_report_boundary(fake), "sigh has absorbed its variance")

  ## An interior fit says nothing at all.
  fake$sigr_at_bound <- FALSE
  expect_silent(sfa:::.sfa_report_boundary(fake))

  ## The cross-sectional analogue, and it cites the test that gives a p-value.
  cs <- structure(list(sigma_u_at_bound = TRUE, wrong_skew = TRUE), class = "sfareg")
  expect_output(sfa:::.sfa_report_boundary(cs), "wrong skew")
  expect_output(sfa:::.sfa_report_boundary(cs), "Waldman 1982")
  expect_output(sfa:::.sfa_report_boundary(cs), "skewness_test")

  ## tHN reports the same condition under its own field name and for a
  ## different reason, so it must be covered AND must not cite wrong skewness.
  th <- structure(list(thn_sigma_u_at_bound = TRUE), class = "sfareg")
  expect_output(sfa:::.sfa_report_boundary(th), "sigma_u is on the zero boundary")
  expect_output(sfa:::.sfa_report_boundary(th), "Heavy-")
  expect_failure(expect_output(sfa:::.sfa_report_boundary(th), "Waldman"))

  ## A model with none of these fields must not error on the missing ones.
  expect_silent(sfa:::.sfa_report_boundary(structure(list(), class = "sfareg")))
})
