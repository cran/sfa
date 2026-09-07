## Out-of-domain penalties in the likelihood closures.
##
## Every branch guards against a non-positive scale, but the penalty has to be
## a large FINITE number. optim() differences the objective to get a gradient,
## and differencing .Machine$double.xmax overflows to a non-finite value, which
## aborts the fit with "non-finite finite-difference value" rather than steering
## the search away. NGE used xmax and died on 3 of 45 fits at N = 150 -- including
## at sigma_u = 1, sigma_v = 0.3 -- before this was changed. NU and NE carried
## the same construct.

test_that("no likelihood branch returns .Machine$double.xmax as its penalty", {
  ## Source-level guard, and only meaningful when the sources are on disk --
  ## under R CMD check the package is installed and R/ is gone, so this skips.
  ## readLines() on a missing path ERRORS rather than returning empty, so the
  ## existence check has to come first.
  ## Every file carrying a likelihood closure, not just sfm.R. Checking only
  ## sfm.R is how ttsfm.R kept SEVEN of these until 2026-09-04: TTNE, TTHN and
  ## TTNLS each returned xmax, and TTHN's fits were failing with optim()'s
  ## ABNORMAL_TERMINATION_IN_LNSRCH -- the same non-finite-gradient mechanism
  ## this test exists to prevent, in a file it never read.
  want <- c("sfm.R", "ttsfm.R", "zsfm.R", "lcsfm.R", "psfm.R", "selsfm.R",
    "ivsfm.R", "copsfm.R")
  roots <- c("../../R", "../../../R", "R")
  root <- roots[dir.exists(roots)]
  skip_if(!length(root), "R/ source not reachable (installed package)")
  files <- file.path(root[1], want)
  files <- files[file.exists(files)]
  skip_if(!length(files), "no likelihood sources found")
  src <- unlist(lapply(files, readLines, warn = FALSE))
  ## `return(.Machine$double.xmax)` is the exact construct that overflows
  ## optim()'s finite-difference gradient.
  bad <- grep("return\\(\\.Machine\\$double\\.xmax\\)", src, value = TRUE)
  expect_equal(length(bad), 0)
})

test_that("models with domain guards fit without aborting", {
  skip_on_cran()
  ## These configurations produced "non-finite finite-difference value" before
  ## the penalties were made finite.
  cfgs <- list(c(0.3, 0.4), c(1, 0.3), c(0.2, 0.8))
  for (mn in c("NGE", "NU", "NE")) {
    yv <- switch(mn, NGE = "y_pcs_ge", NU = "y_pcs_u", NE = "y_pcs_e")
    errs <- 0L
    for (cfg in cfgs) {
      for (r in 1:5) {
        d <- data_gen_cs(N = 150, rand = r, sig_u = cfg[1], sig_v = cfg[2], cons = 0.5,
                         beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.1)
        f <- try(suppressWarnings(suppressMessages(
          sfm(as.formula(paste(yv, "~ x1 + x2")), model_name = mn, data = d))), silent = TRUE)
        if (inherits(f, "try-error")) errs <- errs + 1L
      }
    }
    expect_equal(errs, 0L, info = mn)
  }
})
