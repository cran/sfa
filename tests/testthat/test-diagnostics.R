## Optimizer diagnostics: sfa_diagnostics(), plot.sfareg(), and the
## convergence line that print()/summary() now emit.

diag_data <- function(N = 250, seed = 3) {
  data_gen_cs(N = N, rand = seed, sig_u = 1, sig_v = 0.3, cons = 0.5,
              beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.5)
}

test_that("sfa_diagnostics() refuses non-sfareg input", {
  expect_error(sfa_diagnostics(1:5), "class \"sfareg\"")
  expect_error(sfa_diagnostics(list()), "class \"sfareg\"")
})

test_that("a healthy fit reports convergence, a PD Hessian and no flags", {
  skip_on_cran()
  f <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = diag_data())
  d <- sfa_diagnostics(f)

  expect_s3_class(d, "sfadiag")
  expect_equal(d$convergence$code, 0)
  expect_equal(d$convergence$meaning, "converged")
  expect_true(is.finite(d$convergence$logLik))
  expect_equal(d$convergence$logLik, -f$opt$value)

  expect_true(d$hessian$pos_def)
  expect_length(d$hessian$eigenvalues, length(f$coefficients))
  expect_true(is.finite(d$hessian$condition) && d$hessian$condition > 0)

  ## without keep_objective there is no gradient, and that is the only flag
  expect_null(d$gradient)
  expect_false(d$has_objective)
  expect_true(any(grepl("keep_objective", d$flags)))
})

test_that("keep_objective stores the true objective and enables the gradient", {
  skip_on_cran()
  f <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = diag_data(200),
           keep_objective = TRUE)
  expect_true(is.function(f$objective))
  ## the stored closure must reproduce the optimizer's own value
  expect_equal(f$objective(f$opt$par), f$opt$value, tolerance = 1e-8)

  d <- sfa_diagnostics(f)
  expect_true(d$has_objective)
  expect_length(d$gradient$gradient, length(f$coefficients))
  ## a converged fit has a gradient indistinguishable from zero
  expect_lt(d$gradient$max_rel, 1e-3)
  expect_false(any(grepl("keep_objective", d$flags)))
})

test_that("keep_objective defaults to off", {
  skip_on_cran()
  f <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = diag_data(150))
  expect_null(f$objective)
})

test_that("weak identification is detected and named", {
  skip_on_cran()
  ## NNAK's sigma_u and shape are the documented ridge: the convergence work
  ## measured cor = 0.997 across replications. The fit is deliberately
  ## ill-conditioned, so how much of the vcov survives is platform arithmetic,
  ## not a property of the package -- assert the contract, and the ridge itself
  ## whenever the pair carrying it survives.
  f <- suppressWarnings(sfm(y_pcs_nak ~ x1 + x2, model_name = "NNAK",
                            data = diag_data(300)))
  d <- sfa_diagnostics(f)
  skip_if(is.null(d$correlation),
          "fewer than two parameters have a usable variance on this platform")

  R <- d$correlation$correlation
  expect_true(is.matrix(R) && nrow(R) == ncol(R))
  expect_identical(rownames(R), colnames(R))
  expect_true(all(is.finite(R)))
  ## Deliberately NOT asserting abs(R) <= 1. vcov.sfareg() inverts the stored
  ## Hessian, and a near-singular Hessian -- the whole reason this model is
  ## here -- can invert to a matrix with a positive diagonal that is still
  ## indefinite, which yields |r| > 1. That is information about the fit, not
  ## a bug in the diagnostic, and it is platform arithmetic whether it happens.
  expect_equal(unname(diag(R)), rep(1, nrow(R)), tolerance = 1e-8)
  expect_length(d$correlation$worst_pair, 2L)
  expect_true(all(d$correlation$worst_pair %in% rownames(R)))
  expect_true(is.finite(d$correlation$worst_value))
  ## dropped parameters are named, and never overlap with those reported
  expect_length(intersect(d$correlation$dropped, rownames(R)), 0L)

  if (all(c("mu", "sigu") %in% rownames(R))) {
    expect_setequal(d$correlation$worst_pair, c("mu", "sigu"))
    expect_gt(abs(d$correlation$worst_value), 0.95)
    expect_true(any(grepl("correlated at", d$flags)))
  }
})

test_that("the correlation diagnostic survives a partly unusable vcov", {
  ## The CI-visible case, made deterministic: one parameter with no usable
  ## variance must cost only that parameter, not the whole diagnostic.
  V <- diag(4)
  V[2, 2] <- NaN
  pn <- c("sigv", "sigu", "x1", "x2")
  dimnames(V) <- list(pn, pn)
  V[1, 3] <- V[3, 1] <- 0.99
  ## .sfa_corr_diag() only ever asks stats::vcov() for the matrix, so a stub
  ## class with its own method pins the input exactly. Registering is needed
  ## because dispatch happens inside the sfa namespace, which cannot see a
  ## method defined in this test's environment.
  registerS3method("vcov", "sfa_fake_fit", function(object, ...) object$V)
  fake <- structure(list(V = V), class = "sfa_fake_fit")
  cd <- sfa:::.sfa_corr_diag(fake, pn)

  expect_false(is.null(cd))
  expect_identical(rownames(cd$correlation), c("sigv", "x1", "x2"))
  expect_identical(cd$dropped, "sigu")
  expect_setequal(cd$worst_pair, c("sigv", "x1"))
  expect_equal(cd$worst_value, 0.99, tolerance = 1e-12)
})

test_that("a non-zero convergence code is surfaced rather than swallowed", {
  ## Build a fit object whose optimizer failed, without paying for a real one.
  f <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = diag_data(150))
  f$opt$convergence <- 1L
  f$opt$message <- NULL
  d <- sfa_diagnostics(f)
  expect_match(d$convergence$meaning, "ITERATION LIMIT")
  expect_true(any(grepl("code 1", d$flags)))
  ## the iteration limit is never excused, even with a fine Hessian
  expect_false(isTRUE(d$convergence$benign_nonzero))
  expect_false(isTRUE(d$convergence$unverified_nonzero))
  expect_true(any(grepl("not at an optimum", d$flags)))
  ## and print() must say so, where before it said nothing
  out <- paste(capture.output(print(f)), collapse = "\n")
  expect_match(out, "convergence: 1")
  expect_match(out, "not a converged optimum")
})

test_that("a non-zero code with a clean gradient is called benign, not a failure", {
  skip_on_cran()
  ## Measured behaviour, not a guess: L-BFGS-B routinely returns code 52 from
  ## the final polish stage because the staged minimizer already converged and
  ## it cannot take a step. Such a fit must NOT be reported as failed.
  f <- suppressWarnings(sfm(y_pcs ~ x1 + x2, model_name = "NHN",
                            data = diag_data(150), keep_objective = TRUE))
  d <- sfa_diagnostics(f)
  if (identical(as.integer(d$convergence$code), 0L)) skip("this draw converged cleanly")
  expect_lt(d$gradient$max_rel, 1e-3)
  expect_true(d$hessian$pos_def)
  expect_true(isTRUE(d$convergence$benign_nonzero))
  expect_true(any(grepl("treat the code as noise", d$flags)))
  expect_false(any(grepl("is not at an optimum", d$flags)))
})

test_that("a line-search code without the objective is left unverified", {
  skip_on_cran()
  ## No gradient means no positive evidence either way; the tool must say so
  ## rather than quietly resolving it in the fit's favour.
  f <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = diag_data(150))
  f$opt$convergence <- 52L
  d <- sfa_diagnostics(f)
  expect_false(isTRUE(d$convergence$benign_nonzero))
  expect_true(isTRUE(d$convergence$unverified_nonzero))
  expect_true(any(grepl("cannot be confirmed", d$flags)))
})

test_that("a non-zero code WITH a bad gradient is called a failure", {
  skip_on_cran()
  f <- suppressWarnings(sfm(y_pcs ~ x1 + x2, model_name = "NHN",
                            data = diag_data(150), keep_objective = TRUE))
  f$opt$convergence <- 52L
  f$opt$par <- f$opt$par + 5      ## move well off the optimum
  d <- sfa_diagnostics(f)
  expect_false(isTRUE(d$convergence$benign_nonzero))
  expect_true(any(grepl("is not at an optimum", d$flags)))
})

test_that("a singular Hessian is reported, not silently ignored", {
  f <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = diag_data(150))
  H <- f$opt$hessian
  H[, 1] <- 0; H[1, ] <- 0          ## exactly flat in one direction
  f$opt$hessian <- H
  d <- suppressWarnings(sfa_diagnostics(f))
  expect_false(d$hessian$pos_def)
  expect_true(any(grepl("NOT positive definite|singular", d$flags)))
})

test_that("plot() draws, and degrades rather than misleading", {
  skip_on_cran()
  f <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = diag_data(150),
           keep_objective = TRUE)
  tmp <- tempfile(fileext = ".png")
  grDevices::png(tmp, width = 900, height = 700)
  expect_silent(invisible(plot(f)))
  grDevices::dev.off()
  expect_gt(file.size(tmp), 1000)
  unlink(tmp)

  ## no objective: panels 1-2 still work, 3-4 are dropped
  g <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = diag_data(150))
  tmp2 <- tempfile(fileext = ".png")
  grDevices::png(tmp2, width = 700, height = 400)
  expect_silent(invisible(plot(g)))
  grDevices::dev.off()
  unlink(tmp2)
  ## asking only for an impossible panel is an error, not an empty device
  expect_error(plot(g, which = 3), "Nothing to plot")
})
