## opt.optim()'s third-stage guard (I12).
##
## Stage 3 is a polish pass over a point stages 1-2 already produced. Before
## this guard it could throw "non-finite value supplied by optim" straight out
## of a user's fit: NG's default start values did it on 2 of 12 seeds of the
## test-start-from.R DGP, and NNAK did it on 2 of 30 replications at N = 100
## (data_gen_cs, y_pcs_nak) while optim formed the numerical Hessian. Those
## rates are measurements, not assertions -- which seeds fail is a
## platform-specific optimizer path and is deliberately not tested. What IS
## asserted is the property the guard provides: opt.optim() always returns an
## optim-shaped result, never an error.

## Two bounds wide enough that L-BFGS-B has room to move off the start.
lo <- c(-5, -5)
hi <- c(5, 5)
p0 <- c(1, 1)

call_stage3 <- function(fn, opt.TF = TRUE, optHessian = TRUE) {
  opt.optim(
    fn = fn, start_v = p0, lower.optim = lo, upper.optim = hi,
    maxit.optim = 100, opt.TF = opt.TF, method = "L-BFGS-B",
    optHessian = optHessian, trace = 0, verbose = FALSE
  )
}

test_that("an objective that throws away from the start does not abort the stage", {
  fn <- function(p) if (identical(p, p0)) 10 else stop("objective unusable here")

  res <- expect_warning(call_stage3(fn), "final optim\\(\\) stage")

  ## Falls back to the point stage 2 handed in, not to nothing.
  expect_identical(res$start_v, p0)
  expect_equal(res$start_feval, 10)
  expect_equal(res$opt$par, p0)
  expect_equal(res$opt$value, 10)
  expect_identical(res$opt$convergence, 99L)
  ## Hessian is present and correctly shaped so downstream SE code can index
  ## it; it is all-NA so those standard errors come back NA rather than wrong.
  expect_true(is.matrix(res$opt$hessian))
  expect_identical(dim(res$opt$hessian), c(2L, 2L))
  expect_true(all(is.na(res$opt$hessian)))
})

test_that("an objective that goes non-finite away from the start does not abort the stage", {
  fn <- function(p) if (identical(p, p0)) 10 else NaN

  res <- expect_warning(call_stage3(fn), "final optim\\(\\) stage")

  expect_identical(res$opt$convergence, 99L)
  expect_true(is.finite(res$opt$value))
  expect_equal(res$opt$par, p0)
})

test_that("a well-behaved objective is untouched by the guard", {
  fn <- function(p) sum((p - 0.25)^2)

  res <- expect_silent(call_stage3(fn))

  ## No claim about WHICH optimum is reached -- only that stage 3 ran, reported
  ## an optim code rather than the fallback one, and did not make things worse.
  expect_false(identical(res$opt$convergence, 99L))
  expect_true(is.finite(res$opt$value))
  expect_lte(res$opt$value, fn(p0) + 1e-8)
  expect_true(all(is.finite(res$opt$hessian)))
})

test_that("opt.TF = FALSE still means stage 3 is skipped entirely", {
  fn <- function(p) sum((p - 0.25)^2)

  res <- expect_silent(call_stage3(fn, opt.TF = FALSE))

  expect_null(res$opt)
  expect_identical(res$start_v, p0)
  expect_equal(res$start_feval, fn(p0))
})

test_that("with optHessian = FALSE the fallback carries no Hessian", {
  fn <- function(p) if (identical(p, p0)) 10 else stop("objective unusable here")

  res <- expect_warning(call_stage3(fn, optHessian = FALSE), "final optim\\(\\) stage")

  expect_identical(res$opt$convergence, 99L)
  expect_null(res$opt$hessian)
})

test_that("code 99 is reported by name, not as an unknown optim code", {
  fit <- list(
    opt = list(convergence = 99L, message = "test"),
    coefficients = c(a = 1)
  )
  out <- utils::capture.output(sfa:::.sfa_report_convergence(fit))
  expect_true(any(grepl("final optim\\(\\) polish stage", out)))
  expect_false(any(grepl("see \\?optim", out)))
})
