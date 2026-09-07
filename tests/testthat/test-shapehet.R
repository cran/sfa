## Covariate-dependent SHAPE for NG/NNAK (G5). Experimental -- see ?sfm.
##
## VALIDATION IS DELIBERATELY NOT DONE BY FITTING HERE. `NG` runtimes are both
## long and erratic -- measured at 114s for a plain fit at n = 300 and 7.7s at
## n = 600 on the same machine, with the shape block reversing that ordering --
## so any fitting test would be slow and flaky, and a flaky test in a suite
## that gates CI is worse than no test. The mechanism was checked by hand and
## the evidence, including the recovery FAILURE, is recorded in ?sfm and in
## horserace/FUNCTIONALITY_GAPS.md rather than asserted here.
##
## What is tested is the part that is fast, deterministic, and most likely to
## rot: the argument's contract.

test_that("shapehet is refused where there is no shape parameter to vary", {
  set.seed(1)
  n <- 200
  x1 <- rnorm(n); x2 <- rnorm(n); z <- rnorm(n)
  y <- 1 + 0.5 * x1 + 0.5 * x2 + rnorm(n, 0, 0.3) - abs(rnorm(n))
  d <- data.frame(y, x1, x2, z)

  ## Only the gamma and Nakagami families have a shape at all. For a
  ## half-normal or exponential u, scale heterogeneity is what is meant and
  ## `uhet` already does it -- so this errors rather than silently accepting an
  ## argument it cannot honour.
  for (m in c("NHN", "NE", "NTN", "NR")) {
    expect_error(
      sfm(y ~ x1 + x2, model_name = m, data = d, shapehet = ~z),
      "parameterizes the SHAPE",
      info = m
    )
  }
  expect_error(
    sfm(y ~ x1 + x2, model_name = "NG", data = d, shapehet = "z"),
    "must be a one-sided formula"
  )
  ## A constant factor on the shape IS the shape, so an intercept-only formula
  ## is a specification error rather than a no-op.
  expect_error(
    sfm(y ~ x1 + x2, model_name = "NG", data = d, shapehet = ~1),
    "at least one covariate"
  )
})
