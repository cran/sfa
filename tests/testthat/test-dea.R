## Output-oriented DEA envelopment, .dea_out().
##
## Solved in the package so that npsfm(method = "SZ") does not need an outside
## frontier package installed. During development the implementation was
## checked against an independent reference over all four RTS settings at
## n = 30/80/150 with one to three inputs and one to two outputs, agreeing to
## 3e-12 or better; see horserace/FUNCTIONALITY_GAPS.md entry B5 for the
## recipe. That comparison is deliberately NOT kept as a test -- it would
## reintroduce the dependency it exists to remove -- so what is pinned here is
## the analytic structure instead: a hand-solvable case, the returns-to-scale
## ordering, and the invariants any correct output-oriented envelopment must
## satisfy.

test_that(".dea_out solves a case that can be worked out by hand", {
  skip_if_not_installed("lpSolve")
  ## One input, one output. Points 1-3 lie on the frontier; point 4 uses the
  ## same input as point 2 but produces half as much.
  X <- matrix(c(1, 2, 3, 2), ncol = 1)
  Y <- matrix(c(1, 2, 2.5, 1), ncol = 1)
  th <- sfa:::.dea_out(X, Y, "vrs")

  ## Under VRS the best output reachable with input <= 2 is 2 (point 2
  ## itself; any convex mix with point 3 needs more input). So point 4 scores
  ## 2/1 = 2 and the three frontier points score 1.
  expect_equal(th, c(1, 1, 1, 2), tolerance = 1e-6)
})

test_that("returns to scale order the scores as they must", {
  skip_if_not_installed("lpSolve")
  set.seed(11)
  n <- 40
  X <- matrix(runif(n * 2, 1, 10), n, 2)
  Y <- matrix(runif(n, 1, 10), n, 1)
  th <- lapply(c("vrs", "drs", "irs", "crs"),
               function(r) sfa:::.dea_out(X, Y, r))
  names(th) <- c("vrs", "drs", "irs", "crs")

  ## VRS is the tightest envelope and CRS the loosest, so VRS gives the
  ## smallest expansion factor and CRS the largest, with DRS and IRS between.
  tol <- 1e-6
  expect_true(all(th$vrs <= th$drs + tol))
  expect_true(all(th$vrs <= th$irs + tol))
  expect_true(all(th$drs <= th$crs + tol))
  expect_true(all(th$irs <= th$crs + tol))
})

test_that("every score is a feasible expansion factor", {
  skip_if_not_installed("lpSolve")
  set.seed(12)
  n <- 30
  X <- matrix(runif(n * 3, 1, 10), n, 3)
  Y <- matrix(runif(n * 2, 1, 10), n, 2)
  for (r in c("vrs", "crs", "drs", "irs")) {
    th <- sfa:::.dea_out(X, Y, r)
    expect_length(th, n)
    expect_true(all(is.finite(th)), info = r)
    ## Output orientation: no unit can be contracted, so theta >= 1.
    expect_true(all(th >= 1 - 1e-6), info = r)
    ## At least one unit must define the frontier and score exactly 1.
    expect_true(any(abs(th - 1) < 1e-6), info = r)
  }
})

test_that("a unit that dominates every other is efficient under all RTS", {
  skip_if_not_installed("lpSolve")
  ## Unit 1 uses the least input and makes the most output, so nothing can
  ## dominate it and its score is 1 whatever the scale assumption.
  X <- matrix(c(1, 5, 6, 7), ncol = 1)
  Y <- matrix(c(9, 3, 2, 1), ncol = 1)
  for (r in c("vrs", "crs", "drs", "irs")) {
    expect_equal(sfa:::.dea_out(X, Y, r)[1], 1, tolerance = 1e-6, info = r)
  }
})

test_that("npsfm(method = \"SZ\") runs on the internal solver", {
  skip_on_cran()
  skip_if_not_installed("np")
  skip_if_not_installed("lpSolve")
  d <- cs_small(N = 120)
  fit <- suppressWarnings(npsfm(y_pcs ~ x1 + x2, data = d, method = "SZ"))
  expect_s3_class(fit, "npsfareg")
  expect_length(fit$frontier, nrow(d))
  expect_true(all(is.finite(fit$frontier)))
  ## The DEA step can only push the frontier up, never down.
  expect_true(all(fit$dea.efficiency >= 1 - 1e-6))
})
