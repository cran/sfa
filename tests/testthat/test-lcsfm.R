## The latent class stochastic frontier: lcsfm(model_name = "LCM" | "LCM_Z").

## A two-class mixture with well separated frontiers. Separation is deliberate:
## the point of these tests is that the machinery recovers a mixture it can
## see, not that it can resolve one it cannot.
lcm_data <- function(seed, n = 500) {
  set.seed(seed)
  x1 <- rnorm(n)
  x2 <- rnorm(n)
  cls <- rbinom(n, 1, 0.4) + 1L
  b <- rbind(c(1, 0.5, 0.5), c(3, 1.0, 0.2))
  sv <- c(0.2, 0.2)
  su <- c(0.5, 1.0)
  y <- b[cls, 1] + b[cls, 2] * x1 + b[cls, 3] * x2 +
    rnorm(n, 0, sv[cls]) - abs(rnorm(n, 0, su[cls]))
  list(d = data.frame(y, x1, x2), cls = cls)
}

test_that("LCM returns a well-formed sfareg with the right parameter layout", {
  dat <- lcm_data(1, n = 300)
  f <- lcsfm(y ~ x1 + x2, "LCM", dat$d, n_class = 2)

  expect_s3_class(f, "sfareg")
  expect_identical(f$model_name, "LCM")
  expect_identical(f$n_class, 2L)
  ## J blocks of (sigv, sigu, beta) plus (J-1) logit intercepts.
  expect_identical(nrow(f$out), 2L * (2L + 3L) + 1L)
  expect_identical(colnames(f$out), c("par", "st_err", "t-val"))
  expect_identical(
    rownames(f$out),
    c("sigv_class1", "sigu_class1", "(Intercept)_class1", "x1_class1", "x2_class1",
      "sigv_class2", "sigu_class2", "(Intercept)_class2", "x1_class2", "x2_class2",
      "logit_(Intercept)_class1")
  )
  ## The scales enter through abs(), so the sign is not identified; they are
  ## reported as magnitudes rather than as whatever sign the optimizer left.
  expect_true(all(f$out[c("sigv_class1", "sigu_class1",
                          "sigv_class2", "sigu_class2"), "par"] > 0))
})

test_that("posterior class probabilities are a proper n x J matrix", {
  dat <- lcm_data(2, n = 300)
  f <- lcsfm(y ~ x1 + x2, "LCM", dat$d, n_class = 2)

  expect_identical(dim(f$post.prob), c(300L, 2L))
  expect_equal(unname(rowSums(f$post.prob)), rep(1, 300), tolerance = 1e-8)
  expect_true(all(f$post.prob >= 0 & f$post.prob <= 1))
  expect_identical(dim(f$jlms_class), c(300L, 2L))
  expect_length(f$jlms, 300L)
  expect_length(f$class, 300L)
  expect_true(all(f$class %in% 1:2))
  ## jlms is the posterior-weighted average of the class-conditional scores,
  ## not the modal class's score.
  expect_equal(f$jlms, rowSums(f$post.prob * f$jlms_class), tolerance = 1e-10)
  expect_equal(sum(f$class_prob), 1, tolerance = 1e-8)
})

test_that("LCM recovers a two-class mixture it can see", {
  skip_on_cran()
  dat <- lcm_data(42, n = 800)
  f <- lcsfm(y ~ x1 + x2, "LCM", dat$d, n_class = 2)
  p <- f$out[, "par"]

  ## Label switching: identify the classes by their intercepts rather than by
  ## position, since the labelling is only identified up to permutation.
  ints <- c(p[["(Intercept)_class1"]], p[["(Intercept)_class2"]])
  lo <- which.min(ints)
  hi <- which.max(ints)
  g <- function(nm, k) p[[paste0(nm, "_class", k)]]

  expect_equal(g("(Intercept)", lo), 1.0, tolerance = 0.25)
  expect_equal(g("x1", lo), 0.5, tolerance = 0.15)
  expect_equal(g("x2", lo), 0.5, tolerance = 0.15)
  expect_equal(g("sigu", lo), 0.5, tolerance = 0.25)

  expect_equal(g("(Intercept)", hi), 3.0, tolerance = 0.30)
  expect_equal(g("x1", hi), 1.0, tolerance = 0.20)
  expect_equal(g("x2", hi), 0.2, tolerance = 0.20)
  expect_equal(g("sigu", hi), 1.0, tolerance = 0.35)

  ## Classification, up to the same permutation.
  agree <- max(mean(f$class == dat$cls), mean(f$class != dat$cls))
  expect_gt(agree, 0.80)
})

test_that("LCM_Z lets covariates drive class membership", {
  skip_on_cran()
  set.seed(7)
  n <- 700
  x1 <- rnorm(n); x2 <- rnorm(n); z <- rnorm(n)
  cls <- ifelse(runif(n) < plogis(0.5 + 1.2 * z), 1L, 2L)
  b <- rbind(c(1, 0.5, 0.5), c(3, 1.0, 0.2))
  sv <- c(0.2, 0.2); su <- c(0.5, 1.0)
  y <- b[cls, 1] + b[cls, 2] * x1 + b[cls, 3] * x2 +
    rnorm(n, 0, sv[cls]) - abs(rnorm(n, 0, su[cls]))
  f <- lcsfm(y ~ x1 + x2 | z, "LCM_Z", data.frame(y, x1, x2, z), n_class = 2)

  expect_identical(f$model_name, "LCM_Z")
  expect_true("logit_z_class1" %in% rownames(f$out))
  ## The z coefficient must have the sign that raises the probability of the
  ## class it is generated to favour. Which label that class carries is not
  ## fixed, so read the sign off the intercepts.
  lo <- which.min(c(f$out[["(Intercept)_class1", "par"]],
                    f$out[["(Intercept)_class2", "par"]]))
  zc <- f$out[["logit_z_class1", "par"]]
  expect_gt(if (lo == 1L) zc else -zc, 0.4)
  expect_gt(max(mean(f$class == cls), mean(f$class != cls)), 0.80)
})

test_that("three classes fit and stay well-formed", {
  skip_on_cran()
  set.seed(11)
  n <- 900
  x1 <- rnorm(n); x2 <- rnorm(n)
  cls <- sample(1:3, n, TRUE, prob = c(0.4, 0.35, 0.25))
  b <- rbind(c(0.5, 0.5, 0.5), c(2.5, 1.0, 0.2), c(5.0, 0.2, 1.0))
  y <- b[cls, 1] + b[cls, 2] * x1 + b[cls, 3] * x2 +
    rnorm(n, 0, 0.2) - abs(rnorm(n, 0, c(0.4, 0.8, 0.6)[cls]))
  f <- lcsfm(y ~ x1 + x2, "LCM", data.frame(y, x1, x2), n_class = 3)

  expect_identical(f$n_class, 3L)
  expect_identical(nrow(f$out), 3L * (2L + 3L) + 2L)
  expect_identical(dim(f$post.prob), c(as.integer(n), 3L))
  expect_equal(unname(rowSums(f$post.prob)), rep(1, n), tolerance = 1e-8)
  ## The three fitted intercepts must be distinct -- a collapsed mixture would
  ## put two of them on top of each other.
  ints <- sort(f$out[grepl("^\\(Intercept\\)_class", rownames(f$out)), "par"])
  expect_gt(min(diff(ints)), 0.5)
})

test_that("n_class is validated, and ignored where it means nothing", {
  dat <- lcm_data(3, n = 120)
  expect_error(lcsfm(y ~ x1 + x2, "LCM", dat$d, n_class = 1),
    "must be a single integer >= 2")
  expect_error(lcsfm(y ~ x1 + x2, "LCM", dat$d, n_class = 2.5),
    "must be a single integer >= 2")
  expect_error(lcsfm(y ~ x1 + x2, "LCM", dat$d, n_class = c(2, 3)),
    "must be a single integer >= 2")
})

test_that("the two mixture entry points stay separate", {
  dat <- lcm_data(4, n = 120)
  ## zsfm() is the zero-inefficiency model and does not answer to a latent
  ## class name; lcsfm() does not answer to ZISF. Keeping the boundary sharp
  ## is the point of having two functions rather than one with a mode switch.
  expect_error(zsfm(y ~ x1 + x2, "LCM", dat$d), "not a recognized choice")
  expect_error(lcsfm(y ~ x1 + x2, "ZISF", dat$d), "not a recognized choice")
  ## n_class belongs to lcsfm() alone -- ZISF's two components are not free.
  expect_false("n_class" %in% names(formals(zsfm)))
  expect_true("n_class" %in% names(formals(lcsfm)))
})

test_that("the row-wise log-sum-exp helper stays on the log scale", {
  M <- matrix(c(-1, -2, -3, -4), nrow = 2)
  expect_equal(sfa:::.log_row_sum_exp(M), log(rowSums(exp(M))), tolerance = 1e-12)
  ## Values far below the underflow threshold: exp() of any of these is 0 in
  ## double precision, so the naive log(rowSums(exp(.))) gives -Inf, while the
  ## shifted computation gives the right answer.
  big <- matrix(c(-800, -801, -900, -902), nrow = 2, byrow = TRUE)
  expect_true(all(is.finite(sfa:::.log_row_sum_exp(big))))
  expect_equal(sfa:::.log_row_sum_exp(big)[1], -800 + log1p(exp(-1)),
    tolerance = 1e-10)
  expect_equal(sfa:::.log_row_sum_exp(big)[2], -900 + log1p(exp(-2)),
    tolerance = 1e-10)
  expect_true(all(is.infinite(log(rowSums(exp(big))))))
  ## An all -Inf row is a genuine zero density and must stay -Inf, not NaN.
  expect_identical(sfa:::.log_row_sum_exp(matrix(-Inf, 1, 3))[1], -Inf)
})
