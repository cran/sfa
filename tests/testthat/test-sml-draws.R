## Draws for simulated maximum likelihood, .sml_draws() and sfm()'s options.
##
## The properties pinned here are the ones Train (2002, ch. 9) requires and the
## ones that were previously established inline in sfm.R and would be easy to
## lose in a refactor:
##   - each unit gets its own CONTIGUOUS block, not a strided subsequence;
##   - leading elements are discarded;
##   - draws never touch 0 or 1;
##   - antithetics are exact mirrors on the uniform scale;
##   - and the default reproduces the pre-existing construction EXACTLY, so no
##     fitted result moves.

test_that("the default reproduces the previous construction byte for byte", {
  for (cfg in list(c(500, 100), c(377, 213), c(120, 64))) {
    n <- cfg[1]; Ns <- cfg[2]
    hseq <- randtoolbox::halton(n * Ns + 1000, 1, start = 1, normal = FALSE)[-c(1:1000)]
    old <- matrix(pmin(pmax(hseq, 1e-6), 1 - 1e-6), nrow = n, ncol = Ns, byrow = TRUE)
    new <- sfa:::.sml_draws(n, Ns, dim = 1L, sim_type = "halton",
                            antithetics = FALSE, burn = 1000, clamp = 1e-6)[[1]]
    expect_identical(old, new, info = paste(n, Ns))
  }
})

test_that("every unit gets a contiguous block, not a strided subsequence", {
  ## The failure this guards against: a column-major fill gives unit 1 the
  ## elements h_1, h_{1+n}, h_{1+2n}, ... of a van der Corput sequence, which
  ## for n = 500, n_draws = 100 spans only [0.50, 0.75]. A contiguous block
  ## spans essentially the whole interval.
  M <- sfa:::.sml_draws(500, 100, dim = 1L, sim_type = "halton")[[1]]
  expect_equal(dim(M), c(500L, 100L))
  spans <- apply(M, 1, function(r) diff(range(r)))
  expect_gt(min(spans), 0.9)
  ## And rows must differ from one another -- sharing one block across units is
  ## exactly the defect this constructor exists to prevent.
  expect_false(isTRUE(all.equal(M[1, ], M[2, ])))
  expect_false(isTRUE(all.equal(M[1, ], M[nrow(M), ])))
})

test_that("burn-in is applied and draws stay strictly inside (0, 1)", {
  no_burn <- sfa:::.sml_draws(10, 20, sim_type = "halton", burn = 0)[[1]]
  burned  <- sfa:::.sml_draws(10, 20, sim_type = "halton", burn = 1000)[[1]]
  expect_false(isTRUE(all.equal(no_burn, burned)))
  for (M in list(no_burn, burned)) {
    expect_true(all(M > 0 & M < 1))
    expect_true(all(is.finite(qnorm(M))))
  }
})

test_that("antithetics are exact mirrors and preserve the requested count", {
  M <- sfa:::.sml_draws(8, 20, sim_type = "halton", antithetics = TRUE)[[1]]
  expect_equal(dim(M), c(8L, 20L))
  ## First half and second half must sum to 1 elementwise.
  expect_equal(M[, 1:10] + M[, 11:20], matrix(1, 8, 10), tolerance = 1e-12)
  ## A symmetric transform then gives exact sign mirrors, which is the point:
  ## qnorm(u) = -qnorm(1 - u).
  Z <- qnorm(M)
  expect_equal(Z[, 1:10], -Z[, 11:20], tolerance = 1e-10)
  ## An odd count must still return exactly what was asked for.
  Modd <- sfa:::.sml_draws(4, 7, sim_type = "halton", antithetics = TRUE)[[1]]
  expect_equal(ncol(Modd), 7L)
})

test_that("each sequence type runs and gives distinct, valid draws", {
  got <- lapply(c("halton", "sobol", "torus", "uniform"), function(st)
    sfa:::.sml_draws(30, 40, sim_type = st, seed = 7)[[1]])
  names(got) <- c("halton", "sobol", "torus", "uniform")
  for (nm in names(got)) {
    expect_equal(dim(got[[nm]]), c(30L, 40L), info = nm)
    expect_true(all(got[[nm]] > 0 & got[[nm]] < 1), info = nm)
  }
  ## They must actually be different sequences.
  expect_false(isTRUE(all.equal(got$halton, got$sobol)))
  expect_false(isTRUE(all.equal(got$halton, got$torus)))
  expect_false(isTRUE(all.equal(got$halton, got$uniform)))
})

test_that("multiple dimensions come back as separate, uncorrelated-enough blocks", {
  ## Portability contract for psfm(), which integrates in 2 dimensions: the
  ## constructor must return one matrix per dimension, each blocked by unit.
  D <- sfa:::.sml_draws(50, 60, dim = 2L, sim_type = "halton")
  expect_length(D, 2L)
  for (M in D) expect_equal(dim(M), c(50L, 60L))
  expect_false(isTRUE(all.equal(D[[1]], D[[2]])))
  ## Primes 2 and 3 with a burn-in leave the dimensions near-uncorrelated
  ## already, which is why the old brute-force decorrelation was unnecessary.
  expect_lt(abs(cor(as.vector(D[[1]]), as.vector(D[[2]]))), 0.05)
})

test_that("a seed makes the draws reproducible and restores the RNG stream", {
  set.seed(123); before <- runif(1)
  set.seed(123)
  a <- sfa:::.sml_draws(10, 10, sim_type = "uniform", seed = 42)[[1]]
  after <- runif(1)
  b <- sfa:::.sml_draws(10, 10, sim_type = "uniform", seed = 42)[[1]]
  expect_equal(a, b)
  ## The caller's stream must be where it would have been.
  expect_equal(after, before)
})

test_that("sfm exposes the options and its default is unchanged", {
  skip_on_cran()
  d <- cs_small(N = 250)
  base <- suppressWarnings(sfm(y_pcs_ln ~ x1 + x2, model_name = "NLN", data = d))
  same <- suppressWarnings(sfm(y_pcs_ln ~ x1 + x2, model_name = "NLN", data = d,
                               sim_type = "halton", antithetics = FALSE))
  expect_equal(coef(base), coef(same), tolerance = 1e-10)

  ## A different scheme must actually change the fit, or the argument is inert.
  alt <- suppressWarnings(sfm(y_pcs_ln ~ x1 + x2, model_name = "NLN", data = d,
                              sim_type = "sobol", sim_scrambling = 1L, sim_seed = 5))
  expect_false(isTRUE(all.equal(unname(coef(base)), unname(coef(alt)))))
  expect_true(all(is.finite(coef(alt))))
})

test_that("bad simulation arguments are rejected", {
  expect_error(sfa:::.sml_draws(0, 10), "at least 1")
  ## A NULL or NA size must say so, not fall through to a comparison error.
  expect_error(sfa:::.sml_draws(10, 10, dim = NULL), "single positive")
  expect_error(sfa:::.sml_draws(NA, 10), "single positive")
  expect_error(sfa:::.sml_draws(10, 10, sim_type = "sausage"), "should be one of")
  skip_on_cran()
  d <- cs_small(N = 120)
  expect_error(sfm(y_pcs_ln ~ x1 + x2, model_name = "NLN", data = d, antithetics = "yes"),
               "TRUE or FALSE")
  expect_error(sfm(y_pcs_ln ~ x1 + x2, model_name = "NLN", data = d, sim_burn = -1),
               "non-negative")
})

## A trustworthy reference for the composed density. integrate(f, 0, Inf) is
## NOT one: the normal kernel is a spike of width sigma_v at u = -eps, and the
## default quadrature misses it entirely for small sigma_v -- during
## development it was wrong at 5 of 21 grid points, by up to 46 log units, and
## briefly made a correct simulator look catastrophically broken. Splitting the
## range at the spike and factoring out the peak of the log-integrand fixes it.
ref_composed <- function(e, sv, ldens, half = 18) {
  u0 <- -e
  lo <- max(0, u0 - half * sv)
  hi <- max(u0, 0) + half * sv
  lg <- function(u) ldens(u) + dnorm((e + u) / sv, log = TRUE) - log(sv)
  lv <- lg(seq(lo, hi, length.out = 20001))
  lv <- lv[is.finite(lv)]
  if (!length(lv)) {
    return(-Inf)
  }
  M <- max(lv)
  M + log(integrate(function(u) exp(lg(u) - M), lo, hi,
    rel.tol = 1e-12, subdivisions = 4000L, stop.on.error = FALSE
  )$value)
}

test_that("the reference itself is sound where naive quadrature is not", {
  ## Guards the guard. At sigma_v = 0.05 the two disagree by tens of log units
  ## and the reference is the one that matches a dense trapezoid sum.
  ld <- function(u) dweibull(u, 1.5, 1, log = TRUE)
  trap <- function(e, sv) {
    u <- seq(max(0, -e - 18 * sv), max(-e, 0) + 18 * sv, length.out = 2e5)
    lg <- ld(u) + dnorm((e + u) / sv, log = TRUE) - log(sv)
    M <- max(lg[is.finite(lg)])
    w <- exp(lg - M)
    w[!is.finite(w)] <- 0
    M + log(diff(range(u)) / (length(u) - 1) * (sum(w) - 0.5 * (w[1] + w[length(u)])))
  }
  naive <- function(e, sv) {
    log(integrate(function(u) dweibull(u, 1.5, 1) * dnorm((e + u) / sv) / sv,
      0, Inf, rel.tol = 1e-10, subdivisions = 2000L
    )$value)
  }
  for (e in c(-4.7, -3.5)) {
    expect_equal(ref_composed(e, 0.05, ld), trap(e, 0.05), tolerance = 1e-5)
    expect_gt(abs(ref_composed(e, 0.05, ld) - naive(e, 0.05)), 40)
  }
  ## Where naive quadrature is trustworthy the two must of course agree.
  for (e in c(-2, -0.4, 0.5)) {
    expect_equal(ref_composed(e, 0.3, ld), naive(e, 0.3), tolerance = 1e-6)
  }
})

test_that("the density/quantile helpers select the right distribution", {
  ## The failure this guards against is a swap: NLN silently getting the
  ## Weibull, or the shape and scale arguments transposed. Both would still
  ## produce a plausible-looking fit.
  u <- c(1e-8, 0.1, 0.5, 1, 3, 10, 200)
  p <- c(1e-6, 0.01, 0.25, 0.5, 0.75, 0.99, 1 - 1e-6)
  expect_equal(sfa:::.nsml_ldens("NW", 0.8, 1.5)(u), dweibull(u, shape = 1.5, scale = 0.8, log = TRUE))
  expect_equal(sfa:::.nsml_qdens("NW", 0.8, 1.5)(p), qweibull(p, shape = 1.5, scale = 0.8))
  expect_equal(sfa:::.nsml_ldens("NLN", 0.8, -0.5)(u), dlnorm(u, meanlog = -0.5, sdlog = 0.8, log = TRUE))
  expect_equal(sfa:::.nsml_qdens("NLN", 0.8, -0.5)(p), qlnorm(p, meanlog = -0.5, sdlog = 0.8))
  ## Third argument is the Weibull SHAPE and the lognormal MEANLOG, which are
  ## not interchangeable: transposing them must change the answer.
  expect_false(isTRUE(all.equal(
    sfa:::.nsml_ldens("NW", 0.8, 1.5)(u), dweibull(u, shape = 0.8, scale = 1.5, log = TRUE)
  )))
})

test_that(".sml_mis stays accurate whichever density is the narrow one", {
  ## Drawing only the noise fails where the inefficiency density is the narrow
  ## one; drawing only the inefficiency fails where the noise is. The cells
  ## below straddle both regimes.
  set.seed(11)
  n <- 200L
  eps <- rnorm(n, 0, 0.3) - (-log(1 - runif(n)))^(1 / 1.5)
  Fi <- sfa:::.sml_draws(n, 200L, 1L, sim_type = "halton", burn = 1000, clamp = 1e-6)[[1]]

  cells <- list(
    truth    = c(sv = 0.3, su = 1.0, k = 1.5),
    narrow_u = c(sv = 0.3, su = 0.3, k = 1.5), ## inefficiency is the spike
    peaked   = c(sv = 0.3, su = 1.0, k = 4.0), ## and again, via the shape
    tiny_v   = c(sv = 0.05, su = 1.0, k = 1.5), ## noise is the spike
    wide_v   = c(sv = 1.0, su = 1.0, k = 1.5)
  )
  for (nm in names(cells)) {
    p <- cells[[nm]]
    ld <- sfa:::.nsml_ldens("NW", p[["su"]], p[["k"]])
    exact <- vapply(eps, ref_composed, numeric(1), sv = p[["sv"]], ldens = ld)
    got <- sfa:::.sml_mis(eps, p[["sv"]], Fi, ld,
      sfa:::.nsml_qdens("NW", p[["su"]], p[["k"]])
    )
    expect_true(got$ok, info = nm)
    ## Total error over the sample, which is what the optimizer actually sees.
    ## Measured range at this configuration is 0.12 to 7.05.
    expect_lt(abs(sum(got$ldens - exact)), 10, label = paste("total error,", nm))
    ## And no single observation may dominate. Measured range 0.008 to 2.37.
    expect_lt(max(abs(got$ldens - exact)), 4, label = paste("worst obs,", nm))
  }
})

test_that("the two-proposal estimator is never badly wrong, unlike either alone", {
  ## The property is robustness, NOT universal dominance: in a cell that suits
  ## one proposal, that proposal alone is a little better than the hedge. What
  ## must never happen is the large error a single proposal shows off its own
  ## ground -- 38 to 46 log-likelihood units in the cells below.
  set.seed(11)
  n <- 200L
  eps <- rnorm(n, 0, 0.3) - (-log(1 - runif(n)))^(1 / 1.5)
  Fi <- sfa:::.sml_draws(n, 200L, 1L, sim_type = "halton", burn = 1000, clamp = 1e-6)[[1]]
  lmean <- function(lw) {
    m <- max(lw)
    if (!is.finite(m)) return(-Inf)
    m + log(mean(exp(lw - m)))
  }
  only_noise <- function(e, sv, ld, Fi) {
    ph <- pnorm(e / sv, lower.tail = FALSE)
    log(ph) + lmean(ld(sv * qnorm(ph * (1 - Fi), lower.tail = FALSE) - e))
  }
  only_ineff <- function(e, sv, qd, Fi) {
    lmean(dnorm((e + qd(Fi)) / sv, log = TRUE) - log(sv))
  }
  cells <- list(
    narrow_u = c(sv = 0.3, su = 0.3, k = 1.5),
    peaked   = c(sv = 0.3, su = 1.0, k = 4.0),
    tiny_v   = c(sv = 0.05, su = 1.0, k = 1.5),
    truth    = c(sv = 0.3, su = 1.0, k = 1.5)
  )
  for (nm in names(cells)) {
    p <- cells[[nm]]
    ld <- sfa:::.nsml_ldens("NW", p[["su"]], p[["k"]])
    qd <- sfa:::.nsml_qdens("NW", p[["su"]], p[["k"]])
    exact <- vapply(eps, ref_composed, numeric(1), sv = p[["sv"]], ldens = ld)
    e_mis <- abs(sum(sfa:::.sml_mis(eps, p[["sv"]], Fi, ld, qd)$ldens - exact))
    e_n <- abs(sum(vapply(seq_len(n), function(j) only_noise(eps[j], p[["sv"]], ld, Fi[j, ]), numeric(1)) - exact))
    e_i <- abs(sum(vapply(seq_len(n), function(j) only_ineff(eps[j], p[["sv"]], qd, Fi[j, ]), numeric(1)) - exact))
    ## Never more than a fraction of a log-likelihood unit behind whichever
    ## proposal suits the cell. Measured worst gap 0.204, at the truth.
    expect_lt(e_mis, min(e_n, e_i) + 0.5, label = paste("hedging cost,", nm))
    ## And where a single proposal is badly wrong, well clear of it.
    worse <- max(e_n, e_i)
    if (worse > 5) {
      expect_lt(e_mis, worse / 3, label = paste("vs the failing proposal,", nm))
    }
  }
  ## The cells must actually contain a failure, or this tests nothing: each of
  ## the two proposals has to be the badly wrong one somewhere. Note it is
  ## "peaked" and not "narrow_u" that breaks the noise draw worst -- at
  ## sigma_u = 0.3 BOTH proposals struggle and the inefficiency draw is the
  ## worse of the two, which is exactly why a selector between them fails.
  fails <- vapply(c("peaked", "tiny_v"), function(nm) {
    p <- cells[[nm]]
    ld <- sfa:::.nsml_ldens("NW", p[["su"]], p[["k"]])
    qd <- sfa:::.nsml_qdens("NW", p[["su"]], p[["k"]])
    exact <- vapply(eps, ref_composed, numeric(1), sv = p[["sv"]], ldens = ld)
    e_n <- abs(sum(vapply(seq_len(n), function(j) only_noise(eps[j], p[["sv"]], ld, Fi[j, ]), numeric(1)) - exact))
    e_i <- abs(sum(vapply(seq_len(n), function(j) only_ineff(eps[j], p[["sv"]], qd, Fi[j, ]), numeric(1)) - exact))
    which.max(c(e_n, e_i))
  }, numeric(1))
  expect_equal(unname(fails), c(1, 2))
})

test_that("a non-finite draw is survivable, not fatal", {
  ## The likelihood steers the optimizer away from such a region, but the
  ## predictor still has to be formable, so weights must always come back.
  set.seed(5)
  n <- 40L
  eps <- rnorm(n, 0, 0.3) - rweibull(n, 1.5, 1)
  Fi <- sfa:::.sml_draws(n, 40L, 1L, sim_type = "halton", burn = 1000, clamp = 1e-6)[[1]]
  bad <- sfa:::.sml_mis(eps, 0.3, Fi,
    sfa:::.nsml_ldens("NW", 1, 1.5),
    function(p) { q <- qweibull(p, 1.5, 1); q[1] <- Inf; q }
  )
  expect_false(bad$ok)
  expect_false(is.null(bad$w))
  expect_true(all(is.finite(sfa:::.sml_mis_mean(bad, bad$u)[-1])))
  ## A clean call reports ok.
  good <- sfa:::.sml_mis(eps, 0.3, Fi, sfa:::.nsml_ldens("NW", 1, 1.5),
    sfa:::.nsml_qdens("NW", 1, 1.5)
  )
  expect_true(good$ok)
})

test_that("the posterior mean uses the same weights as the density", {
  ## If predictor and likelihood ever drift apart, u_hat stops being the Bayes
  ## rule under the fitted model. E[1|eps] = 1 is the cheap invariant.
  set.seed(3)
  n <- 120L
  eps <- rnorm(n, 0, 0.3) - rweibull(n, 1.5, 1)
  Fi <- sfa:::.sml_draws(n, 200L, 1L, sim_type = "halton", burn = 1000, clamp = 1e-6)[[1]]
  mis <- sfa:::.sml_mis(eps, 0.3, Fi, sfa:::.nsml_ldens("NW", 1, 1.5),
    sfa:::.nsml_qdens("NW", 1, 1.5)
  )
  expect_equal(sfa:::.sml_mis_mean(mis, matrix(1, n, ncol(mis$u))), rep(1, n), tolerance = 1e-12)
  ## E[u|eps] must be positive and E[exp(-u)|eps] must be a valid efficiency.
  eu <- sfa:::.sml_mis_mean(mis, mis$u)
  ee <- sfa:::.sml_mis_mean(mis, exp(-mis$u))
  expect_true(all(eu > 0))
  expect_true(all(ee > 0 & ee < 1))
  ## Jensen: E[exp(-u)] >= exp(-E[u]).
  expect_true(all(ee >= exp(-eu) - 1e-10))
})

test_that("NW and NLN share one draw rule, and NW's is far below its old one", {
  ## Both use the same two-proposal estimator, so both take
  ## max(200, ceiling(3*sqrt(n))), half the count to each proposal.
  ## NW's rule before the change of variable was max(400, ceiling(8*sqrt(n))).
  skip_on_cran()
  d <- cs_small(N = 400)
  fw <- suppressWarnings(sfm(y_pcs_wb ~ x1 + x2, model_name = "NW", data = d))
  fl <- suppressWarnings(sfm(y_pcs_ln ~ x1 + x2, model_name = "NLN", data = d))
  expect_true(all(is.finite(coef(fw))))
  expect_true(all(is.finite(coef(fl))))
  ## Efficiency predictions stay in range under the substitution.
  expect_true(all(fw$exp_u_hat >= 0 & fw$exp_u_hat <= 1))
  expect_true(all(fw$u_hat >= 0))
})

test_that("a seed randomizes the deterministic sequences by shifting", {
  ## halton() and torus() ignore set.seed() -- they are fixed lattices -- so
  ## before the shift was added, passing `seed` with the default sim_type did
  ## nothing at all while appearing to control the draws. That also left the
  ## cross-sectional simulated-ML models with UNrandomized draws, which the
  ## panel path (.gtre_halton_draws) had not been since 1.2.0.
  a <- .sml_draws(50, 200, 1, sim_type = "halton", seed = 1)[[1]]
  b <- .sml_draws(50, 200, 1, sim_type = "halton", seed = 999)[[1]]
  expect_false(isTRUE(all.equal(a, b)))

  ## Same seed, same draws.
  expect_equal(a, .sml_draws(50, 200, 1, sim_type = "halton", seed = 1)[[1]])

  ## The default is untouched: no seed, no shift.
  expect_equal(
    .sml_draws(50, 200, 1, sim_type = "halton", seed = NULL)[[1]],
    .sml_draws(50, 200, 1, sim_type = "halton", seed = NULL)[[1]]
  )

  ## A shift preserves equidistribution, which is the reason for shifting
  ## rather than permuting: the lattice survives.
  expect_true(all(b > 0 & b < 1))
  expect_equal(mean(b), 0.5, tolerance = 0.01)
  expect_equal(sd(b), 1 / sqrt(12), tolerance = 0.01)
})

test_that("a seeded shift does not leak into the caller's RNG stream", {
  ## The `seed` block takes one snapshot and registers one on.exit. A SECOND
  ## snapshot/restore pair around the shift would run after it and restore the
  ## post-set.seed state instead of the caller's -- which a first attempt at
  ## the shift did, silently reseeding the caller.
  set.seed(1234)
  before <- runif(3)
  set.seed(1234)
  invisible(.sml_draws(20, 50, 1, sim_type = "halton", seed = 77))
  expect_equal(runif(3), before)

  set.seed(4321)
  before2 <- runif(3)
  set.seed(4321)
  invisible(.sml_draws(20, 50, 2, sim_type = "torus", seed = 5))
  expect_equal(runif(3), before2)
})
