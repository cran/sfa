## Copula stochastic frontier: copsfm(), gap K5.

cop_gen <- function(seed, n, rho, su = 1, sv = 0.4, b = c(0.5, 0.8, -0.4)) {
  set.seed(seed)
  x1 <- rnorm(n); x2 <- rnorm(n)
  ## Gaussian copula between the marginals: draw correlated normals, push each
  ## through its own marginal quantile function. w1 = Phi(z1) is F_V(v) and
  ## w2 = Phi(z2) is F_U(u), so (w1, w2) carry exactly the dependence rho.
  z1 <- rnorm(n); z2 <- rho * z1 + sqrt(1 - rho^2) * rnorm(n)
  v <- sv * z1
  u <- su * qnorm((1 + pnorm(z2)) / 2)
  data.frame(y = b[1] + b[2] * x1 + b[3] * x2 + v - u, x1 = x1, x2 = x2)
}

test_that("each copula density integrates to 1 and is 1 at independence", {
  ## The check that catches a mistranscribed density. Both are asserted for
  ## every family, which is why families whose density could not be verified
  ## are deliberately absent from .cop_logc().
  gl <- .gauss_legendre_01(80L)
  G <- expand.grid(a = gl$nodes, b = gl$nodes)
  W <- as.vector(outer(gl$weights, gl$weights))
  for (fam in c("gaussian", "fgm")) {
    for (p in c(-0.6, -0.2, 0, 0.3, 0.7)) {
      I <- sum(W * exp(.cop_logc(G$a, G$b, p, fam)))
      expect_equal(I, 1, tolerance = 1e-4, info = paste(fam, p))
    }
    ## At the independence parameter the density is exactly 1, so log c is 0.
    expect_equal(max(abs(.cop_logc(G$a, G$b, 0, fam))), 0, info = fam)
  }
})

test_that("the composed density reduces to normal/half-normal at independence", {
  ## With c = 1 the quadrature must reproduce the closed form the rest of the
  ## package uses. This checks the composition and the change of variable
  ## u = t/(1-t) together.
  gl <- .gauss_legendre_01(128L)
  tt <- gl$nodes; wt <- gl$weights
  uu <- tt / (1 - tt); jac <- 1 / (1 - tt)^2
  su <- 1; sv <- 0.4; e <- seq(-4, 2, length.out = 9)

  U <- matrix(uu, length(e), length(uu), byrow = TRUE)
  V <- e + U
  lg <- dnorm(V, 0, sv, log = TRUE) + log(2) + dnorm(U, 0, su, log = TRUE)
  lg <- sweep(lg, 2, log(wt) + log(jac), "+")
  got <- .log_row_sum_exp(lg)

  sig <- sqrt(su^2 + sv^2); lam <- su / sv
  want <- log(2) - log(sig) + dnorm(e / sig, log = TRUE) +
    pnorm(-e * lam / sig, log.p = TRUE)
  expect_equal(got, want, tolerance = 1e-8)
})

test_that("copsfm returns a well-formed sfareg", {
  skip_on_cran()
  d <- cop_gen(4, 600, 0.5)
  f <- copsfm(y ~ x1 + x2, data = d, copula = "gaussian")

  expect_s3_class(f, "sfareg")
  expect_identical(f$model_name, "COP")
  expect_identical(f$copula, "gaussian")
  expect_equal(ncol(f$out), 3L)
  expect_identical(rownames(f$out),
    c("(Intercept)", "x1", "x2", "sigma_u", "sigma_v", "rho"))
  expect_gt(f$out["sigma_u", "par"], 0)
  expect_lt(abs(f$out["rho", "par"]), 1)
  expect_equal(length(f$jlms), nrow(d))
  expect_equal(f$efficiency, exp(-f$jlms))
  expect_equal(nobs(f), nrow(d))
  ## fgm names its parameter differently, and reports it.
  g <- copsfm(y ~ x1 + x2, data = d, copula = "fgm")
  expect_true("theta" %in% rownames(g$out))
})

test_that("the frontier slopes are recovered even where rho is not", {
  skip_on_cran()
  ## The dependence parameter needs a large sample (see ?copsfm); the SLOPES do
  ## not, and separating the two is the point of this test.
  d <- cop_gen(4, 2000, 0.6)
  f <- copsfm(y ~ x1 + x2, data = d, copula = "gaussian")
  expect_equal(unname(f$out["x1", "par"]), 0.8, tolerance = 0.08)
  expect_equal(unname(f$out["x2", "par"]), -0.4, tolerance = 0.08)
})

test_that("copsfm rejects malformed calls", {
  d <- cop_gen(2, 300, 0.3)
  expect_error(copsfm(~ x1 + x2, data = d), "two-sided formula")
  expect_error(copsfm(y ~ x1 | x2, data = d), "must not contain a `\\|` segment")
  expect_error(copsfm(y ~ x1 + x2, data = as.matrix(d)), "must be a data.frame")
  expect_error(copsfm(y ~ x1 + x2, data = d, inefdec = NA), "must be TRUE or FALSE")
  expect_error(copsfm(y ~ x1 + x2, data = d, n_nodes = 4), ">= 16")
  ## "clayton" used to belong here, as a family deliberately not implemented.
  ## It is implemented now, so the guard moves to one that still is not: the
  ## point is that an unrecognised family is refused, not that this particular
  ## one is missing.
  expect_error(copsfm(y ~ x1 + x2, data = d, copula = "bb1"))
  expect_error(copsfm(y ~ x1 + x2, data = d, copula = "tawn"))
})

## ---------------------------------------------------------------------------
## Frank, Clayton, Gumbel, Joe and the rotations (added 2026-09-04)
## ---------------------------------------------------------------------------
##
## The families the header of .cop_logc() used to list as deliberately absent.
## They are admitted now because each is checked against something OTHER than
## itself: its own CDF. A density that integrates to 1 can still be the wrong
## density; a density that equals the second mixed partial of the right CDF at
## forty points cannot be.

test_that("each new density IS the second mixed partial of its own CDF", {
  CDF <- list(
    frank = function(u, v, th)
      -1 / th * log1p(expm1(-th * u) * expm1(-th * v) / expm1(-th)),
    clayton = function(u, v, th) (u^(-th) + v^(-th) - 1)^(-1 / th),
    gumbel = function(u, v, th) exp(-(((-log(u))^th + (-log(v))^th)^(1 / th))),
    joe = function(u, v, th)
      1 - ((1 - u)^th + (1 - v)^th - (1 - u)^th * (1 - v)^th)^(1 / th)
  )
  num <- function(f, u, v, th, h = 1e-5) {
    (f(u + h, v + h, th) - f(u + h, v - h, th) -
       f(u - h, v + h, th) + f(u - h, v - h, th)) / (4 * h^2)
  }
  pars <- list(frank = c(-8, -2, 2, 8), clayton = c(0.3, 1, 3, 8),
    gumbel = c(1.3, 2, 4), joe = c(1.3, 2, 5))
  g <- expand.grid(u = seq(0.1, 0.9, by = 0.2), v = seq(0.1, 0.9, by = 0.2))
  for (fam in names(CDF)) {
    for (th in pars[[fam]]) {
      a <- exp(.cop_logc(g$u, g$v, th, fam))
      b <- mapply(function(u, v) num(CDF[[fam]], u, v, th), g$u, g$v)
      expect_equal(a, b, tolerance = 1e-4, info = paste(fam, th))
    }
  }
})

test_that("the new densities integrate to 1 and are 1 at independence", {
  gl <- .gauss_legendre_01(120L)
  G <- expand.grid(a = gl$nodes, b = gl$nodes)
  W <- as.vector(outer(gl$weights, gl$weights))
  pars <- list(frank = c(-8, -2, 2, 8), clayton = c(0.3, 1), gumbel = c(1.3, 2),
    joe = c(1.3, 2))
  for (fam in names(pars)) {
    for (p in pars[[fam]]) {
      I <- sum(W * exp(.cop_logc(G$a, G$b, p, fam)))
      ## Looser than the Gaussian/FGM tolerance ON PURPOSE: these densities
      ## concentrate in a corner, and what is left is quadrature error rather
      ## than density error -- which is why the mixed-partial test above, and
      ## not this one, is the check that would catch a wrong formula.
      expect_equal(I, 1, tolerance = 5e-3, info = paste(fam, p))
    }
  }
  ind <- c(frank = 0, clayton = 0, gumbel = 1, joe = 1)
  for (fam in names(ind)) {
    expect_equal(max(abs(.cop_logc(G$a, G$b, ind[[fam]], fam))), 0, info = fam)
  }
})

test_that("rotations reverse the sign of dependence, which is why they exist", {
  ## Clayton, Gumbel and Joe carry only POSITIVE dependence. Nothing rules out
  ## a negative association between noise and inefficiency, so without the
  ## rotations those families could only ever report the independence boundary
  ## against negatively dependent data.
  set.seed(1)
  n <- 1e5
  u <- runif(n); v <- runif(n)
  rho_s <- function(fam, th) {
    w <- exp(.cop_logc_rot(u, v, th, fam))
    w[!is.finite(w)] <- 0
    12 * sum(w * (u - 0.5) * (v - 0.5)) / sum(w)
  }
  for (base in c("clayton", "gumbel", "joe")) {
    th <- if (base == "clayton") 3 else 2.5
    pos <- rho_s(base, th)
    expect_gt(pos, 0.3)
    expect_lt(rho_s(paste0(base, "90"), th), -0.3)
    expect_lt(rho_s(paste0(base, "270"), th), -0.3)
    ## The survival rotation preserves the sign rather than flipping it.
    expect_gt(rho_s(paste0(base, "180"), th), 0.3)
  }
})

test_that("every advertised family actually fits", {
  skip_on_cran()
  d <- cop_gen(seed = 4, n = 250, rho = 0.3)
  fams <- eval(formals(copsfm)$copula)
  for (fam in fams) {
    f <- suppressWarnings(copsfm(y ~ x1 + x2, data = d, copula = fam,
      n_nodes = 32))
    expect_s3_class(f, "sfareg")
    expect_true(is.finite(as.numeric(logLik(f))), info = fam)
    sp <- sfa:::.cop_spec(fam)
    expect_gte(f$copula_par, sp$lo)
    expect_lte(f$copula_par, sp$hi)
  }
})

test_that("families that do not recover their parameter warn, and the others do not", {
  skip_on_cran()
  d <- cop_gen(seed = 5, n = 150, rho = 0.3)
  fit <- function(fam) copsfm(y ~ x1 + x2, data = d, copula = fam, n_nodes = 32)

  ## Verified to recover: silent. Asserted with expect_silent-style intent but
  ## via the message text, because copsfm() can emit unrelated optimiser
  ## warnings at small n_nodes and those are not what this test is about.
  for (fam in c("gaussian", "fgm", "frank")) {
    w <- tryCatch({ suppressMessages(fit(fam)); NA_character_ },
      warning = function(x) conditionMessage(x))
    expect_false(isTRUE(grepl("did not recover|has not been measured", w)), info = fam)
  }

  ## Measured to collapse: the warning names the measured rate, so the user is
  ## told how bad it is rather than merely that something is wrong.
  expect_warning(fit("gumbel"), "did not recover.*36%")
  expect_warning(fit("joe"), "did not recover.*40%")
  expect_warning(fit("clayton270"), "did not recover.*56%")
  expect_warning(fit("gumbel90"), "did not recover.*60%")

  ## Not measured: says so, rather than implying either result.
  expect_warning(fit("joe180"), "has not been measured")
  expect_warning(fit("clayton90"), "has not been measured")

  ## A missing name must give NA, not an error -- `[[` on a named vector throws.
  expect_true(is.na(unname(sfa:::.COP_COLLAPSE["gumbel270"])))
  expect_silent(sfa:::.cop_warn_family("frank"))
})
