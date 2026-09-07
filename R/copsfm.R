## Copula stochastic frontier: dependence between the noise and the
## inefficiency. Gap K5. See notes/code_history/copsfm.md.
##
## Every other model in this package assumes v and u are independent. That
## assumption is ubiquitous in the frontier literature and, as Smith (2008)
## argued, mostly an artefact of how easily the likelihood is then written
## down: if a farmer misjudges a seasonal shock and mis-allocates inputs in
## response, the shock and the inefficiency are the same event seen twice.
##
## The package's own taxonomy already states the general form (JSS paper, Eq 6):
##
##   f_{V,U}(v, u) = f_V(v) f_U(u) c(F_V(v), F_U(u); rho)
##
## with every existing specification setting c = 1. This file relaxes that, and
## nothing else: the marginals are the ordinary normal and half-normal.
##
## Composing, with eps = v - S*u so that v = eps + S*u,
##
##   f_eps(eps) = int_0^inf f_V(eps + S u) f_U(u) c(F_V(eps + S u), F_U(u)) du
##
## which is a one-dimensional integral per observation, done by Gauss-Legendre
## on a transformed range rather than by simulation.

## Copula densities on the unit square, in LOGS.
##
## Only families whose density can be written down exactly and checked are
## included. Each is verified in the tests two ways: it integrates to 1 over
## [0,1]^2, and at the independence parameter it returns exactly 1.
.cop_logc <- function(w1, w2, par, family) {
  switch(family,
    independent = rep(0, length(w1)),
    gaussian = {
      ## rho in (-1, 1). a = Phi^-1(w1), b = Phi^-1(w2):
      ##   c = (1-rho^2)^-1/2 exp{ [2 rho a b - rho^2 (a^2 + b^2)] / [2(1-rho^2)] }
      rho <- par
      if (abs(rho) >= 1) return(rep(NA_real_, length(w1)))
      a <- stats::qnorm(w1)
      b <- stats::qnorm(w2)
      om <- 1 - rho^2
      -0.5 * log(om) + (2 * rho * a * b - rho^2 * (a^2 + b^2)) / (2 * om)
    },
    fgm = {
      ## Farlie-Gumbel-Morgenstern: c = 1 + theta (1-2 w1)(1-2 w2), |theta| <= 1.
      ## Bounded dependence (Spearman rho is theta/3), so it cannot represent
      ## strong association -- which is exactly why it is a useful contrast to
      ## the Gaussian rather than a redundant second option.
      th <- par
      if (abs(th) > 1) return(rep(NA_real_, length(w1)))
      v <- 1 + th * (1 - 2 * w1) * (1 - 2 * w2)
      ifelse(v > 0, log(v), NA_real_)
    },
    ## Frank: the only other one-parameter family that spans BOTH signs of
    ## dependence over the full range. theta and (1 - e^-theta) always share a
    ## sign, so their product is positive and can be logged directly; the
    ## independence limit theta -> 0 is taken explicitly because the ratio is
    ## 0/0 there.
    frank = {
      th <- par
      if (!is.finite(th)) return(rep(NA_real_, length(w1)))
      if (abs(th) < 1e-8) return(rep(0, length(w1)))
      em <- -expm1(-th)                      # 1 - exp(-th), signed like th
      D <- em - (-expm1(-th * w1)) * (-expm1(-th * w2))
      ok <- is.finite(D) & D != 0
      out <- rep(NA_real_, length(w1))
      out[ok] <- log(th * em) - th * (w1[ok] + w2[ok]) - 2 * log(abs(D[ok]))
      out
    },
    ## Clayton: lower-tail dependence, theta > 0, independence as theta -> 0.
    clayton = {
      th <- par
      if (!is.finite(th) || th <= 0) {
        return(if (isTRUE(all.equal(th, 0))) rep(0, length(w1)) else rep(NA_real_, length(w1)))
      }
      z <- w1^(-th) + w2^(-th) - 1
      ok <- is.finite(z) & z > 0
      out <- rep(NA_real_, length(w1))
      out[ok] <- log1p(th) - (1 + th) * (log(w1[ok]) + log(w2[ok])) -
        (2 + 1 / th) * log(z[ok])
      out
    },
    ## Gumbel: upper-tail dependence, theta >= 1, independence at theta = 1.
    gumbel = {
      th <- par
      if (!is.finite(th) || th < 1) return(rep(NA_real_, length(w1)))
      if (abs(th - 1) < 1e-12) return(rep(0, length(w1)))
      x <- -log(w1); y <- -log(w2)
      A <- (x^th + y^th)^(1 / th)
      ok <- is.finite(A) & A > 0 & x > 0 & y > 0
      out <- rep(NA_real_, length(w1))
      out[ok] <- -A[ok] + (th - 1) * (log(x[ok]) + log(y[ok])) +
        (1 - 2 * th) * log(A[ok]) + log(A[ok] + th - 1) -
        log(w1[ok]) - log(w2[ok])
      out
    },
    ## Joe: upper-tail dependence, heavier than Gumbel; theta >= 1.
    joe = {
      th <- par
      if (!is.finite(th) || th < 1) return(rep(NA_real_, length(w1)))
      if (abs(th - 1) < 1e-12) return(rep(0, length(w1)))
      a <- (1 - w1)^th; b <- (1 - w2)^th
      z <- a + b - a * b
      ok <- is.finite(z) & z > 0
      out <- rep(NA_real_, length(w1))
      out[ok] <- (1 / th - 2) * log(z[ok]) +
        (th - 1) * (log1p(-w1[ok]) + log1p(-w2[ok])) +
        log(th - 1 + z[ok])
      out
    },
    NA_real_
  )
}

## Rotations. Clayton, Gumbel and Joe carry only POSITIVE dependence, which is
## a real restriction here: nothing rules out a negative association between
## noise and inefficiency, and a family that cannot express one will report the
## independence boundary instead of a negative estimate. The 90 and 270 degree
## rotations are exact reflections of the density, so they need no new algebra:
##   c_90(u,v)  = c(1-u, v)      c_180(u,v) = c(1-u, 1-v)      c_270 = c(u, 1-v)
.COP_ROT <- c(gaussian = NA, fgm = NA, frank = NA, independent = NA,
  clayton = 0, gumbel = 0, joe = 0,
  clayton90 = 90, clayton180 = 180, clayton270 = 270,
  gumbel90 = 90, gumbel180 = 180, gumbel270 = 270,
  joe90 = 90, joe180 = 180, joe270 = 270)

.cop_base <- function(family) sub("(90|180|270)$", "", family)

.cop_logc_rot <- function(w1, w2, par, family) {
  rot <- .COP_ROT[[family]]
  base <- .cop_base(family)
  if (is.null(rot) || is.na(rot) || rot == 0) return(.cop_logc(w1, w2, par, base))
  if (rot == 90) return(.cop_logc(1 - w1, w2, par, base))
  if (rot == 180) return(.cop_logc(1 - w1, 1 - w2, par, base))
  .cop_logc(w1, 1 - w2, par, base)
}

## Independence value of each family's parameter, and its admissible range.
.cop_spec <- function(family) {
  ## Upper limits are where the density stops being computable in double
  ## precision, not where the family stops being defined: Clayton at theta = 28
  ## already has a Kendall tau of 0.93, and past that u^-theta overflows.
  ## Independence sits ON the lower bound for Clayton/Gumbel/Joe, so par0 is set
  ## just inside it -- starting the optimiser exactly at a bound is how a search
  ## reports the boundary back as an estimate.
  base <- .cop_base(family)
  switch(base,
    gaussian = list(par0 = 0, lo = -0.95, hi = 0.95, name = "rho"),
    fgm      = list(par0 = 0, lo = -0.99, hi = 0.99, name = "theta"),
    frank    = list(par0 = 1e-4, lo = -35, hi = 35, name = "theta"),
    clayton  = list(par0 = 0.05, lo = 1e-6, hi = 28, name = "theta"),
    gumbel   = list(par0 = 1.0001, lo = 1, hi = 17, name = "theta"),
    joe      = list(par0 = 1.0001, lo = 1, hi = 30, name = "theta"),
    list(par0 = 0, lo = -0.95, hi = 0.95, name = "par")
  )
}

## Which families have been shown to recover their own dependence parameter.
##
## Measured 2026-09-04: 25 samples per family generated FROM that family at
## n = 400, refitted with it, counting how often the estimate came back on the
## independence boundary. A correct density does not imply a recoverable
## parameter, and for most of these families it is not recoverable: on Gumbel
## data with Spearman 0.685 at n = 2000, every family INCLUDING the true one
## returns the boundary and their log-likelihoods differ by less than 0.04.
## The likelihood is flat in the dependence parameter. That is the model, not
## the code -- each density is verified against the second mixed partial of its
## own CDF in test-copsfm.R.
.COP_COLLAPSE <- c(gumbel = 36, joe = 40, clayton270 = 56, gumbel90 = 60)
.COP_VERIFIED <- c("gaussian", "fgm", "frank", "clayton")

.cop_warn_family <- function(family) {
  if (family %in% .COP_VERIFIED) return(invisible(NULL))
  ## `[` not `[[`: a name absent from the vector must give NA, not an error.
  rate <- unname(.COP_COLLAPSE[family])
  if (!is.na(rate)) {
    warning("copsfm(copula = \"", family, "\"): this family did not recover its ",
      "own dependence parameter in testing -- on 25 samples generated from it ",
      "at n = 400, ", rate, "% of fits returned the independence boundary. The ",
      "density is correct; the parameter is weakly identified, and the ",
      "likelihood is close to flat in it. Prefer copula = \"frank\" or ",
      "\"clayton\", which recovered with no boundary collapses, and treat any ",
      "estimate from this family as exploratory. See ?copsfm.",
      call. = FALSE)
  } else {
    warning("copsfm(copula = \"", family, "\"): this family's density is ",
      "verified but its RECOVERY has not been measured. The families that were ",
      "measured collapsed to the independence boundary on 36-60% of samples ",
      "generated from themselves, so do not assume this one behaves better. ",
      "Prefer copula = \"frank\" or \"clayton\". See ?copsfm.",
      call. = FALSE)
  }
  invisible(NULL)
}

copsfm <- function(formula,
                   data,
                   copula = c("gaussian", "fgm", "frank",
                              "clayton", "clayton90", "clayton180", "clayton270",
                              "gumbel", "gumbel90", "gumbel180", "gumbel270",
                              "joe", "joe90", "joe180", "joe270"),
                   inefdec = TRUE,
                   n_nodes = 128,
                   maxit.bobyqa = 10000,
                   maxit.psoptim = 1000,
                   maxit.optim = 1000,
                   start_val = FALSE,
                   PSopt = FALSE,
                   optHessian = TRUE,
                   Method = "L-BFGS-B",
                   verbose = FALSE,
                   rand.psoptim = NULL) {
  call <- match.call()
  copula <- match.arg(copula)
  .cop_warn_family(copula)
  Start.Time <- Sys.time()
  cz <- .SFA_CONSTANTS

  if (missing(formula) || !inherits(formula, "formula") || length(formula) != 3L) {
    stop("copsfm(): `formula` must be a two-sided formula, e.g. y ~ x1 + x2.",
      call. = FALSE
    )
  }
  if (any(grepl("|", deparse(formula), fixed = TRUE))) {
    stop("copsfm(): `formula` must not contain a `|` segment. Determinants of ",
      "the variance are not yet supported alongside a copula.",
      call. = FALSE
    )
  }
  if (missing(data) || !is.data.frame(data)) {
    stop("copsfm(): `data` must be a data.frame.", call. = FALSE)
  }
  if (length(inefdec) != 1L || !is.logical(inefdec) || is.na(inefdec)) {
    stop("copsfm(): `inefdec` must be TRUE or FALSE.", call. = FALSE)
  }
  if (length(n_nodes) != 1L || !is.numeric(n_nodes) || n_nodes < 16) {
    stop("copsfm(): `n_nodes` must be a single number >= 16.", call. = FALSE)
  }

  mf <- stats::model.frame(formula, data = data, na.action = stats::na.omit)
  y <- as.numeric(stats::model.response(mf))
  X <- stats::model.matrix(stats::terms(mf), mf)
  n <- length(y)
  k <- ncol(X)
  if (n <= k + 3L) {
    stop("copsfm(): ", n, " observations cannot identify ", k + 3L,
      " parameters.",
      call. = FALSE
    )
  }
  S <- if (isTRUE(inefdec)) 1 else -1
  cs <- .cop_spec(copula)

  ## Gauss-Legendre nodes on (0, 1), mapped to u in (0, Inf) by u = t/(1-t).
  ## The map is chosen so the integrand's mass -- which sits at moderate u --
  ## lands in the middle of the node range rather than in a tail.
  ##
  ## 128 nodes, not 64. Measured against the closed-form normal/half-normal
  ## density at independence, the maximum absolute error in log f is 1.3e-2 at
  ## 32 nodes, 2.5e-5 at 64 and 8.2e-14 at 128. 64 looks adequate and is not:
  ## an error of 1e-5 per observation is an error of 0.02 in a log-likelihood
  ## over 2000 points, which is the scale at which modes are being compared.
  gl <- .gauss_legendre_01(as.integer(n_nodes))
  tt <- gl$nodes
  wt <- gl$weights
  uu <- tt / (1 - tt)
  jac <- 1 / (1 - tt)^2

  ## log f_eps(eps_i), the composed density, by quadrature over u.
  .log_dens <- function(eps, su, sv, cpar) {
    ## Rows are observations, columns quadrature nodes.
    U <- matrix(uu, nrow = length(eps), ncol = length(uu), byrow = TRUE)
    V <- eps + S * U
    lfv <- stats::dnorm(V, 0, sv, log = TRUE)
    lfu <- log(2) + stats::dnorm(U, 0, su, log = TRUE)
    lg <- lfv + lfu
    if (!identical(copula, "independent") && !is.null(cpar)) {
      ## Marginal CDFs, clamped off the endpoints: qnorm(0) is -Inf.
      w1 <- pmin(pmax(stats::pnorm(V / sv), 1e-12), 1 - 1e-12)
      w2 <- pmin(pmax(2 * stats::pnorm(U / su) - 1, 1e-12), 1 - 1e-12)
      lc <- .cop_logc_rot(as.numeric(w1), as.numeric(w2), cpar, copula)
      if (any(!is.finite(lc))) return(rep(NA_real_, length(eps)))
      lg <- lg + matrix(lc, nrow = nrow(lg))
    }
    lg <- sweep(lg, 2, log(wt) + log(jac), "+")
    .log_row_sum_exp(lg)
  }

  ## Starting values: OLS plus Olson's moments, independence for the copula.
  ols <- stats::lm.fit(x = X, y = y)
  e0 <- as.numeric(ols$residuals)
  m2 <- mean(e0^2); m3 <- mean(e0^3)
  k3 <- sqrt(2 / pi) * (1 - 4 / pi)
  su0 <- if (S * m3 < 0) (S * m3 / k3)^(1 / 3) else 0.5 * stats::sd(e0)
  su0 <- if (is.finite(su0)) max(su0, 1e-3) else max(0.5 * stats::sd(e0), 1e-3)
  sv0 <- sqrt(max(m2 - (1 - 2 / pi) * su0^2, 1e-6))
  b0 <- as.numeric(ols$coefficients)
  if (attr(stats::terms(mf), "intercept") == 1L) {
    b0[1L] <- b0[1L] + S * su0 * sqrt(2 / pi)
  }
  start_v <- c(b0, log(su0), log(sv0), cs$par0)
  par_names <- c(colnames(X), "sigma_u", "sigma_v", cs$name)
  if (isTRUE(start_val)) names(start_v) <- par_names

  i_b <- seq_len(k); i_su <- k + 1L; i_sv <- k + 2L; i_c <- k + 3L

  like.fn <- function(th) {
    if (!all(is.finite(th))) return(cz$MAX_VALUE)
    su <- exp(pmin(th[i_su], 12)); sv <- exp(pmin(th[i_sv], 12))
    cpar <- th[i_c]
    if (cpar < cs$lo || cpar > cs$hi) return(cz$MAX_VALUE)
    eps <- S * as.numeric(y - X %*% th[i_b])
    ll <- .log_dens(eps, su, sv, cpar)
    if (any(!is.finite(ll))) return(cz$MAX_VALUE)
    ## The scaffold MINIMIZES; every likelihood here returns the negative sum.
    -sum(ll)
  }

  span <- pmax(10 * abs(b0), 10)
  lower1 <- c(b0 - span, log(1e-6), log(1e-6), cs$lo)
  upper1 <- c(b0 + span, log(1e4), log(1e4), cs$hi)

  Opt.Bobyqa <- opt.bobyqa(fn = like.fn, start_v = start_v,
    lower.bobyqa = lower1, upper.bobyqa = upper1,
    maxit.bobyqa = maxit.bobyqa, bob.TF = TRUE, rhobeg = NA, rhoend = NA,
    verbose = verbose
  )
  start_v <- Opt.Bobyqa$start_v; bob1 <- Opt.Bobyqa$bob1

  Opt.Psoptim <- opt.psoptim(fn = like.fn, start_v,
    lower.psoptim = lower1, rand.psoptim = rand.psoptim,
    upper.psoptim = upper1, maxit.psoptim, psopt.TF = PSopt,
    rand.order = FALSE, verbose = verbose
  )
  start_v <- Opt.Psoptim$start_v; opt00 <- Opt.Psoptim$opt00

  Opt.Optim <- opt.optim(fn = like.fn, start_v = start_v,
    lower.optim = lower1, upper.optim = upper1, maxit.optim = maxit.optim,
    opt.TF = optHessian, method = Method, optHessian = TRUE, verbose = verbose
  )
  start_v <- Opt.Optim$start_v; opt <- Opt.Optim$opt
  End.Time <- end.time(Start.Time)

  if (optHessian == FALSE & PSopt == FALSE) { opt <- bob1; st_err <- rep(NA, length(opt$par)) }
  if (optHessian == FALSE & PSopt == TRUE)  { opt <- opt00; st_err <- rep(NA, length(opt$par)) }
  if (optHessian == TRUE) {
    st_err <- if (isTRUE(as.numeric(sum(colMeans(opt$hessian))) == 0)) {
      rep(NA, length(opt$par))
    } else {
      suppressWarnings(sqrt(diag(solve(opt$hessian))))
    }
  }

  th <- opt$par
  su <- exp(th[i_su]); sv <- exp(th[i_sv])
  par <- c(th[i_b], su, sv, th[i_c])
  se <- c(st_err[i_b], su * st_err[i_su], sv * st_err[i_sv], st_err[i_c])

  out <- matrix(NA_real_, 3L, length(par))
  rownames(out) <- c("par", "st_err", "t-val")
  colnames(out) <- par_names
  out[1, ] <- par; out[2, ] <- se; out[3, ] <- par / se

  ## E[u | eps] by the same quadrature the likelihood used, so the predictor and
  ## the objective cannot disagree about the model.
  eps <- S * as.numeric(y - X %*% th[i_b])
  Um <- matrix(uu, nrow = n, ncol = length(uu), byrow = TRUE)
  V <- eps + S * Um
  lg <- stats::dnorm(V, 0, sv, log = TRUE) + log(2) + stats::dnorm(Um, 0, su, log = TRUE)
  w1 <- pmin(pmax(stats::pnorm(V / sv), 1e-12), 1 - 1e-12)
  w2 <- pmin(pmax(2 * stats::pnorm(Um / su) - 1, 1e-12), 1 - 1e-12)
  lg <- lg + matrix(.cop_logc_rot(as.numeric(w1), as.numeric(w2), th[i_c], copula), nrow = n)
  lg <- sweep(lg, 2, log(wt) + log(jac), "+")
  wts <- exp(lg - .log_row_sum_exp(lg))
  jlms <- as.numeric(rowSums(wts * Um))

  results <- list(
    t(out), c(opt), End.Time, start_v, "COP", formula, copula, th[i_c],
    jlms, exp(-jlms), S, n, as.integer(n_nodes),
    out["par", ], out["st_err", ], out["t-val", ], call
  )
  class(results) <- "sfareg"
  names(results) <- c(
    "out", "opt", "total_time", "start_v", "model_name", "formula", "copula",
    "copula_par", "jlms", "efficiency", "S", "nobs", "n_nodes",
    "coefficients", "std.errors", "t.values", "call"
  )
  results
}
