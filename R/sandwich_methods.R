## bread() and estfun() for "sfareg", which is all the `sandwich` package needs
## to work against these fits: vcovHC() for heteroskedasticity-consistent
## standard errors, vcovCL() for clustered ones, and coeftest() to print them.
## See notes/code_history/sandwich_methods.md.
##
## Clustered standard errors are routine in applied work on firm-level panels,
## and vcov.sfareg() returns the inverse Hessian and nothing else -- which is
## valid only if the likelihood is correctly specified and the observations are
## independent. Neither is safe to assume in this literature.

## The "bread" of the sandwich: n times the inverse of the negative Hessian of
## the summed log-likelihood, i.e. n * vcov.
bread.sfareg <- function(x, ...) {
  V <- stats::vcov(x)
  n <- stats::nobs(x)
  if (!is.finite(n)) {
    stop("bread(): the number of observations is unavailable for this fit.",
      call. = FALSE
    )
  }
  V * n
}

## The per-observation score contributions, d loglik_i / d theta, as an n x p
## matrix.
##
## There is no analytic gradient for most of these models, so the scores are
## taken by CENTRAL differences of the per-observation log-likelihood, which
## `sfm(keep_objective = TRUE)` retains. Central rather than forward because
## the error is O(h^2) instead of O(h), and these likelihoods are built from
## pnorm/dnorm calls whose own relative error is around 1e-15 -- a forward
## difference loses roughly half the available digits.
##
## The step is scaled per parameter: h_j = eps * max(|theta_j|, 1). A single
## absolute step cannot serve a vector holding both a variance on the log scale
## and a frontier coefficient in levels.
estfun.sfareg <- function(x, ...) {
  if (is.null(x$objective)) {
    stop("estfun(): this fit does not retain its likelihood, so the score ",
      "matrix cannot be built. Refit with keep_objective = TRUE (accepted by ",
      "sfm() and psfm()).",
      call. = FALSE
    )
  }
  if (!is.null(x$robust) && !identical(x$robust, "mle")) {
    stop("estfun(): the robust divergence estimators (robust = ",
      dQuote(x$robust), ") do not maximise a log-likelihood, so a score ",
      "matrix is not defined for them. sandwich-based standard errors do not ",
      "apply; sfm() already reports sandwich-form errors for those fits.",
      call. = FALSE
    )
  }
  par <- x$opt$par
  if (is.null(par)) {
    stop("estfun(): this fit has no parameter vector (was it produced by a ",
      "non-likelihood estimator?).",
      call. = FALSE
    )
  }
  nm <- names(x$coefficients)
  p <- length(par)

  ll_i <- function(theta) {
    v <- tryCatch(x$objective(theta, per_obs = TRUE), error = function(e) NULL)
    if (is.null(v)) {
      ## psfm()'s closures were never given the `per_obs` branch that sfm()'s
      ## gained in 1.2.0, so keep_objective = TRUE stores an objective the
      ## score matrix cannot use. Say which case this is rather than blaming a
      ## stale fit for both.
      stop("estfun(): the stored likelihood does not return per-observation ",
        "contributions, so the score matrix cannot be built.\n",
        if (!is.null(x$model_name)) {
          paste0("  This fit is model_name = ", dQuote(x$model_name), ". ")
        } else {
          "  "
        },
        "psfm()'s panel likelihoods do not yet support `per_obs`, so ",
        "sandwich-form standard errors, influence_sfa() and the other ",
        "score-based diagnostics are currently available for sfm() fits only. ",
        "If this IS an sfm() fit, it was retained by a version older than ",
        "1.2.0 -- refit.",
        call. = FALSE
      )
    }
    as.numeric(v)
  }

  base <- ll_i(par)
  n <- length(base)
  eps <- .Machine$double.eps^(1 / 3)
  sc <- matrix(NA_real_, n, p, dimnames = list(NULL, nm))
  for (j in seq_len(p)) {
    h <- eps * max(abs(par[j]), 1)
    up <- dn <- par
    up[j] <- par[j] + h
    dn[j] <- par[j] - h
    sc[, j] <- (ll_i(up) - ll_i(dn)) / (2 * h)
  }
  ## A non-finite score would silently poison the whole meat matrix.
  sc[!is.finite(sc)] <- 0
  sc
}
