## Hafner, Manner and Simar (2018), Econometric Reviews 37(4):380-400.
## The "wrong skewness" problem: a model that tolerates it.
## See notes/code_history/esfm.md.
##
## The classical frontier implies a NEGATIVELY skewed composed error. When the
## sample skewness comes out positive -- which happens often at small n or when
## sigma_v/sigma_u is large, and is a small-sample accident rather than
## evidence against the model -- the MLE collapses to OLS: sigma_u is exactly
## zero, every firm is fully efficient, and the fit says nothing. That is the
## Type I failure skewness_test() reports.
##
## Hafner, Manner and Simar keep the model but widen the distribution of u.
## One parameter gamma does two jobs: its SIZE is the scale of inefficiency,
## and its SIGN chooses the skewness direction.
##
##   gamma > 0  u has the classical density (half-normal or exponential),
##              scale |gamma| -- exactly the usual model
##   gamma < 0  take that density, mirror it about zero, shift right by
##              B = a0|gamma| and truncate to [0, B]
##   gamma = 0  u is degenerate at zero and w is just noise
##
## The truncation point a0 is FIXED, not estimated: it is chosen so that
## E[u] = k1|gamma| whichever sign gamma takes. So the mean of u depends only
## on |gamma|, which is what keeps this a one-parameter family and avoids the
## identification trouble that bounded-inefficiency models run into.
##
## The payoff is that the classical model is NESTED (gamma > 0), so a wrongly
## skewed sample no longer forces a degenerate answer -- it just estimates a
## negative gamma. And because gamma = 0 is an interior point rather than a
## boundary, the likelihood ratio test of "no inefficiency" is an ordinary
## chi-square(1) rather than the chi-bar-square mixture inefficiency_test()
## needs for the classical model.

## a0 solves E[u | h-] = E[u | h+]; both are quoted to ten digits in the paper
## for the half-normal, and derived here for the exponential from
## a0 = 2(1 - exp(-a0)).
.HMS <- list(
  halfnormal = list(a0 = 1.38920329287, k1 = sqrt(2 / pi)),
  exponential = list(a0 = 1.59362426004, k1 = 1)
)

## log(exp(x) - 1) for x >= 0, without forming exp(x).
.log_expm1 <- function(x) log(expm1(pmax(x, .Machine$double.eps)))

## log[Phi(hi) - Phi(lo)] for hi > lo, stable when both are far in a tail.
.log_phi_diff <- function(lo, hi) {
  llo <- stats::pnorm(lo, log.p = TRUE)
  lhi <- stats::pnorm(hi, log.p = TRUE)
  llo + .log_expm1(lhi - llo)
}

## Per-observation log density of w = v - u, for any real gamma.
.esfm_logdens <- function(w, gamma, sigma_v, dist) {
  if (!is.finite(gamma) || !is.finite(sigma_v) || sigma_v <= 0) {
    return(rep(NA_real_, length(w)))
  }
  a0 <- .HMS[[dist]]$a0
  if (gamma == 0) {
    return(stats::dnorm(w, 0, sigma_v, log = TRUE))
  }
  g <- abs(gamma)
  if (identical(dist, "halfnormal")) {
    s <- sqrt(g^2 + sigma_v^2)
    if (gamma > 0) {
      log(2) - log(s) + stats::dnorm(w / s, log = TRUE) +
        stats::pnorm(-(w / s) * (g / sigma_v), log.p = TRUE)
    } else {
      ## gamma < 0: a normal centred at B = a0|gamma|, truncated to [0, B].
      z <- (w - a0 * gamma) / s
      A <- z * (gamma / sigma_v)
      -log(s) - log(stats::pnorm(a0) - stats::pnorm(0)) +
        stats::dnorm(z, log = TRUE) + .log_phi_diff(A, A + a0 * s / sigma_v)
    }
  } else {
    if (gamma > 0) {
      ## The exponential tilt and log Phi both diverge as gamma -> 0 and cancel;
      ## .log_phi_tilt() does that cancellation in closed form, as sfm()'s "NE"
      ## branch does.
      -log(g) - w^2 / (2 * sigma_v^2) + .log_phi_tilt(-w / sigma_v - sigma_v / g)
    } else {
      B <- a0 * g
      z2 <- w / sigma_v - sigma_v / g
      z1 <- z2 + B / sigma_v
      ## Written through .log_phi_tilt(z2) so the exp(-w/g) tilt never forms.
      -log(g) - a0 - log1p(-exp(-a0)) - w^2 / (2 * sigma_v^2) +
        .log_phi_tilt(z2) +
        .log_expm1(stats::pnorm(z1, log.p = TRUE) - stats::pnorm(z2, log.p = TRUE))
    }
  }
}

#' Stochastic frontier tolerating the wrong skewness
#'
#' @param formula,data Model and data, as for [sfm()].
#' @param dist `"halfnormal"` or `"exponential"` for the base inefficiency
#'   density.
#' @param inefdec `TRUE` for a production frontier, `FALSE` for a cost frontier.
#' @param start Optional numeric start, `c(gamma, sigma_v, beta)`.
#' @return An object of class `"esfm"` (also `"sfareg"`).
#' @export
esfm <- function(formula, data, dist = c("halfnormal", "exponential"),
                 inefdec = TRUE, start = NULL) {
  dist <- match.arg(dist)
  cl <- match.call()
  mf <- stats::model.frame(formula, data = as.data.frame(data))
  Y <- as.numeric(stats::model.response(mf))
  X <- stats::model.matrix(formula, data = as.data.frame(data))
  n <- length(Y)
  k <- ncol(X)
  sgn <- if (isTRUE(inefdec)) 1 else -1
  a0 <- .HMS[[dist]]$a0
  k1 <- .HMS[[dist]]$k1

  ols <- stats::lm.fit(X, Y)
  b0 <- ols$coefficients
  e0 <- sgn * as.numeric(ols$residuals)
  m2 <- mean(e0^2)
  m3 <- mean(e0^3)
  ## Third central moment of w is -a3 gamma^3 for either sign, so its SIGN
  ## picks the branch and its size the scale. This is the whole point: a
  ## wrongly skewed sample starts on the gamma < 0 side instead of collapsing.
  a3p <- if (dist == "halfnormal") sqrt(2 / pi) * (4 - pi) / pi else 2
  a3m <- if (dist == "halfnormal") 0.0167414748 else NA_real_
  g_start <- if (m3 < 0) {
    (-m3 / a3p)^(1 / 3)
  } else if (is.finite(a3m) && a3m > 0) {
    -(m3 / a3m)^(1 / 3)
  } else {
    -(m3 / a3p)^(1 / 3)
  }
  sv_start <- sqrt(max(m2 - (if (dist == "halfnormal") (1 - 2 / pi) else 1) *
    g_start^2, 0.1 * m2))

  nll <- function(p) {
    gam <- p[1]
    sv <- p[2]
    if (!is.finite(gam) || !is.finite(sv) || sv <= 0) {
      return(1e12)
    }
    w <- sgn * as.numeric(Y - X %*% p[3:(k + 2)])
    ll <- .esfm_logdens(w, gam, sv, dist)
    if (anyNA(ll) || any(!is.finite(ll))) {
      return(1e12)
    }
    -sum(ll)
  }

  ## Optimise from BOTH sides of zero. The likelihood is continuous at gamma=0
  ## but not differentiable there, so an optimiser started on one side will not
  ## cross to the other; the two basins have to be visited on purpose.
  starts <- list(
    c(abs(g_start), sv_start, b0),
    c(-abs(g_start), sv_start, b0),
    c(0.5 * abs(g_start), sqrt(m2), b0),
    c(-0.5 * abs(g_start), sqrt(m2), b0)
  )
  if (!is.null(start)) starts <- c(list(as.numeric(start)), starts)
  best <- NULL
  for (s0 in starts) {
    if (any(!is.finite(s0))) next
    op <- tryCatch(stats::optim(s0, nll,
      method = "Nelder-Mead",
      control = list(maxit = 4000, reltol = 1e-12)
    ), error = function(e) NULL)
    if (is.null(op)) next
    op <- tryCatch(stats::optim(op$par, nll,
      method = "BFGS",
      control = list(maxit = 800, reltol = 1e-14)
    ), error = function(e) op)
    if (is.null(best) || op$value < best$value) best <- op
  }
  if (is.null(best)) {
    stop("esfm(): every optimiser start failed.", call. = FALSE)
  }

  par <- best$par
  names(par) <- c("gamma", "sigv", colnames(X))
  gam <- par[["gamma"]]
  sv <- par[["sigv"]]
  H <- tryCatch(numDeriv::hessian(nll, par), error = function(e) NULL)
  se <- rep(NA_real_, length(par))
  if (!is.null(H)) {
    V <- tryCatch(solve(H), error = function(e) NULL)
    if (!is.null(V)) se <- sqrt(pmax(diag(V), 0))
  }
  out <- cbind(par = par, st_err = se, "t-val" = par / se)

  ## The restricted (gamma = 0) fit is OLS with the ML variance.
  ll0 <- sum(stats::dnorm(e0, 0, sqrt(mean(e0^2)), log = TRUE))
  LR <- 2 * (-best$value - ll0)

  res <- list(
    out = out, opt = best, model_name = paste0("EXT-", toupper(dist)),
    formula = formula, dist = dist, inefdec = inefdec, nobs = n,
    coefficients = par, std.errors = se, t.values = par / se, call = cl,
    gamma = gam, sigma_v = sv, wrong_skew_sample = m3 >= 0,
    residual_m3 = m3, ols_residuals = sgn * as.numeric(ols$residuals),
    logLik = -best$value, logLik_restricted = ll0,
    lr_symmetry = LR,
    p_symmetry = stats::pchisq(LR, df = 1, lower.tail = FALSE),
    Efficiency = NULL
  )
  res$Efficiency <- .esfm_te(res, Y, X)
  class(res) <- c("esfm", "sfareg")
  res
}

## E[exp(-u) | w]. For gamma > 0 this is Battese-Coelli; for gamma < 0 the
## posterior is a normal truncated to [0, B], and the expectation is taken by
## numerical integration as the paper prescribes.
.esfm_te <- function(fit, Y, X) {
  gam <- fit$gamma
  sv <- fit$sigma_v
  dist <- fit$dist
  a0 <- .HMS[[dist]]$a0
  sgn <- if (isTRUE(fit$inefdec)) 1 else -1
  b <- fit$out[colnames(X), "par"]
  w <- sgn * as.numeric(Y - X %*% b)
  if (gam == 0) {
    return(rep(1, length(w)))
  }
  if (gam > 0 && identical(dist, "halfnormal")) {
    s2 <- gam^2 + sv^2
    mu_s <- -w * gam^2 / s2
    sig_s <- sqrt(gam^2 * sv^2 / s2)
    return(.te_battese_coelli(mu_s, sig_s))
  }
  if (gam < 0 && identical(dist, "halfnormal")) {
    ## Almanidis et al. (2014), Table 1 row 1, with mu = B and scale |gamma|.
    s <- sqrt(gam^2 + sv^2)
    B <- a0 * abs(gam)
    sig_s <- abs(gam) * sv / s
    mu_s <- -(a0 * gam * sv^2 + w * gam^2) / s^2
    return(vapply(seq_along(w), function(i) {
      lo <- (0 - mu_s[i]) / sig_s
      hi <- (B - mu_s[i]) / sig_s
      den <- stats::pnorm(hi) - stats::pnorm(lo)
      if (!is.finite(den) || den <= 0) {
        return(NA_real_)
      }
      num <- tryCatch(stats::integrate(function(u) {
        exp(-u) * stats::dnorm((u - mu_s[i]) / sig_s) / sig_s
      }, 0, B, rel.tol = 1e-10)$value, error = function(e) NA_real_)
      num / den
    }, numeric(1)))
  }
  ## Exponential: integrate the posterior directly on its own support.
  g <- abs(gam)
  B <- a0 * g
  up <- if (gam > 0) Inf else B
  vapply(seq_along(w), function(i) {
    ker <- function(u) exp(-u / g) * stats::dnorm((w[i] + u) / sv)
    num <- tryCatch(stats::integrate(function(u) exp(-u) * ker(u), 0, up,
      rel.tol = 1e-9
    )$value, error = function(e) NA_real_)
    den <- tryCatch(stats::integrate(ker, 0, up, rel.tol = 1e-9)$value,
      error = function(e) NA_real_
    )
    if (!is.finite(den) || den <= 0) NA_real_ else num / den
  }, numeric(1))
}

#' Likelihood ratio test that the composed error is symmetric
#'
#' @param object An `"esfm"` fit.
#' @return A one-row data frame.
#' @export
symmetry_test <- function(object) {
  if (!inherits(object, "esfm")) {
    stop("`object` must come from esfm().", call. = FALSE)
  }
  data.frame(
    test = "LR, H0: gamma = 0 (no inefficiency)",
    statistic = object$lr_symmetry, df = 1L,
    null = "chi2(1)", p.value = object$p_symmetry,
    stringsAsFactors = FALSE
  )
}

#' @export
print.esfm <- function(x, ...) {
  cat("Extended stochastic frontier (Hafner, Manner and Simar 2018)\n")
  cat("  base density: ", x$dist, "   n = ", x$nobs, "\n", sep = "")
  cat("  gamma = ", format(x$gamma, digits = 5),
    if (x$gamma < 0) "  (negative: the sample's skewness is the 'wrong' way)" else "",
    "\n",
    sep = ""
  )
  print(round(x$out, 4))
  cat("log-likelihood: ", format(x$logLik, digits = 7),
    "   LR(gamma = 0) = ", format(x$lr_symmetry, digits = 5),
    ", p = ", format.pval(x$p_symmetry, digits = 4), "\n",
    sep = ""
  )
  invisible(x)
}
