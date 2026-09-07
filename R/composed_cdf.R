## A composed-error CDF for EVERY cross-sectional model, not just half-normal.
##
## pcomposed() computes P(eps <= q) as an expectation over u,
##
##   F(q) = E_u[ F_v( (q + s*u) / sigma_v ) ],   s = +1 production, -1 cost,
##
## which is the right idea but is hard-wired to a half-normal u. Goodness-of-fit
## testing (Wang, Amsler and Schmidt 2011) needs the same quantity for whichever
## u the model actually assumes, so this file supplies the missing piece: a
## registry of the u-density behind each `model_name`, and a quadrature that
## works for all of them. See notes/code_history/composed_cdf.md.
##
## Conditioning on u leaves the NOISE cdf in closed form, so v is never
## integrated -- the tails inherit pnorm()'s (or pt()'s) own accuracy.

## Half-normal scale from the (lambda, sigma) pair NHN and NTN report.
.lam_sig_to_u <- function(lambda, sigma) sigma * lambda / sqrt(1 + lambda^2)
.lam_sig_to_v <- function(lambda, sigma) sigma / sqrt(1 + lambda^2)

## The registry. `par` is the NAMED parameter vector from fit$out[, "par"].
##
## Every parameterization below is taken from the likelihood in sfm.R and
## cross-checked against DATA_GENERATION_REFERENCE.md, because two of them do
## not read the way the name suggests:
##   THT lists sigu BEFORE sigv, the only model that does.
##   NG's `sigu` is the gamma SCALE and its `mu` is the gamma SHAPE.
##   NNAK's `sigu` is the Nakagami spread: Omega = sigu^2, the RMS of u.
##   NLN's `mu` is a meanlog, not a mean.
.composed_u_spec <- function(model_name, par) {
  g <- function(nm) {
    if (!nm %in% names(par)) {
      stop(".composed_u_spec(): parameter ", dQuote(nm), " not found on this fit for ",
        "model_name = ", dQuote(model_name), ".",
        call. = FALSE
      )
    }
    unname(par[[nm]])
  }
  norm <- function(sv) list(type = "normal", sigma_v = sv, df = NA_real_)
  tnoise <- function(sv, df) list(type = "t", sigma_v = sv, df = df)

  switch(model_name,
    "NHN" = {
      lam <- g("lambda"); sg <- g("sigma")
      su <- .lam_sig_to_u(lam, sg)
      list(ldens = function(u) log(2) + stats::dnorm(u, 0, su, log = TRUE),
        scale = su, upper = Inf, noise = norm(.lam_sig_to_v(lam, sg)))
    },
    "NTN" = {
      lam <- g("lambda"); sg <- g("sigma"); mu <- g("mu")
      su <- .lam_sig_to_u(lam, sg)
      ## u ~ N(mu, su^2) truncated to u >= 0 (Stevenson 1980).
      lden <- stats::pnorm(mu / su, log.p = TRUE)
      list(ldens = function(u) stats::dnorm(u, mu, su, log = TRUE) - lden,
        scale = max(su, abs(mu)), upper = Inf, noise = norm(.lam_sig_to_v(lam, sg)))
    },
    "NE" = {
      su <- g("sigu")
      list(ldens = function(u) stats::dexp(u, rate = 1 / su, log = TRUE),
        scale = su, upper = Inf, noise = norm(g("sigv")))
    },
    "NR" = {
      ## Rayleigh with scale sigu/sqrt(2), so E[u^2] = sigu^2.
      b <- g("sigu") / sqrt(2)
      list(ldens = function(u) log(u) - 2 * log(b) - u^2 / (2 * b^2),
        scale = b, upper = Inf, noise = norm(g("sigv")))
    },
    "NU" = {
      th <- g("theta")
      list(ldens = function(u) ifelse(u >= 0 & u <= th, -log(th), -Inf),
        scale = th, upper = th, noise = norm(g("sigv")))
    },
    "NGE" = {
      ## Generalized exponential GE(2, 1/sigu): F(u) = (1 - exp(-u/sigu))^2.
      su <- g("sigu")
      list(ldens = function(u) log(2) - log(su) - u / su +
             log(-expm1(-u / su)), scale = su, upper = Inf,
        noise = norm(g("sigv")))
    },
    "NLN" = {
      list(ldens = function(u) stats::dlnorm(u, meanlog = g("mu"),
             sdlog = g("sigu"), log = TRUE),
        scale = exp(g("mu")), upper = Inf, noise = norm(g("sigv")))
    },
    "NW" = {
      list(ldens = function(u) stats::dweibull(u, shape = g("k"),
             scale = g("sigu"), log = TRUE),
        scale = g("sigu"), upper = Inf, noise = norm(g("sigv")))
    },
    "NG" = {
      ## `sigu` is the gamma SCALE and `mu` the SHAPE -- see the file header.
      sh <- g("mu"); sc <- g("sigu")
      list(ldens = function(u) stats::dgamma(u, shape = sh, scale = sc, log = TRUE),
        scale = sh * sc, upper = Inf, noise = norm(g("sigv")))
    },
    "NNAK" = {
      ## u = sqrt(G), G ~ Gamma(shape = m, scale = Omega/m), Omega = sigu^2.
      m <- g("mu"); Om <- g("sigu")^2
      list(ldens = function(u) log(2) + m * log(m) - lgamma(m) - m * log(Om) +
             (2 * m - 1) * log(u) - m * u^2 / Om,
        scale = g("sigu"), upper = Inf, noise = norm(g("sigv")))
    },
    "TSL" = {
      ## Truncated skew-Laplace (Wang 2012). Read off the log-likelihood's own
      ## normalising constant: f(u) = ((1+l)/(2l+1))(1/s)[2e^{-u/s} -
      ## e^{-(1+l)u/s}], which integrates to 1 because 2 - 1/(1+l) =
      ## (2l+1)/(1+l).
      s <- g("sigu"); l <- g("lambda")
      list(ldens = function(u) {
        a <- log(2) - u / s
        b <- -(1 + l) * u / s
        log1p(l) - log(2 * l + 1) - log(s) +
          a + log(-expm1(pmin(b - a, -.Machine$double.eps)))
      }, scale = s, upper = Inf, noise = norm(g("sigv")))
    },
    "THT" = {
      ## Tancredi (2002). The ONLY model that lists sigu before sigv, and the
      ## only one whose composed error is NOT an independent convolution.
      ##
      ## The composed error is skew-t, which is a SCALE MIXTURE: v and u are
      ## divided by the SAME sqrt(V/a), V ~ chi^2_a. Treating it as an
      ## independent t noise plus a half-normal -- which is what tHN actually
      ## is -- gets the density wrong by up to 1.74 in logs, so the two models
      ## must not share a branch here even though they look alike.
      ##
      ##   eps = (v0 - u0)/sqrt(W),  W = V/a ~ Gamma(a/2, rate = a/2)
      ##   F(q) = E_W[ F_NHN(q sqrt(W)) ]
      ##
      ## so the half-normal machinery below is reused, wrapped in one extra
      ## quadrature over W. `mix_df` switches that on.
      su <- g("sigu")
      list(ldens = function(u) log(2) + stats::dnorm(u, 0, su, log = TRUE),
        scale = su, upper = Inf, noise = norm(g("sigv")), mix_df = g("a"))
    },
    "tHN" = {
      su <- g("sigu")
      list(ldens = function(u) log(2) + stats::dnorm(u, 0, su, log = TRUE),
        scale = su, upper = Inf, noise = tnoise(g("sigv"), g("nu")))
    },
    stop("pcomposed_model(): no composed-error CDF for model_name = ",
      dQuote(model_name), ". Supported: NHN, NTN, NE, NR, NU, NGE, NLN, NW, ",
      "NG, NNAK, TSL, THT, tHN.",
      call. = FALSE
    )
  )
}

## log F_v(z) for the two noise families, upper or lower tail.
.noise_lcdf <- function(z, noise, lower) {
  if (identical(noise$type, "normal")) {
    stats::pnorm(z, lower.tail = lower, log.p = TRUE)
  } else {
    stats::pt(z, df = noise$df, lower.tail = lower, log.p = TRUE)
  }
}

#' Composed-error distribution function for any cross-sectional model
#'
#' @param q Numeric vector of quantiles.
#' @param model_name The model whose composed error is wanted.
#' @param par A named parameter vector, as in `fit$out[, "par"]`.
#' @param inefdec `TRUE` for a production frontier, `FALSE` for a cost frontier.
#' @param lower.tail,log.p As elsewhere in R.
#' @param n_nodes Quadrature nodes over the inefficiency term.
#' @return A numeric vector, `P(eps <= q)` (or its log / upper tail).
#' @export
pcomposed_model <- function(q, model_name, par, inefdec = TRUE,
                            lower.tail = TRUE, log.p = FALSE, n_nodes = 128) {
  if (!is.numeric(q)) stop("`q` must be numeric.", call. = FALSE)
  if (is.null(names(par))) {
    stop("pcomposed_model(): `par` must be NAMED, e.g. fit$out[, \"par\"].",
      call. = FALSE
    )
  }
  if (length(n_nodes) != 1L || !is.finite(n_nodes) || n_nodes < 16) {
    stop("pcomposed_model(): `n_nodes` must be a single number >= 16.", call. = FALSE)
  }
  sp <- .composed_u_spec(model_name, par)
  sv <- sp$noise$sigma_v
  if (!is.finite(sv) || sv <= 0) {
    stop("pcomposed_model(): the noise scale is not positive on this fit.", call. = FALSE)
  }
  s <- if (isTRUE(inefdec)) 1 else -1
  sgn <- if (isTRUE(lower.tail)) 1 else -1

  nd <- .composed_u_nodes(sp, n_nodes)
  uu <- nd$u
  lw <- nd$lw
  ltot <- nd$ltot

  base <- function(qq) {
    vapply(qq, function(qi) {
      .logsumexp(lw + .noise_lcdf(sgn * (qi + s * uu) / sv, sp$noise, TRUE)) - ltot
    }, numeric(1))
  }

  out <- if (is.null(sp$mix_df)) {
    base(q)
  } else {
    ## Scale mixture: F(q) = E_W[F_base(q sqrt(W))], W ~ Gamma(a/2, rate a/2),
    ## which has mean 1 -- so the t/(1-t) map at scale 1 sits where the mass is.
    a <- sp$mix_df
    gw <- .gauss_legendre_01(as.integer(n_nodes))
    tw <- gw$nodes
    ww <- tw / (1 - tw)
    lww <- log(gw$weights) - 2 * log1p(-tw) +
      stats::dgamma(ww, shape = a / 2, rate = a / 2, log = TRUE)
    kw <- is.finite(lww)
    lww <- lww[kw]
    ww <- ww[kw]
    lwtot <- .logsumexp(lww)
    sq <- sqrt(ww)
    vapply(q, function(qi) .logsumexp(lww + base(qi * sq)) - lwtot, numeric(1))
  }
  if (isTRUE(log.p)) out else exp(out)
}

#' Composed-error CDF evaluated from a fitted model
#'
#' @param object An `"sfareg"` fit from [sfm()].
#' @param q Quantiles. Defaults to the fit's own composed-error residuals,
#'   which requires `data`.
#' @param data The data frame the model was fitted to. `sfm()` does not retain
#'   it, so it must be supplied (or be findable from the fit's `call`) whenever
#'   `q` is not given.
#' @param ... Passed to [pcomposed_model()].
#' @return A numeric vector of probabilities.
#' @export
composed_cdf <- function(object, q = NULL, data = NULL, ...) {
  if (!inherits(object, "sfareg")) {
    stop("`object` must be an \"sfareg\" fit.", call. = FALSE)
  }
  if (is.null(q)) {
    q <- .sfa_eps_hat(object, data, parent.frame())
  }
  pcomposed_model(q, object$model_name, object$out[, "par"],
    inefdec = .sfa_inefdec(object), ...
  )
}

## The fitted composed error eps = y - x'beta, on the production orientation.
## NOT the OLS residual: the frontier coefficients are the ML ones.
##
## sfm() does not store `data` on the fit, so it has to be supplied or
## recovered from the recorded call. Recovering it is a convenience, not a
## guarantee -- if the object has moved away from where it was fitted, the
## caller passes `data` instead and gets a clear error if they do not.
.sfa_eps_hat <- function(object, data = NULL, env = parent.frame()) {
  if (is.null(data)) {
    data <- tryCatch(eval(object$call$data, envir = env), error = function(e) NULL)
  }
  if (is.null(data)) {
    stop("This needs the data the model was fitted to, and sfm() does not ",
      "retain it. Pass `data =` explicitly.",
      call. = FALSE
    )
  }
  data <- as.data.frame(data)
  fx <- stats::formula(Formula::Formula(object$formula), lhs = 1, rhs = 1)
  mf <- stats::model.frame(fx, data = data)
  y <- as.numeric(stats::model.response(mf))
  X <- stats::model.matrix(fx, data = data)
  cf <- object$out[, "par"]
  bn <- intersect(colnames(X), names(cf))
  if (length(bn) != ncol(X)) {
    stop(".sfa_eps_hat(): could not match every design column to a fitted ",
      "coefficient; missing ",
      paste(setdiff(colnames(X), names(cf)), collapse = ", "), ".",
      call. = FALSE
    )
  }
  as.numeric(y - X[, bn, drop = FALSE] %*% cf[bn])
}


## Draw n composed errors from a fitted model. Needed by the parametric
## bootstrap in gof_test(), which has to generate data from the null exactly as
## the model states it -- a bootstrap that resampled residuals instead would be
## testing something else.
##
## Returned on the PRODUCTION orientation (eps = v - u); the caller flips it.
.composed_rgen <- function(n, model_name, par) {
  g <- function(nm) unname(par[[nm]])
  ## u first, then the noise, because two models make the noise depend on it.
  u <- switch(model_name,
    "NHN" = abs(stats::rnorm(n, 0, .lam_sig_to_u(g("lambda"), g("sigma")))),
    "NTN" = {
      su <- .lam_sig_to_u(g("lambda"), g("sigma"))
      mu <- g("mu")
      ## Inverse-CDF draw from N(mu, su^2) truncated below at 0: exact, and
      ## it cannot loop forever the way rejection sampling can when mu << 0.
      lo <- stats::pnorm(0, mu, su)
      stats::qnorm(lo + stats::runif(n) * (1 - lo), mu, su)
    },
    "NE" = stats::rexp(n, rate = 1 / g("sigu")),
    "NR" = (g("sigu") / sqrt(2)) * sqrt(stats::rnorm(n)^2 + stats::rnorm(n)^2),
    "NU" = stats::runif(n, 0, g("theta")),
    "NGE" = -log1p(-sqrt(stats::runif(n))) * g("sigu"),
    "NLN" = stats::rlnorm(n, meanlog = g("mu"), sdlog = g("sigu")),
    "NW" = stats::rweibull(n, shape = g("k"), scale = g("sigu")),
    "NG" = stats::rgamma(n, shape = g("mu"), scale = g("sigu")),
    "NNAK" = sqrt(stats::rgamma(n, shape = g("mu"), scale = g("sigu")^2 / g("mu"))),
    "TSL" = {
      ## F(u) is closed form; invert it on a grid and interpolate, which is
      ## cheaper and steadier than n calls to uniroot().
      s <- g("sigu"); l <- g("lambda")
      Fu <- function(x) {
        ((1 + l) / (2 * l + 1)) *
          (2 * (1 - exp(-x / s)) - (1 - exp(-(1 + l) * x / s)) / (1 + l))
      }
      grid <- seq(0, s * 60, length.out = 20001)
      Fg <- Fu(grid)
      ## Drop the saturated tail: beyond it F is 1 to machine precision and the
      ## repeated values make approx() collapse ties and warn.
      keep <- !duplicated(Fg) & Fg < 1 - 1e-12
      keep[1] <- TRUE
      stats::approx(Fg[keep], grid[keep], xout = stats::runif(n), rule = 2)$y
    },
    "THT" = abs(stats::rnorm(n, 0, g("sigu"))),
    "tHN" = abs(stats::rnorm(n, 0, g("sigu"))),
    stop(".composed_rgen(): no sampler for model_name = ", dQuote(model_name),
      ".", call. = FALSE)
  )
  sp <- .composed_u_spec(model_name, par)
  sv <- sp$noise$sigma_v
  eps <- if (identical(sp$noise$type, "t")) {
    sv * stats::rt(n, df = sp$noise$df) - u
  } else {
    stats::rnorm(n, 0, sv) - u
  }
  if (!is.null(sp$mix_df)) {
    ## THT: v and u share ONE chi-square scale factor. Dividing the already
    ## formed difference is what makes it skew-t rather than a convolution.
    eps <- eps / sqrt(stats::rgamma(n, shape = sp$mix_df / 2, rate = sp$mix_df / 2))
  }
  eps
}


## Was this fitted as a production or a cost frontier? sfm() takes inefdec as
## TRUE/FALSE and does NOT store it on the fit, so it is recovered from the
## recorded call. psfm() uses the string "cost function" for the same idea, so
## both spellings are accepted rather than assuming one convention holds
## everywhere.
.sfa_inefdec <- function(object) {
  if (!is.null(object$inefdec)) {
    if (is.character(object$inefdec)) {
      return(!grepl("cost", object$inefdec, ignore.case = TRUE))
    }
    if (is.logical(object$inefdec) && length(object$inefdec) == 1L) {
      return(isTRUE(object$inefdec))
    }
  }
  v <- tryCatch(eval(object$call$inefdec), error = function(e) NULL)
  if (is.null(v)) {
    return(TRUE) ## sfm()'s own default
  }
  if (is.character(v)) !grepl("cost", v, ignore.case = TRUE) else isTRUE(v)
}


## Quadrature nodes and log-weights for the inefficiency term, given a spec
## from .composed_u_spec(). Factored out so the composed CDF and the
## Chen-Wang moment machinery integrate over the SAME nodes.
.composed_u_nodes <- function(sp, n_nodes = 128) {
  gl <- .gauss_legendre_01(as.integer(n_nodes))
  tt <- gl$nodes
  if (is.finite(sp$upper)) {
    ## Bounded support: map straight onto (0, b). The t/(1-t) map would put
    ## half its nodes past the edge and meet a kink at u = b, which is exactly
    ## where Gauss-Legendre stops converging quickly.
    uu <- sp$upper * tt
    lw <- log(gl$weights) + log(sp$upper) + sp$ldens(uu)
  } else {
    ## u = scale * t/(1-t) carries the scale, so the nodes follow the
    ## distribution instead of sitting at fixed absolute values.
    zz <- tt / (1 - tt)
    uu <- sp$scale * zz
    lw <- log(gl$weights) + log(sp$scale) - 2 * log1p(-tt) + sp$ldens(uu)
  }
  keep <- is.finite(lw)
  lw <- lw[keep]
  uu <- uu[keep]
  if (!length(lw)) {
    stop("the inefficiency density underflowed at every quadrature node.",
      call. = FALSE
    )
  }
  ## Renormalise on the quadrature's own total mass rather than assuming 1, so
  ## a slightly under-resolved density cannot bias the result either way.
  list(u = uu, lw = lw, ltot = .logsumexp(lw), w = exp(lw - .logsumexp(lw)))
}
