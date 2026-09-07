## meanefficiency(): the model-implied mean efficiency, and mean efficiency
## among the most efficient p of the distribution. Gap G6.
## See notes/code_history/meanefficiency.md.
##
## This is a DIFFERENT quantity from the average of the per-observation
## conditional predictions that print()/summary() report. That average is a
## property of this sample; E[exp(-U)] is a property of the fitted model, which
## is what a reader comparing two studies wants.
##
## Every closed form below is checked against numerical integration of the same
## density in test-meanefficiency.R. That is the point of keeping the numeric
## engine even where a closed form exists: the two must agree, and a wrong
## algebraic simplification cannot survive the comparison.

## Distribution of u implied by a fit, as (density, quantile, closed form).
## Parameterizations are the ones the package REPORTS, which are not always the
## ones the generator uses -- `NHN` and `NTN` report lambda/sigma, so sigma_u
## has to be reconstructed as lambda*sigma/sqrt(1+lambda^2).
.u_spec <- function(object) {
  p <- object$out[, "par"]
  nm <- names(p)
  g <- function(k) if (k %in% nm) unname(p[[k]]) else NA_real_

  su_from_lamsig <- function() {
    lam <- g("lambda"); sg <- g("sigma")
    if (!is.finite(lam) || !is.finite(sg)) return(NA_real_)
    lam * sg / sqrt(1 + lam^2)
  }

  m <- object$model_name
  switch(m,
    NHN = , NHN_Z = {
      su <- if ("sigu" %in% nm) g("sigu") else su_from_lamsig()
      list(dist = "half-normal", par = c(sigma_u = su),
           pdf = function(u) 2 * stats::dnorm(u, 0, su),
           qf  = function(pr) stats::qnorm((1 + pr) / 2, 0, su),
           closed = function() 2 * exp(su^2 / 2) * stats::pnorm(-su))
    },
    NE = , NE_Z = {
      su <- g("sigu")
      list(dist = "exponential", par = c(sigma_u = su),
           pdf = function(u) stats::dexp(u, rate = 1 / su),
           qf  = function(pr) stats::qexp(pr, rate = 1 / su),
           closed = function() 1 / (1 + su))
    },
    NTN = {
      su <- if ("sigu" %in% nm) g("sigu") else su_from_lamsig()
      mu <- g("mu")
      Z <- stats::pnorm(mu / su)
      list(dist = "truncated normal", par = c(sigma_u = su, mu = mu),
           pdf = function(u) stats::dnorm(u, mu, su) / Z,
           qf  = function(pr) stats::qnorm((1 - Z) + pr * Z, mu, su),
           closed = function() {
             exp(-mu + su^2 / 2) * stats::pnorm(mu / su - su) / Z
           })
    },
    NR = {
      su <- g("sigu")
      list(dist = "Rayleigh", par = c(sigma_u = su),
           pdf = function(u) (u / su^2) * exp(-u^2 / (2 * su^2)),
           qf  = function(pr) su * sqrt(-2 * log(1 - pr)),
           ## erfc(x) = 2 * pnorm(-x * sqrt(2)), so erfc(su/sqrt(2)) = 2*pnorm(-su).
           closed = function() {
             1 - su * sqrt(pi / 2) * exp(su^2 / 2) * 2 * stats::pnorm(-su)
           })
    },
    NG = {
      sh <- g("mu"); sc <- g("sigu")   # sfm reports mu as SHAPE, sigu as SCALE
      list(dist = "gamma", par = c(shape = sh, scale = sc),
           pdf = function(u) stats::dgamma(u, shape = sh, scale = sc),
           qf  = function(pr) stats::qgamma(pr, shape = sh, scale = sc),
           closed = function() (1 + sc)^(-sh))
    },
    NU = {
      th <- g("theta")
      list(dist = "uniform", par = c(theta = th),
           pdf = function(u) ifelse(u >= 0 & u <= th, 1 / th, 0),
           qf  = function(pr) pr * th,
           closed = function() (1 - exp(-th)) / th)
    },
    NGE = {
      su <- g("sigu")
      list(dist = "generalized exponential", par = c(sigma_u = su),
           pdf = function(u) (2 / su) * (1 - exp(-u / su)) * exp(-u / su),
           qf  = function(pr) -su * log(1 - sqrt(pr)),
           ## F(u) = (1 - e^{-u/s})^2, so the density is a DIFFERENCE of two
           ## exponentials and the transform is 2/((s+1)(s+2)).
           closed = function() 2 / ((su + 1) * (su + 2)))
    },
    NLN = {
      sd <- g("sigu"); ml <- g("mu")
      list(dist = "lognormal", par = c(meanlog = ml, sdlog = sd),
           pdf = function(u) stats::dlnorm(u, meanlog = ml, sdlog = sd),
           qf  = function(pr) stats::qlnorm(pr, meanlog = ml, sdlog = sd),
           closed = NULL)   # no closed form
    },
    NW = {
      sc <- g("sigu"); k <- g("k")
      list(dist = "Weibull", par = c(shape = k, scale = sc),
           pdf = function(u) stats::dweibull(u, shape = k, scale = sc),
           qf  = function(pr) stats::qweibull(pr, shape = k, scale = sc),
           closed = NULL)   # no closed form
    },
    NNAK = {
      mm <- g("mu"); om <- g("sigu")^2
      ## u = sqrt(G) with G ~ Gamma(m, scale = Omega/m), so f(u) = 2u * g(u^2).
      list(dist = "Nakagami", par = c(m = mm, Omega = om),
           pdf = function(u) 2 * u * stats::dgamma(u^2, shape = mm, scale = om / mm),
           qf  = function(pr) sqrt(stats::qgamma(pr, shape = mm, scale = om / mm)),
           closed = NULL)   # parabolic cylinder function; integrated instead
    },
    NULL
  )
}

meanefficiency <- function(object,
                           p = c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99),
                           use_closed_form = TRUE) {
  if (!inherits(object, "sfareg")) {
    stop("`object` must be an \"sfareg\" fit.", call. = FALSE)
  }
  spec <- .u_spec(object)
  if (is.null(spec)) {
    stop("meanefficiency(): no distribution of u is implemented for ",
      "model_name \"",
      if (is.null(object$model_name)) "unknown" else object$model_name,
      "\". Supported: NHN, NHN_Z, NE, NE_Z, NTN, NR, NG, NU, NGE, NLN, NW, NNAK.",
      call. = FALSE
    )
  }
  if (any(!is.finite(spec$par))) {
    stop("meanefficiency(): the fitted parameters of the ", spec$dist,
      " distribution are not all finite, so its mean efficiency is undefined.",
      call. = FALSE
    )
  }
  if (!is.numeric(p) || any(!is.finite(p)) || any(p <= 0) || any(p > 1)) {
    stop("`p` must be numbers in (0, 1].", call. = FALSE)
  }

  integrand <- function(u) exp(-u) * spec$pdf(u)
  num_mean <- function() {
    stats::integrate(integrand, 0, Inf, rel.tol = 1e-10,
      stop.on.error = FALSE
    )$value
  }

  closed_ok <- isTRUE(use_closed_form) && !is.null(spec$closed)
  te <- if (closed_ok) spec$closed() else num_mean()
  method <- if (closed_ok) "closed form" else "numerical integration"

  ## Supra-percentile: mean efficiency among the most efficient p of the
  ## distribution -- that is, conditioning on the SMALLEST p of u, since
  ## efficiency is decreasing in u. Always integrated: only the unconditional
  ## mean has a closed form worth having.
  supra <- vapply(p, function(pr) {
    if (pr >= 1) return(te)
    hi <- spec$qf(pr)
    if (!is.finite(hi) || hi <= 0) return(NA_real_)
    v <- stats::integrate(integrand, 0, hi, rel.tol = 1e-10,
      stop.on.error = FALSE
    )$value
    v / pr
  }, numeric(1))

  list(
    model = object$model_name,
    distribution = spec$dist,
    parameters = spec$par,
    mean_efficiency = te,
    method = method,
    supra = data.frame(p = p, mean_efficiency = supra)
  )
}
