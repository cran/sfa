## pcomposed()/dcomposed(): the distribution of the composed error itself.
## Gap L5, from Amsler, Schmidt and Tsay (2019), J Prod Anal 52:29-35.
## See notes/code_history/pcomposed.md.
##
## Every likelihood in this package evaluates the composed error's DENSITY, and
## none of them needs its CDF -- which is why there has never been one. Copula
## models do need it: a Gaussian copula on the composed error contains
## qnorm(F(eps)), so an F that saturates at 0 or 1 does not merely lose
## precision, it returns -Inf or +Inf and takes the copula density with it.
## That is the failure Amsler et al. set out to avoid, and it is the reason
## this is a prerequisite for extending copsfm() rather than a convenience.
##
## There is no closed form for the skew-normal CDF except at lambda = 0 and
## lambda = 1. Both are implemented as exact shortcuts AND reachable by the
## general path, so the two can be compared against each other.

## The composed density, which is closed form: a skew-normal.
##
##   eps = v - u  (production, inefdec = TRUE)   or   v + u  (cost)
##   v ~ N(0, sigma_v^2),  u ~ N+(0, sigma_u^2),  independent
##
## With sigma^2 = sigma_u^2 + sigma_v^2 and lambda = sigma_u/sigma_v, the
## density of the COST case is (2/sigma) phi(e/sigma) Phi(lambda e/sigma); the
## production case is its mirror image, which is why only the sign of `e`
## changes below rather than the formula.
dcomposed <- function(x, sigma_u, sigma_v, inefdec = TRUE, log = FALSE) {
  .check_scales(sigma_u, sigma_v, "dcomposed")
  s <- sqrt(sigma_u^2 + sigma_v^2)
  lam <- sigma_u / sigma_v
  e <- if (isTRUE(inefdec)) -x else x
  out <- log(2) - log(s) + stats::dnorm(e / s, log = TRUE) +
    stats::pnorm(lam * e / s, log.p = TRUE)
  if (isTRUE(log)) out else exp(out)
}

.check_scales <- function(sigma_u, sigma_v, what) {
  if (length(sigma_u) != 1L || !is.finite(sigma_u) || sigma_u < 0) {
    stop(what, "(): `sigma_u` must be a single finite number >= 0.",
      call. = FALSE
    )
  }
  if (length(sigma_v) != 1L || !is.finite(sigma_v) || sigma_v <= 0) {
    stop(what, "(): `sigma_v` must be a single finite positive number.",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

## P(eps <= q), computed as an expectation over u rather than as a double
## integral:
##
##   F(q) = E_u[ Phi( (q + s*u) / sigma_v ) ],   s = +1 production, -1 cost
##
## Conditioning on u leaves a normal CDF in closed form, so v is never drawn or
## integrated. That is the whole trick in Amsler et al.: it removes one source
## of error entirely, and it means the extreme tails inherit pnorm's own
## accuracy, which is reliable twenty standard deviations out.
##
## The UPPER tail is computed by negating that argument, not as 1 - F. At
## q = -12 with lambda = 1 the lower tail is around 3e-66; forming it as a
## complement would return exactly 1 and the copula that needed it would see
## qnorm(1) = Inf. This is the precision the gap was opened for.
pcomposed <- function(q, sigma_u, sigma_v, inefdec = TRUE,
                      lower.tail = TRUE, log.p = FALSE,
                      method = c("quadrature", "simulate"),
                      n_nodes = 128, R = 1e6, seed = NULL) {
  method <- match.arg(method)
  .check_scales(sigma_u, sigma_v, "pcomposed")
  if (!is.numeric(q)) stop("`q` must be numeric.", call. = FALSE)
  if (length(n_nodes) != 1L || !is.finite(n_nodes) || n_nodes < 16) {
    stop("pcomposed(): `n_nodes` must be a single number >= 16.", call. = FALSE)
  }
  ## Validated up here rather than inside the simulate branch, which the
  ## sigma_u = 0 and lambda = 1 shortcuts below can skip entirely -- an
  ## argument check that a fast path can jump over is not a check.
  if (identical(method, "simulate") &&
    (length(R) != 1L || !is.finite(R) || R < 1000)) {
    stop("pcomposed(): `R` must be a single number >= 1000.", call. = FALSE)
  }

  ## Production has eps = v - u, so a LARGER u makes eps smaller: u enters the
  ## lower-tail argument with a plus sign. Cost is the mirror.
  s <- if (isTRUE(inefdec)) 1 else -1
  sgn <- if (isTRUE(lower.tail)) 1 else -1

  ## sigma_u = 0 is the no-inefficiency boundary, where the composed error is
  ## just the noise. Worth handling exactly rather than sending a degenerate
  ## half-normal through the quadrature: sfm() reports fits on this boundary.
  if (sigma_u == 0) {
    return(stats::pnorm(q / sigma_v,
      lower.tail = lower.tail, log.p = log.p
    ))
  }

  ## NOTE: lambda = 1 has the only non-trivial closed form there is,
  ## P(Q) = Phi(Q/(sqrt(2) sigma_u))^2, and it is deliberately NOT used as a
  ## fast path. Three reasons, in order of importance:
  ##
  ##  1. It is stated for the COST case. Production needs its complement,
  ##     1 - Phi(-Q/(sqrt(2) sigma_u))^2 -- and a complement saturates to
  ##     exactly 1 in the far tail, which is the precise failure this whole
  ##     file exists to avoid.
  ##  2. The quadrature already reproduces it to 1e-13 relative at values down
  ##     to 1e-115, so the shortcut buys speed and no accuracy.
  ##  3. Two paths through one function is two chances to get the
  ##     production/cost sign backwards, and the first draft did.
  ##
  ## It earns its keep as a test oracle instead: test-pcomposed.R checks the
  ## quadrature against it, which is what an exact standard is for.

  if (method == "simulate") {
    if (!is.null(seed)) {
      st <- .rng_snapshot()
      on.exit(.rng_restore(st), add = TRUE)
      set.seed(seed)
    }
    u <- abs(stats::rnorm(as.integer(R))) * sigma_u
    out <- vapply(q, function(qi) {
      lg <- stats::pnorm(sgn * (qi + s * u) / sigma_v, log.p = TRUE)
      .logmeanexp(lg)
    }, numeric(1))
    return(if (isTRUE(log.p)) out else exp(out))
  }

  ## Gauss-Legendre on (0, 1) mapped to u in (0, Inf) by u = t/(1-t), the same
  ## map copsfm() uses, so the two agree about where the mass of a half-normal
  ## sits.
  ## The map carries sigma_u: u = sigma_u * t/(1-t), not t/(1-t). Without it
  ## the nodes sit at a fixed set of absolute values while the half-normal
  ## shrinks, so at sigma_u = 1e-8 every node falls in the far tail, every
  ## density underflows, and the function returns 0 for every q. Scaling makes
  ## the quadrature invariant in sigma_u -- and the sigma_u in the Jacobian
  ## then cancels the one in the density, which is why neither appears below.
  gl <- .gauss_legendre_01(as.integer(n_nodes))
  tt <- gl$nodes
  zz <- tt / (1 - tt)
  uu <- sigma_u * zz
  ## log of weight * jacobian * half-normal density at each node.
  lw <- log(gl$weights) - 2 * log1p(-tt) +
    log(2) + stats::dnorm(zz, log = TRUE)

  out <- vapply(q, function(qi) {
    lg <- stats::pnorm(sgn * (qi + s * uu) / sigma_v, log.p = TRUE)
    .logsumexp(lw + lg)
  }, numeric(1))
  if (isTRUE(log.p)) out else exp(out)
}

## Stable log(sum(exp(x))) and log(mean(exp(x))). The whole point of this file
## is the far tail, where every term underflows to zero if exponentiated first.
.logsumexp <- function(x) {
  x <- x[is.finite(x)]
  if (!length(x)) return(-Inf)
  m <- max(x)
  m + log(sum(exp(x - m)))
}

.logmeanexp <- function(x) {
  n <- length(x)
  x <- x[is.finite(x)]
  if (!length(x)) return(-Inf)
  m <- max(x)
  m + log(sum(exp(x - m))) - log(n)
}
