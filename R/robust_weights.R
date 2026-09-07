## Density-power weights: what the robust criteria actually do to each
## observation, and what they cannot see. Gap L17, from Bernstein, Parmeter
## and Wright (2026). See notes/code_history/robust_tuning.md.
##
## The weight responds to an observation far from the fitted SURFACE. It has
## no purchase on one that is wrong in a regressor and, because the surface
## bends toward it, ends up close to that surface -- so a mis-recorded input
## can be highly influential and still weigh nearly one. That is a property of
## these estimators, not of this code: their bounded influence is bounded with
## respect to RESPONSE contamination, conditional on the design. ?density_weights
## says so, and points at influence_sfa() for the other half.

density_weights <- function(object, sigma_v = NULL, sigma_u = NULL, c = NULL,
                            normalize = TRUE) {
  if (is.numeric(object)) {
    e <- object
    if (is.null(sigma_v) || is.null(sigma_u) || is.null(c))
      stop("When 'object' is a residual vector, sigma_v, sigma_u and c are required.",
           call. = FALSE)
  } else {
    e       <- .robust_residuals(object)
    sigma_v <- sigma_v %||% .robust_get(object, "sigma_v")
    sigma_u <- sigma_u %||% .robust_get(object, "sigma_u")
    c       <- c       %||% .robust_get(object, "c")
    if (is.null(c)) stop("No tuning parameter found on the fit; supply 'c'.",
                         call. = FALSE)
  }
  if (!is.finite(c) || c < 0) stop("c must be non-negative.", call. = FALSE)
  if (c <= 1e-10) return(rep(1, length(e)))

  f <- .dens_nhn(e, sigma_v, sigma_u)
  w <- f^c
  if (isTRUE(normalize)) w <- w / max(w)
  w
}
