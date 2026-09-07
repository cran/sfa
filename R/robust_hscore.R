## Hyvarinen-score evaluation for the robust divergence estimators. Gap L17,
## from Bernstein, Parmeter and Wright (2026).
## See notes/code_history/robust_tuning.md.
##
## The H-score of Hyvarinen (2005), in the form Sugasawa and Yonekura (2021)
## use to tune robust divergences, is
##
##     H(c) = n^-1 sum_i { 2 D''(e_i) + D'(e_i)^2 },   D = (f^c - 1)/c,
##
## with derivatives taken with respect to the composed residual. Evaluated
## directly, the second-derivative term carries f^(c-2). Because c < 1 the
## exponent sits near -2, so any observation whose fitted density underflows
## sends that term to infinity and the whole tuning candidate becomes
## unevaluable -- silently, and in the direction of too little robustness,
## because the candidates that survive are the ones closest to maximum
## likelihood.
##
## Expanding 2D'' + (D')^2 analytically leaves the density raised only to the
## positive powers c and 2c, so extreme residuals drive each term to zero
## instead. .hscore_nhn() evaluates that expansion by signed log-sum-exp.

## Signed log-sum-exp: given log|x_j| and sign(x_j), return log|sum x_j| and its
## sign, without forming any x_j on the natural scale.
.signed_logsumexp <- function(logs, signs) {
  m <- do.call(pmax, c(logs, list(na.rm = TRUE)))
  m[!is.finite(m)] <- 0
  s <- 0
  for (j in seq_along(logs)) s <- s + signs[[j]] * exp(logs[[j]] - m)
  list(sign = sign(s), log = m + log(abs(s)))
}

## Natural-scale evaluation, kept for validation only. This is the form that
## overflows; .hscore_nhn() falls back to it at c = 0, where D = log f and the
## density cancels out of D'/f and f''/f so the expression is already stable.
.hscore_nhn_natural <- function(e, sigma_v, sigma_u, c) {
  sigma  <- sqrt(sigma_v^2 + sigma_u^2)
  lambda <- sigma_u / sigma_v
  z  <- e / sigma
  f  <- .dens_nhn(e, sigma_v, sigma_u)
  A  <- z * pnorm(-lambda * z) + lambda * dnorm(lambda * z)
  B  <- (1 - z^2) * pnorm(-lambda * z) - lambda * z * (2 + lambda^2) * dnorm(lambda * z)
  fp  <- -(2 / sigma^2) * dnorm(z) * A
  fpp <- -(2 / sigma^3) * dnorm(z) * B
  if (c <= 1e-10) {
    lp  <- fp / f
    lpp <- fpp / f - lp^2
    return(mean(2 * lpp + lp^2))
  }
  mean(2 * (c - 1) * f^(c - 2) * fp^2 + 2 * f^(c - 1) * fpp + f^(2 * c - 2) * fp^2)
}

## Log-space evaluation. Same criterion, evaluable everywhere.
.hscore_nhn <- function(e, sigma_v, sigma_u, c) {
  if (c <= 1e-10) return(.hscore_nhn_natural(e, sigma_v, sigma_u, 0))

  sigma  <- sqrt(sigma_v^2 + sigma_u^2)
  lambda <- sigma_u / sigma_v
  z      <- e / sigma
  lphi   <- dnorm(z, log = TRUE)
  lPhi   <- pnorm(-lambda * z, log.p = TRUE)
  lphil  <- dnorm(lambda * z, log = TRUE)
  ls <- log(sigma); l2 <- log(2)

  A <- .signed_logsumexp(
    list(log(abs(z)) + lPhi, log(lambda) + lphil),
    list(sign(z), rep(1, length(z))))
  B <- .signed_logsumexp(
    list(log(abs(1 - z^2)) + lPhi,
         log(lambda * (2 + lambda^2)) + log(abs(z)) + lphil),
    list(sign(1 - z^2), -sign(z)))

  k2 <- l2 - ls
  lT1 <- log(2 * abs(c - 1)) + (c - 2) * k2 + log(4) - 4 * ls +
    c * lphi + (c - 2) * lPhi + 2 * A$log
  lT2 <- l2 + (c - 1) * k2 + l2 - 3 * ls + c * lphi + (c - 1) * lPhi + B$log
  lT3 <- (2 * c - 2) * k2 + log(4) - 4 * ls + 2 * c * lphi +
    (2 * c - 2) * lPhi + 2 * A$log

  tot <- .signed_logsumexp(
    list(lT1, lT2, lT3),
    list(rep(sign(c - 1), length(z)), -B$sign, rep(1, length(z))))
  mean(tot$sign * exp(tot$log))
}

## The user-facing evaluation. Lower is better; see ?hscore.
hscore <- function(e, sigma_v, sigma_u, c, stable = TRUE) {
  stopifnot(is.numeric(e), length(e) > 0L)
  stopifnot(length(sigma_v) == 1L, length(sigma_u) == 1L, length(c) == 1L)
  if (!is.finite(sigma_v) || sigma_v <= 0) stop("sigma_v must be positive.")
  if (!is.finite(sigma_u) || sigma_u <= 0) stop("sigma_u must be positive.")
  if (!is.finite(c) || c < 0) stop("c must be non-negative.")
  if (isTRUE(stable)) .hscore_nhn(e, sigma_v, sigma_u, c)
  else .hscore_nhn_natural(e, sigma_v, sigma_u, c)
}
