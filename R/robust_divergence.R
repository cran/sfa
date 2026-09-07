## Robust divergence estimation (MLqE, Psi-likelihood, MDPD), wrapped
## generically around a model's own per-observation density.

## NHN composed-error density.
.dens_nhn <- function(e, sigma_v, sigma_u, log = FALSE) {
  sigma <- sqrt(sigma_v^2 + sigma_u^2)
  lambda <- sigma_u / sigma_v
  logdens <- log(2) - log(sigma) + dnorm(e / sigma, log = TRUE) + pnorm(-lambda * e / sigma, log.p = TRUE)
  if (log) logdens else pmax(exp(logdens), .Machine$double.xmin)
}

## Integral I_p(theta) = int f(e;theta)^p de for the NHN density, via the
## substitution t = e/sigma: f(e) = sigma^{-1} g(e/sigma; lambda).
.nhn_power_integral <- function(sigma_v, sigma_u, p, rel.tol = 1e-8) {
  if (p <= 1) {
    return(1)
  }
  sigma <- sqrt(sigma_v^2 + sigma_u^2)
  lambda <- sigma_u / sigma_v
  integrand <- function(t) (2 * dnorm(t) * pnorm(-lambda * t))^p
  val <- integrate(integrand, lower = -Inf, upper = Inf, rel.tol = rel.tol, subdivisions = 400)$value
  sigma^(1 - p) * val
}

## Per-observation NEGATIVE objective contribution (NOT summed).
.robust_objective_vec <- function(method = c("mle", "mlqe", "psi", "mdpd"),
                                  loglik, c = NULL, power_integral_fn = NULL) {
  method <- match.arg(method)

  if (method == "mle" || is.null(c) || c <= 1e-10) {
    return(-loglik)
  }

  f <- exp(loglik)

  if (method == "mlqe") {
    return(-(f^c - 1) / c)
  }

  if (is.null(power_integral_fn)) {
    stop(".robust_objective_vec(): method '", method, "' requires power_integral_fn.", call. = FALSE)
  }
  I <- power_integral_fn(1 + c)

  if (method == "psi") {
    return(-(f^c / c - I / (1 + c)))
  }
  if (method == "mdpd") {
    return(-(((1 + c) / c) * f^c - I))
  }
}

## Generic robust-divergence objective, on the same "negative summed" scale
## every likelihood in this package's like.fn()/fn() closures uses.
.robust_objective <- function(method = c("mle", "mlqe", "psi", "mdpd"),
                              loglik, c = NULL, power_integral_fn = NULL) {
  method <- match.arg(method)
  vec <- .robust_objective_vec(method, loglik, c, power_integral_fn)
  if (method == "mle" || is.null(c) || c <= 1e-10) {
    return(sum(vec[is.finite(vec)]))
  }
  sum(vec)
}

## Sandwich standard errors for the NHN robust M-estimators.
.sandwich_se_nhn <- function(par_hat, hessian, per_obs_fn) {
  k <- length(par_hat)
  na_result <- rep(NA_real_, k)

  A_inv <- tryCatch(solve(hessian), error = function(e) NULL)
  if (is.null(A_inv)) {
    return(na_result)
  }

  G <- tryCatch(numDeriv::jacobian(per_obs_fn, par_hat), error = function(e) NULL)
  if (is.null(G) || !all(is.finite(G))) {
    return(na_result)
  }

  B <- t(G) %*% G
  V <- tryCatch(A_inv %*% B %*% A_inv, error = function(e) NULL)
  if (is.null(V)) {
    return(na_result)
  }

  sqrt(pmax(diag(V), 0))
}
