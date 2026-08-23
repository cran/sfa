## ============================================================================
## Robust divergence-based estimation (MLqE, Psi-likelihood, MDPD) as a
## generic wrapper on top of a model's own per-observation composed-error
## density, for use inside sfm()'s like.fn() closures. Exposed via sfm()'s
## `robust=` argument, which wraps a model's existing likelihood rather than
## forking into parallel model_name branches.
##
## Currently implemented for NHN only. MLqE's objective needs nothing
## model-specific beyond a density (`mean((f^c-1)/c)`, no integral
## correction), so it generalizes to sfm()'s other models for free once each
## one's density is isolated the same way `.dens_nhn()` is below. Psi and
## MDPD additionally need `integrate(f^(1+c), -Inf, Inf)`, computed at the
## model's current parameter values on every optimizer evaluation -- more
## per-model derivation and real per-evaluation optimizer overhead.
##
## Mathematical note worth keeping in mind when interpreting results:
## `obj_dpd(c) = (1+c) * obj_psi(c)` exactly (a positive rescaling of the
## same objective), so Psi and MDPD at the same c always produce IDENTICAL
## parameter point estimates -- MLqE is the only of the three that's a
## genuinely distinct estimator here. Both are exposed anyway
## (robust = "psi" / "mdpd") for compatibility with existing three-method
## naming conventions in the literature.
##
## Standard errors: for robust = "mlqe"/"psi"/"mdpd", the naive Hessian-
## inverse SE that ordinary MLE fits use is NOT valid -- it assumes the
## information-matrix equality (Hessian = negative expected outer-product of
## the score), which these M-estimators don't satisfy in general. The
## correct sandwich-form asymptotic variance is implemented instead -- see
## .sandwich_se_nhn() below.
## ============================================================================

## NHN composed-error density. Matches sfm.R's own NHN branch of like.fn()
## exactly (same formula, factored out as a standalone function of
## (e, sigma_v, sigma_u) instead of inlined against the optimizer's raw
## (lambda, sigma) parameter vector).
.dens_nhn <- function(e, sigma_v, sigma_u, log = FALSE) {
  sigma <- sqrt(sigma_v^2 + sigma_u^2)
  lambda <- sigma_u / sigma_v
  logdens <- log(2) - log(sigma) + dnorm(e / sigma, log = TRUE) + pnorm(-lambda * e / sigma, log.p = TRUE)
  if (log) logdens else pmax(exp(logdens), .Machine$double.xmin)
}

## Integral I_p(theta) = int f(e;theta)^p de for the NHN density, via the
## substitution t = e/sigma: f(e) = sigma^{-1} g(e/sigma; lambda),
## g(t) = 2*dnorm(t)*pnorm(-lambda*t) => int f^p de = sigma^{1-p} int g(t)^p dt.
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

## Per-observation NEGATIVE objective contribution (NOT summed). Exact
## algebraic decomposition of .robust_objective() below -- for psi/mdpd the
## power-integral correction term doesn't depend on i, so allocating it
## evenly across observations (dividing by n) reproduces the summed
## objective exactly: sum_i(this) == .robust_objective(...). This
## per-observation form is what M-estimation asymptotic theory needs (a
## genuine sum of per-observation loss contributions) and is used both here
## (via .robust_objective(), a thin sum() wrapper) and by
## .sandwich_vcov_nhn() below, which differentiates it observation-by-
## observation to build the sandwich variance's "meat" matrix.
##
##   mle  : -loglik_i
##   mlqe : -(f_i^c - 1)/c
##   psi  : -( f_i^c/c - I_(1+c)/(1+c) )
##   mdpd : -( ((1+c)/c)*f_i^c - I_(1+c) )
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
## every likelihood in this package's like.fn()/fn() closures uses (opts.R's
## bobyqa/psoptim/optim wrappers all MINIMIZE; see ttsfm.R's header comment
## on this convention, and the real bug it caused there when violated).
##
## `loglik` is the model's own per-observation log-density vector at the
## CURRENT parameter values (exactly what every like.fn() already computes
## before summing for the "mle" case). `power_integral_fn(p)` computes
## int f^p de at those SAME current parameter values -- only needed for
## psi/mdpd, supplied per model (NHN's is .nhn_power_integral above).
##
##   mle  : standard MLE.                    -sum(loglik)
##   mlqe : q = 1-c.                         -sum((f^c - 1)/c)
##   psi  : eta = c.                         -( sum(f^c)/c - n*I_(1+c)/(1+c) )
##   mdpd : alpha = c.                       -( ((1+c)/c)*sum(f^c) - n*I_(1+c) )
## c <= 0 (or NULL) falls back to plain MLE for any method, so robust="mlqe"
## with c=0 (say) behaves identically to robust="mle" rather than hitting a
## 0/0 in (f^c-1)/c.
.robust_objective <- function(method = c("mle", "mlqe", "psi", "mdpd"),
                              loglik, c = NULL, power_integral_fn = NULL) {
  method <- match.arg(method)
  vec <- .robust_objective_vec(method, loglik, c, power_integral_fn)
  if (method == "mle" || is.null(c) || c <= 1e-10) {
    return(sum(vec[is.finite(vec)]))
  }
  sum(vec)
}

## ============================================================================
## Sandwich standard errors for the NHN robust M-estimators.
##
## The naive Hessian-inverse SE (valid for ordinary MLE, which satisfies the
## information-matrix equality) is NOT valid for MLqE/Psi/MDPD in general.
## The correct asymptotic variance for an M-estimator defined by minimizing
## sum_i rho_i(theta) is the sandwich form Var(theta_hat) = A^-1 B A^-1,
## where:
##   A ("bread") = Hessian of the SUMMED objective at theta_hat -- exactly
##     opt$hessian, already computed by optim(hessian=TRUE) for the point
##     estimate's own optimization, reused here rather than recomputed.
##   B ("meat")  = sum_i g_i g_i', g_i = gradient of the PER-OBSERVATION
##     objective rho_i (.robust_objective_vec()'s output) w.r.t. theta at
##     theta_hat. Obtained via numDeriv::jacobian() (numerical, since no
##     analytic per-model score is implemented) applied to the vector-valued
##     per-observation objective -- jacobian() differentiates a
##     vector-valued function w.r.t. a parameter vector, giving exactly the
##     n x k matrix of per-observation gradients needed, one call.
##
## No sign correction is needed between A and B: both are computed from the
## SAME rho_i = .robust_objective_vec(...) (the per-observation NEGATIVE
## objective, i.e. the thing actually minimized) that opt$hessian already
## used, so this is internally consistent by construction.
##
## Sanity property this was checked against (including a Monte Carlo check
## against the empirical spread of repeated fits): as c -> 0, MLqE/Psi/MDPD
## all converge to ordinary MLE, where the information-matrix equality holds
## (A ~= B asymptotically), so this sandwich formula reduces to
## (approximately) the ordinary Hessian-inverse SE in that limit.
##
## Returns a vector of SEs (sqrt of the diagonal), or a vector of NA (silent
## -- the caller is expected to have already decided SEs are wanted) if
## either matrix inversion fails.
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
