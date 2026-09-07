## Heteroskedastic cross-sectional frontiers. sigma_u, sigma_v and (for the
## truncated-normal family) the pre-truncation mean may each be driven by their
## own covariates. This is the engine behind sfm()'s `vhet` and `muhet`
## arguments; see notes/code_history/het.md for the derivations.

## Build one variance-determinant design matrix. NULL means intercept-only,
## i.e. the homoskedastic case, which keeps a single code path for both.
.het_design <- function(f, data, label) {
  n <- nrow(data)
  if (is.null(f)) {
    return(matrix(1, nrow = n, ncol = 1L, dimnames = list(NULL, "(Intercept)")))
  }
  tt <- tryCatch(stats::delete.response(stats::terms(stats::as.formula(f), data = data)),
    error = function(e) {
      stop("sfm(): `", label, "` could not be interpreted as a formula: ",
        conditionMessage(e),
        call. = FALSE
      )
    }
  )
  miss <- setdiff(all.vars(tt), names(data))
  if (length(miss)) {
    stop("sfm(): `", label, "` refers to ",
      if (length(miss) > 1L) "variables " else "variable ",
      paste(sQuote(miss), collapse = ", "), " not found in `data`.",
      call. = FALSE
    )
  }
  mf <- stats::model.frame(tt, data = data, na.action = stats::na.pass)
  X <- stats::model.matrix(tt, mf)
  if (nrow(X) != n || anyNA(X)) {
    stop("sfm(): `", label, "` has missing values on rows that the frontier ",
      "formula keeps. Drop those rows before calling sfm().",
      call. = FALSE
    )
  }
  X
}

## Complete-case filter over the het formulas, applied BEFORE data_proc() so
## that its row bookkeeping and the het design matrices cannot fall out of
## step.
.het_prefilter <- function(data, forms) {
  vs <- unique(unlist(lapply(Filter(Negate(is.null), forms), all.vars)))
  vs <- intersect(vs, names(data))
  if (!length(vs)) {
    return(data)
  }
  keep <- stats::complete.cases(data[, vs, drop = FALSE])
  data[keep, , drop = FALSE]
}

## Method-of-moments starting values from the OLS residuals, per family.
.het_start_moments <- function(e, family) {
  m2 <- mean(e^2)
  m3 <- mean(e^3)
  sd_e <- sqrt(max(m2, .Machine$double.eps))
  ## A production frontier implies m3 < 0. Under wrong skew the moment
  ## equations have no admissible root (Olson, Schmidt and Waldman 1980), so
  ## fall back to a neutral split rather than failing here -- the likelihood
  ## path is entitled to find the boundary solution on its own.
  s_u <- if (identical(family, "exponential")) {
    if (m3 < 0) (-m3 / 2)^(1 / 3) else NA_real_
  } else {
    k <- sqrt(2 / pi) * (1 - 4 / pi)
    if (m3 < 0) (m3 / k)^(1 / 3) else NA_real_
  }
  if (!is.finite(s_u) || s_u <= 0) s_u <- 0.5 * sd_e
  v_share <- if (identical(family, "exponential")) s_u^2 else (1 - 2 / pi) * s_u^2
  s_v <- sqrt(max(m2 - v_share, 0.1 * m2))
  e_u <- if (identical(family, "exponential")) s_u else s_u * sqrt(2 / pi)
  list(sigma_u = s_u, sigma_v = s_v, E_u = e_u)
}

## The per-observation log-likelihood, given the three fitted vectors.
##
## The half-normal case is the truncated-normal case at mu = 0: with mu = 0 the
## last term is -log Phi(0) = log 2, which is exactly the half-normal
## normalizer. Both therefore share one expression rather than two that have to
## be kept in step.
.het_loglik <- function(eps, s_u, s_v, mu, family) {
  if (identical(family, "exponential")) {
    z <- -(eps / s_v) - (s_v / s_u)
    return(-log(s_u) - (eps^2 / (2 * s_v^2)) + .log_phi_tilt(z))
  }
  s2 <- s_u^2 + s_v^2
  s <- sqrt(s2)
  lam <- s_u / s_v
  -log(s) + stats::dnorm((eps + mu) / s, log = TRUE) +
    stats::pnorm((mu / lam - eps * lam) / s, log.p = TRUE) -
    stats::pnorm(mu / s_u, log.p = TRUE)
}

## The posterior of u given the composed residual. Truncated normal in all
## three families; only (mu_star, sigma_star) differ.
.het_posterior <- function(eps, s_u, s_v, mu, family) {
  if (identical(family, "exponential")) {
    return(list(mu_star = -eps - s_v^2 / s_u, sigma_star = s_v))
  }
  s2 <- s_u^2 + s_v^2
  list(
    mu_star = (s_v^2 * mu - s_u^2 * eps) / s2,
    sigma_star = s_u * s_v / sqrt(s2)
  )
}

## Fit one heteroskedastic cross-sectional frontier.
##
## Everything that must stay positive is carried on the log scale, so the
## parameter vector is unconstrained and the optimizer needs no bounds. That is
## also why this path uses nlminb followed by an unbounded BFGS rather than
## sfm()'s bounded four-stage scaffold: the scaffold exists to keep sigmas off
## the zero boundary, and there is no boundary here to keep off.
.sfm_het_fit <- function(family, Y, X, Zu, Zv, Zmu, inefdec_n, z_sigma, z_link,
                         x_names, maxit.nlminb = 500, maxit.optim = 1000,
                         optHessian = TRUE, verbose = FALSE, Zs = NULL,
                         wts = NULL) {
  n <- length(as.numeric(Y))
  X <- as.matrix(X)
  n_b <- ncol(X)
  n_v <- ncol(Zv)
  n_u <- ncol(Zu)
  n_m <- if (identical(family, "truncnormal")) ncol(Zmu) else 0L
  ## The SCALING PROPERTY (Wang and Schmidt 2002; Alvarez, Amsler, Orea and
  ## Schmidt 2006): u_i = h(z_i, delta) * u*_i, so the covariates scale a
  ## common draw rather than entering its variance. Both sigma_u and the
  ## pre-truncation mean move by the SAME factor, which is what fixes the shape
  ## of the distribution and leaves only its scale free -- the identification
  ## argument those papers make. Zs NULL is the ordinary heteroskedastic path,
  ## unchanged.
  scaling <- !is.null(Zs)
  n_s <- if (scaling) ncol(Zs) else 0L

  ib <- seq_len(n_b)
  iv <- n_b + seq_len(n_v)
  iu <- n_b + n_v + seq_len(n_u)
  im <- if (n_m > 0L) n_b + n_v + n_u + seq_len(n_m) else integer(0)
  is_ <- if (n_s > 0L) n_b + n_v + n_u + n_m + seq_len(n_s) else integer(0)

  ## exp() of an unbounded linear predictor is the one place this
  ## parameterization can overflow, so clip the predictor rather than repairing
  ## the result.
  cl <- .SFA_CONSTANTS$EXP_CLIP_UPPER
  .eta <- function(Z, d) pmin(pmax(as.numeric(Z %*% d), -cl), cl)

  .parts <- function(p) {
    ## Under the scaling property one factor h multiplies BOTH the scale and
    ## the pre-truncation mean. Zu and Zmu are intercept-only there, so sigma_u
    ## and mu are the scalars of the common draw u* and h carries all the
    ## covariate dependence.
    h <- if (scaling) exp(.eta(Zs, p[is_])) else 1
    list(
      eps = as.numeric(inefdec_n * (Y - X %*% p[ib])),
      s_v = z_sigma(.eta(Zv, p[iv])),
      s_u = z_sigma(.eta(Zu, p[iu])) * h,
      mu = (if (n_m > 0L) as.numeric(Zmu %*% p[im]) else 0) * h
    )
  }

  nll <- function(p) {
    if (!all(is.finite(p))) {
      return(1e12)
    }
    q <- .parts(p)
    if (!all(is.finite(q$s_v)) || !all(is.finite(q$s_u)) ||
      any(q$s_v <= 0) || any(q$s_u <= 0) || !all(is.finite(q$mu))) {
      return(1e12)
    }
    ll <- .het_loglik(q$eps, q$s_u, q$s_v, q$mu, family)
    ## A finite floor, not -Inf and not -.Machine$double.xmax: optim()
    ## differences the objective for its gradient, and differencing xmax
    ## overflows to a non-finite value that aborts the fit outright.
    ll[!is.finite(ll)] <- -1e12 / n
    if (!is.null(wts)) ll <- wts * ll
    v <- -sum(ll)
    if (is.finite(v)) v else 1e12
  }

  ## Starting values.
  ols <- stats::lm.fit(X, as.numeric(Y))
  b0 <- unname(ols$coefficients)
  b0[!is.finite(b0)] <- 0
  e0 <- as.numeric(inefdec_n * ols$residuals)
  mom <- .het_start_moments(e0, family)
  ## OLS absorbs -E[u] into the intercept; shift it back so the frontier starts
  ## at the frontier rather than through the middle of the cloud.
  icol <- match("(Intercept)", x_names)
  if (!is.na(icol)) b0[icol] <- b0[icol] + inefdec_n * mom$E_u
  ## The link decides what the intercept of each variance block means: on the
  ## "sd" link the predictor is log sigma, on "var" it is log sigma^2.
  .d0 <- function(s) if (identical(z_link, "sd")) log(s) else log(s^2)
  ## Put the moment estimate on the block's INTERCEPT, which is not always the
  ## first column -- `vhet = ~ 0 + z` has none, and a user may order terms
  ## however they like. Falling back to column one keeps a sane start either
  ## way; this is a starting value, not a constraint.
  .block0 <- function(Z, s) {
    d <- rep(0, ncol(Z))
    j <- match("(Intercept)", colnames(Z), nomatch = 1L)
    d[j] <- .d0(s)
    d
  }
  dv0 <- .block0(Zv, mom$sigma_v)
  du0 <- .block0(Zu, mom$sigma_u)
  dm0 <- if (n_m > 0L) rep(0, n_m) else numeric(0)
  ## delta starts at zero, i.e. h = 1 and the fit begins from the homoskedastic
  ## model. Any other start would assert a direction for the scaling before
  ## the data has been asked.
  ds0 <- if (n_s > 0L) rep(0, n_s) else numeric(0)
  start_v <- c(b0, dv0, du0, dm0, ds0)

  nms <- c(
    x_names,
    paste0("Zv.", colnames(Zv)),
    paste0("Zu.", colnames(Zu)),
    if (n_m > 0L) paste0("Zmu.", colnames(Zmu)) else NULL,
    if (n_s > 0L) paste0("scale.", colnames(Zs)) else NULL
  )
  names(start_v) <- nms

  Start.Time <- start.time()
  fit1 <- tryCatch(
    suppressWarnings(stats::nlminb(start_v, nll,
      control = list(iter.max = maxit.nlminb, trace = if (verbose) 10L else 0L)
    )),
    error = function(e) NULL
  )
  if (!is.null(fit1) && all(is.finite(fit1$par)) && is.finite(fit1$objective) &&
    fit1$objective <= nll(start_v)) {
    start_v <- fit1$par
  }

  opt <- tryCatch(
    suppressWarnings(stats::optim(start_v, nll,
      method = "BFGS", hessian = isTRUE(optHessian),
      control = list(maxit = maxit.optim, trace = if (verbose) 1L else 0L)
    )),
    error = function(e) NULL
  )
  if (is.null(opt) || !is.finite(opt$value)) {
    stop("sfm(): the heteroskedastic fit did not converge. Try simplifying ",
      "`vhet`/`muhet`, or rescaling their covariates -- the variance blocks ",
      "are on a log scale, so a covariate spanning several orders of ",
      "magnitude drives the linear predictor into the clipping bound.",
      call. = FALSE
    )
  }
  End.Time <- end.time(Start.Time)
  names(opt$par) <- nms

  st_err <- rep(NA_real_, length(opt$par))
  if (isTRUE(optHessian) && !is.null(opt$hessian)) {
    vc <- tryCatch(suppressWarnings(solve(opt$hessian)), error = function(e) NULL)
    if (!is.null(vc)) st_err <- suppressWarnings(sqrt(pmax(diag(vc), 0)))
  }

  out <- matrix(NA_real_, nrow = 3L, ncol = length(opt$par),
    dimnames = list(c("par", "st_err", "t-val"), nms)
  )
  out[1, ] <- opt$par
  out[2, ] <- st_err
  out[3, ] <- opt$par / st_err

  q <- .parts(opt$par)
  post <- .het_posterior(q$eps, q$s_u, q$s_v, q$mu, family)
  sig_star <- rep_len(post$sigma_star, length(post$mu_star))

  list(
    out = out, opt = opt, start_v = start_v, total_time = End.Time,
    objective = nll, n_blocks = c(beta = n_b, v = n_v, u = n_u, mu = n_m),
    exp_u_hat = .te_battese_coelli(post$mu_star, sig_star),
    u_hat = .jlms_u(post$mu_star, sig_star),
    u_posterior = list(mu_star = post$mu_star, sigma_star = sig_star),
    sigma_u = q$s_u, sigma_v = q$s_v,
    mu = rep_len(q$mu, n), eps = q$eps
  )
}
