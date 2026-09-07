## npsfm(): nonparametric stochastic frontier models

.npsfm_need_np <- function() {
  if (!requireNamespace("np", quietly = TRUE)) {
    stop("npsfm() requires the 'np' package for kernel regression and ",
      "bandwidth selection, which is not installed. ",
      'Install it with install.packages("np").',
      call. = FALSE
    )
  }
  invisible(TRUE)
}

## Fan, Li and Weersink's concentrated log-likelihood, their equation (14) and
## the discussion beneath (15).
.flw_neg_cll <- function(lambda, e) {
  if (!is.finite(lambda) || lambda <= 0) {
    return(.SFA_CONSTANTS$MAX_VALUE)
  }
  sc <- sqrt(2 * lambda^2 / (pi * (1 + lambda^2)))
  s2 <- mean(e^2) / (1 - sc^2)
  if (!is.finite(s2) || s2 <= 0) {
    return(.SFA_CONSTANTS$MAX_VALUE)
  }
  si <- sqrt(s2)
  ## Mean-correct the residuals: E[eps] = -E[u] = -sigma*sc
  ep <- e - si * sc
  ll <- -length(ep) * log(si) +
    sum(pnorm(-ep * lambda / si, log.p = TRUE)) -
    0.5 * sum(ep^2) / s2
  if (!is.finite(ll)) {
    return(.SFA_CONSTANTS$MAX_VALUE)
  }
  -ll
}

## Kernel regression of ydat on xdat, returning fitted values, gradients and
## residuals.
.npsfm_kreg <- function(xdat, ydat, regtype, bw.sel, bw = NULL, nmulti = 1L) {
  if (is.null(bw)) {
    bws <- np::npregbw(
      xdat = xdat, ydat = ydat, regtype = regtype,
      bwmethod = bw.sel, nmulti = nmulti
    )
  } else {
    bws <- np::npregbw(
      xdat = xdat, ydat = ydat, regtype = regtype,
      bandwidth.compute = FALSE
    )
    bws$bw <- bw
  }
  fit <- np::npreg(bws = bws, gradients = TRUE, residuals = TRUE)
  ## Namespace-qualify everything: np is only in Suggests, so it is not
  ## attached when npsfm() runs.
  list(
    bws = bws,
    fitted = as.numeric(stats::fitted(fit)),
    grad = as.matrix(np::gradients(fit)),
    resid = as.numeric(stats::residuals(fit))
  )
}


## FLW: Fan, Li and Weersink (1996)
.npsfm_flw <- function(xdat, ydat, dist, regtype, bw.sel, bw, cost) {
  cost_sc <- if (cost) -1 else 1
  k <- .npsfm_kreg(
    xdat = xdat, ydat = cost_sc * ydat, regtype = regtype,
    bw.sel = bw.sel, bw = bw
  )
  e <- k$resid
  m2 <- mean(e^2)
  m3 <- mean(e^3)
  m4 <- mean(e^4)

  sigma_u <- sigma_v <- lambda <- sigma <- NA_real_
  theta <- bhat <- NA_real_

  if (dist == "hn") {
    ## Identification of sigma_u rests on the residuals being negatively
    ## skewed.
    if (m3 >= 0) {
      warning("FLW: residuals are not negatively skewed (third moment ",
        signif(m3, 3), "), so sigma_u is not identified and will be ",
        "driven to zero. Check the `cost` orientation.",
        call. = FALSE
      )
    }
    ## Search lambda on a wide bounded interval; the concentrated likelihood
    ## is one-dimensional so this is both cheap and reliable.
    opt <- stats::optimize(.flw_neg_cll, interval = c(1e-4, 50), e = e)
    lambda <- opt$minimum
    sc <- sqrt(2 * lambda^2 / (pi * (1 + lambda^2)))
    sigma <- sqrt(m2 / (1 - sc^2))
    sigma_v <- sigma / sqrt(1 + lambda^2)
    sigma_u <- sigma_v * lambda
    mu <- sigma * sc
  } else if (dist == "exp") {
    ## u ~ Exp(rate theta): mu3(eps) = -2/theta^3
    if (m3 >= 0) {
      warning("FLW: third residual moment is non-negative (wrong skew); ",
        "the exponential scale is not identified from these data.",
        call. = FALSE
      )
      theta <- NA_real_
      mu <- 0
    } else {
      theta <- (-2 / m3)^(1 / 3)
      sigma_u <- 1 / theta
      sigma_v <- sqrt(max(m2 - 1 / theta^2, 0))
      mu <- 1 / theta
    }
  } else if (dist == "gamma") {
    ## u ~ Gamma(P, theta); uses the third and fourth central moments
    m4star <- m4 - 3 * m2^2
    theta <- -3 * m3 / m4star
    p_hat <- -theta^3 * m3 / 2
    sigma_v <- sqrt(max(m2 - p_hat / theta^2, 0))
    sigma_u <- if (is.finite(p_hat) && p_hat > 0) sqrt(p_hat) / theta else NA_real_
    mu <- p_hat / theta
  } else if (dist == "unif") {
    ## u ~ U(0, b)
    bhat2 <- (3 * m2^2 - m4) * 120
    if (!is.finite(bhat2) || bhat2 < 0) {
      warning("FLW: implied uniform upper bound is negative; setting it to zero.",
        call. = FALSE
      )
      bhat <- 0
    } else {
      bhat <- bhat2^(1 / 4)
    }
    sigma_v <- sqrt(max(m2 - bhat^2 / 12, 0))
    sigma_u <- bhat / sqrt(12)
    mu <- bhat / 2
  }

  ## The kernel fit estimates E[y|x] = m(x) - E[u]; shift it back up to the
  ## frontier. Gradients are unaffected because the shift is constant.
  frontier <- cost_sc * (k$fitted + mu)
  frontier_grad <- cost_sc * k$grad
  eps <- cost_sc * e

  list(
    frontier = frontier, frontier_grad = frontier_grad, residuals = eps,
    sigma_u = sigma_u, sigma_v = sigma_v, lambda = lambda, sigma = sigma,
    theta = theta, b = bhat, mean_correction = mu,
    conditional_mean = cost_sc * k$fitted, bws = k$bws
  )
}


## SVKZ: Simar, Van Keilegom and Zelenyuk (2017)
.npsfm_svkz <- function(xdat, ydat, bw.sel, cost) {
  cost_sc <- if (cost) -1 else 1
  yy <- cost_sc * ydat

  ## r1: conditional mean
  r1 <- .npsfm_kreg(xdat = xdat, ydat = yy, regtype = "ll", bw.sel = bw.sel)
  e1 <- r1$resid

  ## r2, r3: conditional second and third moments of the residual
  r2 <- .npsfm_kreg(xdat = xdat, ydat = e1^2, regtype = "ll", bw.sel = bw.sel)
  r3 <- .npsfm_kreg(xdat = xdat, ydat = e1^3, regtype = "ll", bw.sel = bw.sel)

  ## SVKZ (3.5): sigma_u(x)^3 = sqrt(pi/2) * (pi/(pi-4)) * r3(x).
  sigu3 <- sqrt(pi / 2) * (pi / (pi - 4)) * r3$fitted
  wrong_skew <- sigu3 <= 0
  sigma_u <- pmax(sigu3, 0)^(1 / 3)

  ## (3.8): sigma_v(x)^2 = r2(x) - sigma_u(x)^2 (pi-2)/pi, floored at zero
  sigma_v2 <- pmax(r2$fitted - sigma_u^2 * (pi - 2) / pi, 0)

  ## (3.2)/(3.9): mean correction E[u|x] = sqrt(2/pi) sigma_u(x)
  mu <- sqrt(2 / pi) * sigma_u

  ## Gradient of the correction, by the chain rule through (3.5).
  dsigu3 <- sqrt(pi / 2) * (pi / (pi - 4)) * r3$grad
  scale_g <- (1 / 3) * ifelse(wrong_skew, 0, pmax(sigu3, .SFA_CONSTANTS$MIN_POSITIVE)^(-2 / 3))
  sigma_u_grad <- dsigu3 * scale_g
  sigma_u_grad[wrong_skew, ] <- 0

  frontier <- cost_sc * (r1$fitted + mu)
  frontier_grad <- cost_sc * (r1$grad + sqrt(2 / pi) * sigma_u_grad)

  list(
    frontier = frontier, frontier_grad = frontier_grad,
    residuals = cost_sc * e1,
    sigma_u = sigma_u, sigma_v = sqrt(sigma_v2),
    sigma_u_grad = sigma_u_grad, mean_correction = mu,
    conditional_mean = cost_sc * r1$fitted, wrong_skew = wrong_skew,
    bws = list(r1 = r1$bws, r2 = r2$bws, r3 = r3$bws)
  )
}


## User-facing entry point
npsfm <- function(formula, data, method = c("FLW", "SVKZ", "PSZ", "KPST", "MY", "SZ"),
                  dist = c("hn", "exp", "gamma", "unif"),
                  regtype = c("lc", "ll"), bw.sel = c("cv.ls", "cv.aic"),
                  bw = NULL, cost = FALSE, eff = TRUE,
                  maxit = 5000, tol = 1e-3, iter = 25,
                  rts = c("vrs", "crs", "drs", "irs"),
                  prior.fit = NULL, log.form = TRUE, verbose = FALSE) {
  .npsfm_need_np()

  ## np prints "Multistart 1 of 1" progress from every bandwidth search, which
  ## floods the console when npsfm() is called in a loop.
  old_np_msg <- getOption("np.messages")
  options(np.messages = isTRUE(verbose))
  on.exit(options(np.messages = old_np_msg), add = TRUE)

  method <- .match_model_name(method, eval(formals()$method), arg = "method")
  regtype <- match.arg(regtype)
  bw.sel <- match.arg(bw.sel)
  dist <- match.arg(dist)
  rts <- match.arg(rts)

  ## "KPST" is the alternative name for the Park-Simar-Zelenyuk estimator that
  ## the source scripts accepted; keep it working but normalize internally.
  if (method == "KPST") method <- "PSZ"

  if (method != "FLW" && dist != "hn") {
    stop("method = \"", method, "\" is derived for the normal-half normal ",
      "model only; dist = \"", dist, "\" is available for method = \"FLW\" ",
      "alone.",
      call. = FALSE
    )
  }

  cl <- match.call()

  ## The frontier is the first pipe segment; npsfm() has no variance
  ## determinants of its own.
  parsed <- .parse_pipe_formula(formula)
  if (!is.null(parsed$formula_z)) {
    stop("npsfm() takes a single-part formula (y ~ x1 + x2). ",
      "Variance determinants enter nonparametrically through the ",
      "covariates themselves; there is no `| z` segment.",
      call. = FALSE
    )
  }

  mf <- stats::model.frame(parsed$formula_x, data = data, na.action = stats::na.omit)
  ydat <- as.numeric(stats::model.response(mf))
  xdat <- mf[, -1, drop = FALSE]
  if (!ncol(xdat)) {
    stop("npsfm() needs at least one covariate on the right-hand side.", call. = FALSE)
  }
  n_obs <- length(ydat)

  start_time <- Sys.time()
  fit <- switch(method,
    FLW  = .npsfm_flw(xdat, ydat, dist, regtype, bw.sel, bw, cost),
    SVKZ = .npsfm_svkz(xdat, ydat, bw.sel, cost),
    PSZ  = .npsfm_psz(xdat, ydat, bw.sel, bw, cost, maxit, verbose),
    MY   = .npsfm_my(xdat, ydat, bw.sel, bw, cost, tol, iter, maxit, verbose),
    SZ   = .npsfm_sz(xdat, ydat, bw.sel, bw, cost, rts, prior.fit, log.form)
  )
  total_time <- Sys.time() - start_time

  ## Observation-level inefficiency predictions. The composed residual is
  ## measured against the corrected frontier, so eps = y - frontier.
  eps <- ydat - fit$frontier
  u_hat <- exp_u_hat <- NULL
  if (isTRUE(eff)) {
    sc <- if (cost) -1 else 1
    su <- fit$sigma_u
    sv <- fit$sigma_v
    if (dist == "hn" && all(is.finite(c(su, sv))) && all(sv > 0)) {
      s2 <- su^2 + sv^2
      mu_star <- -sc * eps * su^2 / s2
      sig_star <- su * sv / sqrt(s2)
      u_hat <- .jlms_u(mu_star, sig_star)
      exp_u_hat <- .te_battese_coelli(mu_star, sig_star)
    } else if (dist == "exp" && is.finite(fit$theta) && is.finite(sv) && sv > 0) {
      mu_star <- -sc * eps - sv^2 * fit$theta
      u_hat <- .jlms_u(mu_star, rep(sv, length(mu_star)))
      exp_u_hat <- .te_battese_coelli(mu_star, rep(sv, length(mu_star)))
    }
  }

  out <- list(
    method = method, dist = dist, formula = formula, call = cl,
    frontier = fit$frontier, frontier.grad = fit$frontier_grad,
    conditional.mean = fit$conditional_mean,
    residuals = eps, mean.correction = fit$mean_correction,
    sigma.u = fit$sigma_u, sigma.v = fit$sigma_v,
    lambda = fit$lambda, sigma = fit$sigma,
    theta = fit$theta, b = fit$b,
    sigma.u.grad = fit$sigma_u_grad, sigma.v.grad = fit$sigma_v_grad,
    wrong.skew = fit$wrong_skew, convergence = fit$convergence,
    iterations = fit$iterations, converged = fit$converged,
    tol.reached = fit$tol_reached, prior.fit = fit$prior.fit,
    dea.efficiency = fit$dea.efficiency, rts = fit$rts,
    u_hat = u_hat, exp_u_hat = exp_u_hat,
    bws = fit$bws, cost = cost, regtype = regtype, bw.sel = bw.sel,
    nobs = n_obs, total_time = total_time, data = data
  )
  out <- out[!vapply(out, is.null, logical(1))]
  class(out) <- "npsfareg"
  out
}


## Methods.
print.npsfareg <- function(x, ...) {
  cat("Nonparametric stochastic frontier model\n")
  cat("Method:       ", x$method,
    switch(x$method,
      FLW = " (Fan, Li and Weersink 1996)",
      SVKZ = " (Simar, Van Keilegom and Zelenyuk 2017)"
    ), "\n",
    sep = ""
  )
  cat("Distribution: ", x$dist, "\n", sep = "")
  cat("Frontier:     ", if (x$cost) "cost" else "production", "\n", sep = "")
  cat("Observations: ", x$nobs, "\n", sep = "")
  if (length(x$sigma.u) == 1L) {
    cat("\nsigma_u = ", signif(x$sigma.u, 5),
      "   sigma_v = ", signif(x$sigma.v, 5), "\n",
      sep = ""
    )
    if (!is.null(x$lambda) && is.finite(x$lambda)) {
      cat("lambda  = ", signif(x$lambda, 5),
        "   sigma   = ", signif(x$sigma, 5), "\n",
        sep = ""
      )
    }
  } else {
    cat("\nsigma_u(x) and sigma_v(x) vary by observation:\n")
    print(summary(data.frame(sigma_u = x$sigma.u, sigma_v = x$sigma.v)))
    if (!is.null(x$wrong.skew) && any(x$wrong.skew)) {
      cat("\n", sum(x$wrong.skew), " of ", x$nobs,
        " points had the wrong residual skew and were floored at sigma_u = 0.\n",
        sep = ""
      )
    }
  }
  if (!is.null(x$exp_u_hat)) {
    cat("\nMean predicted efficiency E[exp(-u)|e]: ",
      signif(mean(x$exp_u_hat, na.rm = TRUE), 4), "\n",
      sep = ""
    )
  }
  invisible(x)
}

summary.npsfareg <- function(object, ...) {
  print(object)
  cat("\nCall:\n")
  print(object$call)
  cat("\nFrontier fit:\n")
  print(summary(object$frontier))
  cat("\nComposed residuals (y - frontier):\n")
  print(summary(object$residuals))
  invisible(object)
}

fitted.npsfareg <- function(object, ...) object$frontier

residuals.npsfareg <- function(object, ...) object$residuals

nobs.npsfareg <- function(object, ...) object$nobs


## Kernel weights, for the local-likelihood estimators below.
.npsfm_kw <- function(xm, j, bw) {
  z <- sweep(xm, 2L, xm[j, ], "-")
  z <- sweep(z, 2L, bw, "/")
  w <- exp(rowSums(dnorm(z, log = TRUE)))
  ## A degenerate weight vector (every point too far away) would make the local
  ## problem rank deficient; fall back to a flat weighting rather than fail.
  if (!any(is.finite(w)) || sum(w) <= .SFA_CONSTANTS$MIN_POSITIVE) {
    w <- rep(1, nrow(xm))
  }
  w
}

## Log density of the normal-half normal composed error, in (lambda, sigma).
## Constants that do not depend on the parameters are dropped.
.lnhn <- function(e, lambda, sigma) {
  -log(sigma) + pnorm(-e * lambda / sigma, log.p = TRUE) - 0.5 * e^2 / sigma^2
}


## PSZ / KPST: Park, Simar and Zelenyuk. Local maximum likelihood.
.psz_nll <- function(par, yv, wm, K, q) {
  r0 <- par[1]
  r1 <- par[2:q]
  a0 <- par[q + 1]
  a1 <- par[(q + 2):(2 * q)]
  b0 <- par[2 * q + 1]
  b1 <- par[(2 * q + 2):(3 * q)]

  ei <- yv - r0 - as.numeric(wm %*% r1)
  ls2v <- a0 + as.numeric(wm %*% a1)
  ls2u <- b0 + as.numeric(wm %*% b1)

  ## Keep the exponentials inside representable range
  ls2v <- pmin(pmax(ls2v, -.SFA_CONSTANTS$EXP_CLIP_UPPER), .SFA_CONSTANTS$EXP_CLIP_UPPER)
  ls2u <- pmin(pmax(ls2u, -.SFA_CONSTANTS$EXP_CLIP_UPPER), .SFA_CONSTANTS$EXP_CLIP_UPPER)

  s2 <- exp(ls2u) + exp(ls2v)
  if (any(!is.finite(s2)) || any(s2 <= 0)) {
    return(.SFA_CONSTANTS$MAX_VALUE)
  }
  zz <- -ei * exp(ls2u / 2 - ls2v / 2) / sqrt(s2)
  ll <- -log(s2) / 2 - ei^2 / (2 * s2) + pnorm(zz, log.p = TRUE)
  val <- sum(ll * K)
  if (!is.finite(val)) {
    return(.SFA_CONSTANTS$MAX_VALUE)
  }
  -val
}

.npsfm_psz <- function(xdat, ydat, bw.sel, bw, cost, maxit, verbose) {
  cost_sc <- if (cost) -1 else 1
  yv <- cost_sc * ydat
  xm <- as.matrix(xdat)
  n <- length(yv)
  q <- ncol(xm) + 1L

  ## Seed from an FLW fit: it supplies the local intercept/slope start and
  ## global variance components, which is the "prior fit" the paper assumes.
  seed <- .npsfm_flw(xdat, ydat, "hn", "ll", bw.sel, bw, cost)
  bwv <- as.numeric(seed$bws$bw)
  s2v0 <- log(max(seed$sigma_v^2, .SFA_CONSTANTS$MIN_POSITIVE))
  s2u0 <- log(max(seed$sigma_u^2, .SFA_CONSTANTS$MIN_POSITIVE))

  par_store <- matrix(NA_real_, n, 3L * q)
  conv <- integer(n)

  for (j in seq_len(n)) {
    if (verbose && j %% 25L == 0L) {
      message("PSZ: observation ", j, " of ", n)
    }
    K <- .npsfm_kw(xm, j, bwv)
    wm <- sweep(xm, 2L, xm[j, ], "-")
    start <- c(
      cost_sc * seed$frontier[j], cost_sc * seed$frontier_grad[j, ],
      s2v0, rep(0, q - 1L),
      s2u0, rep(0, q - 1L)
    )
    fitj <- tryCatch(
      minqa::bobyqa(start, .psz_nll,
        control = list(maxfun = maxit),
        yv = yv, wm = wm, K = K, q = q
      ),
      error = function(e) NULL
    )
    if (is.null(fitj)) {
      par_store[j, ] <- start
      conv[j] <- -1L
    } else {
      par_store[j, ] <- fitj$par
      conv[j] <- fitj$ierr
    }
  }

  sigma_v <- sqrt(exp(pmin(par_store[, q + 1L], .SFA_CONSTANTS$EXP_CLIP_UPPER)))
  sigma_u <- sqrt(exp(pmin(par_store[, 2L * q + 1L], .SFA_CONSTANTS$EXP_CLIP_UPPER)))
  mu <- sqrt(2 / pi) * sigma_u

  ## NO mean correction here.
  frontier <- cost_sc * par_store[, 1L]
  frontier_grad <- cost_sc * par_store[, 2:q, drop = FALSE]

  list(
    frontier = frontier, frontier_grad = frontier_grad,
    residuals = ydat - frontier,
    sigma_u = sigma_u, sigma_v = sigma_v, mean_correction = mu,
    conditional_mean = frontier - cost_sc * mu,
    sigma_v_grad = par_store[, (q + 2L):(2L * q), drop = FALSE],
    sigma_u_grad = par_store[, (2L * q + 2L):(3L * q), drop = FALSE],
    convergence = conv, bws = seed$bws
  )
}


## MY: Martins-Filho and Yao. Iterative local likelihood.
.my_local_nll <- function(par, yv, wm, K, lambda, sigma) {
  e <- yv - par[1] - as.numeric(wm %*% par[-1])
  val <- sum(.lnhn(e, lambda, sigma) * K)
  if (!is.finite(val)) {
    return(.SFA_CONSTANTS$MAX_VALUE)
  }
  -val
}

.my_global_nll <- function(par, e) {
  lambda <- par[1]
  sigma <- par[2]
  if (!is.finite(lambda) || !is.finite(sigma) || lambda <= 0 || sigma <= 0) {
    return(.SFA_CONSTANTS$MAX_VALUE)
  }
  val <- sum(.lnhn(e, lambda, sigma))
  if (!is.finite(val)) {
    return(.SFA_CONSTANTS$MAX_VALUE)
  }
  -val
}

.npsfm_my <- function(xdat, ydat, bw.sel, bw, cost, tol, iter, maxit, verbose) {
  cost_sc <- if (cost) -1 else 1
  yv <- cost_sc * ydat
  xm <- as.matrix(xdat)
  n <- length(yv)
  q <- ncol(xm) + 1L

  ## Step 1 of the paper: FLW supplies the starting frontier and scale.
  seed <- .npsfm_flw(xdat, ydat, "hn", "ll", bw.sel, bw, cost)
  bwv <- as.numeric(seed$bws$bw)
  theta <- c(
    max(seed$lambda, 1e-3),
    max(sqrt(seed$sigma_u^2 + seed$sigma_v^2), 1e-3)
  )
  ab <- cbind(cost_sc * seed$frontier, cost_sc * seed$frontier_grad)

  ## Precompute the kernel weights and centred designs once; they do not
  ## change across iterations.
  Klist <- vector("list", n)
  Wlist <- vector("list", n)
  for (j in seq_len(n)) {
    Klist[[j]] <- .npsfm_kw(xm, j, bwv)
    Wlist[[j]] <- sweep(xm, 2L, xm[j, ], "-")
  }

  step <- 1L
  dif <- Inf
  while (step <= iter && dif > tol) {
    for (j in seq_len(n)) {
      fitj <- tryCatch(
        stats::optim(ab[j, ], .my_local_nll,
          method = "BFGS", control = list(maxit = maxit),
          yv = yv, wm = Wlist[[j]], K = Klist[[j]],
          lambda = theta[1], sigma = theta[2]
        )$par,
        error = function(e) ab[j, ]
      )
      ab[j, ] <- fitj
    }
    ## Step 3: global update of (lambda, sigma) from the local intercepts
    e_opt <- yv - ab[, 1L]   ## composed errors eps = v - u, not centered
    new_theta <- tryCatch(
      stats::optim(theta, .my_global_nll,
        method = "L-BFGS-B",
        lower = c(1e-4, 1e-4), upper = c(50, 50), e = e_opt
      )$par,
      error = function(e) theta
    )
    dif <- sum((theta - new_theta)^2)
    theta <- new_theta
    if (verbose) {
      message("MY: iteration ", step, ", change ", signif(dif, 4))
    }
    step <- step + 1L
  }
  converged <- dif <= tol

  lambda <- theta[1]
  sigma <- theta[2]
  sigma_v <- sigma / sqrt(1 + lambda^2)
  sigma_u <- sigma_v * lambda
  mu <- sqrt(2 / pi) * sigma_u

  ## As in PSZ: the local objective is the composed-error likelihood, so the
  ## local intercept is the frontier itself and takes no mean correction.
  frontier <- cost_sc * ab[, 1L]
  frontier_grad <- cost_sc * ab[, 2:q, drop = FALSE]

  list(
    frontier = frontier, frontier_grad = frontier_grad,
    residuals = ydat - frontier,
    sigma_u = sigma_u, sigma_v = sigma_v, lambda = lambda, sigma = sigma,
    mean_correction = mu, conditional_mean = frontier - cost_sc * mu,
    iterations = step - 1L, tol_reached = dif, converged = converged,
    bws = seed$bws
  )
}


## Output-oriented DEA envelopment, solved one linear program per unit.
.dea_out <- function(X, Y, rts = c("vrs", "crs", "drs", "irs")) {
  rts <- match.arg(rts)
  if (!requireNamespace("lpSolve", quietly = TRUE)) {
    stop("method = \"SZ\" needs the 'lpSolve' package for its DEA step, ",
      "which is not installed. ",
      'Install it with install.packages("lpSolve").',
      call. = FALSE
    )
  }
  X <- as.matrix(X)
  Y <- as.matrix(Y)
  n <- nrow(X)
  k <- ncol(X)
  m <- ncol(Y)

  obj <- c(1, rep(0, n))
  vapply(seq_len(n), function(o) {
    con <- rbind(cbind(0, t(X)), cbind(-Y[o, ], t(Y)))
    dir <- c(rep("<=", k), rep(">=", m))
    rhs <- c(X[o, ], rep(0, m))
    if (rts != "crs") {
      con <- rbind(con, c(0, rep(1, n)))
      dir <- c(dir, switch(rts, vrs = "=", drs = "<=", irs = ">="))
      rhs <- c(rhs, 1)
    }
    sol <- lpSolve::lp("max", obj, con, dir, rhs)
    if (!identical(sol$status, 0L)) NA_real_ else sol$solution[1]
  }, numeric(1))
}


## SZ: Simar and Zelenyuk (2011). DEA monotonization of a prior fit.
.npsfm_sz <- function(xdat, ydat, bw.sel, bw, cost, rts, prior.fit, log.form) {
  if (cost) {
    stop("method = \"SZ\" is implemented for production frontiers only ",
      "(cost = FALSE).",
      call. = FALSE
    )
  }

  ## Stage 1: a smooth frontier to monotonize. Default to SVKZ, as the paper does.
  base <- NULL
  if (is.null(prior.fit)) {
    base <- .npsfm_svkz(xdat, ydat, bw.sel, cost = FALSE)
    prior.fit <- base$frontier
  }
  if (length(prior.fit) != length(ydat)) {
    stop("`prior.fit` must have one fitted value per observation (",
      length(ydat), "); got ", length(prior.fit), ".",
      call. = FALSE
    )
  }

  ## DEA operates on levels. The frontier is in logs when the model was fitted
  ## to logged data, which is the usual convention in this literature.
  Y <- if (log.form) exp(prior.fit) else prior.fit
  X <- if (log.form) exp(as.matrix(xdat)) else as.matrix(xdat)

  eff <- .dea_out(X = X, Y = matrix(Y, ncol = 1L), rts = rts)
  dea_front <- Y * eff

  frontier <- if (log.form) log(dea_front) else dea_front

  list(
    frontier = frontier,
    frontier_grad = if (!is.null(base)) base$frontier_grad else NULL,
    residuals = ydat - frontier,
    sigma_u = if (!is.null(base)) base$sigma_u else NA_real_,
    sigma_v = if (!is.null(base)) base$sigma_v else NA_real_,
    mean_correction = if (!is.null(base)) base$mean_correction else NA_real_,
    conditional_mean = if (!is.null(base)) base$conditional_mean else prior.fit,
    prior.fit = prior.fit, dea.efficiency = eff, rts = rts,
    bws = if (!is.null(base)) base$bws else NULL
  )
}
