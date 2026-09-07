## Two classical time-varying fixed-effects frontier estimators. Neither is
## maximum likelihood and neither assumes a distribution for u -- they buy
## time-varying inefficiency with a parametric restriction on how it moves
## instead. See notes/code_history/panel_classical.md.

## Shared setup: the response, the regressors with any intercept removed, and
## a firm/period index. The intercept goes because both estimators absorb it
## into the firm effect, exactly as psfm()'s SSFE does.
.classical_panel_setup <- function(formula_x, data, individual, time, y_var,
                                   model_name) {
  X <- stats::model.matrix(formula_x, data = data)
  icol <- match("(Intercept)", colnames(X))
  if (!is.na(icol)) X <- X[, -icol, drop = FALSE]
  if (!ncol(X)) {
    stop("psfm(model_name = \"", model_name, "\"): the frontier has no ",
      "regressors once the intercept is absorbed into the firm effect.",
      call. = FALSE
    )
  }
  y <- as.numeric(data[[y_var]])
  id <- as.factor(data[[individual]])
  tv <- if (is.null(time)) {
    stats::ave(seq_len(nrow(data)), id, FUN = seq_along)
  } else {
    z <- data[[time]]
    if (is.factor(z)) as.numeric(as.character(z)) else as.numeric(z)
  }
  list(X = X, y = y, id = id, t = tv, idx = split(seq_len(nrow(data)), id, drop = TRUE))
}

## Cornwell, Schmidt and Sickles (1990): a firm-specific QUADRATIC IN TIME,
## alpha_it = theta_i1 + theta_i2 t + theta_i3 t^2, swept out firm by firm.
## Every firm gets its own three-parameter trajectory, so inefficiency may
## cross over between firms -- the point of the estimator, and what
## time-invariant models like SSFE and PL80 rule out by assumption.
.fit_css90 <- function(setup, inefdec) {
  X <- setup$X
  y <- setup$y
  idx <- setup$idx
  n <- length(y)
  K <- ncol(X)

  ## Fewer than four periods and the firm's own quadratic is saturated: it
  ## contributes nothing to beta and its residuals are identically zero.
  short <- vapply(idx, length, integer(1)) < 4L
  if (all(short)) {
    stop("psfm(model_name = \"CSS\"): every firm has fewer than 4 periods. ",
      "The firm-specific quadratic has three parameters, so T_i >= 4 is ",
      "needed before a firm can contribute to beta. Use SSFE or PL80 for ",
      "short panels.",
      call. = FALSE
    )
  }
  if (any(short)) {
    warning("psfm(model_name = \"CSS\"): ", sum(short), " of ", length(idx),
      " firms have fewer than 4 periods and are absorbed exactly by their own ",
      "quadratic; they contribute nothing to beta and their fitted ",
      "inefficiency is zero by construction.",
      call. = FALSE
    )
  }

  ## Residualize y and X on each firm's own (1, t, t^2) rather than building
  ## 3N dummy columns.
  Wl <- lapply(idx, function(ii) cbind(1, setup$t[ii], setup$t[ii]^2))
  My <- numeric(n)
  MX <- matrix(0, n, K)
  for (k in seq_along(idx)) {
    ii <- idx[[k]]
    W <- Wl[[k]]
    Q <- qr(W)
    My[ii] <- qr.resid(Q, y[ii])
    MX[ii, ] <- qr.resid(Q, X[ii, , drop = FALSE])
  }

  ## A regressor entirely spanned by the firm quadratics -- a pure time trend,
  ## say -- leaves a zero column here, which is a rank failure, not a small
  ## coefficient.
  keep <- apply(MX, 2, function(z) max(abs(z)) > 1e-10)
  if (!all(keep)) {
    stop("psfm(model_name = \"CSS\"): regressor(s) ",
      paste(sQuote(colnames(X)[!keep]), collapse = ", "),
      " are spanned by the firm-specific quadratics in time and are not ",
      "separately identified. Drop them; their effect is absorbed into the ",
      "firm trajectories.",
      call. = FALSE
    )
  }

  fit <- stats::lm.fit(MX, My)
  beta <- fit$coefficients
  names(beta) <- colnames(X)
  ## 3 parameters per firm, plus K.
  df <- n - K - 3L * length(idx)
  s2 <- sum(fit$residuals^2) / max(df, 1L)
  XtX <- crossprod(MX)
  vc <- tryCatch(s2 * solve(XtX), error = function(e) NULL)
  se <- if (is.null(vc)) rep(NA_real_, K) else sqrt(pmax(diag(vc), 0))

  ## The firm trajectories, recovered from the residual of the frontier.
  r <- y - as.numeric(X %*% beta)
  alpha_it <- numeric(n)
  theta <- matrix(NA_real_, nrow = length(idx), ncol = 3L,
    dimnames = list(names(idx), c("const", "t", "t2"))
  )
  for (k in seq_along(idx)) {
    ii <- idx[[k]]
    th <- qr.coef(qr(Wl[[k]]), r[ii])
    th[!is.finite(th)] <- 0
    theta[k, ] <- th
    alpha_it[ii] <- as.numeric(Wl[[k]] %*% th)
  }

  list(beta = beta, se = se, df = df, alpha_it = alpha_it, theta = theta,
       sigma2 = s2, residuals = r - alpha_it)
}

## Lee and Schmidt (1993): one COMMON temporal pattern, scaled per firm --
## alpha_it = delta_t * alpha_i. Far more parsimonious than CSS (N + T - 1
## parameters rather than 3N), and it forbids crossover: every firm's
## inefficiency moves in the same shape, only the amplitude differs.
##
## alpha_i * delta_t is a rank-one matrix, so this is a rank-one factor model
## and is fitted by alternating least squares, which handles an unbalanced
## panel without special-casing.
.fit_ls93 <- function(setup, inefdec, maxit = 500L, tol = 1e-10) {
  X <- setup$X
  y <- setup$y
  id <- setup$id
  n <- length(y)
  K <- ncol(X)
  tf <- as.factor(setup$t)
  Tn <- nlevels(tf)
  if (Tn < 2L) {
    stop("psfm(model_name = \"LS\"): needs at least two distinct periods.",
      call. = FALSE
    )
  }

  ## Start from the time-invariant solution, delta_t == 1, which is exactly
  ## Schmidt and Sickles (1984); the iteration then relaxes it.
  delta <- rep(1, Tn)
  names(delta) <- levels(tf)
  b <- stats::lm.fit(X, y)$coefficients
  b[!is.finite(b)] <- 0
  r <- y - as.numeric(X %*% b)
  alpha <- as.numeric(tapply(r, id, mean))
  names(alpha) <- levels(id)

  ssr_old <- Inf
  it <- 0L
  for (it in seq_len(maxit)) {
    d_i <- delta[as.integer(tf)]
    a_i <- alpha[as.integer(id)]

    ## beta, given the current factor.
    fb <- stats::lm.fit(X, y - d_i * a_i)
    b <- fb$coefficients
    b[!is.finite(b)] <- 0
    r <- y - as.numeric(X %*% b)

    ## alpha, given delta: a per-firm regression of the residual on delta_t.
    num <- as.numeric(tapply(r * d_i, id, sum))
    den <- as.numeric(tapply(d_i^2, id, sum))
    alpha <- num / pmax(den, .SFA_CONSTANTS$MIN_POSITIVE)
    a_i <- alpha[as.integer(id)]

    ## delta, given alpha: a per-period regression of the residual on alpha_i.
    numd <- as.numeric(tapply(r * a_i, tf, sum))
    dend <- as.numeric(tapply(a_i^2, tf, sum))
    delta <- numd / pmax(dend, .SFA_CONSTANTS$MIN_POSITIVE)

    ## alpha_i*delta_t is unchanged by rescaling one against the other, so fix
    ## the scale. delta_1 = 1 is Lee and Schmidt's own normalization; fall back
    ## to a norm constraint if the first period's effect is near zero, where
    ## dividing by it would be arbitrarily badly conditioned.
    sc <- delta[1L]
    if (!is.finite(sc) || abs(sc) < 1e-8) sc <- sqrt(sum(delta^2) / Tn)
    if (is.finite(sc) && abs(sc) > 1e-12) {
      delta <- delta / sc
      alpha <- alpha * sc
    }
    names(delta) <- levels(tf)
    names(alpha) <- levels(id)

    ssr <- sum((r - delta[as.integer(tf)] * alpha[as.integer(id)])^2)
    if (is.finite(ssr) && abs(ssr_old - ssr) <= tol * (abs(ssr) + tol)) {
      ssr_old <- ssr
      break
    }
    ssr_old <- ssr
  }
  converged <- it < maxit
  if (!converged) {
    warning("psfm(model_name = \"LS\"): the alternating least-squares ",
      "iteration hit maxit = ", maxit, " without meeting its tolerance. ",
      "Treat delta and the inefficiency scores as provisional.",
      call. = FALSE
    )
  }

  alpha_it <- delta[as.integer(tf)] * alpha[as.integer(id)]
  resid <- y - as.numeric(X %*% b) - alpha_it
  ## N + T - 1 free factor parameters after the scale normalization.
  df <- n - K - (nlevels(id) + Tn - 1L)
  s2 <- sum(resid^2) / max(df, 1L)
  vc <- tryCatch(s2 * solve(crossprod(X)), error = function(e) NULL)
  se <- if (is.null(vc)) rep(NA_real_, K) else sqrt(pmax(diag(vc), 0))
  names(b) <- colnames(X)

  list(beta = b, se = se, df = df, alpha_it = as.numeric(alpha_it),
       alpha = alpha, delta = delta, sigma2 = s2, residuals = resid,
       iterations = it, converged = converged, ssr = ssr_old)
}

## Both estimators identify inefficiency only RELATIVE to the best firm in
## each period, so the frontier is renormalized period by period. That is what
## makes them time-varying: with a time-invariant effect the comparison would
## be the same in every period.
.classical_panel_scores <- function(alpha_it, tvec, inefdec) {
  best <- if (isTRUE(inefdec)) {
    stats::ave(alpha_it, tvec, FUN = max)
  } else {
    stats::ave(alpha_it, tvec, FUN = min)
  }
  u_hat <- if (isTRUE(inefdec)) best - alpha_it else alpha_it - best
  u_hat <- pmax(u_hat, 0)
  list(u_hat = u_hat, exp_u_hat = exp(-u_hat), frontier = best)
}

## Schmidt and Sickles (1984), the random-effects side. `SSFE` uses their
## within/LSDV estimator; these two complete the family with the GLS estimator
## and its correlated-random-effects correction. All three read inefficiency
## off the firm effect as distance from the best firm, and none assumes a
## distribution for u.
##
## The trade is the classic one. RE is more efficient than FE and identifies
## coefficients on TIME-INVARIANT regressors, which the within estimator sweeps
## away -- but only if the effects are uncorrelated with the regressors, and in
## a frontier model that assumption says inefficiency is unrelated to input
## choice, which is exactly what a manager would not do. CRE is the answer:
## Mundlak (1978) adds the within-firm means of the time-varying regressors, so
## the correlation is modelled rather than assumed away, and the remaining
## slopes equal the within estimates while time-invariant regressors stay
## identified.
.fit_ss_re <- function(formula_x, data, individual, model_name, inefdec,
                       x_vars_vec) {
  cre <- identical(model_name, "SSCRE")
  form <- formula_x
  mund <- character(0)

  if (cre) {
    ## Group means of the time-varying regressors. A regressor already constant
    ## within firm has a mean identical to itself, so adding it would be exactly
    ## collinear -- those are precisely the ones RE exists to identify, and they
    ## are left alone.
    X <- stats::model.matrix(formula_x, data = data)
    icol <- match("(Intercept)", colnames(X))
    if (!is.na(icol)) X <- X[, -icol, drop = FALSE]
    id <- as.factor(data[[individual]])
    varies <- vapply(seq_len(ncol(X)), function(j) {
      any(tapply(X[, j], id, function(z) stats::var(z) > 1e-12), na.rm = TRUE)
    }, logical(1))
    if (!any(varies)) {
      stop("psfm(model_name = \"SSCRE\"): no regressor varies within firm, so ",
        "there is nothing for the Mundlak means to correct. Use ",
        "model_name = \"SSRE\".",
        call. = FALSE
      )
    }
    Xm <- X[, varies, drop = FALSE]
    mund <- paste0(".mean_", make.names(colnames(Xm)))
    gm <- apply(Xm, 2, function(z) as.numeric(stats::ave(z, id, FUN = mean)))
    gm <- matrix(gm, nrow = nrow(X), dimnames = list(NULL, mund))
    data <- cbind(data, as.data.frame(gm))
    form <- stats::update(stats::formula(formula_x),
      stats::reformulate(c(".", mund))
    )
  }

  fit <- plm::plm(form, data,
    effect = "individual", model = "random", index = individual
  )
  cf <- stats::coef(fit)
  se <- summary(fit)$coefficients[, "Std. Error"]

  ## The BLUPs of the firm effects. Under CRE these are the effects NET of the
  ## part the Mundlak means explain, which is the point: what is left is the
  ## component not predictable from input use.
  alpha_hat <- plm::ranef(fit)
  alpha_hat <- stats::setNames(as.numeric(alpha_hat), names(alpha_hat))

  u_i <- if (isTRUE(inefdec)) {
    max(alpha_hat) - alpha_hat
  } else {
    alpha_hat - min(alpha_hat)
  }
  ## One value per firm, expanded to one per observation so the return shape
  ## matches CSS/LS and every other panel model here.
  key <- as.character(data[[individual]])
  u_hat <- unname(u_i[match(key, names(u_i))])

  list(beta = cf, se = se, alpha_i = alpha_hat,
       u_hat = u_hat, exp_u_hat = exp(-u_hat),
       mundlak = mund, theta = plm::ercomp(fit)$theta,
       ercomp = plm::ercomp(fit), plm_fit = fit)
}
