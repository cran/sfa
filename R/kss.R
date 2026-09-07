## Kneip, Sickles and Song (2012): the firm effect is a smooth function of
## time lying in an L-dimensional space whose basis is estimated FROM THE DATA
## rather than assumed. See notes/code_history/kss.md.
##
## This nests the two estimators in panel_classical.R. Writing the effect as
## alpha_it = sum_r theta_ir g_r(t):
##
##   SSFE  L = 1 with g_1(t) constant           -- effects fixed over time
##   LS    L = 1 with g_1 free                  -- one common shape, free
##   CSS   L = 3 with g fixed to {1, t, t^2}    -- basis assumed, not estimated
##   KSS   L and g both estimated               -- neither assumed
##
## What is bought is that no functional form is imposed on the temporal
## pattern. What is paid is two tuning choices (the rank L and the smoothing)
## and a balanced panel.

## Discrete penalized smoother over the time grid. A second-difference penalty
## rather than a full B-spline basis: T is typically 6-20 in a productivity
## panel, and at that length the two coincide for practical purposes while this
## one is exact, needs no knot placement, and stays defined at T = 3.
.kss_smoother <- function(Tn, kappa) {
  if (Tn < 3L || kappa <= 0) {
    return(diag(Tn))
  }
  D <- diff(diag(Tn), differences = 2L)
  solve(diag(Tn) + kappa * crossprod(D))
}

## Generalized cross-validation for the smoothing parameter, pooled over firms.
## The curves share a smoother, so the GCV score sums over them.
.kss_gcv <- function(R, kappa_grid) {
  Tn <- ncol(R)
  n <- nrow(R)
  sc <- vapply(kappa_grid, function(k) {
    S <- .kss_smoother(Tn, k)
    tr <- sum(diag(S))
    if (Tn - tr <= .Machine$double.eps) {
      return(Inf)
    }
    rss <- sum((R - R %*% t(S))^2) / n
    rss / (1 - tr / Tn)^2
  }, numeric(1))
  if (all(!is.finite(sc))) {
    return(0)
  }
  kappa_grid[which.min(sc)]
}

## Bai and Ng (2002) IC_p2 for the number of factors, evaluated on the raw
## residual matrix against a basis estimated from the smoothed one. Not KSS's own selection
## rule -- theirs is a threshold on the eigenvalue sequence -- but it is the
## standard criterion for exactly this problem, and the book that motivates
## this entry cites Bai and Ng (2002) and Bai (2009) as the comparable
## estimators. `kss_L` overrides it.
.kss_select_L <- function(R, G_full, L_max) {
  n <- nrow(R)
  Tn <- ncol(R)
  nt <- n * Tn
  pen <- ((n + Tn) / nt) * log(min(n, Tn))
  ic <- vapply(seq_len(L_max), function(L) {
    G <- G_full[, seq_len(L), drop = FALSE]
    V <- sum((R - R %*% G %*% t(G))^2) / nt
    if (!is.finite(V) || V <= 0) Inf else log(V) + L * pen
  }, numeric(1))
  if (all(!is.finite(ic))) 1L else which.min(ic)
}

## Fit the KSS model. Follows the book's equations 74-76: cross-sectional
## centering by period first, which imposes sum_i u_i(t) = 0 and removes the
## common function beta_0(t), then the factor structure on what is left.
.fit_kss <- function(y, X, id, tv, kss_L, kss_smooth, kss_L_max,
                     maxit = 200L, tol = 1e-10) {
  idf <- as.factor(id)
  tf <- as.factor(tv)
  n <- nlevels(idf)
  Tn <- nlevels(tf)
  K <- ncol(X)

  ## KSS is defined on a balanced panel, as is the reference Stata
  ## implementation. Silently unbalancing it would change the estimator, so
  ## this refuses instead.
  if (length(y) != n * Tn || any(table(idf, tf) != 1L)) {
    stop("psfm(model_name = \"KSS\"): needs a BALANCED panel -- every firm ",
      "observed in every period. Found ", length(y), " rows against ",
      n, " firms x ", Tn, " periods. Use CSS or LS, which do not require it.",
      call. = FALSE
    )
  }
  if (Tn < 3L) {
    stop("psfm(model_name = \"KSS\"): needs at least 3 periods to estimate a ",
      "temporal basis. Use SSFE or SSRE.",
      call. = FALSE
    )
  }
  ## Two caps, because they answer different questions. L_hard is what the
  ## algebra permits; L_auto_cap is what the SELECTION CRITERION can be trusted
  ## with. Bai-Ng needs T large relative to L: allowing the criterion up to
  ## T - 1 lets the factors span nearly the whole time dimension, V(L)
  ## collapses toward zero, and no penalty can compete. Measured at n = 100
  ## over 5 seeds and two true ranks: with the cap at T - 1 the criterion
  ## returned the cap on 10 of 10 designs at T = 6 and again at T = 8; with
  ## floor(T/2) it recovers the true rank exactly from T = 8 upward. An
  ## explicit `kss_L` is the user's own call and is allowed up to L_hard.
  L_hard <- max(1L, min(as.integer(kss_L_max), Tn - 1L, n - 1L))
  L_auto_cap <- max(1L, min(L_hard, Tn %/% 2L))

  ## Equation 74: centre by PERIOD, across firms.
  ymat <- matrix(0, n, Tn)
  ymat[cbind(as.integer(idf), as.integer(tf))] <- y
  Xa <- array(0, c(n, Tn, K))
  for (k in seq_len(K)) {
    Xa[cbind(as.integer(idf), as.integer(tf), k)] <- X[, k]
  }
  yc <- ymat - matrix(colMeans(ymat), n, Tn, byrow = TRUE)
  Xc <- Xa
  for (k in seq_len(K)) {
    Xc[, , k] <- Xa[, , k] - matrix(colMeans(Xa[, , k]), n, Tn, byrow = TRUE)
  }
  Xflat <- matrix(0, n * Tn, K, dimnames = list(NULL, colnames(X)))
  for (k in seq_len(K)) Xflat[, k] <- as.numeric(Xc[, , k])
  yflat <- as.numeric(yc)

  ## Step 1: an initial beta, then alternate.
  b <- stats::lm.fit(Xflat, yflat)$coefficients
  b[!is.finite(b)] <- 0
  kappa_grid <- c(0, 10^seq(-3, 2, length.out = 40))
  kappa <- if (identical(kss_smooth, "auto")) NA_real_ else as.numeric(kss_smooth)
  L <- if (identical(kss_L, "auto")) NA_integer_ else as.integer(kss_L)
  alpha <- matrix(0, n, Tn)
  ssr_old <- Inf
  it <- 0L

  for (it in seq_len(maxit)) {
    R <- yc - matrix(Xflat %*% b, n, Tn)

    ## Step 2: smooth the firm curves, then take the eigenfunctions of their
    ## empirical covariance. Smoothing FIRST is what makes the estimated basis
    ## a set of smooth functions rather than T arbitrary period effects -- it
    ## is the whole difference between this and an unrestricted factor model.
    if (identical(kss_smooth, "auto")) kappa <- .kss_gcv(R, kappa_grid)
    S <- .kss_smoother(Tn, kappa)
    Rs <- R %*% t(S)

    ev <- eigen(crossprod(Rs) / n, symmetric = TRUE)
    G_full <- ev$vectors[, seq_len(L_hard), drop = FALSE]

    ## Score the criterion on the RAW residuals, not the smoothed ones, even
    ## though the basis comes from the smoothed. Bai-Ng works because V(L)
    ## flattens out at the noise variance once the real factors are in; the
    ## smoother has already removed most of that noise, so scoring on Rs makes
    ## every extra factor keep paying and the rank runs to L_max. Measured on a
    ## 2-factor design with eigenvalues 3.23, 0.50, 0.045, ...: scoring on Rs
    ## selected 7, scoring on R selects 2.
    if (identical(kss_L, "auto")) L <- .kss_select_L(R, G_full, L_auto_cap)
    L <- max(1L, min(as.integer(L), L_hard))
    G <- G_full[, seq_len(L), drop = FALSE]

    ## Step 2b/3: loadings by least squares on the estimated basis, then the
    ## effects rebuilt from them. G is orthonormal, so the projection is the
    ## least-squares fit.
    theta <- R %*% G
    alpha <- theta %*% t(G)

    ## beta given the effects.
    fb <- stats::lm.fit(Xflat, yflat - as.numeric(alpha))
    b <- fb$coefficients
    b[!is.finite(b)] <- 0

    ssr <- sum((yflat - Xflat %*% b - as.numeric(alpha))^2)
    if (is.finite(ssr) && abs(ssr_old - ssr) <= tol * (abs(ssr) + tol)) {
      ssr_old <- ssr
      break
    }
    ssr_old <- ssr
  }
  ## A selection that lands ON the cap did not choose freely, so say so rather
  ## than reporting it as though the criterion settled on it.
  if (identical(kss_L, "auto") && L >= L_auto_cap && L_auto_cap < L_hard) {
    warning("psfm(model_name = \"KSS\"): the dimension criterion selected ",
      L, ", which is the largest it is allowed to consider at T = ", Tn,
      " (floor(T/2)). The factor space may be richer than this fit reports. ",
      "Inspect $kss$eigenvalues, and set kss_L explicitly if a larger ",
      "dimension is warranted; a longer panel is the real remedy.",
      call. = FALSE
    )
  }

  converged <- it < maxit
  if (!converged) {
    warning("psfm(model_name = \"KSS\"): the iteration hit maxit = ", maxit,
      " without meeting its tolerance. Treat the basis and the scores as ",
      "provisional.",
      call. = FALSE
    )
  }
  names(b) <- colnames(X)

  resid <- yflat - as.numeric(Xflat %*% b) - as.numeric(alpha)
  ## K slopes, plus L loadings per firm and L basis functions of length T,
  ## less L^2 for the orthonormality normalization.
  df <- n * Tn - K - (n * L + Tn * L - L^2)
  s2 <- sum(resid^2) / max(df, 1L)
  vc <- tryCatch(s2 * solve(crossprod(Xflat)), error = function(e) NULL)
  se <- if (is.null(vc)) rep(NA_real_, K) else sqrt(pmax(diag(vc), 0))

  ## Back to the original row order.
  ## Everything above works on the (firm x period) matrix; `ord` maps it back
  ## to the order the rows arrived in.
  ord <- cbind(as.integer(idf), as.integer(tf))
  list(beta = b, se = se, df = df, sigma2 = s2,
       alpha_it = alpha[ord],
       residuals = matrix(resid, n, Tn)[ord],
       basis = G, loadings = theta, L = L, kappa = kappa,
       eigenvalues = ev$values, iterations = it, converged = converged)
}
