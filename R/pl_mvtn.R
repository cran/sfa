## Pitt and Lee (1981, JoE) Model III: the multivariate truncated normal panel
## likelihood from their Appendix 2.

## Pieces of the equicorrelated system. Returns NULL on an inadmissible
## parameter vector, which every caller treats as "return the penalty".
.pl_mvtn_parts <- function(sig_u, rho, sig_v, BigT) {
  if (any(!is.finite(c(sig_u, rho, sig_v))) || sig_u <= 0 || sig_v <= 0) {
    return(NULL)
  }
  ## rho must keep Sigma positive definite: -1/(T-1) < rho < 1
  if (rho <= -1 / (BigT - 1) + 1e-8 || rho >= 1 - 1e-8) {
    return(NULL)
  }

  a <- sig_u^2 * (1 - rho) ## Sigma = a I + b J
  b <- sig_u^2 * rho
  if (a <= 0 || a + BigT * b <= 0) {
    return(NULL)
  }
  logdet_Sig <- (BigT - 1) * log(a) + log(a + BigT * b)

  ## Sigma^{-1} = (1/a) I - cc J
  cc <- b / (a * (a + BigT * b))

  ## Q^{-1} = Sigma^{-1} + I/sigma_v^2 = a2 I - cc J
  a2 <- 1 / a + 1 / sig_v^2
  den <- a2 - BigT * cc
  if (!is.finite(a2) || a2 <= 0 || !is.finite(den) || den <= 0) {
    return(NULL)
  }
  logdet_Q <- -((BigT - 1) * log(a2) + log(den))

  ## Q = (1/a2) I + qq J
  qq <- cc / (a2 * den)
  Q <- diag(1 / a2, BigT) + qq
  if (any(!is.finite(Q))) {
    return(NULL)
  }

  list(Q = Q, logdet_Sig = logdet_Sig, logdet_Q = logdet_Q,
       a = a, b = b, Sigma = diag(a, BigT) + b)
}

## Orthant probability Pr(w <= 0), w ~ N(mean, varcov), on the log scale.
.pl_mvtn_log_orthant <- function(mean, varcov) {
  p <- tryCatch(
    mnormt::sadmvn(lower = rep(-Inf, length(mean)), upper = rep(0, length(mean)),
                   mean = mean, varcov = varcov),
    error = function(e) NA_real_
  )
  if (!is.finite(p) || p <= 0) {
    return(NA_real_)
  }
  log(p)
}

## Negative summed log-likelihood. `par` is (sig_v, sig_u, rho, beta...);
## eps_list holds one T-vector of y - X beta per firm, rebuilt each call.
.pl_mvtn_nll <- function(par, Y_list, X_list, BigT, inefdec_n = 1) {
  sig_v <- par[1]
  sig_u <- par[2]
  rho <- par[3]
  beta <- par[-(1:3)]

  pr <- .pl_mvtn_parts(sig_u, rho, sig_v, BigT)
  if (is.null(pr)) {
    return(.SFA_CONSTANTS$MAX_VALUE)
  }
  logP0 <- .pl_mvtn_log_orthant(rep(0, BigT), pr$Sigma)
  if (!is.finite(logP0)) {
    return(.SFA_CONSTANTS$MAX_VALUE)
  }

  Qi <- pr$Q
  s2 <- sig_v^2
  const <- -(BigT / 2) * log(2 * pi) - BigT * log(sig_v) -
    0.5 * pr$logdet_Sig + 0.5 * pr$logdet_Q - logP0

  tot <- 0
  for (i in seq_along(Y_list)) {
    ## u <= 0 is PL's sign convention, i.e. inefficiency SUBTRACTS from a
    ## production frontier. inefdec_n flips it for a cost frontier.
    eps <- inefdec_n * (Y_list[[i]] - as.numeric(X_list[[i]] %*% beta))
    Qe <- Qi %*% eps
    mu <- as.numeric(Qe) / s2
    quad <- sum(eps^2) / s2 - as.numeric(crossprod(eps, Qe)) / s2^2
    lP <- .pl_mvtn_log_orthant(mu, Qi)
    if (!is.finite(lP)) {
      return(.SFA_CONSTANTS$MAX_VALUE)
    }
    tot <- tot + const - 0.5 * quad + lP
  }
  if (!is.finite(tot)) {
    return(.SFA_CONSTANTS$MAX_VALUE)
  }
  -tot
}


## Self-contained fitting routine, called from psfm() before the shared
## start_panel() machinery.
.psfm_pl_mvtn <- function(data, y_var, x_vars_vec, individual, inefdec_n,
                          maxit.optim, Method, optHessian, start_val,
                          verbose, call, formula) {
  id_chr <- as.character(data[, individual])
  indiv <- unique(id_chr)
  gid <- match(id_chr, indiv)
  t_i <- as.numeric(table(gid))
  if (length(unique(t_i)) != 1L) {
    stop("psfm(model_name = \"PL80_MVTN\") requires a BALANCED panel: Sigma is ",
      "a single T x T matrix shared by every firm. Found T ranging from ",
      min(t_i), " to ", max(t_i), ". Use model_name = \"PL80\" for the ",
      "time-invariant model, which handles unbalanced panels.",
      call. = FALSE
    )
  }
  BigT <- t_i[1]
  if (BigT < 2L) {
    stop("psfm(model_name = \"PL80_MVTN\") needs at least two periods per firm; ",
      "with T = 1 there is no cross-period dependence for Sigma to describe.",
      call. = FALSE
    )
  }

  ## "(Intercept)" is a NAME in x_vars_vec, not a column of `data` -- the same
  ## thing the GTRE_FML branch works around by setdiff()-ing it out.
  xn <- x_vars_vec
  xcols <- setdiff(xn, "(Intercept)")
  Xall <- as.matrix(as.data.frame(data)[, xcols, drop = FALSE])
  if ("(Intercept)" %in% xn) {
    Xall <- cbind(`(Intercept)` = 1, Xall)
    xn <- c("(Intercept)", xcols)
  } else {
    xn <- xcols
  }
  Yall <- as.numeric(as.data.frame(data)[, y_var])
  Y_list <- split(Yall, gid)
  X_list <- lapply(split(seq_along(gid), gid), function(ix) Xall[ix, , drop = FALSE])

  ## Start values: pooled OLS for beta, and the usual moment split of the
  ## residual variance between the two components.
  ols <- stats::lm.fit(Xall, Yall)
  b0 <- as.numeric(ols$coefficients)
  b0[!is.finite(b0)] <- 0
  r <- as.numeric(ols$residuals)
  s2 <- stats::var(r)
  start_v <- if (isTRUE(is.numeric(start_val)) && length(start_val) == length(b0) + 3L) {
    as.numeric(start_val)
  } else {
    c(sqrt(max(0.5 * s2, 1e-3)), sqrt(max(0.5 * s2, 1e-3)), 0.3, b0)
  }

  Start.Time <- start.time()
  np <- length(start_v)
  lower <- c(1e-3, 1e-3, -1 / (BigT - 1) + 0.02, rep(-Inf, np - 3))
  upper <- c(Inf, Inf, 0.95, rep(Inf, np - 3))

  opt <- stats::optim(
    par = start_v, fn = .pl_mvtn_nll, method = Method,
    lower = lower, upper = upper, hessian = isTRUE(optHessian),
    control = list(maxit = maxit.optim, trace = if (isTRUE(verbose)) 1 else 0),
    Y_list = Y_list, X_list = X_list, BigT = BigT, inefdec_n = inefdec_n
  )
  total_time <- start.time() - Start.Time

  par_names <- c("sigv", "sigu", "rho", xn)
  out <- matrix(NA_real_, 3, np, dimnames = list(c("par", "st_err", "t-val"), par_names))
  out[1, ] <- opt$par
  se <- rep(NA_real_, np)
  if (isTRUE(optHessian) && !is.null(opt$hessian)) {
    V <- tryCatch(solve(opt$hessian), error = function(e) NULL)
    if (!is.null(V) && all(is.finite(diag(V))) && all(diag(V) >= 0)) {
      se <- sqrt(diag(V))
    }
  }
  out[2, ] <- se
  out[3, ] <- out[1, ] / se

  ## E[u_it | eps_i] has no closed form here -- u is truncated multivariate
  ## normal and the conditional mean is another T-dimensional integral.
  pr <- .pl_mvtn_parts(opt$par[2], opt$par[3], opt$par[1], BigT)
  u_hat <- rep(NA_real_, length(Yall))
  if (!is.null(pr)) {
    s2v <- opt$par[1]^2
    for (i in seq_along(Y_list)) {
      eps <- inefdec_n * (Y_list[[i]] - as.numeric(X_list[[i]] %*% opt$par[-(1:3)]))
      mu <- as.numeric(pr$Q %*% eps) / s2v
      u_hat[gid == i] <- pmax(-pmin(mu, 0), 0)
    }
  }

  results <- list(
    out = out, opt = opt, total_time = total_time, start_v = start_v,
    model_name = "PL80_MVTN", formula = formula,
    coefficients = out[1, ], std.errors = out[2, ], t.values = out[3, ],
    u_hat = u_hat, exp_u_hat = pmin(pmax(exp(-u_hat), 0), 1),
    BigT = BigT, call = call, data = data
  )
  class(results) <- "sfareg"
  results
}
