## Amsler, Prokhorov and Schmidt (2016) JoE 190:280-288, and (2017) JoE
## 199:131-140. Endogenous regressors in a stochastic frontier: correlation
## between the regressors and the NOISE v. See notes/code_history/ivsfm.md for
## the derivation, the parameterization that keeps sigma_c positive, and why
## the two papers share one likelihood.
ivsfm <- function(formula,
                  endogenous,
                  instruments,
                  data,
                  model_name = c("IVLIML", "IVCF", "C2SLS"),
                  uhet = NULL,
                  inefdec = TRUE,
                  maxit.bobyqa = 10000,
                  maxit.psoptim = 1000,
                  maxit.optim = 1000,
                  start_val = FALSE,
                  PSopt = FALSE,
                  optHessian = TRUE,
                  Method = "L-BFGS-B",
                  verbose = FALSE,
                  rand.psoptim = NULL) {
  call <- match.call()
  model_name <- .match_model_name(model_name, eval(formals()$model_name))
  Start.Time <- Sys.time()

  for (nm in c("formula", "endogenous", "instruments")) {
    if (eval(call("missing", as.name(nm))) || !inherits(get(nm), "formula")) {
      stop("ivsfm(): `", nm, "` must be a formula.", call. = FALSE)
    }
  }
  if (length(formula) != 3L) {
    stop("ivsfm(): `formula` must have a left-hand side, e.g. y ~ x1 + x2.",
      call. = FALSE
    )
  }
  if (missing(data) || !is.data.frame(data)) {
    stop("ivsfm(): `data` must be a data.frame.", call. = FALSE)
  }
  if (length(inefdec) != 1L || !is.logical(inefdec) || is.na(inefdec)) {
    stop("ivsfm(): `inefdec` must be TRUE or FALSE.", call. = FALSE)
  }
  ## Every other entry point reads `|` as a variance segment. Here the equations
  ## are separate arguments, so a pipe is an error rather than a silent
  ## reinterpretation of the user's formula.
  for (nm in c("formula", "endogenous", "instruments")) {
    if (any(grepl("|", deparse(get(nm)), fixed = TRUE))) {
      stop("ivsfm(): `", nm, "` must not contain a `|` segment. The frontier, ",
        "the endogenous regressors and the instruments are separate arguments; ",
        "environmental variables go in `uhet`.",
        call. = FALSE
      )
    }
  }

  S <- if (isTRUE(inefdec)) 1 else -1
  cz <- .SFA_CONSTANTS

  ## ---- design matrices ----------------------------------------------------
  mf <- stats::model.frame(formula, data = data, na.action = stats::na.omit)
  y <- as.numeric(stats::model.response(mf))
  X <- stats::model.matrix(stats::terms(mf), mf)
  keep <- as.integer(rownames(mf))
  dsub <- data[keep, , drop = FALSE]

  endog_v <- all.vars(endogenous)
  inst_v <- all.vars(instruments)
  q_v <- if (is.null(uhet)) character(0) else all.vars(uhet)

  miss <- setdiff(c(endog_v, inst_v, q_v), names(data))
  if (length(miss)) {
    stop("ivsfm(): variable(s) not found in `data`: ",
      paste(miss, collapse = ", "), ".",
      call. = FALSE
    )
  }
  bad <- setdiff(endog_v, c(colnames(X), q_v))
  if (length(bad)) {
    stop("ivsfm(): endogenous variable(s) ", paste(bad, collapse = ", "),
      " appear in neither the frontier nor `uhet`. An endogenous variable must ",
      "be a regressor somewhere in the model.",
      call. = FALSE
    )
  }
  if (any(inst_v %in% endog_v)) {
    stop("ivsfm(): ", paste(intersect(inst_v, endog_v), collapse = ", "),
      " listed as both endogenous and an instrument.",
      call. = FALSE
    )
  }

  ## P: the endogenous block. Both papers put the endogenous regressors of the
  ## frontier AND of the variance function into one vector p = (x2, q2), which
  ## is exactly why one likelihood serves both.
  P <- as.matrix(dsub[, endog_v, drop = FALSE])
  m <- ncol(P)

  ## Q: environmental variables, giving sigma_u,i = sigma_u * exp(q_i'delta).
  ## Without them this is the 2016 model; with them, the 2017 one.
  Q <- if (length(q_v)) {
    stats::model.matrix(stats::reformulate(q_v), dsub)[, -1L, drop = FALSE]
  } else {
    NULL
  }
  n_q <- if (is.null(Q)) 0L else ncol(Q)

  ## Z: instruments = the INCLUDED exogenous regressors plus the EXCLUDED ones.
  X_exog <- X[, setdiff(colnames(X), endog_v), drop = FALSE]
  W <- stats::model.matrix(stats::reformulate(inst_v), dsub)[, -1L, drop = FALSE]
  Z <- cbind(X_exog, W)
  Z <- Z[, !duplicated(t(Z)), drop = FALSE]
  L <- ncol(Z)
  n <- length(y)

  if (ncol(W) < m) {
    stop("ivsfm(): ", ncol(W), " excluded instrument(s) for ", m,
      " endogenous variable(s). The order condition needs at least as many ",
      "excluded instruments as endogenous variables.",
      call. = FALSE
    )
  }
  if (n <= ncol(X) + L * m + m + 2L) {
    stop("ivsfm(): ", n, " observations cannot identify this specification.",
      call. = FALSE
    )
  }

  ## ---- reduced form, and 2SLS --------------------------------------------
  ## .safe_inverse() returns a list; $value is the inverse.
  ZtZi <- .safe_inverse(crossprod(Z))$value
  Pi_ols <- ZtZi %*% crossprod(Z, P)          # L x m
  Xi_ols <- P - Z %*% Pi_ols                  # reduced-form residuals
  Sxx_ols <- crossprod(Xi_ols) / n

  Xhat <- Z %*% (ZtZi %*% crossprod(Z, X))
  b_2sls <- as.numeric(.safe_inverse(crossprod(Xhat, X))$value %*% crossprod(Xhat, y))

  ## APS (2016) Eq (11): the residual uses the ACTUAL endogenous regressors,
  ## not their fitted values. Using x-hat here is the classic mistake and
  ## produces the wrong second and third moments.
  e_2sls <- as.numeric(y - X %*% b_2sls)

  ## Olson's moment estimators. k3 is negative (1 - 4/pi = -0.273), so the
  ## wrong-skew case is S * m3 >= 0.
  m2 <- mean(e_2sls^2)
  m3 <- mean(e_2sls^3)
  k3 <- sqrt(2 / pi) * (1 - 4 / pi)
  wrong_skew <- !(S * m3 < 0)
  su0 <- if (!wrong_skew) (S * m3 / k3)^(1 / 3) else 0.5 * stats::sd(e_2sls)
  su0 <- if (is.finite(su0)) max(su0, 1e-3) else max(0.5 * stats::sd(e_2sls), 1e-3)
  sv0 <- sqrt(max(m2 - (1 - 2 / pi) * su0^2, 1e-6))

  has_int <- attr(stats::terms(mf), "intercept") == 1L

  ## ---- C2SLS: 2SLS, then correct the intercept ---------------------------
  if (model_name == "C2SLS") {
    if (wrong_skew) {
      warning("ivsfm(): the 2SLS residuals are skewed the wrong way for a ",
        model_name, " fit, so the moment estimator of sigma_u has no real ",
        "solution and a fallback was used. Treat sigma_u and the corrected ",
        "intercept as uninformative.",
        call. = FALSE
      )
    }
    par <- b_2sls
    if (has_int) par[1L] <- par[1L] + S * su0 * sqrt(2 / pi)
    par <- c(su0, sv0, par)
    nm_par <- c("sigma_u", "sigma_v", colnames(X))

    out <- matrix(NA_real_, 3L, length(par))
    rownames(out) <- c("par", "st_err", "t-val")
    colnames(out) <- nm_par
    out[1, ] <- par

    eps <- S * (y - X %*% par[-(1:2)])
    ssq <- su0^2 + sv0^2
    mus <- -as.numeric(eps) * su0^2 / ssq
    sst <- sqrt(su0^2 * sv0^2 / ssq)
    zz <- mus / sst
    jlms <- mus + sst * stats::dnorm(zz) / pmax(stats::pnorm(zz), cz$MIN_POSITIVE)

    results <- list(
      t(out), NULL, end.time(Start.Time), par, model_name, formula, endogenous,
      instruments, jlms, exp(-jlms), Pi_ols, Sxx_ols, b_2sls, wrong_skew, S, n,
      out["par", ], out["st_err", ], out["t-val", ], call
    )
    class(results) <- "sfareg"
    names(results) <- c(
      "out", "opt", "total_time", "start_v", "model_name", "formula",
      "endogenous", "instruments", "jlms", "efficiency", "Pi", "Sigma_xi",
      "b_2sls", "wrong_skew", "S", "nobs",
      "coefficients", "std.errors", "t.values", "call"
    )
    return(results)
  }

  ## ---- the LIML / control-function likelihood ----------------------------
  ## Parameter layout:
  ##   beta (k) | log sigma_u | log sigma_v | delta (n_q) | t (m) |
  ##   [ Pi (L*m) | chol(Sigma_xi) (m(m+1)/2)  -- IVLIML only ]
  ##
  ## `t` carries the v-to-xi correlation through r = t/sqrt(1+||t||^2), which
  ## lands in the open unit ball for ANY real t. That matters: sigma_c^2 =
  ## sigma_v^2 (1 - ||r||^2), so an unconstrained rho could drive the
  ## conditional noise variance negative and the likelihood undefined. This
  ## transform makes that unreachable rather than merely penalized.
  k <- ncol(X)
  i_b <- 1:k
  i_su <- k + 1L
  i_sv <- k + 2L
  i_d <- if (n_q) (k + 3L):(k + 2L + n_q) else integer(0)
  i_t <- (k + 3L + n_q):(k + 2L + n_q + m)
  n_base <- k + 2L + n_q + m
  full <- identical(model_name, "IVLIML")
  i_pi <- if (full) (n_base + 1L):(n_base + L * m) else integer(0)
  i_ch <- if (full) (n_base + L * m + 1L):(n_base + L * m + m * (m + 1L) / 2L) else integer(0)

  .chol_from <- function(v) {
    A <- matrix(0, m, m)
    A[lower.tri(A, diag = TRUE)] <- v
    diag(A) <- exp(pmin(diag(A), 12))
    A
  }
  Lc_ols <- t(chol(Sxx_ols + diag(1e-8, m)))
  ch0 <- Lc_ols
  diag(ch0) <- log(pmax(diag(ch0), 1e-6))
  ch0 <- ch0[lower.tri(ch0, diag = TRUE)]

  like.fn <- function(th) {
    if (!all(is.finite(th))) return(cz$MAX_VALUE)
    beta <- th[i_b]
    su <- exp(pmin(th[i_su], 12))
    sv <- exp(pmin(th[i_sv], 12))
    tt <- th[i_t]
    r <- tt / sqrt(1 + sum(tt^2))

    Pim <- if (full) matrix(th[i_pi], L, m) else Pi_ols
    Lc <- if (full) .chol_from(th[i_ch]) else Lc_ols
    Xi <- P - Z %*% Pim

    ## L^{-1} xi, by triangular solve rather than an explicit inverse.
    Zi <- tryCatch(forwardsolve(Lc, t(Xi)), error = function(e) NULL)
    if (is.null(Zi)) return(cz$MAX_VALUE)
    Zi <- t(Zi)                                  # n x m

    ## mu_c,i = Sigma_v,xi Sigma_xixi^{-1} xi_i  =  sigma_v * r' L^{-1} xi_i
    mu_c <- sv * as.numeric(Zi %*% r)
    sc <- sv * sqrt(max(1 - sum(r^2), cz$MIN_POSITIVE))

    su_i <- if (n_q) su * exp(pmin(as.numeric(Q %*% th[i_d]), cz$EXP_CLIP_UPPER / 4)) else su
    s_i <- sqrt(su_i^2 + sc^2)
    lam_i <- su_i / sc

    eps <- S * as.numeric(y - X %*% beta) - mu_c
    z1 <- pmin(pmax(-lam_i * eps / s_i, cz$CLIP_Z1_LOWER), cz$CLIP_Z1_UPPER)

    ll1 <- log(2) - log(s_i) - cz$LOG_SQRT_2PI - 0.5 * (eps / s_i)^2 +
      stats::pnorm(z1, log.p = TRUE)
    ## Reduced form: log |Sigma_xixi|^{-1/2} = -sum(log(diag(L))).
    ll2 <- -m * cz$LOG_SQRT_2PI - sum(log(diag(Lc))) - 0.5 * rowSums(Zi^2)

    tot <- ll1 + ll2
    if (any(!is.finite(tot))) return(cz$MAX_VALUE)
    -sum(tot)
  }

  start_v <- c(b_2sls, log(su0), log(sv0), rep(0, n_q), rep(0, m))
  if (has_int) start_v[1L] <- start_v[1L] + S * su0 * sqrt(2 / pi)
  if (full) start_v <- c(start_v, as.numeric(Pi_ols), ch0)

  nm_par <- c(
    colnames(X), "log_sigma_u", "log_sigma_v",
    if (n_q) paste0("delta_", colnames(Q)) else NULL,
    paste0("t_", seq_len(m))
  )
  if (full) {
    nm_par <- c(
      nm_par,
      paste0("Pi_", rep(seq_len(L), m), "_", rep(seq_len(m), each = L)),
      paste0("chol_", seq_along(ch0))
    )
  }
  if (isTRUE(start_val)) names(start_v) <- nm_par

  span <- pmax(10 * abs(start_v), 10)
  lower1 <- start_v - span
  upper1 <- start_v + span
  lower1[c(i_su, i_sv)] <- log(1e-6)
  upper1[c(i_su, i_sv)] <- log(1e4)

  Opt.Bobyqa <- opt.bobyqa(
    fn = like.fn, start_v = start_v, lower.bobyqa = lower1,
    upper.bobyqa = upper1, maxit.bobyqa = maxit.bobyqa, bob.TF = TRUE,
    rhobeg = NA, rhoend = NA, verbose = verbose
  )
  start_v <- Opt.Bobyqa$start_v
  bob1 <- Opt.Bobyqa$bob1

  Opt.Psoptim <- opt.psoptim(
    fn = like.fn, start_v, lower.psoptim = lower1,
    rand.psoptim = rand.psoptim, upper.psoptim = upper1,
    maxit.psoptim, psopt.TF = PSopt, rand.order = FALSE, verbose = verbose
  )
  start_v <- Opt.Psoptim$start_v
  opt00 <- Opt.Psoptim$opt00

  Opt.Optim <- opt.optim(
    fn = like.fn, start_v = start_v, lower.optim = lower1,
    upper.optim = upper1, maxit.optim = maxit.optim,
    opt.TF = optHessian, method = Method, optHessian = TRUE, verbose = verbose
  )
  start_v <- Opt.Optim$start_v
  opt <- Opt.Optim$opt
  End.Time <- end.time(Start.Time)

  if (optHessian == FALSE & PSopt == FALSE) {
    opt <- bob1
    st_err <- rep(NA, length(opt$par))
  }
  if (optHessian == FALSE & PSopt == TRUE) {
    opt <- opt00
    st_err <- rep(NA, length(opt$par))
  }
  if (optHessian == TRUE) {
    st_err <- if (isTRUE(as.numeric(sum(colMeans(opt$hessian))) == 0)) {
      rep(NA, length(opt$par))
    } else {
      suppressWarnings(sqrt(diag(solve(opt$hessian))))
    }
  }

  ## Report on the scale users think in: sigma_u, sigma_v and the correlations
  ## rho, not their unconstrained transforms. Delta-method throughout.
  ##
  ## rho's standard errors used to be NA here, on the ground that r = t/s with
  ## s = sqrt(1 + t't) is a vector map whose Jacobian is not diagonal. It is
  ## not diagonal, but it is short:
  ##
  ##   dr_i/dt_j = delta_ij/s - t_i t_j/s^3   =>   J = (I - r r') / s,
  ##
  ## so Var(rho) = J V_tt J'. At the null rho = 0 the map is the identity, which
  ## is what makes endogeneity_test() well behaved there. See
  ## notes/code_history/ivsfm.md.
  th <- opt$par
  su <- exp(th[i_su])
  sv <- exp(th[i_sv])
  tt <- th[i_t]
  ss <- sqrt(1 + sum(tt^2))
  rho <- tt / ss

  V <- if (optHessian == TRUE) {
    suppressWarnings(tryCatch(solve(opt$hessian), error = function(e) NULL))
  } else {
    NULL
  }
  J_rho <- (diag(m) - tcrossprod(rho)) / ss
  V_rho <- if (!is.null(V) && all(is.finite(V[i_t, i_t, drop = FALSE]))) {
    J_rho %*% V[i_t, i_t, drop = FALSE] %*% t(J_rho)
  } else {
    matrix(NA_real_, m, m)
  }
  se_rho <- suppressWarnings(sqrt(diag(V_rho)))

  par <- c(th[i_b], su, sv, if (n_q) th[i_d] else NULL, rho)
  se <- c(
    st_err[i_b], su * st_err[i_su], sv * st_err[i_sv],
    if (n_q) st_err[i_d] else NULL, se_rho
  )
  rep_names <- c(
    colnames(X), "sigma_u", "sigma_v",
    if (n_q) paste0("delta_", colnames(Q)) else NULL,
    paste0("rho_", endog_v)
  )

  out <- matrix(NA_real_, 3L, length(par))
  rownames(out) <- c("par", "st_err", "t-val")
  colnames(out) <- rep_names
  out[1, ] <- par
  out[2, ] <- se
  out[3, ] <- par / se

  ## ---- inefficiency, APS (2016) Eqs (14)-(15) -----------------------------
  ## The predictor conditions on eta as well as epsilon: eta is correlated with
  ## v and so carries information about u even though u is independent of eta.
  ## That is what sigma_c < sigma_v buys, and it is a strictly better predictor
  ## than the plain JLMS.
  Pim <- if (full) matrix(th[i_pi], L, m) else Pi_ols
  Lc <- if (full) .chol_from(th[i_ch]) else Lc_ols
  Xi <- P - Z %*% Pim
  Zi <- t(forwardsolve(Lc, t(Xi)))
  mu_c <- sv * as.numeric(Zi %*% rho)
  sc <- sv * sqrt(max(1 - sum(rho^2), cz$MIN_POSITIVE))
  su_i <- if (n_q) su * exp(as.numeric(Q %*% th[i_d])) else rep(su, n)
  s_i <- sqrt(su_i^2 + sc^2)
  eps_t <- S * as.numeric(y - X %*% th[i_b]) - mu_c

  mus <- -eps_t * su_i^2 / s_i^2
  sst <- su_i * sc / s_i
  zz <- mus / sst
  jlms <- mus + sst * stats::dnorm(zz) / pmax(stats::pnorm(zz), cz$MIN_POSITIVE)

  ## vcov_rho travels with the fit so endogeneity_test() can form the JOINT
  ## Wald statistic; the diagonal alone would only give one-at-a-time t-ratios.
  dimnames(V_rho) <- list(paste0("rho_", endog_v), paste0("rho_", endog_v))

  results <- list(
    t(out), c(opt), End.Time, start_v, model_name, formula, endogenous,
    instruments, jlms, exp(-jlms), Pim, tcrossprod(Lc), rho, V_rho, sc, b_2sls,
    wrong_skew, S, n, uhet,
    out["par", ], out["st_err", ], out["t-val", ], call
  )
  class(results) <- "sfareg"
  names(results) <- c(
    "out", "opt", "total_time", "start_v", "model_name", "formula",
    "endogenous", "instruments", "jlms", "efficiency", "Pi", "Sigma_xi",
    "rho", "vcov_rho", "sigma_c", "b_2sls", "wrong_skew", "S", "nobs", "uhet",
    "coefficients", "std.errors", "t.values", "call"
  )
  results
}
