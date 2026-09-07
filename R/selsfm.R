## Greene (2010), J Prod Anal 34:15-24. Two-step: probit for the selection
## equation, then maximum SIMULATED likelihood for the frontier on the selected
## subsample. See notes/code_history/selsfm.md for the derivation, the sigma_v
## typo in the paper, and why the standard errors are conditional on alpha-hat.
selsfm <- function(selection,
                   frontier,
                   data,
                   model_name = c("greene", "kts"),
                   n_nodes = 64,
                   Nsim = "auto",
                   sim_type = c("halton", "sobol", "torus", "uniform"),
                   antithetics = FALSE,
                   seed = NULL,
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
  model_name <- match.arg(model_name)
  sim_type <- match.arg(sim_type)
  Start.Time <- Sys.time()

  if (missing(selection) || !inherits(selection, "formula")) {
    stop("selsfm(): `selection` must be a two-sided formula, e.g. d ~ z1 + z2.",
      call. = FALSE
    )
  }
  if (missing(frontier) || !inherits(frontier, "formula")) {
    stop("selsfm(): `frontier` must be a two-sided formula, e.g. y ~ x1 + x2.",
      call. = FALSE
    )
  }
  if (length(selection) != 3L || length(frontier) != 3L) {
    stop("selsfm(): both `selection` and `frontier` must have a left-hand side.",
      call. = FALSE
    )
  }
  if (missing(data) || !is.data.frame(data)) {
    stop("selsfm(): `data` must be a data.frame.", call. = FALSE)
  }
  if (length(inefdec) != 1L || !is.logical(inefdec) || is.na(inefdec)) {
    stop("selsfm(): `inefdec` must be TRUE or FALSE.", call. = FALSE)
  }
  for (nm in c("maxit.bobyqa", "maxit.psoptim", "maxit.optim")) {
    v <- get(nm)
    if (length(v) != 1L || !is.numeric(v) || !is.finite(v) || v < 1) {
      stop("selsfm(): `", nm, "` must be a single positive number.", call. = FALSE)
    }
  }

  ## The pipe syntax parameterizes variance in the rest of the package; here the
  ## two equations are separate arguments, so a pipe is an error rather than
  ## something quietly ignored.
  for (nm in c("selection", "frontier")) {
    if (any(grepl("|", deparse(get(nm)), fixed = TRUE))) {
      stop("selsfm(): `", nm, "` must not contain a `|` segment. The selection ",
        "and frontier equations are supplied as separate arguments.",
        call. = FALSE
      )
    }
  }

  ## Kumbhakar, Tsionas and Sipilainen (2009) is a different model, not a
  ## variant of Greene's, so it gets its own likelihood rather than a branch
  ## inside this one. See R/selsfm_kts.R.
  if (identical(model_name, "kts")) {
    return(.selsfm_kts_fit(selection, frontier, data, n_nodes, inefdec,
      maxit.bobyqa, maxit.psoptim, maxit.optim, start_val, PSopt, optHessian,
      Method, verbose, rand.psoptim, call))
  }

  ## ---- selection equation, estimated on the FULL sample -------------------
  mf_s <- stats::model.frame(selection, data = data, na.action = stats::na.omit)
  d_raw <- stats::model.response(mf_s)
  Z <- stats::model.matrix(stats::terms(mf_s), mf_s)

  d_num <- if (is.logical(d_raw)) {
    as.numeric(d_raw)
  } else if (is.factor(d_raw)) {
    if (nlevels(d_raw) != 2L) {
      stop("selsfm(): the selection response is a factor with ", nlevels(d_raw),
        " levels; it must be binary.",
        call. = FALSE
      )
    }
    as.numeric(d_raw) - 1
  } else {
    as.numeric(d_raw)
  }
  if (!all(d_num %in% c(0, 1))) {
    stop("selsfm(): the selection response must be binary (0/1, TRUE/FALSE, ",
      "or a two-level factor). Found values: ",
      paste(utils::head(sort(unique(d_num)), 5), collapse = ", "), ".",
      call. = FALSE
    )
  }
  if (sum(d_num) < 2L || sum(1 - d_num) < 1L) {
    stop("selsfm(): the selection indicator must contain both selected and ",
      "unselected observations; got ", sum(d_num), " selected of ",
      length(d_num), ".",
      call. = FALSE
    )
  }

  probit <- stats::glm.fit(
    x = Z, y = d_num,
    family = stats::binomial(link = "probit")
  )
  if (!isTRUE(probit$converged)) {
    warning("selsfm(): the first-stage probit did not converge; the selection ",
      "index it supplies to the frontier stage is unreliable.",
      call. = FALSE
    )
  }
  alpha_hat <- probit$coefficients
  a_full <- as.numeric(Z %*% alpha_hat)

  ## ---- frontier equation, on the SELECTED subsample only ------------------
  ## The frontier's own rows are found by matching against the selection model
  ## frame, so that a row dropped for a missing z is dropped here too.
  keep_s <- as.integer(rownames(mf_s))
  sel_rows <- keep_s[d_num == 1]

  mf_f <- stats::model.frame(frontier,
    data = data[sel_rows, , drop = FALSE],
    na.action = stats::na.omit
  )
  y <- as.numeric(stats::model.response(mf_f))
  X <- stats::model.matrix(stats::terms(mf_f), mf_f)

  ## Rows the frontier lost to missing x line up back to the selection index.
  f_rows <- as.integer(rownames(mf_f))
  a_i <- a_full[match(f_rows, keep_s)]

  n_sel <- length(y)
  n_beta <- ncol(X)
  if (n_sel <= n_beta + 3L) {
    stop("selsfm(): ", n_sel, " selected observations cannot identify ",
      n_beta + 3L, " parameters.",
      call. = FALSE
    )
  }
  if (anyNA(a_i)) {
    stop("selsfm(): could not align the selection index with the frontier ",
      "rows. Check that `data` has unique row names.",
      call. = FALSE
    )
  }

  ## Production frontier subtracts inefficiency, cost frontier adds it, so
  ## v = (y - x'b) + S * sigma_u * |U|.
  S <- if (isTRUE(inefdec)) 1 else -1

  ## ---- simulation draws ---------------------------------------------------
  .auto_nsim <- function(n) max(200L, as.integer(ceiling(3 * sqrt(n))))
  Nsim_floor <- .auto_nsim(n_sel)
  R <- if (identical(Nsim, "auto")) {
    Nsim_floor
  } else {
    if (length(Nsim) != 1L || !is.numeric(Nsim) || !is.finite(Nsim)) {
      stop("selsfm(): `Nsim` must be \"auto\" or a single number.", call. = FALSE)
    }
    max(as.integer(Nsim), 10L)
  }
  if (R < Nsim_floor) {
    warning("selsfm(): Nsim = ", R, " is below the ", Nsim_floor,
      " that this sample size calls for. Simulated ML is consistent only if ",
      "the draw count grows with n.",
      call. = FALSE
    )
  }

  absU <- abs(stats::qnorm(.sml_draws(
    n_units = n_sel, n_draws = R, dim = 1L, sim_type = sim_type,
    antithetics = antithetics, seed = seed
  )[[1L]]))

  ## ---- starting values: OLS plus Olson's moment estimators ----------------
  ols <- stats::lm.fit(x = X, y = y)
  e0 <- as.numeric(ols$residuals)
  m2 <- mean(e0^2)
  m3 <- mean(e0^3)
  ## Olson inverts the second and third central moments of the composed error.
  ## k3 is itself negative (1 - 4/pi = -0.273), so S * m3 >= 0 is the wrong-skew
  ## case, where the cube root has no positive solution; fall back rather than
  ## fail, and let the fit warn.
  k3 <- sqrt(2 / pi) * (1 - 4 / pi)
  su0 <- if (S * m3 < 0) (S * m3 / k3)^(1 / 3) else 0.5 * stats::sd(e0)
  su0 <- if (is.finite(su0)) max(su0, 1e-3) else max(0.5 * stats::sd(e0), 1e-3)
  sv0 <- sqrt(max(m2 - (1 - 2 / pi) * su0^2, 1e-6))
  b0 <- as.numeric(ols$coefficients)
  if (attr(stats::terms(mf_f), "intercept") == 1L) {
    b0[1L] <- b0[1L] + S * su0 * sqrt(2 / pi)
  }
  start_v <- c(su0, sv0, 0, b0)

  par_names <- c("sigma_u", "sigma_v", "rho", colnames(X))
  if (isTRUE(start_val)) names(start_v) <- par_names

  ## ---- simulated log-likelihood -------------------------------------------
  ## Greene's (15). Everything is kept in logs: the summand underflows to zero
  ## in the tails long before the log-sum-exp does.
  cz <- .SFA_CONSTANTS
  like.fn <- function(theta) {
    su <- abs(theta[1L])
    sv <- abs(theta[2L])
    rho <- theta[3L]
    beta <- theta[4L:length(theta)]

    if (!all(is.finite(theta))) return(cz$MAX_VALUE)
    su <- max(su, cz$MIN_POSITIVE)
    sv <- max(sv, cz$MIN_POSITIVE)
    if (abs(rho) >= 1) return(cz$MAX_VALUE)

    r0 <- as.numeric(y - X %*% beta)
    ## v_ir, one column per draw.
    v <- r0 + S * su * absU
    root <- sqrt(1 - rho^2)

    log_dens <- -0.5 * (v / sv)^2 - log(sv) - cz$LOG_SQRT_2PI
    zarg <- (rho * (v / sv) + a_i) / root
    zarg <- pmin(pmax(zarg, cz$CLIP_Z1_LOWER), cz$CLIP_Z1_UPPER)
    log_phi <- stats::pnorm(zarg, log.p = TRUE)

    ll_i <- .log_row_sum_exp(log_dens + log_phi) - log(R)
    if (any(!is.finite(ll_i))) return(cz$MAX_VALUE)
    ## The scaffold MINIMIZES; every likelihood in this package returns the
    ## negative summed log-likelihood.
    -sum(ll_i)
  }

  ## ---- three-stage minimizer ----------------------------------------------
  ## Bounds are built here rather than through lower.start(), whose layout is
  ## keyed to model_name; rho is the only parameter with a natural box.
  b_span <- pmax(10 * abs(b0), 10)
  lower1 <- c(1e-6, 1e-6, -0.995, b0 - b_span)
  upper1 <- c(max(10 * su0, 10), max(10 * sv0, 10), 0.995, b0 + b_span)

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

  ## Both scales enter through abs(), so their sign is not identified.
  opt$par[1:2] <- abs(opt$par[1:2])

  out <- matrix(NA_real_, 3L, length(opt$par))
  rownames(out) <- c("par", "st_err", "t-val")
  colnames(out) <- par_names
  out[1, ] <- opt$par
  out[2, ] <- st_err
  out[3, ] <- opt$par / st_err

  ## ---- observation-specific inefficiency, Greene's (21) -------------------
  ## Simulated Bayes over the same draws the likelihood used: E[u_i | eps_i] is
  ## a weighted mean of the draws, the weights being each draw's share of the
  ## simulated density.
  su <- opt$par[1L]
  sv <- opt$par[2L]
  rho <- opt$par[3L]
  beta <- opt$par[4L:length(opt$par)]
  r0 <- as.numeric(y - X %*% beta)
  v <- r0 + S * su * absU
  root <- sqrt(1 - rho^2)
  lg <- -0.5 * (v / sv)^2 - log(sv) - cz$LOG_SQRT_2PI
  zarg <- pmin(pmax((rho * (v / sv) + a_i) / root, cz$CLIP_Z1_LOWER), cz$CLIP_Z1_UPPER)
  lg <- lg + stats::pnorm(zarg, log.p = TRUE)
  w <- exp(lg - .log_row_sum_exp(lg))
  jlms <- su * rowSums(w * absU)
  eff <- exp(-jlms)

  results <- list(
    t(out), c(opt), End.Time, start_v, "SEL", frontier, selection, jlms, eff,
    alpha_hat, probit, R, sim_type, n_sel, length(d_num), S, n_sel,
    out["par", ], out["st_err", ], out["t-val", ], call
  )
  class(results) <- "sfareg"
  names(results) <- c(
    "out", "opt", "total_time", "start_v", "model_name", "formula",
    "selection_formula", "jlms", "efficiency", "alpha", "probit", "Nsim",
    "sim_type", "n_selected", "n_total", "S", "nobs",
    "coefficients", "std.errors", "t.values", "call"
  )
  ## `nobs` is the SELECTED count, not the rows supplied. The second-stage
  ## likelihood is a sum over the selected observations only, so that is the n
  ## that BIC() must divide by; without this nobs.sfareg() falls through to
  ## re-evaluating the call's `data` and reports the whole sample.
  results
}
