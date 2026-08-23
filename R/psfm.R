psfm <- function(formula,
                 model_name = c(
                   "TRE_Z", "GTRE_Z", "TRE", "GTRE", "GTRE_FML", "TFE", "TFE_WMLE",
                   "FD", "GTRE_SEQ1", "GTRE_SEQ2", "SSFE", "PL80", "PL80_MVTN",
                   "BC92", "K1990", "K1990modified"
                 ),
                 data,
                 maxit.bobyqa = 5000,
                 maxit.nlminb = 500,
                 maxit.psoptim = 100,
                 maxit.optim = 1000,
                 REPORT = 1,
                 trace = 3,
                 pgtol = 0,
                 individual,
                 halton_num = NULL,
                 start_val = FALSE,
                 gamma = FALSE,
                 PSopt = FALSE,
                 optHessian = TRUE,
                 inefdec = TRUE,
                 Method = "L-BFGS-B",
                 verbose = FALSE,
                 rand.gtre = NULL,
                 rand.psoptim = NULL,
                 OPG_calc = FALSE,
                 estimator = c("fiml", "sml", "seq1", "seq2"),
                 collinear_action = c("start_only", "error", "warn_drop"),
                 time = NULL,
                 tfe_lambda_max = 100) {
  ## call/model_name resolution moved ahead of .check_model_formula_pipes() --
  ## see sfm.R's identical fix for why (calling the pipe check on the raw,
  ## unresolved multi-choice model_name default errored "the condition has
  ## length > 1" for any caller relying on the default rather than specifying
  ## model_name explicitly).
  call <- match.call()
  model_name <- .match_model_name(model_name, eval(formals()$model_name))
  collinear_action <- match.arg(collinear_action)

  ## ------------------------------------------------------------------
  ## `estimator` selects HOW the four-component GTRE model is estimated,
  ## the way sfm()'s `estimator` chooses between maximum likelihood and
  ## corrected OLS for one and the same model. The four ways of fitting
  ## GTRE were previously four separate model_name values, which made
  ## them look like four different models rather than four routes to the
  ## same one:
  ##
  ##   "fiml" (default)  full information ML through the closed-skew-normal
  ##                     representation. Deterministic -- nothing is
  ##                     integrated by Monte Carlo, so repeated fits to the
  ##                     same data agree exactly. Requires a BALANCED panel.
  ##   "sml"             simulated ML over Halton draws. Handles unbalanced
  ##                     panels. This is what "GTRE" meant through 1.1.3.
  ##   "seq1" / "seq2"   the two-step moment-based estimators. Not maximum
  ##                     likelihood, so they carry no logLik().
  ##
  ## The older names still work and are unchanged: GTRE_FML, GTRE_SEQ1 and
  ## GTRE_SEQ2 name the same three routes directly.
  ## ------------------------------------------------------------------
  estimator_supplied <- !missing(estimator)
  estimator <- .match_model_name(estimator, eval(formals()$estimator), arg = "estimator")

  if (identical(model_name, "GTRE")) {
    if (!estimator_supplied) {
      ## Same situation as the TFE rename above: an existing script asking
      ## for "GTRE" now gets a different estimator than it did, so say so
      ## rather than quietly changing the answer.
      warning("model_name = \"GTRE\" now defaults to estimator = \"fiml\" ",
        "(full information ML via the closed-skew-normal representation). ",
        "Through sfa 1.1.3 it fitted the simulated-ML estimator, which is now ",
        "estimator = \"sml\". Pass `estimator` explicitly to silence this warning.",
        call. = FALSE
      )
    }

    ## FIML is built at a single T, so it cannot fit an unbalanced panel. If
    ## the user merely took the default, fall back to simulated ML and say so
    ## -- erroring would make GTRE unusable by default on a whole class of
    ## data. If they asked for "fiml" deliberately, honour that by telling
    ## them it cannot be done rather than silently fitting something else.
    if (identical(estimator, "fiml")) {
      t_i <- tryCatch(
        as.numeric(table(as.character(data[, individual]))),
        error = function(e) NULL
      )
      if (!is.null(t_i) && length(t_i) && length(unique(t_i)) != 1L) {
        if (estimator_supplied) {
          stop("estimator = \"fiml\" requires a BALANCED panel: the ",
            "closed-skew-normal representation is built at a single T. Found T ",
            "ranging from ", min(t_i), " to ", max(t_i), ". Use estimator = ",
            "\"sml\", which handles unbalanced panels.",
            call. = FALSE
          )
        }
        warning("Panel is unbalanced (T ranges from ", min(t_i), " to ", max(t_i),
          "), which the default estimator = \"fiml\" cannot fit. Falling back to ",
          "estimator = \"sml\" (simulated ML). Pass `estimator` explicitly to ",
          "choose deliberately.",
          call. = FALSE
        )
        estimator <- "sml"
      }
    }

    model_name <- switch(estimator,
      fiml = "GTRE_FML",
      sml  = "GTRE",
      seq1 = "GTRE_SEQ1",
      seq2 = "GTRE_SEQ2"
    )
  } else if (estimator_supplied) {
    warning("`estimator` only applies to model_name = \"GTRE\"; ignored for \"",
      model_name, "\".",
      call. = FALSE
    )
  }

  ## The meaning of model_name = "TFE" CHANGED in sfa 1.1.3. Through 1.1.2 it
  ## was Chen, Schmidt and Wang's (2014) within maximum-likelihood estimator,
  ## which is not a "true fixed effects" estimator at all in Greene's (2005)
  ## sense -- it is specifically the estimator that AVOIDS estimating the firm
  ## effects. Now that Greene's TFE is also available, "TFE" names Greene's
  ## estimator (the standard meaning in the literature) and the CSW estimator
  ## is "TFE_WMLE". Warn rather than silently returning a different estimator
  ## than an existing script asked for; this warning is intended to stay for
  ## one release cycle.
  if (identical(model_name, "TFE")) {
    warning("model_name = \"TFE\" now fits Greene's (2005) true fixed effects MLE. ",
      "Through sfa 1.1.2 this name selected Chen, Schmidt and Wang's (2014) ",
      "within MLE, which is now model_name = \"TFE_WMLE\". Specify the name ",
      "explicitly to silence this warning.",
      call. = FALSE
    )
  }

  .check_model_formula_pipes(formula, model_name)

  ## Accept an ordinary data.frame (or tibble/data.table) as well as a
  ## plm::pdata.frame. Everything downstream -- data_proc(), start_panel()'s
  ## plm() calls, the per-firm split in each likelihood -- assumes the panel
  ## index is already attached, and handed a plain data.frame the fit used to
  ## die with the opaque message "empty model" (preceded by an unrelated
  ## -looking "'-' not meaningful for factors" warning) that gave no hint the
  ## real problem was a missing index. psfm() already receives `individual`
  ## (and `time`), so it has everything needed to build the index itself.
  data <- .as_panel_data(data, individual, time)

  ## Pre-estimation collinearity check for the models whose starting values
  ## come from plm(..., model = "random"). Run here, BEFORE data_proc(), so
  ## that collinear_action = "warn_drop" can actually reshape the model the
  ## rest of the function estimates. See .check_collinearity() for why the
  ## between-individual design can be singular while the pooled one is fine.
  .hr_uses_re_start <- !(model_name %in% c("TFE", "TFE_WMLE", "FD", "SSFE", "PL80", "BC92", "K1990", "K1990modified")) &&
    isFALSE(is.numeric(start_val))
  collinear_chk <- NULL
  if (.hr_uses_re_start) {
    fx_probe <- tryCatch(stats::formula(Formula::Formula(formula), lhs = 1, rhs = 1),
      error = function(e) NULL
    )
    if (!is.null(fx_probe)) {
      collinear_chk <- tryCatch(.check_collinearity(fx_probe, data, individual),
        error = function(e) NULL
      )
    }

    if (!is.null(collinear_chk) && length(collinear_chk$between_drop)) {
      if (identical(collinear_action, "error")) {
        stop(.collinearity_message(collinear_chk, "error"), call. = FALSE)
      }

      if (identical(collinear_action, "warn_drop")) {
        warning(.collinearity_message(collinear_chk, "warn_drop"), call. = FALSE)
        drop_terms <- .terms_for_columns(fx_probe, data, collinear_chk$between_drop)
        if (length(drop_terms)) {
          keep <- setdiff(attr(stats::terms(fx_probe), "term.labels"), drop_terms)
          formula <- stats::reformulate(if (length(keep)) keep else "1",
            response = all.vars(fx_probe)[1]
          )
        }
        collinear_chk <- NULL ## handled here; nothing left for start_panel
      }
    }
  }

  DR1 <- data_proc(formula, data, model_name, individual, inefdec)

  if (1 == 1) {
    formula <- DR1$formula
    data_orig <- DR1$data_orig
    form_parts <- DR1$form_parts
    formula_x <- DR1$formula_x
    y_var <- DR1$y_var
    model_name <- DR1$model_name
    data_x <- DR1$data_x
    intercept <- DR1$intercept
    inefdec_n <- DR1$inefdec_n
    inefdec_TF <- DR1$inefdec_TF
    x_vars_vec <- DR1$x_vars_vec
    n_x_vars <- DR1$n_x_vars
    x_vars <- DR1$x_vars
    x_x_vec <- DR1$x_x_vec
    fancy_vars <- DR1$fancy_vars
    fancy_vars_z <- DR1$fancy_vars_z
    n_z_vars <- DR1$n_z_vars
    N <- DR1$N
    data_z <- DR1$data_z
    if (length(unlist(form_parts)) > 3) {
      formula_z <- DR1$formula_z
      intercept_z <- DR1$intercept_z
      n_z_vars <- DR1$n_z_vars
      z_vars <- DR1$z_vars
      z_vars_vec <- DR1$z_vars_vec
      z_z_vec <- DR1$z_z_vec
    }
    if (length(unlist(form_parts)) > 4) { ## might need to incorporate this above
      formula_zp <- DR1$formula_zp
      intercept_zp <- DR1$intercept_zp
      n_zp_vars <- DR1$n_zp_vars
      zp_vars <- DR1$zp_vars
      zp_vars_vec <- DR1$zp_vars_vec
      zp_zp_vec <- DR1$zp_zp_vec
    }

    ## Pitt and Lee (1981) Model III -- the multivariate truncated normal
    ## likelihood from their Appendix 2. Handled here, before the shared
    ## start_panel() machinery: its parameters are (sigma_v, sigma_u, rho,
    ## beta) and a pooled OLS supplies every start value it needs, so there is
    ## nothing for start_panel()'s model-specific dispatch to contribute.
    if (identical(model_name, "PL80_MVTN")) {
      return(.psfm_pl_mvtn(
        data = data, y_var = y_var, x_vars_vec = x_vars_vec,
        individual = individual, inefdec_n = inefdec_n,
        maxit.optim = maxit.optim, Method = Method, optHessian = optHessian,
        start_val = start_val, verbose = verbose, call = call,
        formula = formula
      ))
    }

    Start_Panel <- start_panel(formula_x, data, model_name, start_val, intercept, x_vars_vec,
      individual = individual, collinear_chk = collinear_chk
    )

    alpha_hat <- Start_Panel$alpha_hat
    beta_0 <- Start_Panel$beta_0
    beta_0_st <- Start_Panel$beta_0_st
    beta_hat <- Start_Panel$beta_hat
    beta_se <- Start_Panel$beta_se
    epsilon_hat <- Start_Panel$epsilon_hat
    exp_eta <- Start_Panel$exp_eta
    exp_u <- Start_Panel$exp_u
    lambda <- Start_Panel$lambda
    out <- Start_Panel$out
    plm_gtre <- Start_Panel$plm_gtre
    sfa_alp <- Start_Panel$sfa_alp
    sfa_eps <- Start_Panel$sfa_eps
    sigma <- Start_Panel$sigma
    sigma_h <- Start_Panel$sigma_h
    sigma_r <- Start_Panel$sigma_r
    sigma_u <- Start_Panel$sigma_u
    sigma_v <- Start_Panel$sigma_v
    start_v <- Start_Panel$start_v
    plm_tfe <- Start_Panel$plm_tfe
    plm_fd <- Start_Panel$plm_fd

    DR2 <- data_proc2(
      data, data_x, fancy_vars, fancy_vars_z, data_z, y_var,
      x_vars_vec, halton_num, individual, N, model_name, rand.gtre
    )

    data <- DR2$data
    Y <- DR2$Y
    data_i_vars <- DR2$data_i_vars
  }

  if (model_name %in% c("GTRE", "TRE")) {
    R <- DR2$R
    R_H <- DR2$R_H
    indiv <- DR2$indiv
    t <- DR2$t
    data_i <- DR2$data_i
    eps <- DR2$eps
    R_h1 <- DR2$R_h1
    R_h2 <- DR2$R_h2
    data_x <- DR2$data_x

    ## ---- stacked copies for the vectorized likelihood ---------------------
    ## Built once per fit from the same per-firm lists, so row order matches
    ## exactly. R_h1/R_h2 have identical rows within a firm (the draws are
    ## firm-level and shared across firms), so only their first row is needed
    ## and the composed error is rank one: eps[i, j] = (y_i - x_i'beta) + c_j.
    ## That lets outer() build the n x R matrix without materializing per-firm
    ## blocks.
    .yv <- unlist(lapply(Y, function(m) m[, 1L]))
    .Xall <- do.call(rbind, data_i_vars)
    .gid <- rep(seq_len(N), times = vapply(Y, nrow, integer(1)))
    .h1 <- R_h1[[1]][1, ]
    .h2 <- R_h2[[1]][1, ]

    fn_1_loop <- function(x) {
      if (model_name == "GTRE") {
        x_x_vec <- x[5:as.numeric(n_x_vars + 4)]
      }
      if (model_name == "TRE") {
        x_x_vec <- x[4:as.numeric(n_x_vars + 3)]
      }
      fn1 <- function(ii) {
        if (model_name == "GTRE") {
          eps_neg <- eps
          eps_neg[[ii]] <- Y[[ii]] + x[3] * R_h1[[ii]] + x[4] * R_h2[[ii]] * inefdec_n
        }

        if (model_name == "GTRE") {
          eps[[ii]] <- Y[[ii]] - x[3] * R_h1[[ii]] + x[4] * R_h2[[ii]] * inefdec_n
        }
        if (model_name == "TRE") {
          eps[[ii]] <- Y[[ii]] - x[3] * R_h1[[ii]]
        }

        ## Single matrix multiply instead of n_x_vars separate Ti x R
        ## temporary matrices/subtractions; the resulting length-Ti vector
        ## recycles column-wise against the Ti x R matrix, same as before.
        frontier_mean <- as.vector(data_i_vars[[ii]] %*% x_x_vec)
        eps[[ii]] <- eps[[ii]] - frontier_mean

        if (model_name == "GTRE") {
          eps_neg[[ii]] <- eps_neg[[ii]] - frontier_mean

          eps_neg[[ii]] <- inefdec_n * eps_neg[[ii]]
          z1_neg <- eps_neg[[ii]] / x[2]
          z2_neg <- -eps_neg[[ii]] * x[1] / x[2]
        }

        eps[[ii]] <- inefdec_n * eps[[ii]]
        z0 <- 2 / x[2]
        z1 <- eps[[ii]] / x[2]
        z2 <- -eps[[ii]] * x[1] / x[2]
        z1[z1 > .SFA_CONSTANTS$CLIP_Z1_UPPER] <- .SFA_CONSTANTS$CLIP_Z1_UPPER
        z1[z1 < .SFA_CONSTANTS$CLIP_Z1_LOWER] <- .SFA_CONSTANTS$CLIP_Z1_LOWER
        z2[z2 > .SFA_CONSTANTS$CLIP_Z2_UPPER] <- .SFA_CONSTANTS$CLIP_Z2_UPPER
        z2[z2 < .SFA_CONSTANTS$CLIP_Z2_LOWER] <- .SFA_CONSTANTS$CLIP_Z2_LOWER

        prod_vec_n0 <- log(max(mean(.col_prods(
          z0 * dnorm(z1) * pmax(pnorm(z2), eps[[ii]] * 0 + .Machine$double.xmin)
        )), .Machine$double.xmin))
        if (model_name == "TRE") {
          prod_vec_n <- prod_vec_n0
        }
        if (model_name == "GTRE") {
          prod_vec_n1 <- log(max(mean(.col_prods(
            z0 * dnorm(z1_neg) * pmax(pnorm(z2_neg), eps_neg[[ii]] * 0 + .Machine$double.xmin)
          )), .Machine$double.xmin))

          prod_vec_n <- 0.5 * (prod_vec_n0 + prod_vec_n1)
        }

        return(-prod_vec_n)
      }

      fn1_apply <- unlist(lapply(1:N, fn1))
      fn1_apply[is.nan(fn1_apply)] <- sqrt(.SFA_CONSTANTS$MAX_VALUE / length(x))
      fn1_apply[is.infinite(fn1_apply)] <- sqrt(.SFA_CONSTANTS$MAX_VALUE / length(x))

      return(sum(fn1_apply[is.finite(fn1_apply)]))
    }

    ## Vectorized equivalent: one dnorm/pnorm pair over the whole stacked
    ## (n x R) matrix instead of N pairs over Ti x R blocks. OFF by default --
    ## measured, it is not worth the change:
    ##
    ##   model   loop    vectorized   speedup   agreement
    ##   TRE     1.9s    1.4s         1.35x     coefficients match to 5.7e-05
    ##   GTRE    7.0s    6.8s         1.03x     coefficients differ by 7.6e-02
    ##
    ## The speedup is small because R (the number of Halton draws) defaults to
    ## ceiling(sqrt(n)) + 100, so each firm block is already a few thousand
    ## elements -- big enough that R's per-call overhead is not what dominates.
    ## Stacking removes that overhead but adds an outer() and a rowsum() of the
    ## same total size, so most of the gain is given straight back.
    ##
    ## GTRE additionally DISAGREES, and the reason is a genuine asymmetry in the
    ## loop version rather than an error in this one: z1_neg/z2_neg (the
    ## negative half of the +/- r mixture) are computed just above the clipping
    ## block and are never clipped, while z1/z2 are. The two halves of the same
    ## mixture are therefore treated differently. This version clips both,
    ## consistently, which is why it lands on a slightly different optimum.
    ## Worth resolving deliberately -- it is not something to change silently as
    ## a side effect of a speed edit.
    fn_1_vec <- function(x) {
      x_x_vec <- if (model_name == "GTRE") {
        x[5:as.numeric(n_x_vars + 4)]
      } else {
        x[4:as.numeric(n_x_vars + 3)]
      }
      lam <- x[1]
      sig <- x[2]
      if (!is.finite(lam) || !is.finite(sig) || sig <= 0) {
        return(sqrt(.SFA_CONSTANTS$MAX_VALUE))
      }

      base <- .yv - as.vector(.Xall %*% x_x_vec)
      c_pos <- if (model_name == "GTRE") -x[3] * .h1 + x[4] * .h2 * inefdec_n else -x[3] * .h1
      ll <- .gtre_sim_logdens(outer(base, c_pos, "+") * inefdec_n, lam, sig, .gid, N)

      if (model_name == "GTRE") {
        ## Second half of the +/- r mixture: r is symmetric, so the simulated
        ## density averages the draw and its reflection.
        c_neg <- x[3] * .h1 + x[4] * .h2 * inefdec_n
        ll <- 0.5 * (ll + .gtre_sim_logdens(outer(base, c_neg, "+") * inefdec_n, lam, sig, .gid, N))
      }
      ll[!is.finite(ll)] <- -sqrt(.SFA_CONSTANTS$MAX_VALUE / length(x))
      -sum(ll[is.finite(ll)])
    }

    fn_1 <- function(x) {
      if (isTRUE(getOption("sfa.gtre_vectorized", FALSE))) fn_1_vec(x) else fn_1_loop(x)
    }

    Start.Time <- start.time()

    ## ---- phased start: efficiency block first, beta held fixed -------------
    ## These likelihoods are dominated by the cost of a single evaluation
    ## (simulated ML over Halton draws), so the way to make them fast is to
    ## need fewer evaluations rather than to make each one cheaper. The
    ## frontier coefficients arrive from a random-effects regression already
    ## close to their MLE and carry little of the composed-error information,
    ## whereas the variance components are where the curvature and the bad
    ## local optima are. .phased_start() therefore grid-searches the variance
    ## block with beta pinned, polishes just that block, and hands the result
    ## to the full stack below -- which still does the complete maximization.
    ##
    ## The block is the leading 3 (TRE: lambda, sigma, sigma_r) or 4 (GTRE:
    ## plus sigma_h) entries, matching lower.start()'s own per-model layout.
    ## OFF by default -- see .phased_start()'s header for the benchmark that
    ## rejected it. Enable with options(sfa.phased_start = TRUE) to experiment.
    if (isTRUE(getOption("sfa.phased_start", FALSE))) {
      ## Efficiency block, then the frontier intercept: for GTRE the layout is
      ## (lambda, sigma, sigma_r, sigma_h, beta_0, x...), for TRE the same
      ## without sigma_h. intercept = NULL when the formula has none.
      .eff_idx <- if (model_name == "GTRE") 1:4 else 1:3
      .int_idx <- if (isTRUE(intercept == 1)) max(.eff_idx) + 1L else NULL
      start_v <- .phased_start(
        fn = fn_1, start_v = start_v, idx = .eff_idx,
        int_idx = .int_idx, lower = .SFA_CONSTANTS$MIN_POSITIVE,
        grid = c(0.5, 0.75, 1.5, 2), verbose = verbose
      )
    }

    Lower.Start <- lower.start(start_v, model_name, differ = 3)
    Opt.Bobyqa <- opt.bobyqa(fn = fn_1, start_v = start_v, lower.bobyqa = Lower.Start$lower1, maxit.bobyqa = maxit.bobyqa, bob.TF = TRUE, verbose = verbose)
    start_v <- Opt.Bobyqa$start_v
    start_feval <- Opt.Bobyqa$start_feval
    bob1 <- Opt.Bobyqa$bob1

    Lower.Start <- lower.start(start_v, model_name, differ = 2)
    Opt.Psoptim <- opt.psoptim(
      fn = fn_1, start_v = start_v, lower.psoptim = Lower.Start$lower1,
      rand.psoptim = rand.psoptim, upper.psoptim = Lower.Start$upper1, maxit.psoptim, psopt.TF = PSopt,
      verbose = verbose
    )
    start_v <- Opt.Psoptim$start_v
    start_feval <- Opt.Psoptim$start_feval
    opt00 <- Opt.Psoptim$opt00

    Lower.Start <- lower.start(start_v, model_name, differ = 0.5)
    Opt.Optim <- opt.optim(
      fn = fn_1, start_v = start_v, lower.optim = Lower.Start$lower1, upper.optim = Lower.Start$upper1_open,
      maxit.optim = maxit.optim, opt.TF = optHessian, method = Method, optHessian = optHessian, verbose = verbose
    )
    start_v <- Opt.Optim$start_v
    start_feval <- Opt.Optim$start_feval
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
    out <- matrix(0, nrow = 3, ncol = length(start_v))
    rownames(out) <- c("par", "st_err", "t-val")
    colnames(out) <- if (model_name == "GTRE") {
      c("lambda", "sigma", "sigr", "sigh", colnames(data_x))
    } else {
      c("lambda", "sigma", "sigr", colnames(data_x))
    }

    st_err <- if (isTRUE(any(opt$hessian == 0)) | optHessian == FALSE) {
      rep(NA, length(opt$par))
    } else {
      suppressWarnings(sqrt(pmax(diag(solve(opt$hessian)), 0)))
    }
    t_val <- opt$par / st_err
    out[1, ] <- opt$par
    out[2, ] <- st_err
    out[3, ] <- t_val

    ## TE Measurements : GTRE
    if (model_name == "GTRE") {
      beta <- opt$par[-c(1:4)]
      lamb <- opt$par[1]
      sig <- opt$par[2]
      sig_u <- (lamb * sig) / sqrt(1 + lamb^2)
      sig_v <- sig_u / lamb
      sig_r <- opt$par[3]
      sig_h <- opt$par[4]

      e_i <- as.list(rep(0, N))
      A <- as.list(rep(0, N))
      SIG <- as.list(rep(0, N))
      VEE <- as.list(rep(0, N))
      LAM <- as.list(rep(0, N))
      ARR <- as.list(rep(0, N))
      n <- sum(t)
      U <- rep(0, n)

      ## Lambda_i/ARR_i via .safe_linear_combo() (matrix_utils.R), which
      ## ridges/Cholesky-inverts VEE_i/SIG_i/K robustly and guards against
      ## solve()'s result carrying spurious/inconsistent dimnames (a
      ## cbind()-of-unnamed-vector-and-matrix artifact that can otherwise
      ## make isSymmetric() reject a numerically-symmetric matrix -- see
      ## .safe_symmetrize()'s header in matrix_utils.R). GTRE_Z's own
      ## .gtre_te() already uses this same helper.
      for (i in seq_len(N)) {
        e_i[[i]] <- pmin(inefdec_n * (Y[[i]][, 1] - rowSums(t(t(data_i_vars[[i]]) * beta))), Y[[i]][, 1] * 0)
        A[[i]] <- -cbind(rep(1, t[i]), diag(t[i]))
        SIG[[i]] <- sig_v^2 * diag(t[i]) + sig_r^2 * rep(1, t[i]) %*% t(rep(1, t[i]))
        VEE[[i]] <- rbind(c(sig_h^2, rep(0, t[i])), cbind(rep(0, t[i]), sig_u^2 * diag(t[i])))
        lc <- .safe_linear_combo(VEE_i = VEE[[i]], A_i = A[[i]], SIG_i = SIG[[i]])
        LAM[[i]] <- lc$LAM
        ARR[[i]] <- lc$ARR
      }

      res_d_fn <- function(i) {
        ptmvnorm(
          lowerx = rep(0, t[i] + 1), upperx = rep(Inf, t[i] + 1),
          mean = as.numeric(ARR[[i]] %*% e_i[[i]]),
          sigma = LAM[[i]]
        )[1]
      }

      res_n_fn <- function(i) {
        ptmvnorm(
          lowerx = rep(0, t[i] + 1), upperx = rep(Inf, t[i] + 1),
          mean = as.numeric(ARR[[i]] %*% e_i[[i]] + LAM[[i]] %*% c(-1, rep(0, t[i]))),
          sigma = LAM[[i]]
        )[1]
      }

      res_d <- lapply(seq(1, N, 1), res_d_fn)
      res_n <- lapply(seq(1, N, 1), res_n_fn)

      H_fn <- function(i) {
        (max(res_n[[i]], .SFA_CONSTANTS$MIN_POSITIVE) / max(res_d[[i]], .SFA_CONSTANTS$MIN_POSITIVE)) *
          exp(t(c(-1, rep(0, t[i]))) %*% ARR[[i]] %*% e_i[[i]] + 0.5 * t(c(-1, rep(0, t[i]))) %*% LAM[[i]] %*% c(-1, rep(0, t[i])))
      }

      H <- unlist(lapply(seq(1, N, 1), H_fn))
      H <- pmin(H, rep(1, length(H)))

      new_t_exp <- as.list(rep(0, n))

      for (i in seq_len(N)) {
        for (j in seq_len(t[i])) {
          h <- cumsum(t)[i] - t[i] + j
          new_t <- rep(0, t[i] + 1)
          new_t[j + 1] <- -1
          new_t_exp[[h]] <- new_t
        }
      }


      t_cum <- c(cumsum(t))
      t_exp <- e_i_exp <- A_exp <- SIG_exp <- VEE_exp <- LAM_exp <- ARR_exp <- res_d_exp <- as.list(rep(0, n))

      for (m in seq_len(N)) {
        B <- t_cum[m]
        A <- B + 1 - t[m]
        t_exp[A:B] <- rep(t[m], t[m])
        e_i_exp[A:B] <- rep(e_i[m], t[m])
        A_exp[A:B] <- rep(A[m], t[m])
        SIG_exp[A:B] <- rep(SIG[m], t[m])
        VEE_exp[A:B] <- rep(VEE[m], t[m])
        LAM_exp[A:B] <- rep(LAM[m], t[m])
        ARR_exp[A:B] <- rep(ARR[m], t[m])
        res_d_exp[A:B] <- rep(res_d[m], t[m])
      }

      res_n_t_fn <- function(i) {
        ptmvnorm(
          lowerx = rep(0, t_exp[[i]] + 1), upperx = rep(Inf, t_exp[[i]] + 1),
          mean = as.numeric(ARR_exp[[i]] %*% e_i_exp[[i]] + LAM_exp[[i]] %*% new_t_exp[[i]]),
          sigma = LAM_exp[[i]]
        )[1]
      }

      res_n_t <- lapply(seq(1, n, 1), res_n_t_fn)


      U_fn <- function(i) {
        (max(res_n_t[[i]], .SFA_CONSTANTS$MIN_POSITIVE) / max(res_d_exp[[i]], .SFA_CONSTANTS$MIN_POSITIVE)) *
          exp(t(new_t_exp[[i]]) %*% ARR_exp[[i]] %*% e_i_exp[[i]] + 0.5 * t(new_t_exp[[i]]) %*% LAM_exp[[i]] %*% new_t_exp[[i]])
      }

      U <- unlist(lapply(seq(1, n, 1), U_fn))
      U <- pmin(U, rep(1, length(U)))
    }
    ## TE Measurements : TRE
    if (model_name == "TRE") {
      beta <- opt$par[-c(1:3)]
      lamb <- opt$par[1]
      sig <- opt$par[2]
      sig_u <- (lamb * sig) / sqrt(1 + lamb^2)
      sig_v <- sig_u / lamb

      Y_mean <- rep(0, N)
      X_mean <- matrix(0, N, ncol = length(x_vars_vec))
      colnames(X_mean) <- x_vars_vec

      for (ii in seq_len(N)) {
        data_i[[ii]] <- data[which(data[, c(individual)] == indiv[ii]), ]
        Y_mean[ii] <- mean(as.numeric(data_i[[ii]][, y_var]))
        X_mean[ii, ] <- colMeans(data.frame(data_i[[ii]][, c(x_vars_vec)]))
      }

      r_hat_m <- Y_mean - rowSums(t(beta * t(X_mean))) + inefdec_n * sqrt(2 / pi) * sig_u
      r_hat_m_exp <- rep(0, sum(t))
      t_cum <- c(cumsum(t))

      for (m in 1:length(t)) {
        B <- t_cum[m]
        A <- B + 1 - t[m]
        r_hat_m_exp[A:B] <- rep(r_hat_m[m], t[m])
      }

      eps_hat <- pmin(inefdec_n * (data[, y_var] - rowSums(t(t(data[, c(x_vars_vec)]) * beta)) - r_hat_m_exp), data[, y_var] * 0)
      sig_star <- sig_u * sig_v / sig
      inner <- (lamb * eps_hat) / sig
      exp_u_hat <- ((1 - pnorm((sig_u * sig_v / sig) + inner)) / (1 - pnorm(inner))) * exp((sig_u^2 / sig^2) * (eps_hat + 0.5 * sig_v^2))
      U <- exp_u_hat
    }

    # cor(U,exp(-p_data_trial$u))
    # cor(H,exp(-unique(p_data_trial$h) ))

    if (model_name == "GTRE") {
      results <- list(t(out), c(opt), data, End.Time, start_v, model_name, formula, U, H, out["par", ], out["st_err", ], out["t-val", ], call)
      class(results) <- "sfareg"
      names(results) <- c("out", "opt", "data", "total_time", "start_v", "model_name", "formula", "U", "H", "coefficients", "std.errors", "t.values", "call")
    }
    if (model_name == "TRE") {
      results <- list(t(out), c(opt), data, End.Time, start_v, model_name, formula, U, out["par", ], out["st_err", ], out["t-val", ], call)
      class(results) <- "sfareg"
      names(results) <- c("out", "opt", "data", "total_time", "start_v", "model_name", "formula", "U", "coefficients", "std.errors", "t.values", "call")
    }

    return(results)
  }
  if (model_name == "GTRE_Z") {
    ## Default starting values for variance equations
    delta <- rep(0.1, length(z_vars))
    delta_p <- rep(0.1, length(zp_vars))

    ## Starting vector
    if (isTRUE(is.numeric(start_val))) {
      start_v <- start_val
    } else {
      start_v <- if (is.na(beta_0_st)) {
        unname(c(sigma_v, sigma_r, beta_hat, delta, delta_p))
      } else {
        unname(c(sigma_v, sigma_r, beta_0, beta_hat, delta, delta_p))
      }
    }

    ## Output label matrix
    out <- matrix(0, nrow = 3, ncol = length(start_v))
    rownames(out) <- c("par", "st_err", "t-val")
    colnames(out) <- c("sigv", "sigr", colnames(data[, x_vars_vec, drop = FALSE]), z_vars, zp_vars)

    ## Number of simulation draws
    R <- if (isTRUE(is.numeric(halton_num))) halton_num else ceiling(sqrt(nrow(data))) + 100

    ## Halton draws
    R_H <- randtoolbox::halton(
      R + .SFA_CONSTANTS$HALTON_DISCARD,
      2,
      start = 1,
      normal = FALSE
    )[-c(1:.SFA_CONSTANTS$HALTON_DISCARD), c(1:2)]

    ## First column: standard normal via qnorm
    ## Second column: half-normal via inverse error function
    R_H <- cbind(qnorm(R_H[, 1]), sqrt(2) * pracma::erfinv(R_H[, 2]))

    if (!is.null(rand.gtre)) {
      .rng_state <- .rng_snapshot()
      on.exit(.rng_restore(.rng_state), add = TRUE)
      set.seed(rand.gtre)
    }

    ## Optional decorrelation step retained from current code
    ## (sampling order unchanged so the RNG stream -- and therefore the
    ## chosen draw for a given rand.gtre seed -- is identical to before;
    ## only the correlation search is vectorized)
    mat <- matrix(0, nrow = R, ncol = 9999)
    for (v in 1:9999) {
      mat[, v] <- sample(R_H[, 1])
    }

    cor_vec <- abs(cor(mat, R_H[, 2]))

    R_H <- cbind(mat[, which.min(cor_vec)], R_H[, 2])
    rm(mat, cor_vec, v)

    ## Build firm-level lists
    indiv <- noquote(as.vector(unique(data[, individual])))
    t_vec <- rep(0, N)

    data_i <- vector("list", N)
    Y <- vector("list", N)
    data_i_vars <- vector("list", N)
    data_z_vars <- vector("list", N)
    data_zp_vars <- vector("list", N)
    R_h1 <- vector("list", N)
    R_h2 <- vector("list", N)

    for (ii in seq_len(N)) {
      data_i[[ii]] <- data[which(data[, individual] == indiv[ii]), , drop = FALSE]
      t_vec[ii] <- nrow(data_i[[ii]])
      R_h1[[ii]] <- t(matrix(rep(R_H[, 1], t_vec[[ii]]), R, t_vec[[ii]]))
      R_h2[[ii]] <- abs(t(matrix(rep(R_H[, 2], t_vec[[ii]]), R, t_vec[[ii]])))
      Y[[ii]] <- matrix(rep(data_i[[ii]][, y_var], R), t_vec[[ii]], R)
      ## Stored as plain numeric matrices (not data.frames): fn() is called
      ## thousands of times during optimization, and data.frame subsetting +
      ## as.matrix() coercion on every call was a meaningful chunk of runtime.
      data_i_vars[[ii]] <- as.matrix(data_i[[ii]][, x_vars_vec, drop = FALSE])
      data_z_vars[[ii]] <- as.matrix(data_i[[ii]][, z_vars, drop = FALSE])
      data_zp_vars[[ii]] <- as.matrix(data_i[[ii]][, zp_vars, drop = FALSE])
    }

    prep <- list(
      data         = data,                  n_z_vars     = length(z_vars),   t            = t_vec,
      individual   = individual,            n_zp_vars    = length(zp_vars),  data_i       = data_i,
      y_var        = y_var,                 N            = N,                Y            = Y,
      x_vars_vec   = x_vars_vec,            R            = R,                data_i_vars  = data_i_vars,
      z_vars       = z_vars,                start_v      = start_v,          data_z_vars  = data_z_vars,
      zp_vars      = zp_vars,               out_template = out,              data_zp_vars = data_zp_vars,
      n_x_vars     = n_x_vars,              indiv        = indiv,            R_h1         = R_h1,
      R_h2         = R_h2,
      formula      = formula
    )

    t <- t_vec
    n_z_vars <- length(z_vars)
    n_zp_vars <- length(zp_vars)
    out_template <- out

    ## GTRE simulated negative log-likelihood
    fn <- function(x) {
      beta_start <- 3
      beta_end <- beta_start + n_x_vars - 1

      delta_start <- beta_end + 1
      delta_end <- delta_start + n_z_vars - 1

      deltap_start <- delta_end + 1
      deltap_end <- deltap_start + n_zp_vars - 1

      beta <- x[beta_start:beta_end]

      delta <- x[delta_start:delta_end]
      delta_p <- x[deltap_start:deltap_end]

      sig_v <- x[1]
      sig_r <- x[2]

      ## Enforce positivity on sigma_v and sigma_r
      sig_v <- max(sig_v, .SFA_CONSTANTS$MIN_POSITIVE)
      sig_r <- max(sig_r, .SFA_CONSTANTS$MIN_POSITIVE)

      ## Bound the variance linear predictors before exponentiating.
      eta_bound <- 40

      ## One firm-level log-likelihood contribution
      ll_i <- function(ii) {
        Ti <- t[ii]

        ## Persistent inefficiency scale for firm i
        ## IMPORTANT:
        ## This is computed inside the firm-specific contribution,
        ## not outside, so it is truly unit-specific.
        sigma_h_fun <- if (n_zp_vars > 1) {
          eta_h_i <- as.numeric(data_zp_vars[[ii]] %*% delta_p)
          if (any(!is.finite(eta_h_i))) {
            return((.SFA_CONSTANTS$MAX_VALUE)^0.1)
          }
          eta_h_i <- pmin(pmax(eta_h_i, -eta_bound), eta_bound)
          mean(sqrt(exp(eta_h_i)))
        } else {
          eta_h_i <- x[deltap_end]
          if (!is.finite(eta_h_i)) {
            return((.SFA_CONSTANTS$MAX_VALUE)^0.1)
          }
          eta_h_i <- pmin(pmax(eta_h_i, -eta_bound), eta_bound)
          sqrt(exp(eta_h_i))
        }


        ## Construct epsilon_it draw-by-draw
        eps_ii <- Y[[ii]] - sig_r * R_h1[[ii]] + sigma_h_fun * R_h2[[ii]] * inefdec_n

        ## Remove frontier mean: one matrix multiply (Ti x k times k) instead
        ## of n_x_vars separate Ti x R temporary matrices/subtractions. The
        ## resulting length-Ti vector is recycled column-wise against the
        ## Ti x R matrix by R's normal recycling rules, which is exactly
        ## equivalent to matrix(rep(...), Ti, R) but without allocating it.
        eps_ii <- eps_ii - as.vector(data_i_vars[[ii]] %*% beta)

        eps_ii <- inefdec_n * eps_ii

        ## Transient inefficiency scale
        sigma_u_fun <- if (n_z_vars > 1) {
          eta_u_i <- as.numeric(data_z_vars[[ii]] %*% delta)
          if (any(!is.finite(eta_u_i))) {
            return((.SFA_CONSTANTS$MAX_VALUE)^0.1)
          }
          eta_u_i <- pmin(pmax(eta_u_i, -eta_bound), eta_bound)
          sqrt(exp(eta_u_i))
        } else {
          eta_u_i <- x[delta_end]
          if (!is.finite(eta_u_i)) {
            return((.SFA_CONSTANTS$MAX_VALUE)^0.1)
          }
          eta_u_i <- pmin(pmax(eta_u_i, -eta_bound), eta_bound)
          sqrt(exp(eta_u_i))
        }

        sigma_fun <- sqrt(sig_v^2 + sigma_u_fun^2)
        lamb_fun <- sigma_u_fun / sig_v

        ## sigma_fun/lamb_fun are length-Ti (or length-1, which recycles
        ## trivially). Dividing/multiplying the Ti x R matrix by them relies
        ## on the same column-wise recycling as above, so the explicit
        ## matrix(rep(sigma_fun, R), Ti, R) / matrix(rep(lamb_fun, R), Ti, R)
        ## allocations are no longer needed.
        sim_terms <- (2 / sigma_fun) *
          dnorm(eps_ii / sigma_fun) *
          pmax(
            pnorm(-eps_ii * lamb_fun / sigma_fun),
            .SFA_CONSTANTS$MIN_POSITIVE
          )
        prod_vec_n <- log(mean(.col_prods(sim_terms)))
        -prod_vec_n
      }

      ll_vec <- unlist(lapply(seq_len(N), ll_i))

      ll_vec[which(ll_vec == Inf)] <- (.SFA_CONSTANTS$MAX_VALUE)^0.1
      ll_vec[which(ll_vec == -Inf)] <- -(.SFA_CONSTANTS$MAX_VALUE)^0.1

      sum(ll_vec[is.finite(ll_vec)])
    }

    ## Staged optimization for GTRE
    Start.Time <- start.time()

    lower.BOB <- .generate_sfa_bounds(formula, prep) # default to -Inf

    ## ---- Stage 1: bobyqa
    Opt.Bobyqa <- opt.bobyqa(
      fn = fn,
      start_v = start_v,
      lower.bobyqa = lower.BOB,
      maxit.bobyqa = maxit.bobyqa,
      bob.TF = TRUE,
      verbose = verbose
    )

    start_v <- Opt.Bobyqa$start_v
    start_feval <- Opt.Bobyqa$start_feval
    bob1 <- Opt.Bobyqa$bob1

    ## ---- Stage 2: psoptim
    differ <- 10
    # lower1 <- c(rep(.SFA_CONSTANTS$MIN_POSITIVE, 2), start_v[-c(1:2)] - differ)
    lower1 <- .generate_sfa_bounds(formula, prep, inf_sub = min(start_v[-c(1:2)]) - differ)
    Opt.Psoptim <- opt.psoptim(
      fn = fn,
      start_v,
      lower.psoptim = lower1,
      upper.psoptim = c(start_v + differ),
      rand.psoptim = rand.psoptim,
      maxit.psoptim = maxit.psoptim,
      psopt.TF = PSopt,
      verbose = verbose
    )

    start_v <- Opt.Psoptim$start_v
    start_feval <- Opt.Psoptim$start_feval
    opt00 <- Opt.Psoptim$opt00

    ## ---- Stage 3: optim
    differ <- 1
    # lower1 <- c(rep(.SFA_CONSTANTS$MIN_POSITIVE, 2), start_v[-c(1:2)] - differ)
    lower1 <- .generate_sfa_bounds(formula, prep, inf_sub = min(start_v[-c(1:2)]) - differ)

    Opt.Optim <- opt.optim(
      fn = fn,
      start_v = start_v,
      lower.optim = lower1,
      upper.optim = c(start_v + differ),
      maxit.optim = maxit.optim,
      opt.TF = optHessian,
      method = Method,
      optHessian = optHessian,
      verbose = verbose
    )

    start_v <- Opt.Optim$start_v
    start_feval <- Opt.Optim$start_feval
    opt <- Opt.Optim$opt

    End.Time <- end.time(Start.Time)

    ## Preserve current fallback logic
    if (optHessian == FALSE && PSopt == FALSE) {
      opt <- bob1
    }

    if (optHessian == FALSE && PSopt == TRUE) {
      opt <- opt00
    }


    ## ------------------------------------------------------------
    ## Post-estimation GTRE technical efficiency recovery
    ##
    ## This uses the refactored TE code we developed earlier.
    ## ------------------------------------------------------------
    .gtre_te <- function(opt, prep, inefdec_n) {
      ## Basic dimensions and indexing
      n <- sum(t)
      id_obs <- rep(seq_len(N), t)
      t_cum <- cumsum(t)
      t_start <- c(1, head(t_cum, -1) + 1)

      ## Build variance-design matrices
      Z_mat <- .make_var_design(data, z_vars, rows = NULL, int_name = "int_u")
      Zp_mat <- .make_var_design(data, zp_vars, rows = t_start, int_name = "int_h")

      # n_z_eff  <- ncol(Z_mat)
      # n_zp_eff <- ncol(Zp_mat)

      ## Parameter indexing
      beta_start <- 3
      beta_end <- beta_start + n_x_vars - 1

      delta_start <- beta_end + 1
      delta_end <- delta_start + n_z_vars - 1

      deltap_start <- delta_end + 1
      deltap_end <- deltap_start + n_zp_vars - 1

      beta <- opt$par[beta_start:beta_end]
      delta <- opt$par[delta_start:delta_end]
      delta_p <- opt$par[deltap_start:deltap_end]

      sig_v <- max(opt$par[1], 1e-8)
      sig_r <- max(opt$par[2], 1e-8)

      min_sd <- 1e-8
      eta_bound <- 40

      ## Compute and check the variance linear predictors before exp().
      ## The earlier pmax() guard only protected against near-zero standard
      ## deviations.  It did not protect against overflow from exp(Z %*% delta),
      ## which can produce Inf values and cause the posterior covariance system
      ## to fail later in .safe_inverse().
      eta_u <- as.numeric(Z_mat %*% delta)
      eta_h <- as.numeric(Zp_mat %*% delta_p)

      if (any(!is.finite(eta_u)) || any(!is.finite(eta_h))) {
        stop(
          "Non-finite variance linear predictor in .gtre_te(). Check Z*delta and Zp*delta_p.",
          call. = FALSE
        )
      }

      ## Bound before exponentiating to prevent numerical overflow.
      eta_u <- pmin(pmax(eta_u, -eta_bound), eta_bound)
      eta_h <- pmin(pmax(eta_h, -eta_bound), eta_bound)

      sig_u_all <- pmax(sqrt(exp(eta_u)), min_sd)
      sig_h_all <- pmax(sqrt(exp(eta_h)), min_sd)

      if (any(!is.finite(sig_u_all)) || any(!is.finite(sig_h_all))) {
        stop(
          "Non-finite sigma_u or sigma_h values in .gtre_te(). Check bounded variance predictors.",
          call. = FALSE
        )
      }

      sig_u_split <- split(sig_u_all, id_obs)

      ## Residual vectors epsilon_i
      e_i <- Map(
        f = function(Yi, Xi) {
          pmin(
            inefdec_n * (Yi[, 1] - rowSums(t(t(Xi) * beta))),
            Yi[, 1] * 0
          )
        },
        Yi = Y,
        Xi = data_i_vars
      )

      ## Firm-level matrices
      A_i <- lapply(t, function(Ti) -cbind(rep(1, Ti), diag(Ti)))

      SIG <- lapply(
        t,
        function(Ti) sig_v^2 * diag(Ti) + sig_r^2 * tcrossprod(rep(1, Ti))
      )

      VEE <- Map(
        f = function(sig_h_i, sig_u_i) {
          Ti <- length(sig_u_i)
          rbind(
            c(sig_h_i^2, rep(0, Ti)),
            cbind(rep(0, Ti), diag(sig_u_i^2, nrow = Ti, ncol = Ti))
          )
        },
        sig_h_i = sig_h_all,
        sig_u_i = sig_u_split
      )

      ## Build posterior covariance objects one firm at a time so that
      ## any numerical failure reports the offending firm and panel length.
      ## This makes failures much easier to debug than a generic Map/mapply
      ## traceback.
      post_obj <- lapply(seq_len(N), function(ii) {
        if (any(!is.finite(VEE[[ii]])) || any(!is.finite(A_i[[ii]])) || any(!is.finite(SIG[[ii]]))) {
          stop(
            sprintf(
              "Non-finite input matrix in .gtre_te() for firm %s (Ti = %s).",
              ii, t[ii]
            ),
            call. = FALSE
          )
        }

        tryCatch(
          .safe_linear_combo(
            VEE_i      = VEE[[ii]],
            A_i        = A_i[[ii]],
            SIG_i      = SIG[[ii]],
            base_ridge = 1e-10,
            ridge_mult = 10,
            max_tries  = 8
          ),
          error = function(e) {
            stop(
              sprintf(
                "Failed to compute GTRE posterior matrices in .gtre_te() for firm %s (Ti = %s): %s",
                ii, t[ii], conditionMessage(e)
              ),
              call. = FALSE
            )
          }
        )
      })

      LAM <- lapply(post_obj, `[[`, "LAM")
      ARR <- lapply(post_obj, `[[`, "ARR")

      ridge_report <- data.frame(
        firm        = seq_len(N),
        ridge_VEE   = vapply(post_obj, `[[`, numeric(1), "ridge_VEE"),
        ridge_SIG   = vapply(post_obj, `[[`, numeric(1), "ridge_SIG"),
        ridge_K     = vapply(post_obj, `[[`, numeric(1), "ridge_K"),
        method_VEE  = vapply(post_obj, `[[`, character(1), "invVEE_method"),
        method_SIG  = vapply(post_obj, `[[`, character(1), "invSIG_method"),
        method_K    = vapply(post_obj, `[[`, character(1), "invK_method")
      )

      ## Persistent TE
      res_d <- mapply(
        FUN = function(Ti, ARR_i, e_i, LAM_i) {
          ptmvnorm(
            lowerx = rep(0, Ti + 1),
            upperx = rep(Inf, Ti + 1),
            mean   = as.numeric(ARR_i %*% e_i),
            sigma  = LAM_i
          )[1]
        },
        Ti = t,
        ARR_i = ARR,
        e_i = e_i,
        LAM_i = LAM,
        SIMPLIFY = FALSE
      )

      res_n <- mapply(
        FUN = function(Ti, ARR_i, e_i, LAM_i) {
          shift_vec <- c(-1, rep(0, Ti))
          ptmvnorm(
            lowerx = rep(0, Ti + 1),
            upperx = rep(Inf, Ti + 1),
            mean   = as.numeric(ARR_i %*% e_i + LAM_i %*% shift_vec),
            sigma  = LAM_i
          )[1]
        },
        Ti = t,
        ARR_i = ARR,
        e_i = e_i,
        LAM_i = LAM,
        SIMPLIFY = FALSE
      )

      H <- mapply(
        FUN = function(Ti, ARR_i, e_i, LAM_i, rd, rn) {
          shift_vec <- c(-1, rep(0, Ti))
          (max(rn, .SFA_CONSTANTS$MIN_POSITIVE) /
            max(rd, .SFA_CONSTANTS$MIN_POSITIVE)) *
            exp(
              t(shift_vec) %*% ARR_i %*% e_i +
                0.5 * t(shift_vec) %*% LAM_i %*% shift_vec
            )
        },
        Ti = t,
        ARR_i = ARR,
        e_i = e_i,
        LAM_i = LAM,
        rd = res_d,
        rn = res_n
      )

      H <- pmin(H, 1)

      ## Transient TE
      U_list <- mapply(
        FUN = function(Ti, ARR_i, e_i, LAM_i, rd) {
          ## vapply: each element is the 1x1 matrix product below, which
          ## satisfies a numeric(1) template (type and length are what vapply
          ## checks). Fixing the type here stops a zero-length Ti or a changed
          ## ptmvnorm() return from silently yielding a list.
          vapply(seq_len(Ti), function(j) {
            shift_vec <- rep(0, Ti + 1)
            shift_vec[j + 1] <- -1

            rn_t <- ptmvnorm(
              lowerx = rep(0, Ti + 1),
              upperx = rep(Inf, Ti + 1),
              mean   = as.numeric(ARR_i %*% e_i + LAM_i %*% shift_vec),
              sigma  = LAM_i
            )[1]

            (max(rn_t, .SFA_CONSTANTS$MIN_POSITIVE) /
              max(rd, .SFA_CONSTANTS$MIN_POSITIVE)) *
              exp(
                t(shift_vec) %*% ARR_i %*% e_i +
                  0.5 * t(shift_vec) %*% LAM_i %*% shift_vec
              )
          }, numeric(1))
        },
        Ti = t,
        ARR_i = ARR,
        e_i = e_i,
        LAM_i = LAM,
        rd = res_d,
        SIMPLIFY = FALSE
      )

      U <- unlist(U_list, use.names = FALSE)
      U <- pmin(U, 1)

      list(
        U = U,
        H = H,
        ridge_report = ridge_report
      )
    }


    ## ------------------------------------------------------------
    ## Finalize GTRE result object
    ##
    ## Handles:
    ##   - standard errors
    ##   - t-values
    ##   - result labeling
    ##   - sfareg return object
    ## ------------------------------------------------------------
    .gtre_finalize <- function(prep,
                               data,
                               formula,
                               model_name,
                               call,
                               U,
                               H) {
      ## Uses `opt` from the enclosing GTRE_Z block scope (assigned via the
      ## bobyqa -> psoptim -> optim scaffold above, same as every other
      ## model branch in this file) rather than taking it as a formal
      ## argument here -- this function is a nested closure defined inside
      ## that same block, so it resolves `opt` lexically. (Previously took
      ## an `opt_obj` parameter that was never actually populated -- the
      ## `.gtre_optimize()` call meant to build it was commented out below
      ## -- and was never referenced in this body anyway, which is what R
      ## CMD check's "no visible binding for global variable 'opt_obj'"
      ## NOTE was flagging. Removed rather than wired up, since `opt` is
      ## already the correct, real value.)
      out <- out_template

      ## Standard errors (Hessian-based)
      if (is.null(opt$hessian) || any(diag(opt$hessian) <= 0)) {
        st_err <- rep(NA, length(opt$par))
      } else {
        st_err <- tryCatch(
          suppressWarnings(sqrt(pmax(diag(solve(opt$hessian)), 0))),
          error = function(e) rep(NA_real_, length(opt$par))
        )
      }

      t_val <- opt$par / st_err
      out[1, ] <- opt$par
      out[2, ] <- st_err
      out[3, ] <- t_val

      ## Rename intercept labels if needed
      if (length(colnames(out)[which(colnames(out) == "(Intercept)")]) > 2) {
        colnames(out)[which(colnames(out) == "(Intercept)")] <-
          c("(Intercept x)", "(Intercept u)", "(Intercept h)")
      }

      ## ------------------------------------------------------------------
      ## OPG (outer product of gradients) standard errors
      ##
      ## The likelihood is a sum over firms:  fn(theta) = sum_i ll_i(theta)
      ## fn is deterministic conditional on the fixed Halton draws (R_h1,
      ## R_h2), so finite-difference Jacobian of the per-firm contributions
      ## is accurate without simulation noise.
      ##
      ## Procedure:
      ##   1. Build fn_vec: identical to fn but returns the length-N vector
      ##      of per-firm negative log-likelihood contributions.
      ##   2. Use numDeriv::jacobian to compute the N x p score matrix S
      ##      where S[i,j] = d(ll_i)/d(theta_j)  (negated because fn is
      ##      the *negative* log-likelihood).
      ##   3. OPG meat = t(S) %*% S
      ##   4. OPG vcov = solve(OPG_meat)
      ##   5. Sandwich vcov = H^{-1} %*% OPG_meat %*% H^{-1} (most robust;
      ##      computed only when the Hessian is available and invertible).
      ##
      ## out_opg has the same row/column layout as out:
      ##   row 1: par  (identical to out)
      ##   row 2: OPG SE
      ##   row 3: OPG t-value  (par / OPG SE)
      ## out_opg_sandwich uses sandwich SEs when available, else OPG SEs.
      ## ------------------------------------------------------------------

      out_opg <- out_sandwich <- NA

      if (OPG_calc == TRUE) {
        out_opg <- out_sandwich <- out
        ## Declared here so the sandwich block below can reuse the same meat
        ## without a second jacobian() call. Assign to it with plain `<-`:
        ## tryCatch({...}) evaluates its expression in THIS frame rather than a
        ## new one, so `<-` reaches this binding. `<<-` does not -- it starts
        ## its search one frame out, skips this binding entirely, and ends up
        ## writing to globalenv() while the local OPG_meat stays NULL.
        OPG_meat <- NULL

        opg_se <- tryCatch(
          {
            ## fn_vec: returns per-firm negative log-likelihood vector
            fn_vec <- function(x) {
              beta_start <- 3
              beta_end <- beta_start + n_x_vars - 1
              delta_start <- beta_end + 1
              delta_end <- delta_start + n_z_vars - 1
              deltap_start <- delta_end + 1
              deltap_end <- deltap_start + n_zp_vars - 1
              beta <- x[beta_start:beta_end]
              delta <- x[delta_start:delta_end]
              delta_p <- x[deltap_start:deltap_end]
              sig_v <- max(x[1], .SFA_CONSTANTS$MIN_POSITIVE)
              sig_r <- max(x[2], .SFA_CONSTANTS$MIN_POSITIVE)
              eta_bound <- 40

              ll_i_vec <- function(ii) {
                Ti <- t[ii]
                sigma_h_fun <- if (n_zp_vars > 1) {
                  eta_h_i <- as.numeric(data_zp_vars[[ii]] %*% delta_p)
                  if (any(!is.finite(eta_h_i))) {
                    return((.SFA_CONSTANTS$MAX_VALUE)^0.1)
                  }
                  eta_h_i <- pmin(pmax(eta_h_i, -eta_bound), eta_bound)
                  mean(sqrt(exp(eta_h_i)))
                } else {
                  eta_h_i <- x[deltap_end]
                  if (!is.finite(eta_h_i)) {
                    return((.SFA_CONSTANTS$MAX_VALUE)^0.1)
                  }
                  eta_h_i <- pmin(pmax(eta_h_i, -eta_bound), eta_bound)
                  sqrt(exp(eta_h_i))
                }
                eps_ii <- Y[[ii]] - sig_r * R_h1[[ii]] + sigma_h_fun * R_h2[[ii]] * inefdec_n
                eps_ii <- eps_ii - as.vector(data_i_vars[[ii]] %*% beta)
                eps_ii <- inefdec_n * eps_ii
                sigma_u_fun <- if (n_z_vars > 1) {
                  eta_u_i <- as.numeric(data_z_vars[[ii]] %*% delta)
                  if (any(!is.finite(eta_u_i))) {
                    return((.SFA_CONSTANTS$MAX_VALUE)^0.1)
                  }
                  eta_u_i <- pmin(pmax(eta_u_i, -eta_bound), eta_bound)
                  sqrt(exp(eta_u_i))
                } else {
                  eta_u_i <- x[delta_end]
                  if (!is.finite(eta_u_i)) {
                    return((.SFA_CONSTANTS$MAX_VALUE)^0.1)
                  }
                  eta_u_i <- pmin(pmax(eta_u_i, -eta_bound), eta_bound)
                  sqrt(exp(eta_u_i))
                }
                sigma_fun <- sqrt(sig_v^2 + sigma_u_fun^2)
                lamb_fun <- sigma_u_fun / sig_v
                sim_terms <- (2 / sigma_fun) *
                  dnorm(eps_ii / sigma_fun) *
                  pmax(
                    pnorm(-eps_ii * lamb_fun / sigma_fun),
                    .SFA_CONSTANTS$MIN_POSITIVE
                  )
                val <- -log(mean(.col_prods(sim_terms)))
                if (!is.finite(val)) val <- (.SFA_CONSTANTS$MAX_VALUE)^0.1
                val
              }

              ll_vec <- vapply(seq_len(N), ll_i_vec, numeric(1))
              ll_vec[!is.finite(ll_vec)] <- (.SFA_CONSTANTS$MAX_VALUE)^0.1
              ll_vec ## return the vector, not the sum
            }

            ## Score matrix: N x p, each row is per-firm gradient of the
            ## *negative* log-likelihood. Negate to get score contributions.
            score_mat <- -numDeriv::jacobian(func = fn_vec, x = opt$par)

            ## OPG meat and variance
            OPG_meat <- crossprod(score_mat) ## t(S) %*% S
            OPG_vcov <- tryCatch(
              solve(OPG_meat),
              error = function(e) {
                ## Fall back to pseudoinverse -- handles rank-deficient OPG meat
                warning("OPG matrix singular, using pseudoinverse (Moore-Penrose). ",
                  "Some SEs may be unreliable.",
                  call. = FALSE
                )
                tryCatch(MASS::ginv(OPG_meat), error = function(e2) NULL)
              }
            )

            if (is.null(OPG_vcov)) {
              warning("OPG matrix is singular; OPG standard errors set to NA.", call. = FALSE)
              rep(NA_real_, length(opt$par))
            } else {
              sqrt(pmax(diag(OPG_vcov), 0))
            }
          },
          error = function(e) {
            warning("OPG SE computation failed: ", conditionMessage(e),
              ". OPG standard errors set to NA.",
              call. = FALSE
            )
            rep(NA_real_, length(opt$par))
          }
        )

        out_opg[2, ] <- opg_se
        out_opg[3, ] <- opt$par / opg_se
        colnames(out_opg) <- colnames(out)

        ## Sandwich SE: H^{-1} OPG H^{-1}  (only when Hessian is usable)
        ## Re-uses OPG_meat already computed above -- no second Jacobian call.
        sandwich_se <- tryCatch(
          {
            if (is.null(opt$hessian) || any(diag(opt$hessian) <= 0) || all(is.na(opg_se))) {
              rep(NA_real_, length(opt$par))
            } else {
              H_inv <- tryCatch(solve(opt$hessian), error = function(e) NULL)
              if (is.null(H_inv)) {
                rep(NA_real_, length(opt$par))
              } else {
                sw_vcov <- H_inv %*% OPG_meat %*% H_inv
                sqrt(pmax(diag(sw_vcov), 0))
              }
            }
          },
          error = function(e) {
            warning("Sandwich SE computation failed: ", conditionMessage(e),
              ". Sandwich standard errors set to NA.",
              call. = FALSE
            )
            rep(NA_real_, length(opt$par))
          }
        )

        out_sandwich[2, ] <- sandwich_se
        out_sandwich[3, ] <- opt$par / sandwich_se
        colnames(out_sandwich) <- colnames(out)
      }

      total_time <- End.Time

      results <- list(
        t(out),
        c(opt),
        data,
        total_time,
        start_v,
        model_name,
        formula,
        U,
        H,
        out["par", ],
        out["st_err", ],
        out["t-val", ],
        call,
        t(out_opg),
        t(out_sandwich)
      )

      class(results) <- "sfareg"
      names(results) <- c(
        "out",
        "opt",
        "data",
        "total_time",
        "start_v",
        "model_name",
        "formula",
        "U",
        "H",
        "coefficients",
        "std.errors",
        "t.values",
        "call",
        "out_opg",
        "out_sandwich"
      )

      results
    }


    ## ------------------------------------------------------------
    ## Main internal GTRE wrapper
    ##
    ## This is the single entry point for the GTRE branch.
    ##
    ## Recommended usage:
    ##   results <- .gtre(...)
    ##
    ## Then the higher-level package code can dispatch by model_name.
    ## ------------------------------------------------------------
    .gtre <- function(data,
                      individual,
                      y_var,
                      x_vars_vec,
                      z_vars,
                      zp_vars,
                      n_x_vars,
                      beta_hat,
                      beta_0,
                      beta_0_st,
                      sigma_v,
                      sigma_r,
                      start_val,
                      halton_num,
                      rand.gtre,
                      N,
                      inefdec_n,
                      maxit.bobyqa,
                      rand.psoptim,
                      maxit.psoptim,
                      PSopt,
                      maxit.optim,
                      optHessian,
                      Method,
                      formula,
                      call,
                      verbose = FALSE) {
      ## Optimization already happened earlier in this GTRE_Z block via the
      ## standard bobyqa -> psoptim -> optim scaffold (same pattern as every
      ## other model branch in this file), producing `opt`/`bob1`/`opt00` in
      ## the enclosing scope. `.gtre()` is a nested closure defined inside
      ## that same block, so `opt` below resolves lexically to that result
      ## -- there is no separate `.gtre_optimize()` step to call here.

      ## 3. Post-estimation technical efficiency
      te_obj <- .gtre_te(
        opt       = opt,
        prep      = prep,
        inefdec_n = inefdec_n
      )

      ## 4. Finalize result object
      .gtre_finalize(
        prep       = prep,
        data       = data,
        formula    = formula,
        model_name = "GTRE_Z",
        call       = call,
        U          = te_obj$U,
        H          = te_obj$H
      )
    }

    results <- .gtre(
      data          = data,
      individual    = individual,
      y_var         = y_var,
      x_vars_vec    = x_vars_vec,
      z_vars        = z_vars,
      zp_vars       = zp_vars,
      n_x_vars      = n_x_vars,
      beta_hat      = beta_hat,
      beta_0        = beta_0,
      beta_0_st     = beta_0_st,
      sigma_v       = sigma_v,
      sigma_r       = sigma_r,
      start_val     = start_val,
      halton_num    = halton_num,
      rand.gtre     = rand.gtre,
      N             = N,
      inefdec_n     = inefdec_n,
      maxit.bobyqa  = maxit.bobyqa,
      rand.psoptim  = rand.psoptim,
      maxit.psoptim = maxit.psoptim,
      PSopt         = PSopt,
      maxit.optim   = maxit.optim,
      optHessian    = optHessian,
      Method        = Method,
      formula       = formula,
      call          = call,
      verbose       = verbose
    )
    return(results)
  }
  if (model_name == "TRE_Z") {
    delta <- rep(0.1, length(z_vars))

    if (isTRUE(is.numeric(start_val))) {
      start_v <- start_val
    }
    if (isFALSE(is.numeric(start_val))) {
      start_v <- if (is.na(beta_0_st)) {
        unname(c(sigma_v, sigma_r, beta_hat, delta))
      } else {
        unname(c(sigma_v, sigma_r, beta_0, beta_hat, delta))
      }
    }

    out <- matrix(0, nrow = 3, ncol = length(start_v))
    rownames(out) <- c("par", "st_err", "t-val")
    colnames(out) <- c("sigv", "sigr", colnames(data_x), z_vars)

    if (isTRUE(is.numeric(start_val))) {
      start_v <- start_val
    }
    if (isTRUE(is.numeric(halton_num))) {
      R <- halton_num
    } else {
      R <- ceiling(sqrt(nrow(data))) + 100
    } ## Integral reps

    R_H <- randtoolbox::halton(R + .SFA_CONSTANTS$HALTON_DISCARD, 2, start = 1, normal = FALSE)[-c(1:.SFA_CONSTANTS$HALTON_DISCARD), c(1:2)]
    R_H <- cbind(qnorm(R_H[, 1]), sqrt(2) * pracma::erfinv(R_H[, 2])) ## using inverse error function for R_H2

    if (!is.null(rand.gtre)) {
      .rng_state <- .rng_snapshot()
      on.exit(.rng_restore(.rng_state), add = TRUE)
      set.seed(rand.gtre)
    }

    mat <- matrix(0, nrow = R, ncol = 9999)
    for (v in 1:9999) {
      mat[, v] <- sample(R_H[, 1])
    }

    cor_vec <- abs(cor(mat, R_H[, 2]))

    R_H <- cbind(mat[, which.min(cor_vec)], R_H[, 2])
    rm(cor_vec, v, mat)

    # if(verbose){print(paste( "Primes 2 and 3 are in use, with 1,000 discards.  Correlation between R and H draws is:", round(cor(R_H)[1,2],10), sep = "" ),quote = FALSE) }

    indiv <- noquote(as.vector(unique(data[, c(individual)])))
    t <- rep(0, N)
    data_i <- Y <- eps <- data_i_vars <- data_z_vars <- R_h1 <- R_h2 <- as.list(rep(0, N))

    for (ii in seq_len(N)) {
      data_i[[ii]] <- data[which(data[, c(individual)] == indiv[ii]), ]
      t[ii] <- nrow(data_i[[ii]])
      R_h1[[ii]] <- t(matrix(rep(R_H[, 1], t[[ii]]), R, t[[ii]]))
      R_h2[[ii]] <- abs(t(matrix(rep(R_H[, 2], t[[ii]]), R, t[[ii]])))
      Y[[ii]] <- matrix(rep(data_i[[ii]][, y_var], R), t[[ii]], R)
      data_i_vars[[ii]] <- as.matrix(data_i[[ii]][, c(x_vars_vec), drop = FALSE])
      data_z_vars[[ii]] <- as.matrix(data_i[[ii]][, c(z_vars), drop = FALSE])
    }

    fn <- function(x) {
      x_x_vec <- x[3:as.numeric(n_x_vars + 2)]

      for (qq in seq_len(n_z_vars)) {
        v <- qq + 2 + n_x_vars
        z_z_vec[qq] <- x[v]
      }

      fn1 <- function(ii) {
        eps[[ii]] <- Y[[ii]] - x[2] * R_h1[[ii]]

        ## Single matrix multiply instead of n_x_vars separate Ti x R
        ## temporary matrices/subtractions.
        eps[[ii]] <- eps[[ii]] - as.vector(data_i_vars[[ii]] %*% x_x_vec)
        eps[[ii]] <- inefdec_n * eps[[ii]]

        sigma_u_fun <- sqrt(exp(as.vector(data_z_vars[[ii]] %*% z_z_vec)))
        sigma_v_fun <- x[1]
        sigma_fun <- sqrt(sigma_v_fun^2 + sigma_u_fun^2)
        lamb_fun <- sigma_u_fun / sigma_v_fun

        ## sigma_fun/lamb_fun are length-Ti; dividing/multiplying the
        ## Ti x R matrix by them recycles column-wise, so the
        ## matrix(rep(...), Ti, R) allocations are no longer needed.
        prod_vec_n <- log(mean(.col_prods((2 / sigma_fun) *
          dnorm(eps[[ii]] / sigma_fun) *
          pmax(pnorm(-eps[[ii]] * lamb_fun / sigma_fun), .SFA_CONSTANTS$MIN_POSITIVE))))

        return(-prod_vec_n)
      }

      fn1_apply <- unlist(lapply(1:N, fn1))

      fn1_apply[which(fn1_apply == Inf)] <- (.SFA_CONSTANTS$MAX_VALUE)^.1
      fn1_apply[which(fn1_apply == -Inf)] <- -(.SFA_CONSTANTS$MAX_VALUE)^.1

      return(sum(fn1_apply[is.finite(fn1_apply)]))
    }

    Start.Time <- start.time()

    Opt.Bobyqa <- opt.bobyqa(fn = fn, start_v = start_v, lower.bobyqa = -Inf, maxit.bobyqa = maxit.bobyqa, bob.TF = TRUE, verbose = verbose)
    start_v <- Opt.Bobyqa$start_v
    start_feval <- Opt.Bobyqa$start_feval
    bob1 <- Opt.Bobyqa$bob1

    differ <- 10
    lower1 <- c(rep(.SFA_CONSTANTS$MIN_POSITIVE, 2), start_v[-c(1:2)] - differ)

    Opt.Psoptim <- opt.psoptim(
      fn = fn, start_v, lower.psoptim = lower1, upper.psoptim = c(start_v + differ),
      rand.psoptim = rand.psoptim,
      maxit.psoptim = maxit.psoptim, psopt.TF = PSopt, verbose = verbose
    )
    start_v <- Opt.Psoptim$start_v
    start_feval <- Opt.Psoptim$start_feval
    opt00 <- Opt.Psoptim$opt00

    differ <- 1
    lower1 <- c(rep(.SFA_CONSTANTS$MIN_POSITIVE, 2), start_v[-c(1:2)] - differ)

    Opt.Optim <- opt.optim(
      fn = fn, start_v = start_v, lower.optim = lower1, upper.optim = c(start_v + differ),
      maxit.optim = maxit.optim, opt.TF = optHessian, method = Method, optHessian = optHessian, verbose = verbose
    )
    start_v <- Opt.Optim$start_v
    start_feval <- Opt.Optim$start_feval
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

    st_err <- if (isTRUE(any(opt$hessian == 0) | optHessian == FALSE)) {
      rep(NA, length(opt$par))
    } else {
      suppressWarnings(sqrt(pmax(diag(solve(opt$hessian)), 0)))
    }
    t_val <- opt$par / st_err
    out[1, ] <- opt$par
    out[2, ] <- st_err
    out[3, ] <- t_val
    if (length(colnames(out)[which(colnames(out) == "(Intercept)")]) == 2) {
      colnames(out)[which(colnames(out) == "(Intercept)")] <- c("(Intercept x)", "(Intercept u)")
    }

    ## TE Measurements
    NX <- n_x_vars + 2
    NZ1 <- n_x_vars + 3
    NZ2 <- n_x_vars + n_z_vars + 2

    beta <- opt$par[c(3:NX)]
    delta <- opt$par[c(NZ1:NZ2)]

    sig_v <- opt$par[1]
    sig_u <- sqrt(exp((as.matrix(data.frame(data[, c(z_vars)]))) %*% delta))
    lamb <- sig_u / sig_v
    sig <- sqrt(sig_u^2 + sig_v^2)

    eps_hat <- pmin(inefdec_n * (data[, y_var] - rowSums(t(t(data.frame(data[, c(x_vars_vec)])) * beta))), 5)

    sig_star <- (sig_u * sig_v) / sig
    inner <- (lamb * eps_hat) / sig
    U <- ((1 - pnorm((sig_u * sig_v / sig) + inner)) / pmax(1 - pnorm(inner), .SFA_CONSTANTS$MIN_POSITIVE)) * exp((sig_u^2 / sig^2) * (eps_hat + 0.5 * sig_v^2))
    U <- pmax(U, 0)
    U <- pmin(U, 1)

    ## Results
    results <- list(t(out), c(opt), data, End.Time, start_v, model_name, formula, U, out["par", ], out["st_err", ], out["t-val", ], call)
    class(results) <- "sfareg"
    names(results) <- c("out", "opt", "data", "total_time", "start_v", "model_name", "formula", "U", "coefficients", "std.errors", "t.values", "call")
    return(results)
  }
  if (model_name == "GTRE_FML") {
    ## Four-component GTRE by FULL information maximum likelihood, via the
    ## model's closed-skew-normal representation. See .csn_gtre_loglik() in
    ## matrix_utils.R for the algebra; ported from `base code/FML.R`.
    ##
    ## Contrast with model_name = "GTRE", which estimates the same model by
    ## SIMULATED maximum likelihood over Halton draws. This branch integrates
    ## nothing by Monte Carlo, so the likelihood is deterministic and the
    ## reported log-likelihood is directly comparable across runs. The price is
    ## one (T+1)-dimensional normal CDF per firm per evaluation: cost grows
    ## sharply in T, and this is much slower than "GTRE" for long panels.
    ##
    ## Parameters are reported in the order the likelihood consumes them, which
    ## keeps the mapping to FML.R transparent:
    ##   (Intercept), x-vars, sigr, sigv, sigh, sigu
    ## sigr/sigv are the two-sided components (individual effect, noise) and
    ## sigh/sigu the one-sided ones (time-invariant and time-varying
    ## inefficiency). These are RAW standard deviations, not the lambda/sigma
    ## reparameterization "GTRE" reports.
    Start.Time <- start.time()

    ## The CSN derivation fixes T across firms (A, V and Sigma are built once at
    ## a single dimension), so an unbalanced panel would silently be wrong
    ## rather than merely unsupported.
    id_chr <- as.character(data[, individual])
    indiv <- unique(id_chr)
    gid <- match(id_chr, indiv)
    t_i <- as.numeric(table(gid))
    if (length(unique(t_i)) != 1L) {
      stop("psfm(model_name = \"GTRE_FML\") requires a BALANCED panel: the ",
        "closed-skew-normal representation is built at a single T. Found ",
        "T ranging from ", min(t_i), " to ", max(t_i), ". Use model_name = ",
        "\"GTRE\" (simulated ML), which handles unbalanced panels.",
        call. = FALSE
      )
    }
    BigT <- t_i[1]

    Xmat <- as.matrix(data[, setdiff(x_vars_vec, "(Intercept)"), drop = FALSE])
    Y_vec <- inefdec_n * as.numeric(data[, y_var])
    X_s <- inefdec_n * Xmat
    Kx <- ncol(Xmat)

    ## Starting values. Two candidates are built and the better one is chosen by
    ## the likelihood itself, below, once like.fml exists.
    ##
    ## Candidate 1 is the random-effects decomposition start_panel() produces for
    ## the whole GTRE family. NOTE the intercept: this branch used to seed it at
    ## `beta_0_st`, the RAW panel-regression intercept, while GTRE and TRE seed
    ## theirs at `beta_0 = beta_0_st + exp_u + exp_eta`. In a composed-error
    ## model the regression intercept sits below the frontier by
    ## E[u] + E[h] = (sigma_u + sigma_h)*sqrt(2/pi) -- 1.12 at the standard test
    ## DGP -- so FML was starting a full unit low, in exactly the direction of
    ## the sigma_h = 0 boundary optimum where the intercept absorbs E[h]. Seeding
    ## at `beta_0` puts it where the other GTRE models start.
    ##
    ## Candidate 2 is the two-step moment estimator, which is what Colombi (2010,
    ## sec. 3) and Colombi, Martini and Vittadini (2011, sec. 3.2) recommend as
    ## the starting point for this exact likelihood.
    if (isTRUE(is.numeric(start_val))) {
      start_cands <- list(start_val)
    } else {
      b0_re <- if (is.na(beta_0)) mean(Y_vec) else unname(beta_0)
      cand_re <- unname(c(
        b0_re, beta_hat[seq_len(Kx)],
        max(sigma_r, 0.05), max(sigma_v, 0.05),
        max(sigma_h, 0.05), max(sigma_u, 0.05)
      ))

      cand_ts <- tryCatch(
        {
          ts <- .gtre_two_step(epsilon_hat, alpha_hat, beta_0_st)
          b0 <- if (is.na(ts$beta_0)) mean(Y_vec) else unname(ts$beta_0)
          unname(c(
            b0, beta_hat[seq_len(Kx)],
            sqrt(max(ts$sigmaSq_hr * (1 - ts$gamma_hr), 2.5e-3)),
            sqrt(max(ts$sigmaSq_uv * (1 - ts$gamma_uv), 2.5e-3)),
            sqrt(max(ts$sigmaSq_hr * ts$gamma_hr, 2.5e-3)),
            sqrt(max(ts$sigmaSq_uv * ts$gamma_uv, 2.5e-3))
          ))
        },
        error = function(e) NULL
      )

      start_cands <- Filter(
        function(z) !is.null(z) && all(is.finite(z)),
        list(cand_re, cand_ts)
      )
      if (!length(start_cands)) start_cands <- list(cand_re)
    }
    start_cands <- lapply(start_cands, function(z) {
      z[!is.finite(z)] <- 0.1
      z
    })
    start_v <- start_cands[[1]]

    out <- matrix(0, nrow = 3, ncol = length(start_v))
    rownames(out) <- c("par", "st_err", "t-val")
    colnames(out) <- c(
      "(Intercept)", setdiff(x_vars_vec, "(Intercept)"),
      "sigr", "sigv", "sigh", "sigu"
    )

    ## Quadrature nodes for the rank-one CDF reduction, built once per fit.
    gh_fml <- .gauss_hermite_nodes(64L)

    like.fml <- function(x) {
      .csn_gtre_loglik(x,
        Y = Y_vec, X = X_s, gid = gid,
        ngroups = N, BigT = BigT, gh = gh_fml
      )
    }

    ## Pick between the random-effects and two-step starts on the likelihood
    ## itself. Two evaluations, against a fit costing thousands, and it cannot
    ## do worse than either candidate alone. The candidates agree on most draws
    ## and the choice is then immaterial; it matters on the minority where the
    ## random-effects start would otherwise walk into the sigma_h = 0 corner.
    if (length(start_cands) > 1L) {
      .cand_obj <- vapply(
        start_cands,
        function(z) tryCatch(like.fml(z), error = function(e) Inf),
        numeric(1)
      )
      if (any(is.finite(.cand_obj))) start_v <- start_cands[[which.min(.cand_obj)]]
      fml_starts <- list(
        n_tried = length(start_cands),
        loglik = -.cand_obj,
        chosen = c("random-effects", "two-step")[which.min(.cand_obj)]
      )
    } else {
      fml_starts <- NULL
    }

    ## Phased start, as in the GTRE/TRE branch: the four variance components
    ## are the last four entries here (sigr, sigv, sigh, sigu), the frontier
    ## block the first 1 + Kx. See .phased_start() in matrix_utils.R.
    if (isTRUE(getOption("sfa.phased_start", FALSE))) {
      ## GTRE_FML layout is (beta_0, x..., sigr, sigv, sigh, sigu): variance
      ## block last, intercept first.
      .eff_idx_fml <- (length(start_v) - 3L):length(start_v)
      start_v <- .phased_start(
        fn = like.fml, start_v = start_v, idx = .eff_idx_fml,
        int_idx = 1L, lower = .SFA_CONSTANTS$MIN_POSITIVE,
        grid = c(0.5, 0.75, 1.5, 2), verbose = verbose
      )
    }

    lower_f <- c(rep(-Inf, 1 + Kx), rep(.SFA_CONSTANTS$MIN_POSITIVE, 4))
    upper_f <- rep(Inf, length(start_v))
    start_v <- pmin(pmax(start_v, lower_f), upper_f)

    ## Deterministic likelihood, so a quasi-Newton method is well behaved here;
    ## bobyqa follows as a derivative-free safety net. psoptim is skipped
    ## unless asked for -- each objective evaluation costs N multivariate normal
    ## CDFs, which makes a particle swarm prohibitively expensive.
    Opt.Nlminb <- opt.nlminb(
      fn = like.fml, start_v = start_v, lower.nlminb = lower_f,
      upper.nlminb = upper_f, maxit.nlminb = maxit.nlminb, nlminb.TF = TRUE,
      verbose = verbose
    )
    start_v <- Opt.Nlminb$start_v

    Opt.Bobyqa <- opt.bobyqa(
      fn = like.fml, start_v = start_v, lower.bobyqa = lower_f,
      upper.bobyqa = upper_f, maxit.bobyqa = maxit.bobyqa,
      bob.TF = TRUE, verbose = verbose
    )
    start_v <- Opt.Bobyqa$start_v
    bob1 <- Opt.Bobyqa$bob1

    differ <- 2
    Opt.Psoptim <- opt.psoptim(
      fn = like.fml, start_v,
      lower.psoptim = c(
        start_v[seq_len(1 + Kx)] - differ,
        rep(.SFA_CONSTANTS$MIN_POSITIVE, 4)
      ),
      upper.psoptim = start_v + differ, rand.psoptim = rand.psoptim,
      maxit.psoptim = maxit.psoptim, psopt.TF = PSopt, verbose = verbose
    )
    start_v <- Opt.Psoptim$start_v
    opt00 <- Opt.Psoptim$opt00

    Opt.Optim <- opt.optim(
      fn = like.fml, start_v = start_v, lower.optim = lower_f,
      upper.optim = pmin(start_v + differ, upper_f),
      maxit.optim = maxit.optim, opt.TF = optHessian, method = Method,
      optHessian = optHessian, verbose = verbose
    )
    start_v <- Opt.Optim$start_v
    opt <- Opt.Optim$opt

    End.Time <- end.time(Start.Time)

    if (optHessian == FALSE & PSopt == FALSE) {
      opt <- bob1
    }
    if (optHessian == FALSE & PSopt == TRUE) {
      opt <- opt00
    }

    st_err <- if (isTRUE(any(opt$hessian == 0)) | optHessian == FALSE) {
      rep(NA, length(opt$par))
    } else {
      tryCatch(suppressWarnings(sqrt(pmax(diag(solve(opt$hessian)), 0))),
        error = function(e) rep(NA_real_, length(opt$par))
      )
    }
    out[1, ] <- opt$par
    out[2, ] <- st_err
    out[3, ] <- opt$par / st_err

    ## Transient (time-varying) and persistent (time-invariant) efficiency,
    ## each from the corresponding one-sided component against the total
    ## two-sided scale, following the same Battese-Coelli predictor the rest of
    ## the package uses.
    np <- length(opt$par)
    sig_r <- opt$par[np - 3]
    sig_v <- opt$par[np - 2]
    sig_h <- opt$par[np - 1]
    sig_u <- opt$par[np]
    eps <- as.numeric(Y_vec - opt$par[1] - X_s %*% opt$par[2:(1 + Kx)])
    s2u <- sig_u^2
    s2v <- sig_v^2
    U <- .te_battese_coelli(
      mu_star = -eps * s2u / (s2u + s2v),
      sigma_star = sig_u * sig_v / sqrt(s2u + s2v)
    )
    e_bar <- as.numeric(.gsum(eps, gid, N)) / BigT
    s2h <- sig_h^2
    s2r <- sig_r^2 + s2v / BigT
    H <- .te_battese_coelli(
      mu_star = -e_bar * s2h / (s2h + s2r),
      sigma_star = sig_h * sqrt(s2r) / sqrt(s2h + s2r)
    )
    names(H) <- indiv

    results <- list(
      t(out), c(opt), End.Time, start_v, U, H, model_name, formula, data,
      out["par", ], out["st_err", ], out["t-val", ], call, fml_starts
    )
    class(results) <- "sfareg"
    names(results) <- c(
      "out", "opt", "total_time", "start_v", "U", "H", "model_name", "formula", "data",
      "coefficients", "std.errors", "t.values", "call", "start_search"
    )
    return(results)
  }
  if (model_name == "TFE_WMLE") {
    Start.Tfe <- start.tfe(formula_x, data, model_name, start_val, intercept, x_vars_vec, gamma, individual, N, y_var, n_x_vars)

    data_i <- Start.Tfe$data_i
    data_i_vars <- Start.Tfe$data_i_vars
    data_i_vars_dm <- Start.Tfe$data_i_vars_dm
    eps <- Start.Tfe$eps
    I_t <- Start.Tfe$I_t
    I_t1 <- Start.Tfe$I_t1
    indiv <- Start.Tfe$indiv
    one_t <- Start.Tfe$one_t
    one_t1 <- Start.Tfe$one_t1
    out <- Start.Tfe$out
    start_v <- Start.Tfe$start_v
    t <- Start.Tfe$t
    upper <- Start.Tfe$upper
    Y <- Start.Tfe$Y

    ## GH nodes for the rank-one CDF shortcut below don't depend on x --
    ## computed once per model fit rather than once per likelihood eval.
    gh_tfe <- .gauss_hermite_nodes(64L)

    like.tfe <- function(x) {
      x_x_vec <- x[3:as.numeric(n_x_vars + 2)]
      lambda_eff <- if (gamma == FALSE) x[1] else sqrt(x[1] / (1 - x[1]))

      ## optim()'s numerical Hessian (hessian=TRUE) perturbs par via
      ## unconstrained finite differences that ignore lower/upper entirely
      ## (see stats:::C_optimhess) -- so even with the gamma_cap bound below
      ## on the optimizer's own search, a solution close to that boundary
      ## can still see x[1] > 1 during Hessian evaluation, making lambda_eff
      ## NaN. Guard here rather than only at the bounds: return a large
      ## finite penalty (matching the convention in the CSW reference
      ## implementation this was ported from) instead of letting NaN
      ## propagate into .log_mvn_cdf_rank1()'s comparisons and crash.
      if (!is.finite(lambda_eff)) {
        return(1e12)
      }

      fn1 <- function(i) {
        eps_t <- Y[[i]]

        ## Was demean(as.numeric(data_i_vars[[i]][,qq])) recomputed here on
        ## every likelihood evaluation -- data_i_vars_dm precomputes this
        ## once in start.tfe() since it doesn't depend on x.
        for (qq in seq_len(n_x_vars)) {
          eps_t <- eps_t - x_x_vec[qq] * data_i_vars_dm[[i]][, qq]
        }
        eps_t <- eps_t * inefdec_n
        eps_t1 <- eps_t[1:t[[i]] - 1]

        ## l1/l2 use covariances with special rank-one/equicorrelated
        ## structure (Sigma = I + c*11' for l2, sigma2*(I - (1/T)*11') for
        ## l1). .log_within_mvn_density()/.log_mvn_cdf_rank1() (matrix_utils.R)
        ## are an exact closed form / Gauss-Hermite-quadrature evaluation of
        ## exactly this structure -- ~30x faster than a general
        ## Genz-algorithm multivariate normal CDF (e.g. mnormt::pmnorm())
        ## and (for l2) deterministic rather than stochastic.
        log_l1 <- .log_within_mvn_density(eps_t1, x[2] * x[2])
        log_l2 <- .log_mvn_cdf_rank1(upper = -(lambda_eff / x[2]) * eps_t, c = (lambda_eff * lambda_eff) / t[[i]], gh = gh_tfe)

        prod_vec_n <- log_l1 + log_l2 ## Log likelihood
        return(-prod_vec_n)
      }
      fn1_apply <- unlist(lapply(1:N, fn1))
      return(sum(fn1_apply[is.finite(fn1_apply)]))
    }

    Start.Time <- start.time()

    differ <- 2

    ## x[1] is gamma = sigma_u^2/sigma^2 in (0,1) when gamma==TRUE, vs. an
    ## unbounded-above lambda = sigma_u/sigma_v when gamma==FALSE. Without an
    ## upper bound, bobyqa's search (which otherwise defaults to upper=Inf)
    ## can and does cross gamma=1, making sqrt(x[1]/(1-x[1])) (in like.tfe())
    ## produce NaN and crash the CDF/density evaluation. gamma_cap=Inf when
    ## gamma==FALSE leaves that case unbounded, as intended.
    gamma_cap <- if (gamma) 1 - .SFA_CONSTANTS$MIN_POSITIVE else Inf

    Opt.Bobyqa <- opt.bobyqa(
      fn = like.tfe, start_v = start_v, lower.bobyqa = c(rep(.SFA_CONSTANTS$MIN_POSITIVE, 2), rep(-Inf, n_x_vars)),
      upper.bobyqa = c(gamma_cap, Inf, rep(Inf, n_x_vars)),
      maxit.bobyqa = maxit.bobyqa, bob.TF = TRUE, verbose = verbose
    )
    start_v <- Opt.Bobyqa$start_v
    start_feval <- Opt.Bobyqa$start_feval
    bob1 <- Opt.Bobyqa$bob1

    upper_ps <- pmin(start_v + differ, c(gamma_cap, rep(Inf, 1 + n_x_vars)))
    Opt.Psoptim <- opt.psoptim(
      fn = like.tfe, start_v, lower.psoptim = c(rep(.SFA_CONSTANTS$MIN_POSITIVE, 2), start_v[-c(1:2)] - differ),
      rand.psoptim = rand.psoptim,
      upper.psoptim = upper_ps, maxit.psoptim = maxit.psoptim, psopt.TF = PSopt, verbose = verbose
    )
    start_v <- Opt.Psoptim$start_v
    start_feval <- Opt.Psoptim$start_feval
    opt00 <- Opt.Psoptim$opt00

    upper_opt <- pmin(start_v + differ, c(gamma_cap, rep(Inf, 1 + n_x_vars)))
    Opt.Optim <- opt.optim(
      fn = like.tfe, start_v = start_v, lower.optim = c(rep(.SFA_CONSTANTS$MIN_POSITIVE, 2), rep(-Inf, n_x_vars)),
      upper.optim = upper_opt, maxit.optim = maxit.optim, opt.TF = optHessian, method = Method, optHessian = optHessian, verbose = verbose
    )
    start_v <- Opt.Optim$start_v
    start_feval <- Opt.Optim$start_feval
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

    beta <- opt$par[-c(1:2)] ## estimate the r_i's and then e(exp(u)|eps)'s
    lamb <- opt$par[1]
    sig <- opt$par[2]
    sig_u <- (lamb * sig) / sqrt(1 + lamb^2)
    sig_v <- sig_u / lamb
    Y_mean <- rep(0, N)
    X_mean <- matrix(0, N, ncol = length(x_vars_vec))
    colnames(X_mean) <- x_vars_vec

    ## data_i[[ii]] is already the correct per-individual subset from
    ## start.tfe() (same `indiv` ordering) -- was being recomputed from
    ## scratch here via an identical data[which(...)] subset.
    for (ii in seq_len(N)) {
      Y_mean[ii] <- mean(as.numeric(data_i[[ii]][, y_var]))
      X_mean[ii, ] <- colMeans(data.frame(data_i[[ii]][, c(x_vars_vec)]))
    }

    r_hat_m <- Y_mean - rowSums(t(beta * t(X_mean))) + inefdec_n * sqrt(2 / pi) * sig_u
    r_hat_m_exp <- rep(0, sum(t))
    t_cum <- c(cumsum(t))

    for (m in 1:length(t)) {
      B <- t_cum[m]
      A <- B + 1 - t[m]
      r_hat_m_exp[A:B] <- rep(r_hat_m[m], t[m])
    }

    eps_hat <- pmin(inefdec_n * (data[, y_var] - rowSums(t(t(data[, c(x_vars_vec)]) * beta)) - r_hat_m_exp), data[, y_var] * 0)
    sig_star <- sig_u * sig_v / sig
    inner <- (lamb * eps_hat) / sig
    exp_u_hat <- ((1 - pnorm((sig_u * sig_v / sig) + inner)) / (1 - pnorm(inner))) * exp((sig_u^2 / sig^2) * (eps_hat + 0.5 * sig_v^2))

    st_err <- if (isTRUE(any(opt$hessian == 0)) | optHessian == FALSE) {
      rep(NA, length(opt$par))
    } else {
      suppressWarnings(sqrt(pmax(diag(solve(opt$hessian)), 0)))
    }
    t_val <- opt$par / st_err
    out[1, ] <- opt$par
    out[2, ] <- st_err
    out[3, ] <- t_val
    results <- list(t(out), c(opt), End.Time, start_v, r_hat_m, exp_u_hat, model_name, formula, data, out["par", ], out["st_err", ], out["t-val", ], call)
    class(results) <- "sfareg"
    names(results) <- c("out", "opt", "total_time", "start_v", "r_hat_m", "exp_u_hat", "model_name", "formula", "data", "coefficients", "std.errors", "t.values", "call")
    return(results)
  }
  if (model_name == "TFE") {
    ## Greene's (2005) true fixed effects stochastic frontier: the ordinary
    ## normal/half-normal composed-error likelihood with one intercept per
    ## individual estimated jointly with (lambda, sigma, beta). See
    ## .tfe_alpha_profile() in matrix_utils.R for the likelihood, why the N
    ## firm effects are concentrated out rather than handed to the optimizer,
    ## and how they are solved for.
    ##
    ## Contrast with "TFE_WMLE" (Chen, Schmidt and Wang, 2014), which
    ## estimates the same MODEL but from the within-transformed data, so the
    ## firm effects never enter the likelihood. TFE is the estimator with the
    ## incidental parameters problem -- at small T its sigma_u is biased
    ## upward and its efficiency scores are correspondingly pessimistic --
    ## which is exactly the bias TFE_WMLE was designed to remove. Both are
    ## provided because TFE remains the benchmark that applied work reports.
    ##
    ## Starting values, the parameter ordering and the `out` column names are
    ## shared with TFE_WMLE via start.tfe(), so psfm_bootstrap() and the
    ## sfareg methods need no separate case for this model.
    Start.Tfe <- start.tfe(formula_x, data, model_name, start_val, intercept, x_vars_vec, gamma, individual, N, y_var, n_x_vars)

    out <- Start.Tfe$out
    start_v <- Start.Tfe$start_v

    ## Work in the production orientation throughout: e = inefdec_n*(y - x'b) - alpha
    ## is v - u whichever frontier was requested, so one likelihood covers both.
    ## The firm effect on the original scale is inefdec_n*alpha (recovered below).
    Xmat <- as.matrix(data[, x_vars_vec, drop = FALSE])
    y_s <- inefdec_n * as.numeric(data[, y_var])
    X_s <- inefdec_n * Xmat

    ## Group id built by match() against unique(), NOT by assuming the rows
    ## arrive sorted by individual -- see .gsum()'s note on why the obvious
    ## rowsum()/tapply() shortcuts silently misassign firm effects when the
    ## individual labels are numeric-looking strings.
    id_chr <- as.character(data[, individual])
    indiv <- unique(id_chr)
    gid <- match(id_chr, indiv)

    like.tfe.greene <- function(x) {
      lambda_eff <- if (gamma == FALSE) x[1] else sqrt(x[1] / (1 - x[1]))
      sig <- x[2]
      ## Same guard as like.tfe(): optim()'s numerical Hessian steps outside
      ## the optimizer's own bounds, so gamma can exceed 1 there.
      if (!is.finite(lambda_eff) || !is.finite(sig) || lambda_eff <= 0 || sig <= 0) {
        return(1e12)
      }

      r <- y_s - as.vector(X_s %*% x[3:as.numeric(n_x_vars + 2)])
      alpha <- .tfe_alpha_profile(r, gid, N, lambda_eff, sig)
      if (is.null(alpha)) {
        return(1e12)
      }

      z <- (r - alpha[gid]) / sig
      ## Evaluated entirely in logs: pnorm(log.p = TRUE) stays accurate for
      ## arguments far into the lower tail, where the raw CDF underflows to 0
      ## and log(0) would hand the optimizer an -Inf.
      ll <- sum(log(2) - log(sig) + dnorm(z, log = TRUE) +
        pnorm(-lambda_eff * z, log.p = TRUE))
      if (!is.finite(ll)) {
        return(1e12)
      }
      return(-ll)
    }

    Start.Time <- start.time()

    differ <- 2

    ## ---- the sigma_v -> 0 degeneracy, and why lambda is capped -----------
    ## Because alpha_i is unrestricted, this likelihood ALWAYS has a
    ## supremum on the sigma_v = 0 (lambda -> Inf) boundary: setting
    ## alpha_i = max_t(y_it - x_it'beta) makes every composed error weakly
    ## negative, Phi(-lambda e/sigma) -> 1, and what is left is the
    ## deterministic-frontier likelihood, which is finite and attained.
    ## Verified numerically on four simulated panels (T = 10 and 20,
    ## lambda_true from 0.5 to 3.33): in every one, profiling the likelihood
    ## over lambda at its own best sigma dips through a barrier and then
    ## rises again, ending HIGHER at lambda = 1e7 than at the interior local
    ## maximum. When lambda_true is high and T small the profile can be
    ## monotone increasing throughout, i.e. no interior maximum exists at
    ## all. The interior local mode -- the one reached by a local ascent from
    ## the within estimator, which is what the applied literature reports and
    ## what recovers the true parameters here -- is therefore the target, and
    ## an unbounded search is a liability rather than thoroughness.
    ##
    ## tfe_lambda_max keeps the search in that basin. A fit that pins at the
    ## cap has run into the degeneracy, and says so (warning below) rather
    ## than quietly reporting sigma_v ~ 0 with meaningless standard errors.
    if (!is.numeric(tfe_lambda_max) || length(tfe_lambda_max) != 1 ||
      !is.finite(tfe_lambda_max) || tfe_lambda_max <= 0) {
      stop("tfe_lambda_max must be a single finite positive number.", call. = FALSE)
    }

    ## Same cap expressed in whichever parameterization x[1] carries.
    gamma_cap <- if (gamma) {
      min(
        tfe_lambda_max^2 / (1 + tfe_lambda_max^2),
        1 - .SFA_CONSTANTS$MIN_POSITIVE
      )
    } else {
      tfe_lambda_max
    }
    lower_g <- c(rep(.SFA_CONSTANTS$MIN_POSITIVE, 2), rep(-Inf, n_x_vars))
    upper_g <- c(gamma_cap, Inf, rep(Inf, n_x_vars))

    ## start.tfe() derives its lambda from a within regression that knows
    ## nothing about the cap, so on a low-noise panel (or a deliberately tight
    ## tfe_lambda_max) the starting value can already sit above the bound --
    ## which bobyqa rejects outright with "Starting values violate bounds"
    ## rather than projecting onto the feasible set. Clamp first.
    start_v <- pmin(pmax(start_v, lower_g), upper_g)

    ## nlminb leads the stack here (as in the PL80/BC92 branch): the profile
    ## likelihood is smooth -- alpha is solved by Newton to machine precision
    ## rather than by a nested optimizer stopping at its own tolerance -- so a
    ## quasi-Newton method uses it well, and it gets the derivative-free
    ## stages started from a far better point than the within estimator alone.
    Opt.Nlminb <- opt.nlminb(
      fn = like.tfe.greene, start_v = start_v, lower.nlminb = lower_g,
      upper.nlminb = upper_g, maxit.nlminb = maxit.nlminb, nlminb.TF = TRUE,
      verbose = verbose
    )
    start_v <- Opt.Nlminb$start_v

    Opt.Bobyqa <- opt.bobyqa(
      fn = like.tfe.greene, start_v = start_v, lower.bobyqa = lower_g,
      upper.bobyqa = upper_g, maxit.bobyqa = maxit.bobyqa,
      bob.TF = TRUE, verbose = verbose
    )
    start_v <- Opt.Bobyqa$start_v
    start_feval <- Opt.Bobyqa$start_feval
    bob1 <- Opt.Bobyqa$bob1

    upper_ps <- pmin(start_v + differ, upper_g)
    Opt.Psoptim <- opt.psoptim(
      fn = like.tfe.greene, start_v, lower.psoptim = c(rep(.SFA_CONSTANTS$MIN_POSITIVE, 2), start_v[-c(1:2)] - differ),
      rand.psoptim = rand.psoptim,
      upper.psoptim = upper_ps, maxit.psoptim = maxit.psoptim, psopt.TF = PSopt, verbose = verbose
    )
    start_v <- Opt.Psoptim$start_v
    start_feval <- Opt.Psoptim$start_feval
    opt00 <- Opt.Psoptim$opt00

    upper_opt <- pmin(start_v + differ, upper_g)
    Opt.Optim <- opt.optim(
      fn = like.tfe.greene, start_v = start_v, lower.optim = lower_g,
      upper.optim = upper_opt, maxit.optim = maxit.optim, opt.TF = optHessian,
      method = Method, optHessian = optHessian, verbose = verbose
    )
    start_v <- Opt.Optim$start_v
    start_feval <- Opt.Optim$start_feval
    opt <- Opt.Optim$opt

    End.Time <- end.time(Start.Time)

    if (optHessian == FALSE & PSopt == FALSE) {
      opt <- bob1
    }
    if (optHessian == FALSE & PSopt == TRUE) {
      opt <- opt00
    }

    beta <- opt$par[-c(1:2)]
    lamb <- if (gamma == FALSE) opt$par[1] else sqrt(opt$par[1] / (1 - opt$par[1]))
    sig <- opt$par[2]
    sig_u <- (lamb * sig) / sqrt(1 + lamb^2)
    sig_v <- sig_u / lamb

    ## Pinned at the cap => the search left the interior basin for the
    ## deterministic-frontier boundary described above. Report it: the
    ## returned lambda is then the constraint, not an estimate, and the
    ## numerical Hessian at a bound gives standard errors that mean nothing.
    if (isTRUE(lamb >= tfe_lambda_max * (1 - 1e-6))) {
      warning("psfm(model_name = \"TFE\"): lambda converged to its upper bound (",
        format(tfe_lambda_max), "), i.e. sigma_v -> 0. Greene's true fixed ",
        "effects likelihood always has a supremum on that boundary (the ",
        "deterministic frontier alpha_i = max_t(y_it - x_it'beta)), and for ",
        "these data no interior maximum was found below the bound. Treat the ",
        "reported lambda as the constraint rather than an estimate; consider ",
        "model_name = \"TFE_WMLE\", which is not subject to this degeneracy.",
        call. = FALSE
      )
    }

    ## Firm effects at the final estimates. Unlike TFE_WMLE's r_hat_m these
    ## need no sqrt(2/pi)*sig_u mean adjustment: CSW's within likelihood
    ## cannot identify alpha_i and recovers it from the mean residual (which
    ## carries -E[u] and has to be shifted back), whereas here alpha_i is
    ## estimated by maximum likelihood as the frontier intercept directly.
    r <- y_s - as.vector(X_s %*% beta)
    alpha_hat <- .tfe_alpha_profile(r, gid, N, lamb, sig)
    if (is.null(alpha_hat)) alpha_hat <- rep(NA_real_, N)
    r_hat_m <- inefdec_n * alpha_hat
    names(r_hat_m) <- indiv

    eps_hat <- r - alpha_hat[gid]
    exp_u_hat <- .te_battese_coelli(
      mu_star = -eps_hat * sig_u^2 / sig^2,
      sigma_star = sig_u * sig_v / sig
    )
    u_hat <- .jlms_u(
      mu_star = -eps_hat * sig_u^2 / sig^2,
      sigma_star = sig_u * sig_v / sig
    )

    ## solve() on a singular Hessian errors outright; a fit that converged but
    ## cannot be inverted should still return its point estimates with NA
    ## standard errors rather than aborting the whole call.
    st_err <- if (isTRUE(any(opt$hessian == 0)) | optHessian == FALSE) {
      rep(NA, length(opt$par))
    } else {
      tryCatch(suppressWarnings(sqrt(pmax(diag(solve(opt$hessian)), 0))),
        error = function(e) rep(NA_real_, length(opt$par))
      )
    }
    t_val <- opt$par / st_err
    out[1, ] <- opt$par
    out[2, ] <- st_err
    out[3, ] <- t_val

    results <- list(
      t(out), c(opt), End.Time, start_v, r_hat_m, exp_u_hat, u_hat, model_name,
      formula, data, out["par", ], out["st_err", ], out["t-val", ], call
    )
    class(results) <- "sfareg"
    names(results) <- c(
      "out", "opt", "total_time", "start_v", "r_hat_m", "exp_u_hat", "u_hat", "model_name",
      "formula", "data", "coefficients", "std.errors", "t.values", "call"
    )
    return(results)
  }
  if (model_name == "FD") {
    Start.Time <- start.time()

    if (isTRUE(is.numeric(start_val))) {
      start_v <- start_val
    } else {
      beta_hat <- plm_fd$coefficients[c(x_vars_vec)]
      epsilon_hat <- plm_fd$residuals
      beta_se <- as.data.frame(summary(plm_fd)[1])$coefficients.Std..Error
      sfa_eps <- pcs_c(Y = as.numeric(epsilon_hat))[[1]]$par

      exp_u <- sfa_eps[3]
      sigma_v <- sqrt(sfa_eps[2]^2 / (1 + sfa_eps[1]^2)) # coef(sfa_eps,extraPar=TRUE)[c("sigmaV")]
      sigma_u <- sigma_v * sfa_eps[1] # coef(sfa_eps,extraPar=TRUE)[c("sigmaU")]

      sigmaSq_u <- sigma_u^2
      sigmaSq_v <- sigma_v^2
      mu <- 0.1
      delta <- rep(0.1, length(z_vars))
      start_v <- unname(c(sigmaSq_u, sigmaSq_v, mu, beta_hat, delta))
    }


    out <- matrix(0, nrow = 3, ncol = length(start_v))
    rownames(out) <- c("par", "st_err", "t-val")
    colnames(out) <- c("sig_u2", "sig_v2", "mu", x_vars_vec, z_vars)


    indiv <- noquote(as.vector(unique(data[, c(individual)])))
    t <- rep(0, N)
    data_i <- as.list(rep(0, N))
    Y <- as.list(rep(0, N))
    Y_diff <- as.list(rep(0, N))
    eps <- as.list(rep(0, N))
    data_i_vars <- as.list(rep(0, N))
    data_i_vars_diff <- as.list(rep(0, N))
    SIGMA <- as.list(rep(0, N))
    data_z_vars <- as.list(rep(0, N))

    for (ii in seq_len(N)) {
      data_i[[ii]] <- data[which(data[, c(individual)] == indiv[ii]), ]
      t[ii] <- nrow(data_i[[ii]])
      Y[[ii]] <- as.numeric(data_i[[ii]][, y_var])
      ## Precomputed once (used to be recomputed inside like.fd()'s fn1() on
      ## every likelihood evaluation via diff(Y[[i]], lag=1) /
      ## diff(as.numeric(data_i_vars[[i]][,qq]), lag=1) -- neither depends on
      ## the parameter vector x being optimized, so this was pure repeated
      ## work across every bobyqa/psoptim/optim call.
      Y_diff[[ii]] <- diff(Y[[ii]], lag = 1)
      data_i_vars[[ii]] <- data.frame(data_i[[ii]][, c(x_vars_vec)])
      data_i_vars_diff[[ii]] <- diff(as.matrix(data_i_vars[[ii]]), lag = 1)
      data_z_vars[[ii]] <- data.frame(data_i[[ii]][, c(z_vars)])

      if (t[ii] == 2) {
        SIGMA[[ii]] <- 2
      }
      if (t[ii] == 3) {
        SIGMA[[ii]] <- matrix(c(2, -1, -1, 2), 2, 2)
      }
      if (t[ii] > 3) {
        SIGMA[[ii]] <- matrix(0, nrow = t[ii] - 1, ncol = t[ii] - 1)
        diag(SIGMA[[ii]]) <- 2
        diag(SIGMA[[ii]][-1, ]) <- -1
        diag(SIGMA[[ii]][, -1]) <- -1
      }
    }

    like.fd <- function(x) {
      for (q in seq_len(n_x_vars)) {
        v <- q + 3
        x_x_vec[q] <- x[v]
      }

      for (qq in seq_len(n_z_vars)) {
        m <- v + qq
        z_z_vec[qq] <- x[m]
      }

      fn1 <- function(i) {
        ## Y_diff[[i]]/data_i_vars_diff[[i]] precomputed once above (were
        ## recomputed via diff() here on every likelihood evaluation despite
        ## not depending on x).
        eps_t <- Y_diff[[i]]
        eps_h <- rep(0, t[i])

        for (qq in seq_len(n_x_vars)) {
          eps_t <- eps_t - x_x_vec[qq] * data_i_vars_diff[[i]][, qq]
        }
        for (qq in seq_len(n_z_vars)) {
          eps_h <- eps_h + z_z_vec[qq] * as.numeric(data_z_vars[[i]][, qq])
        }

        eps_h <- diff(exp(eps_h), lag = 1)
        eps_t <- eps_t * inefdec_n ## not exactly sure if this is sufficient for cost
        SIG <- x[2] * SIGMA[[i]]
        sig_star2 <- (t(eps_h) %*% qr.solve(SIG) %*% eps_h + (1 / x[1]))^-1
        mu_star <- ((x[3] / x[1]) - t(eps_t) %*% qr.solve(SIG) %*% eps_h) * sig_star2

        l1 <- -0.5 * (t[i] - 1) * log(2 * pi)
        l2 <- -0.5 * log(t[i])
        l3 <- -0.5 * (t[i] - 1) * log(x[2])
        l4 <- -0.5 * t(eps_t) %*% qr.solve(SIG) %*% eps_t
        l5 <- 0.5 * ((mu_star^2 / sig_star2) - ((x[3]^2) / x[1]))
        l6 <- log(sqrt(sig_star2) * max(pnorm(mu_star / sqrt(sig_star2)), .SFA_CONSTANTS$MIN_POSITIVE))
        l7 <- -log(sqrt(x[1]) * max(pnorm(x[3] / sqrt(x[1])), .SFA_CONSTANTS$MIN_POSITIVE))


        prod_vec_n <- sum(l1, l2, l3, l4, l5, l6, l7)
        return(-prod_vec_n)
      }

      fn1_apply <- unlist(lapply(1:N, fn1))

      return(sum(fn1_apply[is.finite(fn1_apply)]))
    }

    differ <- 2

    Opt.Bobyqa <- opt.bobyqa(
      fn = like.fd, start_v = start_v, lower.bobyqa = c(rep(0.000001, 2), rep(-Inf, n_x_vars + 1), rep(-Inf, length(z_vars))),
      maxit.bobyqa = maxit.bobyqa, bob.TF = TRUE, verbose = verbose
    )
    start_v <- Opt.Bobyqa$start_v
    start_feval <- Opt.Bobyqa$start_feval
    bob1 <- Opt.Bobyqa$bob1

    Opt.Psoptim <- opt.psoptim(
      fn = like.fd, start_v = start_v, lower.psoptim = c(rep(.SFA_CONSTANTS$MIN_POSITIVE, 2), start_v[-c(1:2)] - differ),
      rand.psoptim = rand.psoptim,
      upper.psoptim = c(start_v + differ), maxit.psoptim, psopt.TF = PSopt, verbose = verbose
    )
    start_v <- Opt.Psoptim$start_v
    start_feval <- Opt.Psoptim$start_feval
    opt00 <- Opt.Psoptim$opt00

    Opt.Optim <- opt.optim(
      fn = like.fd, start_v = start_v, lower.optim = c(rep(0.000001, 2), rep(-Inf, n_x_vars + 1), rep(-Inf, length(z_vars))),
      upper.optim = c(start_v + differ), maxit.optim = maxit.optim, opt.TF = optHessian, method = Method, optHessian = optHessian, verbose = verbose
    )
    start_v <- Opt.Optim$start_v
    start_feval <- Opt.Optim$start_feval
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

    ## TE Measurements: ui's
    NX <- n_x_vars + 3
    NZ1 <- n_x_vars + 4
    NZ2 <- n_x_vars + n_z_vars + 3
    beta <- opt$par[c(4:NX)]
    delta <- opt$par[c(NZ1:NZ2)]
    sigu2 <- opt$par[1]
    sigv2 <- opt$par[2]
    mu <- opt$par[3]
    u_hat <- rep(0, sum(t))
    h_hat <- exp(as.matrix(data.frame(subset(data, select = z_vars))) %*% delta)

    ## eps_t/eps_h/SIG/sig_star2/mu_star below do not depend on `tt` at all
    ## (only the final u_hat[num] line does, via h_hat[num]) -- this used to
    ## recompute all of them t[i] times per individual i (t[i] identical
    ## repeats) instead of once. Hoisted the tt-invariant part out of the
    ## inner loop; math is unchanged, including reusing the same
    ## diff(as.numeric(data_z_vars[[i]][,qq]), lag=1) for eps_h that this
    ## block already used (note: like.fd()'s fn1() above uses the
    ## UNdifferenced data_z_vars for eps_h -- a pre-existing discrepancy
    ## between the likelihood and this TE-prediction step left untouched
    ## here since this is a performance pass, not a correctness one).
    for (i in seq_len(N)) {
      eps_t <- Y_diff[[i]]
      eps_h <- rep(0, t[i] - 1)

      for (qq in seq_len(n_x_vars)) {
        m <- 3 + qq
        eps_t <- eps_t - opt$par[m] * data_i_vars_diff[[i]][, qq]
      }

      for (qq in seq_len(n_z_vars)) {
        m <- 3 + length(n_x_vars) + qq
        eps_h <- eps_h + opt$par[m] * diff(as.numeric(data_z_vars[[i]][, qq]), lag = 1)
      }

      eps_t <- eps_t * inefdec_n
      SIG <- sigv2 * SIGMA[[i]]
      sig_star2 <- (t(eps_h) %*% qr.solve(SIG) %*% eps_h + (1 / sigu2))^-1
      mu_star <- ((mu / sigu2) - t(eps_t) %*% qr.solve(SIG) %*% eps_h) * sig_star2

      for (tt in seq_len(t[i])) {
        num <- if (i > 1) {
          cumsum(t)[i - 1] + tt
        } else {
          tt
        }
        u_hat[num] <- h_hat[num] * (mu_star + (sqrt(sig_star2) * dnorm(mu_star / sqrt(sig_star2)) / max(pnorm(mu_star / sqrt(sig_star2)), .SFA_CONSTANTS$MIN_POSITIVE)))
      }
    }

    exp_u_hat <- exp(-u_hat)
    exp_u_hat <- pmax(exp_u_hat, 0)
    exp_u_hat <- pmin(exp_u_hat, 1)

    st_err <- if (isTRUE(any(opt$hessian == 0)) | optHessian == FALSE) {
      rep(NA, length(opt$par))
    } else {
      suppressWarnings(sqrt(diag(qr.solve(opt$hessian))))
    }
    t_val <- opt$par / st_err
    out[1, ] <- opt$par
    out[2, ] <- st_err
    out[3, ] <- t_val
    results <- list(t(out), c(opt), End.Time, start_v, model_name, formula, u_hat, h_hat, exp_u_hat, data, out["par", ], out["st_err", ], out["t-val", ], call)
    class(results) <- "sfareg"
    names(results) <- c("out", "opt", "total_time", "start_v", "model_name", "formula", "u_hat", "h_hat", "exp_u_hat", "data", "coefficients", "std.errors", "t.values", "call")
    return(results)
  }
  if (model_name == "GTRE_SEQ1") {
    Start.Time <- start.time()

    ## Sequential Method -- each stage is an ordinary cross-sectional
    ## intercept-only normal-half-normal fit of one composite residual
    ## (epsilon_hat, alpha_hat) against a constant, via .fit_nhn_intercept()
    ## (matrix_utils.R; see that helper's own header for why sfm() itself
    ## can't be reused directly here). reparam_gs() below delta-methods the
    ## helper's (sigma_v,sigma_u) fit into (gamma,sigmaSq) -- matching this
    ## block's own downstream naming -- via numDeriv::jacobian(), the same
    ## technique used for PL80/BC92's reparameterization.
    fit_eps <- .fit_nhn_intercept(as.numeric(epsilon_hat), inefdec_n = inefdec_n)
    fit_alp <- .fit_nhn_intercept(as.numeric(alpha_hat), inefdec_n = inefdec_n)

    reparam_gs <- function(p) {
      sv <- p[1]
      su <- p[2]
      c(gamma = (su * su) / (sv * sv + su * su), sigmaSq = sv * sv + su * su)
    }
    gs_se <- function(fit) {
      p <- fit$par[1:2]
      if (isTRUE(any(fit$hessian == 0))) {
        return(c(NA_real_, NA_real_))
      }
      jac <- tryCatch(numDeriv::jacobian(reparam_gs, p), error = function(e) NULL)
      vc <- tryCatch(solve(fit$hessian)[1:2, 1:2], error = function(e) NULL)
      if (is.null(jac) || is.null(vc)) {
        c(NA_real_, NA_real_)
      } else {
        suppressWarnings(sqrt(pmax(diag(jac %*% vc %*% t(jac)), 0)))
      }
    }

    sigma_v <- fit_eps$par[1]
    sigma_u <- fit_eps$par[2]
    ## .fit_nhn_intercept() returns (TWO-sided sd, ONE-sided sd, intercept).
    ## alpha_i = r_i - h_i has r two-sided and h one-sided, so par[1] is
    ## sigma_r and par[2] is sigma_h -- verified on synthetic r - h with known
    ## components: par = (0.2075, 0.4062) against truths of 0.20 and 0.40.
    ## These were assigned the other way round.
    ##
    ## The swap was easy to miss on real data because ranef() returns SHRUNKEN
    ## individual effects, which compresses the skewness that identifies the
    ## one-sided part; the fit then attributes most of the variance to par[1],
    ## and the swapped assignment happens to reproduce the expected ordering
    ## sigma_h > sigma_r. Only the reported other_parms were affected --
    ## gamma_hr and sigmaSq_hr come from reparam_gs(par[1:2]), which reads the
    ## positions correctly.
    ##
    ## That shrinkage is also the leading explanation for this model's
    ## convergence failure on sigmaSq_hr; see CONVERGENCE_STATUS.md.
    sigma_r <- fit_alp$par[1]
    sigma_h <- fit_alp$par[2]
    exp_u <- fit_eps$par[3]
    exp_eta <- fit_alp$par[3]
    beta_0 <- beta_0_st + exp_u + exp_eta
    Lambda <- sigma_u / sigma_v
    Sigma <- sqrt(sigma_u^2 + sigma_v^2)
    sigma_se_r <- NA
    sigma_r_se <- NA
    sigma_h_se <- NA
    lambda_se_r <- NA
    beta_0_se <- NA
    gs_eps <- reparam_gs(fit_eps$par[1:2])
    gs_alp <- reparam_gs(fit_alp$par[1:2])
    gamma_uv <- gs_eps[["gamma"]]
    sigmaSq_uv <- gs_eps[["sigmaSq"]]
    gamma_hr <- gs_alp[["gamma"]]
    sigmaSq_hr <- gs_alp[["sigmaSq"]]
    se_eps <- gs_se(fit_eps)
    se_alp <- gs_se(fit_alp)
    gamma_uv_se <- se_eps[1]
    sigmaSq_uv_se <- se_eps[2]
    gamma_hr_se <- se_alp[1]
    sigmaSq_hr_se <- se_alp[2]

    other_parms <- as.matrix(c(Lambda, Sigma, beta_0))
    rownames(other_parms) <- c("lambda", "sigma", "beta_0")
    start_v <- unname(c(gamma_uv, sigmaSq_uv, gamma_hr, sigmaSq_hr, beta_hat))
    out <- matrix(0, nrow = 3, ncol = length(start_v))
    rownames(out) <- c("par", "st_err", "t-val")
    colnames(out) <- if (isTRUE(intercept == 0)) {
      c("gamma_uv", "sigmaSq_uv", "gamma_hr", "sigmaSq_hr", colnames(data_x))
    } else {
      c("gamma_uv", "sigmaSq_uv", "gamma_hr", "sigmaSq_hr", colnames(data_x)[-c(1)])
    }

    End.Time <- end.time(Start.Time)
    st_err <- rep(NA, ncol(out))
    st_err <- if (isTRUE(intercept == 0)) {
      c(gamma_uv_se, sigmaSq_uv_se, gamma_hr_se, sigmaSq_hr_se, summary(plm_gtre)$coefficients[, 2])
    } else {
      c(gamma_uv_se, sigmaSq_uv_se, gamma_hr_se, sigmaSq_hr_se, summary(plm_gtre)$coefficients[-c(1), 2])
    }
    t_val <- start_v / st_err
    out[1, ] <- start_v
    out[2, ] <- st_err
    out[3, ] <- t_val
    results <- list(t(out), End.Time, other_parms, model_name, formula, data, out["par", ], out["st_err", ], out["t-val", ], call)
    class(results) <- "sfareg"
    names(results) <- c("out", "total_time", "other_parms", "model_name", "formula", "data", "coefficients", "std.errors", "t.values", "call")
    return(results)
  }
  if (model_name == "SSFE") {
    ## Schmidt & Sickles (1984, JBES) fixed-effects (LSDV) stochastic
    ## frontier estimator -- a genuinely different, deterministic model
    ## from both "TFE" (Greene's 2005 true fixed effects MLE) and
    ## "TFE_WMLE" (Chen, Schmidt and Wang's 2014 within MLE) above.
    ## Schmidt-Sickles needs no numerical
    ## optimization at all: fit the standard fixed-effects (within) panel
    ## regression, then read technical inefficiency directly off the
    ## estimated firm effects -- the best-performing firm (highest fixed
    ## effect for a production frontier, lowest for a cost frontier) is
    ## treated as fully efficient, and every other firm's inefficiency is
    ## its distance from that benchmark. The estimator itself is entirely
    ## deterministic -- no optimizer to fail to converge.
    Start.Time <- start.time()

    ## index = individual is required: by this point `data` has been through
    ## data_proc2() and is a plain data.frame, so plm() would otherwise treat
    ## the first two columns as the panel index. See start.tfe() in
    ## starting.values.R for the full account of this failure mode.
    plm_ss <- plm(formula_x, data,
      effect = "individual", model = "within",
      index = individual
    )

    beta_hat_ss <- plm_ss$coefficients[x_vars_vec]
    beta_se_ss <- summary(plm_ss)$coefficients[x_vars_vec, "Std. Error"]

    alpha_hat <- as.numeric(fixef(plm_ss))
    names(alpha_hat) <- names(fixef(plm_ss))

    ## Firm-specific inefficiency from the fixed effects (Schmidt &
    ## Sickles, 1984): distance from the best-performing firm.
    if (isTRUE(inefdec)) {
      u_hat <- max(alpha_hat) - alpha_hat
    } else {
      u_hat <- alpha_hat - min(alpha_hat)
    }
    exp_u_hat <- exp(-u_hat)

    End.Time <- end.time(Start.Time)

    out <- matrix(0, nrow = 3, ncol = length(beta_hat_ss))
    rownames(out) <- c("par", "st_err", "t-val")
    colnames(out) <- x_vars_vec
    out[1, ] <- beta_hat_ss
    out[2, ] <- beta_se_ss
    out[3, ] <- out[1, ] / out[2, ]

    ## No log-likelihood/opt object for this estimator (it isn't MLE), so
    ## logLik.sfareg()/AIC()/BIC() will return NA with a warning for SSFE
    ## fits, same as psfm()'s GTRE_SEQ1/GTRE_SEQ2 -- see R/sfareg_methods.R.
    results <- list(
      t(out), End.Time, model_name, formula, data, alpha_hat, u_hat, exp_u_hat,
      out["par", ], out["st_err", ], out["t-val", ], call
    )
    class(results) <- "sfareg"
    names(results) <- c(
      "out", "total_time", "model_name", "formula", "data", "alpha_hat", "u_hat", "exp_u_hat",
      "coefficients", "std.errors", "t.values", "call"
    )
    return(results)
  }
  if (model_name %in% c("PL80", "BC92", "K1990", "K1990modified")) {
    ## Time-invariant panel stochastic frontier (Pitt and Lee, 1980, JoE)
    ## and three time-decay generalizations of it -- the "error components
    ## frontier": y_it = x_it'beta + v_it - B_it*u_i (production; sign
    ## flips for cost via inefdec_n), v_it ~ iid N(0, sigma_v^2), u_i ~ iid
    ## N+(0, sigma_u^2) constant across t within firm i, and B_it a
    ## model-specific decay pattern:
    ##   PL80          : B_it = 1                             (no decay)
    ##   BC92          : B_it = exp(-eta*(t - Tref))           (Battese &
    ##                   Coelli 1992; Tref = the LAST time period in the
    ##                   *whole panel*, not per-firm -- confirmed by
    ##                   matching frontier::sfa()'s log-likelihood on an
    ##                   unbalanced-panel test case; the per-firm-last-
    ##                   period alternative gives a different, non-matching
    ##                   logLik there)
    ##   K1990         : B_it = (1 + exp(b*t + c*t^2))^-1      (Kumbhakar
    ##                   1990; matches npsf::sf()'s "K1990" time pattern)
    ##   K1990modified : B_it = 1 + d*(t-T_i) + e*(t-T_i)^2    (Kumbhakar
    ##                   1990's modified pattern; T_i is each FIRM's OWN
    ##                   last observed period, unlike BC92's whole-panel
    ##                   Tref)
    ## See .build_Bit() (matrix_utils.R) for the actual B_it construction --
    ## centralized there since it's the only thing that differs across
    ## these four models. Native implementation (PL80/BC92 previously
    ## wrapped frontier::sfa() -- removed as a dependency since frontier is
    ## a direct competitor to this package).
    ##
    ## Closed-form log-likelihood, derived by integrating u_i out of the
    ## joint density of (epsilon_i, u_i) (completing the square in u_i,
    ## same technique as the standard cross-sectional NHN derivation, here
    ## applied to a firm's whole T_i-vector sharing one u_i):
    ##   log f(eps_i) = log(2) - (T_i/2)*log(2*pi) - (T_i-1)*log(sigma_v)
    ##     - (1/2)*log(denom_i) - SSE_i/(2*sigma_v^2)
    ##     + sigma_u^2*S_i^2/(2*sigma_v^2*denom_i) + log(Phi(z_i))
    ## where SSE_i = sum(eps_it^2), S_i = sum(B_it*eps_it),
    ## S2_i = sum(B_it^2), denom_i = S2_i*sigma_u^2 + sigma_v^2,
    ## z_i = -sigma_u*S_i/(sigma_v*sqrt(denom_i)). Setting B_it=1 recovers
    ## Pitt & Lee (1980) exactly; setting T_i=1 (cross-sectional) recovers
    ## the ordinary Aigner-Lovell-Schmidt (1977) normal-half-normal
    ## likelihood exactly -- both checked analytically. This derivation
    ## never assumed a specific functional form for B_it, so it applies to
    ## K1990/K1990modified unchanged -- only .build_Bit()'s construction of
    ## B_it itself differs by model.
    ##
    ## Verified against frontier::sfa() (PL80/BC92 only -- K1990/
    ## K1990modified have no frontier::sfa() equivalent to check against,
    ## see the parameter-recovery validation in PROJECT_STATUS.md's "Round
    ## 20" entry instead): computed this closed form at frontier::sfa()'s
    ## own reported MLE (both PL80 and BC92, balanced and unbalanced
    ## panels) and matched its reported logLik to ~1e-5 (its own reporting
    ## precision) in every case, before that dependency was removed.
    Start.Time <- start.time()

    time_var <- if (is.null(time)) {
      ave(seq_len(nrow(data)), data[[individual]], FUN = seq_along)
    } else {
      ## data[[time]] can come back as a factor/pseries here (the shared
      ## data_proc()/data_proc2() setup upstream coerces the working data
      ## via plm for panel handling) -- as.numeric() on a factor would
      ## return level *indices*, not the level values themselves, so go
      ## through as.character() first.
      tv <- data[[time]]
      if (is.factor(tv)) as.numeric(as.character(tv)) else as.numeric(tv)
    }
    Tref_pl <- max(time_var)

    X_pl <- model.matrix(formula_x, data = data)
    y_pl <- as.numeric(data[[y_var]])
    K_pl <- ncol(X_pl)
    idx_by_firm <- split(seq_len(nrow(data)), data[[individual]], drop = TRUE)
    time_by_firm <- lapply(idx_by_firm, function(idx) time_var[idx])

    ## Number and names of this model's extra (decay) parameters beyond
    ## (sigma_v, sigma_u, beta).
    n_decay <- c(PL80 = 0L, BC92 = 1L, K1990 = 2L, K1990modified = 2L)[[model_name]]
    decay_names <- switch(model_name,
      PL80 = character(0),
      BC92 = "time",
      K1990 = c("b", "c"),
      K1990modified = c("d", "e")
    )

    like.pl <- function(x) {
      sigma_v <- x[1]
      sigma_u <- x[2]
      beta <- x[3:(2 + K_pl)]
      decay_par <- if (n_decay > 0) x[(3 + K_pl):(2 + K_pl + n_decay)] else numeric(0)

      if (!all(is.finite(c(sigma_v, sigma_u, decay_par))) || sigma_v <= 0 || sigma_u <= 0) {
        return(1e12)
      }

      resid_all <- inefdec_n * (y_pl - as.vector(X_pl %*% beta))

      ll_i <- function(idx, t_i) {
        e_i <- resid_all[idx]
        Ti <- length(e_i)
        B_it <- .build_Bit(model_name, t_i, Tref_pl, decay_par)
        SSE <- sum(e_i * e_i)
        S1 <- sum(B_it * e_i)
        S2 <- sum(B_it * B_it)
        denom <- pmax(S2 * sigma_u * sigma_u + sigma_v * sigma_v, .SFA_CONSTANTS$MIN_POSITIVE)
        z <- -sigma_u * S1 / (sigma_v * sqrt(denom))
        log(2) - (Ti / 2) * log(2 * pi) - (Ti - 1) * log(sigma_v) - 0.5 * log(denom) -
          SSE / (2 * sigma_v * sigma_v) + sigma_u * sigma_u * S1 * S1 / (2 * sigma_v * sigma_v * denom) +
          pnorm(z, log.p = TRUE)
      }
      ll_vec <- mapply(ll_i, idx_by_firm, time_by_firm)
      -sum(ll_vec[is.finite(ll_vec)])
    }

    ## Starting values: OLS beta, residual variance split evenly between
    ## sigma_u^2/sigma_v^2 (a robust, standard default -- moment/skewness-
    ## based starting values are known to be fragile for SFA models with
    ## weakly-skewed residuals, e.g. the "right-skewed" warning frontier's
    ## own OLS step often emits on simulated data). Decay parameters start
    ## at 0, i.e. B_it starts at PL80's no-decay/K1990modified's no-decay
    ## value (1) for eta/d/e; K1990's b=c=0 starts at the model's own
    ## midpoint B_it=0.5 rather than 1 -- there is no zero-parameter value
    ## of K1990's logistic form that reduces to "no decay", so 0 is simply
    ## a neutral starting point for the search.
    ols_pl <- lm.fit(X_pl, y_pl)
    resid_var <- max(mean((inefdec_n * ols_pl$residuals)^2), .SFA_CONSTANTS$MIN_POSITIVE)
    ## sigma_u starts ABOVE sigma_v rather than at an equal split of the
    ## residual variance. For a composed error the one-sided component
    ## carries most of the residual spread, so sigma_v = sigma_u is a poor
    ## prior -- and for the two-parameter decay patterns it is not merely
    ## inefficient, it changes the answer: with the equal split K1990
    ## converged to a log-likelihood of -429.30 on a test panel whose true
    ## parameters score -423.03, while starting sigma_u at the residual RMSE
    ## and sigma_v at half of it reaches -423.54 (bobyqa and nlminb agree on
    ## that point independently, so it is the start, not the optimizer).
    start_v <- c(
      sqrt(resid_var) / 2, sqrt(resid_var), unname(coef(ols_pl)),
      rep(0, n_decay)
    )
    n_par <- length(start_v)

    lower_pl <- c(.SFA_CONSTANTS$MIN_POSITIVE, .SFA_CONSTANTS$MIN_POSITIVE, rep(-Inf, n_par - 2))

    ## nlminb first. psfm()'s default maxit.bobyqa is 100 -- fine for the
    ## 4-5 parameter models it was tuned on, but far too few for the
    ## 7-parameter K1990/K1990modified fits, which were stopping ~6 log-
    ## likelihood units short of the optimum purely on the iteration cap.
    ## nlminb reaches the same point bobyqa reaches with 3000 evaluations,
    ## so this buys accuracy rather than trading it away. bobyqa still runs
    ## afterwards from wherever nlminb finished.
    Opt.Nlminb <- opt.nlminb(
      fn = like.pl, start_v = start_v, lower.nlminb = lower_pl,
      maxit.nlminb = maxit.nlminb, nlminb.TF = TRUE, verbose = verbose
    )
    start_v <- Opt.Nlminb$start_v

    Opt.Bobyqa <- opt.bobyqa(
      fn = like.pl, start_v = start_v, lower.bobyqa = lower_pl,
      maxit.bobyqa = maxit.bobyqa, bob.TF = TRUE, verbose = verbose
    )
    start_v <- Opt.Bobyqa$start_v
    start_feval <- Opt.Bobyqa$start_feval
    bob1 <- Opt.Bobyqa$bob1

    differ <- 2
    Opt.Psoptim <- opt.psoptim(
      fn = like.pl, start_v, lower.psoptim = c(rep(.SFA_CONSTANTS$MIN_POSITIVE, 2), start_v[-c(1:2)] - differ),
      rand.psoptim = rand.psoptim,
      upper.psoptim = c(start_v + differ), maxit.psoptim = maxit.psoptim, psopt.TF = PSopt, verbose = verbose
    )
    start_v <- Opt.Psoptim$start_v
    start_feval <- Opt.Psoptim$start_feval
    opt00 <- Opt.Psoptim$opt00

    Opt.Optim <- opt.optim(
      fn = like.pl, start_v = start_v, lower.optim = lower_pl, upper.optim = rep(Inf, n_par),
      maxit.optim = maxit.optim, opt.TF = TRUE, method = Method, optHessian = optHessian, verbose = verbose
    )
    start_v <- Opt.Optim$start_v
    start_feval <- Opt.Optim$start_feval
    opt <- Opt.Optim$opt

    if (optHessian == FALSE & PSopt == FALSE) {
      opt <- bob1
    }
    if (optHessian == FALSE & PSopt == TRUE) {
      opt <- opt00
    }

    par_hat <- opt$par
    sigma_v_hat <- par_hat[1]
    sigma_u_hat <- par_hat[2]
    beta_hat_pl <- par_hat[3:(2 + K_pl)]
    decay_hat <- if (n_decay > 0) par_hat[(3 + K_pl):(2 + K_pl + n_decay)] else numeric(0)
    sigmaSq_hat <- sigma_v_hat^2 + sigma_u_hat^2
    gamma_hat <- sigma_u_hat^2 / sigmaSq_hat

    ## Delta-method Jacobian from the optimizer's own (sigma_v,sigma_u,beta,
    ## [decay]) parameterization to the reported (beta,sigmaSq,gamma,[decay])
    ## parameterization -- matches frontier::sfa()'s output convention
    ## (PL80/BC92) for backward compatibility with existing
    ## psfm(model_name="PL80"/"BC92") callers; K1990/K1990modified's decay
    ## parameters are reported raw (untransformed), same treatment BC92's
    ## "time"/eta already got.
    par_names <- c(colnames(X_pl), "sigmaSq", "gamma", decay_names)
    reparam <- function(p) {
      sv <- p[1]
      su <- p[2]
      b <- p[3:(2 + K_pl)]
      dcy <- if (n_decay > 0) p[(3 + K_pl):(2 + K_pl + n_decay)] else NULL
      c(b, sv^2 + su^2, su^2 / (sv^2 + su^2), dcy)
    }
    out <- matrix(0, nrow = 3, ncol = length(par_names))
    rownames(out) <- c("par", "st_err", "t-val")
    colnames(out) <- par_names
    out[1, ] <- reparam(par_hat)

    st_err <- if (isTRUE(any(opt$hessian == 0)) | optHessian == FALSE) {
      rep(NA, length(par_names))
    } else {
      jac <- tryCatch(numDeriv::jacobian(reparam, par_hat), error = function(e) NULL)
      vc <- tryCatch(solve(opt$hessian), error = function(e) NULL)
      if (is.null(jac) || is.null(vc)) {
        rep(NA, length(par_names))
      } else {
        suppressWarnings(sqrt(pmax(diag(jac %*% vc %*% t(jac)), 0)))
      }
    }
    out[2, ] <- st_err
    out[3, ] <- out[1, ] / out[2, ]

    ## Per-observation predicted efficiency, E[exp(-B_it*u_i)|eps_i] --
    ## generalizes Battese & Coelli's (1988) JLMS point estimator (their
    ## formula is this with B_it=1) to a truncated-normal posterior for
    ## u_i implied by the same completed-square derivation above:
    ## u_i|eps_i ~ N+(mu_star_i, sigma_star_i^2).
    resid_hat_all <- inefdec_n * (y_pl - as.vector(X_pl %*% beta_hat_pl))
    exp_u_hat <- rep(NA_real_, nrow(data))
    for (nm in names(idx_by_firm)) {
      idx <- idx_by_firm[[nm]]
      t_i <- time_by_firm[[nm]]
      e_i <- resid_hat_all[idx]
      B_it <- .build_Bit(model_name, t_i, Tref_pl, decay_hat)
      S1 <- sum(B_it * e_i)
      S2 <- sum(B_it * B_it)
      denom <- pmax(S2 * sigma_u_hat^2 + sigma_v_hat^2, .SFA_CONSTANTS$MIN_POSITIVE)
      sigma_star2 <- (sigma_v_hat^2 * sigma_u_hat^2) / denom
      sigma_star <- sqrt(sigma_star2)
      mu_star <- -sigma_star2 * S1 / sigma_v_hat^2
      ratio <- mu_star / sigma_star
      exp_u_hat[idx] <- exp(-B_it * mu_star + (B_it^2) * sigma_star2 / 2) *
        pmax(pnorm(ratio - B_it * sigma_star), .SFA_CONSTANTS$MIN_POSITIVE) /
        pmax(pnorm(ratio), .SFA_CONSTANTS$MIN_POSITIVE)
    }

    End.Time <- end.time(Start.Time)

    results <- list(
      t(out), opt, End.Time, model_name, formula, data, exp_u_hat,
      out["par", ], out["st_err", ], out["t-val", ], call
    )
    class(results) <- "sfareg"
    names(results) <- c(
      "out", "opt", "total_time", "model_name", "formula", "data", "exp_u_hat",
      "coefficients", "std.errors", "t.values", "call"
    )
    return(results)
  }
  if (model_name == "GTRE_SEQ2") {
    Start.Time <- start.time()
    ## Sequential Method following 1995 paper
    ## take second and third moments of alpha_hat and epsilon_hat
    alp_2m <- mean(alpha_hat^2)
    alp_3m <- min(0, mean(alpha_hat^3))
    eps_2m <- mean(epsilon_hat^2)
    eps_3m <- min(0, mean(epsilon_hat^3))

    ## The moment inversion itself now lives in .gtre_two_step() (matrix_utils.R)
    ## so that GTRE_FML can seed its FIML search from the same estimator, which
    ## is what Colombi (2010, sec. 3) and Colombi, Martini and Vittadini (2011,
    ## sec. 3.2) recommend. The delta-method standard errors below still need
    ## the raw moments, so those stay here.
    .ts <- .gtre_two_step(epsilon_hat, alpha_hat, beta_0_st)
    gamma_uv <- .ts$gamma_uv
    gamma_hr <- .ts$gamma_hr
    sigmaSq_uv <- .ts$sigmaSq_uv
    sigmaSq_hr <- .ts$sigmaSq_hr
    beta_0 <- .ts$beta_0

    ## calculate the ten needed central moments
    mu_2_eps <- sigmaSq_uv * ((1 - gamma_uv) + gamma_uv * ((pi - 2) / pi))
    mu_3_eps <- sigmaSq_uv^(3 / 2) * (sqrt(2 / pi) * (1 - (4 / pi)) * gamma_uv^(3 / 2))
    mu_4_eps <- sigmaSq_uv^2 * (3 * (1 - gamma_uv)^2 + ((6 * (pi - 2) * gamma_uv * (1 - gamma_uv)) / pi) + gamma_uv^2 * (3 - (4 / pi) - (12 / pi^2)))
    mu_5_eps <- sigmaSq_uv^(5 / 2) * gamma_uv^(3 / 2) * sqrt(2 / pi) * (10 * (1 - (4 / pi)) * (1 - gamma_uv) + (7 - (20 / pi) - (16 / pi^2)) * gamma_uv)
    mu_6_eps <- sigmaSq_uv^3 * (15 * (1 - gamma_uv)^3 +
      (45 * (pi - 2) * (1 - gamma_uv)^2 * gamma_uv / pi) +
      15 * (3 - (4 / pi) - (12 / pi^2)) * (1 - gamma_uv) * gamma_uv^2 +
      (15 - (6 / pi) - (100 / pi^2) - (40 / pi^3)) * gamma_uv^3)

    mu_2_alp <- sigmaSq_hr * ((1 - gamma_hr) + gamma_hr * ((pi - 2) / pi))
    mu_3_alp <- sigmaSq_hr^(3 / 2) * (sqrt(2 / pi) * (1 - (4 / pi)) * gamma_hr^(3 / 2))
    mu_4_alp <- sigmaSq_hr^2 * (3 * (1 - gamma_hr)^2 + ((6 * (pi - 2) * gamma_hr * (1 - gamma_hr)) / pi) + gamma_hr^2 * (3 - (4 / pi) - (12 / pi^2)))
    mu_5_alp <- sigmaSq_hr^(5 / 2) * gamma_hr^(3 / 2) * sqrt(2 / pi) * (10 * (1 - (4 / pi)) * (1 - gamma_hr) + (7 - (20 / pi) - (16 / pi^2)) * gamma_hr)
    mu_6_alp <- sigmaSq_hr^3 * (15 * (1 - gamma_hr)^3 +
      (45 * (pi - 2) * (1 - gamma_hr)^2 * gamma_hr / pi) +
      15 * (3 - (4 / pi) - (12 / pi^2)) * (1 - gamma_hr) * gamma_hr^2 +
      (15 - (6 / pi) - (100 / pi^2) - (40 / pi^3)) * gamma_hr^3)

    var_2m_eps <- (1 / nrow(data)) * (mu_4_eps - mu_2_eps^2)
    var_3m_eps <- (1 / nrow(data)) * (mu_6_eps - mu_3_eps^3 - 6 * mu_2_eps * mu_4_eps + 9 * mu_2_eps)
    cov_23m_eps <- (1 / nrow(data)) * (mu_5_eps - 4 * mu_2_eps * mu_3_eps)

    var_2m_alp <- (1 / N) * (mu_4_alp - mu_2_alp^2)
    var_3m_alp <- (1 / N) * (mu_6_alp - mu_3_alp^3 - 6 * mu_2_alp * mu_4_alp + 9 * mu_2_alp)
    cov_23m_alp <- (1 / N) * (mu_5_alp - 4 * mu_2_alp * mu_3_alp)

    ## define 8 needed derivatives
    d_beta_0_d_m3_eps <- (pi / (pi - 4)) * (1 / 3) * (sqrt(pi / 2) * (pi / (pi - 4)) * eps_3m)^(-2 / 3)
    d_sigma2_d_m3_eps <- sqrt(2 / pi) * (pi / (pi - 4)) * (2 / 3) * (sqrt(pi / 2) * (pi / (pi - 4)) * eps_3m)^(-1 / 3)
    d_gamma_d_m2_eps <- -(sqrt(pi / 2) * (pi / (pi - 4)) * eps_3m)^(-2 / 3) * (eps_3m * (sqrt(pi / 2) * (pi / (pi - 4)) * eps_3m) + (2 / pi))
    d_gamma_d_m3_eps <- eps_2m * sqrt(pi / 2) * (pi / (pi - 4)) * (2 / 3) * (sqrt(pi / 2) * (pi / (pi - 4)) * eps_3m)^(-5 / 3) * (eps_3m * (sqrt(pi / 2) * (pi / (pi - 4)) * eps_3m)^(-2 / 3) + (2 / pi))^(-2)

    d_beta_0_d_m3_alp <- (pi / (pi - 4)) * (1 / 3) * (sqrt(pi / 2) * (pi / (pi - 4)) * alp_3m)^(-2 / 3)
    d_sigma2_d_m3_alp <- sqrt(2 / pi) * (pi / (pi - 4)) * (2 / 3) * (sqrt(pi / 2) * (pi / (pi - 4)) * alp_3m)^(-1 / 3)
    d_gamma_d_m2_alp <- -(sqrt(pi / 2) * (pi / (pi - 4)) * alp_3m)^(-2 / 3) * (alp_3m * (sqrt(pi / 2) * (pi / (pi - 4)) * alp_3m) + (2 / pi))
    d_gamma_d_m3_alp <- alp_2m * sqrt(pi / 2) * (pi / (pi - 4)) * (2 / 3) * (sqrt(pi / 2) * (pi / (pi - 4)) * alp_3m)^(-5 / 3) * (alp_3m * (sqrt(pi / 2) * (pi / (pi - 4)) * alp_3m)^(-2 / 3) + (2 / pi))^(-2)

    beta_0_se <- sqrt(beta_se[1] + d_beta_0_d_m3_eps^2 * var_3m_eps + d_beta_0_d_m3_alp^2 * var_3m_alp)

    sigmaSq_uv_se <- sqrt(var_2m_eps + d_sigma2_d_m3_eps^2 * var_3m_eps + d_sigma2_d_m3_eps * cov_23m_eps)
    gamma_uv_se <- sqrt(d_gamma_d_m2_eps^2 * var_2m_eps + d_gamma_d_m3_eps^2 * var_3m_eps + d_gamma_d_m2_eps * d_gamma_d_m3_eps * cov_23m_eps)

    sigmaSq_hr_se <- sqrt(var_2m_alp + d_sigma2_d_m3_alp^2 * var_3m_alp + d_sigma2_d_m3_alp * cov_23m_alp)
    gamma_hr_se <- sqrt(d_gamma_d_m2_alp^2 * var_2m_alp + d_gamma_d_m3_alp^2 * var_3m_alp + d_gamma_d_m2_alp * d_gamma_d_m3_alp * cov_23m_alp)

    start_v <- if (is.na(beta_0_st)) {
      unname(c(gamma_uv, sigmaSq_uv, gamma_hr, sigmaSq_hr, beta_hat))
    } else {
      unname(c(gamma_uv, sigmaSq_uv, gamma_hr, sigmaSq_hr, beta_0, beta_hat))
    }
    End.Time <- end.time(Start.Time)
    out <- matrix(0, nrow = 3, ncol = length(start_v))
    rownames(out) <- c("par", "st_err", "t-val")
    colnames(out) <- if (isTRUE(intercept == 0)) {
      c("gamma_uv", "sigmaSq_uv", "gamma_hr", "sigmaSq_hr", colnames(data_x))
    } else {
      c("gamma_uv", "sigmaSq_uv", "gamma_hr", "sigmaSq_hr", colnames(data_x))
    }

    st_err <- rep(NA, ncol(out))
    st_err <- if (is.na(beta_0_st)) {
      c(gamma_uv_se, sigmaSq_uv_se, gamma_hr_se, sigmaSq_hr_se, summary(plm_gtre)$coefficients[c(x_vars_vec), 2])
    } else {
      c(gamma_uv_se, sigmaSq_uv_se, gamma_hr_se, sigmaSq_hr_se, beta_0_se, summary(plm_gtre)$coefficients[-c(1), 2])
    }
    t_val <- start_v / st_err
    out[1, ] <- start_v
    out[2, ] <- st_err
    out[3, ] <- t_val

    ## look at individual sigmas
    Sigma_u <- sqrt(gamma_uv * sigmaSq_uv)
    Sigma_v <- sqrt((1 - gamma_uv) * sigmaSq_uv)
    Sigma_h <- sqrt(gamma_hr * sigmaSq_hr)
    Sigma_r <- sqrt((1 - gamma_hr) * sigmaSq_hr)
    Lambda <- sigma_u / sigma_v
    Sigma <- sqrt(sigmaSq_uv)

    other_parms <- as.matrix(c(Sigma_u, Sigma_v, Sigma_h, Sigma_r, Lambda, Sigma))
    rownames(other_parms) <- c("sigma_u", "sigma_v", "sigma_h", "sigma_r", "lambda", "sigma")

    results <- list(t(out), End.Time, other_parms, model_name, formula, data, out["par", ], out["st_err", ], out["t-val", ], call)
    class(results) <- "sfareg"
    names(results) <- c("out", "total_time", "other_parms", "model_name", "formula", "data", "coefficients", "std.errors", "t.values", "call")
    return(results)
  } else {
    stop(paste0(
      "model_name '", model_name, "' is a recognized choice for psfm() but has no implementation branch. ",
      "Valid choices are: \"TRE_Z\", \"GTRE_Z\", \"TRE\", \"GTRE\", \"GTRE_FML\", \"TFE\", \"TFE_WMLE\", \"FD\", \"GTRE_SEQ1\", \"GTRE_SEQ2\", \"SSFE\", \"PL80\", \"BC92\", \"K1990\", \"K1990modified\"."
    ), call. = FALSE)
  }
}
