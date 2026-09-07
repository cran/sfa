start_cs <- function(formula_x, data_orig, x_vars_vec, intercept, model_name, n_x_vars, start_val, n_z_vars, z_vars, n_class = 1) {
  plm_lm <- lm(formula_x, data_orig)
  beta_hat <- if (isTRUE(intercept == 0)) {
    plm_lm$coefficients[x_vars_vec]
  } else {
    plm_lm$coefficients[x_vars_vec][-1]
  }
  epsilon_hat <- plm_lm$residuals
  beta_0_st <- if (isTRUE(intercept == 0)) {
    NA
  } else {
    plm_lm$coefficients[c(1)]
  }
  sigma_u <- 0.1
  sigma_v <- 0.1
  mu <- 0.1
  beta_0 <- beta_0_st
  lambda <- sigma_u / sigma_v
  sigma <- sqrt(sigma_u^2 + sigma_v^2)
  gam <- 1

  start_v_ntn <- if (is.na(beta_0_st)) {
    unname(c(lambda, sigma, mu, beta_hat))
  } else {
    unname(c(lambda, sigma, mu, beta_0, beta_hat))
  }
  ## THT's third parameter is the degrees of freedom of the skew-t.
  .exkurt <- mean((epsilon_hat - mean(epsilon_hat))^4) / stats::var(epsilon_hat)^2 - 3
  a_start <- if (is.finite(.exkurt) && .exkurt > 0.2) min(max(4 + 6 / .exkurt, 3), 30) else 10
  start_v_t <- if (is.na(beta_0_st)) {
    unname(c(sigma_u, sigma_v, a_start, beta_hat))
  } else {
    unname(c(sigma_u, sigma_v, a_start, beta_0, beta_hat))
  }
  ## tHN uses the SAME kurtosis-based df start as THT, but the conventional
  ## (sigma_v, sigma_u) order rather than THT's inverted one.
  start_v_thn <- if (is.na(beta_0_st)) {
    unname(c(sigma_v, sigma_u, a_start, beta_hat))
  } else {
    unname(c(sigma_v, sigma_u, a_start, beta_0, beta_hat))
  }
  start_v_ng <- if (is.na(beta_0_st)) {
    unname(c(sigma_v, sigma_u, 1, beta_hat))
  } else {
    unname(c(sigma_v, sigma_u, 1, beta_0, beta_hat))
  }
  ## NNAK starts from the method of moments, not from the constants above.
  .moment_start <- local({
    mm <- tryCatch(stats::model.matrix(plm_lm), error = function(e) NULL)
    yy <- tryCatch(stats::model.response(stats::model.frame(plm_lm)),
                   error = function(e) NULL)
    if (is.null(mm) || is.null(yy)) {
      NULL
    } else {
      cf <- tryCatch(.cols_fit(yy, mm, "NHN", intercept_col = NA), error = function(e) NULL)
      if (is.null(cf) || isTRUE(cf$wrong_skew) ||
        !is.finite(cf$sigma_u) || !is.finite(cf$sigma_v) ||
        cf$sigma_u <= 0 || cf$sigma_v <= 0) {
        NULL
      } else {
        cf
      }
    }
  })

  ## TSL carries (sigma_v, sigma_u, lambda), and takes its scales from the
  ## EXPONENTIAL moment inversion rather than the half-normal one used above.
  .tsl_start <- local({
    mm <- tryCatch(stats::model.matrix(plm_lm), error = function(e) NULL)
    yy <- tryCatch(stats::model.response(stats::model.frame(plm_lm)),
                   error = function(e) NULL)
    if (is.null(mm) || is.null(yy)) {
      NULL
    } else {
      cf <- tryCatch(.cols_fit(yy, mm, "NE", intercept_col = NA), error = function(e) NULL)
      if (is.null(cf) || isTRUE(cf$wrong_skew) || !is.finite(cf$sigma_u) || cf$sigma_u <= 0) {
        NULL
      } else {
        sv_tsl <- sqrt(max(abs(unname(cf$moments[["m2"]]) - cf$sigma_u^2),
                           .SFA_CONSTANTS$MIN_POSITIVE))
        list(sigma_u = cf$sigma_u, sigma_v = sv_tsl, eu = cf$eu)
      }
    }
  })

  start_v_tsl <- if (!is.null(.tsl_start)) {
    if (is.na(beta_0_st)) {
      unname(c(.tsl_start$sigma_v, .tsl_start$sigma_u, 1, beta_hat))
    } else {
      unname(c(.tsl_start$sigma_v, .tsl_start$sigma_u, 1,
               beta_0_st + .tsl_start$eu, beta_hat))
    }
  } else if (is.na(beta_0_st)) {
    unname(c(sigma_v, sigma_u, 1, beta_hat))
  } else {
    unname(c(sigma_v, sigma_u, 1, beta_0, beta_hat))
  }

  start_v_nnak <- if (!is.null(.moment_start)) {
    ## Only the intercept moves; the OLS slopes are already consistent.
    if (is.na(beta_0_st)) {
      unname(c(.moment_start$sigma_v, .moment_start$sigma_u, 0.5, beta_hat))
    } else {
      unname(c(.moment_start$sigma_v, .moment_start$sigma_u, 0.5,
               beta_0_st + .moment_start$eu, beta_hat))
    }
  } else if (is.na(beta_0_st)) {
    unname(c(sigma_v, sigma_u, 0.5, beta_hat))
  } else {
    unname(c(sigma_v, sigma_u, 0.5, beta_0, beta_hat))
  }
  ## NE is started from the bias-corrected moment estimator, not the flat 0.1
  ## above; see R/ne_start.R.  Falls back to the flat start only if the
  ## residuals are degenerate.
  start_v_ne <- tryCatch(
    .ne_start(epsilon_hat, beta_0_st, beta_hat, rule = "bc"),
    error = function(e) {
      if (is.na(beta_0_st)) {
        unname(c(sigma_v, sigma_u, beta_hat))
      } else {
        unname(c(sigma_v, sigma_u, beta_0, beta_hat))
      }
    }
  )
  start_v_nhn <- if (is.na(beta_0_st)) {
    unname(c(lambda, sigma, beta_hat))
  } else {
    unname(c(lambda, sigma, beta_0, beta_hat))
  }
  ## NR is started from the Rayleigh moment equations, not the flat 0.1.
  .nr_mom <- .nr_start(epsilon_hat, beta_0_st, beta_hat)
  start_v_nr <- if (!is.null(.nr_mom)) {
    .nr_mom
  } else if (is.na(beta_0_st)) {
    unname(c(sigma_v, sigma_u, beta_hat))
  } else {
    unname(c(sigma_v, sigma_u, beta_0, beta_hat))
  }
  start_v_zisf <- if (is.na(beta_0_st)) {
    unname(c(gam, sigma_v, sigma_u, beta_hat))
  } else {
    unname(c(gam, sigma_v, sigma_u, beta_0, beta_hat))
  }

  if (model_name %in% c("NHN_Z", "NE_Z")) {
    start_v_nez <- if (is.na(beta_0_st)) {
      unname(c(sigma_v, beta_hat, rep(mu, n_z_vars)))
    } else {
      unname(c(sigma_v, beta_0, beta_hat, rep(mu, n_z_vars)))
    }
    start_v <- start_v_nez
    out <- matrix(0, nrow = 3, ncol = length(start_v))
    colnames(out) <- c("sigma_v", c(names(plm_lm$coefficients)), z_vars)
    lower_bob <- c(.Machine$double.eps, rep(-.Machine$double.xmax^.1, length(start_v[-c(1)])))
  }
  ## LATENT CLASS. Unlike every other model here the parameter vector has no
  ## fixed length: J blocks of (sigv, sigu, beta) followed by the (J-1) blocks
  ## of multinomial-logit coefficients, class J being the reference.
  ##
  ## The starting values are NOT J perturbations of one fit. A finite mixture
  ## whose components start identical sits at a saddle point of the likelihood
  ## -- the posterior class probabilities are then equal for every observation
  ## and the score with respect to the class split is zero -- so the classes
  ## have to be separated before the optimizer is handed the problem. Splitting
  ## the OLS residuals at their J-quantiles and refitting within each group is
  ## the usual EM-style initialisation and costs one lm() per class. It also
  ## makes the labelling reproducible rather than arbitrary: class 1 starts on
  ## the lowest-residual group. See the note on label switching in ?zsfm.
  ## LCM_CN: the CONTAMINATED NORMAL frontier. Every parameter is common
  ## across components except the noise scale, so the layout is
  ## [sigv_1..sigv_J, sigu, beta, logit], not one block per class. Start the
  ## noise scales spread around the pooled one -- a contaminating component is
  ## by definition the wider one -- and everything else at its pooled value.
  if (identical(model_name, "LCM_CN")) {
    .J <- max(2L, as.integer(n_class))
    e_lc <- as.numeric(epsilon_hat)
    s0 <- stats::sd(e_lc)
    if (!is.finite(s0) || s0 <= 0) s0 <- 1
    spread <- seq(0.6, 1.8, length.out = .J)
    lam0 <- 1
    sigv0 <- s0 / sqrt(1 + lam0^2) * spread
    sigu0 <- s0 * lam0 / sqrt(1 + lam0^2)
    b0 <- unname(plm_lm$coefficients)
    b0[!is.finite(b0)] <- 0
    start_v <- c(sigv0, sigu0, b0, rep(0, (.J - 1L)))
    .bnames <- names(plm_lm$coefficients)
    out <- matrix(0, nrow = 3, ncol = length(start_v))
    colnames(out) <- c(
      paste0("sigv_class", seq_len(.J)), "sigu", .bnames,
      paste0("logit_(Intercept)_class", seq_len(.J - 1L))
    )
    lower_bob <- c(
      rep(.Machine$double.eps, .J + 1L),
      rep(-Inf, length(b0) + (.J - 1L))
    )
  }

  if (model_name %in% c("LCM", "LCM_Z")) {
    .J <- max(2L, as.integer(n_class))
    X_lc <- stats::model.matrix(plm_lm)
    e_lc <- as.numeric(epsilon_hat)
    qs <- stats::quantile(e_lc, probs = seq(0, 1, length.out = .J + 1L), names = FALSE)
    ## Ties in the residuals can collapse a break; jitter the interior breaks
    ## apart rather than letting cut() drop a class.
    qs <- sort(unique(qs))
    grp <- if (length(qs) < .J + 1L) {
      ## Degenerate residual distribution: fall back to an even split of the
      ## ranks, which always yields J non-empty groups.
      cut(rank(e_lc, ties.method = "first"),
        breaks = .J, labels = FALSE, include.lowest = TRUE)
    } else {
      cut(e_lc, breaks = qs, labels = FALSE, include.lowest = TRUE)
    }
    ## One sigma per class from that class's own residual spread, split into
    ## sigv/sigu by the pooled lambda so the two scales start on the same
    ## footing as the single-class models.
    .lam0 <- 1
    .blocks <- lapply(seq_len(.J), function(j) {
      ok <- which(grp == j)
      b_j <- tryCatch(
        stats::lm.fit(X_lc[ok, , drop = FALSE], plm_lm$model[[1L]][ok])$coefficients,
        error = function(e) plm_lm$coefficients
      )
      b_j[!is.finite(b_j)] <- plm_lm$coefficients[!is.finite(b_j)]
      s_j <- stats::sd(e_lc[ok])
      if (!is.finite(s_j) || s_j <= 0) s_j <- max(stats::sd(e_lc), 1e-3)
      c(s_j / sqrt(1 + .lam0^2), s_j * .lam0 / sqrt(1 + .lam0^2), unname(b_j))
    })
    ## Class J is the reference, so only J-1 blocks of logit coefficients are
    ## free. They start at zero: equal prior class probabilities. The classes
    ## are already separated by the beta blocks above, so this is not the
    ## saddle point the identical-component start would be.
    n_q <- if (model_name == "LCM") 1L else n_z_vars
    start_v <- c(unlist(.blocks), rep(0, (.J - 1L) * n_q))
    .bnames <- names(plm_lm$coefficients)
    .qnames <- if (model_name == "LCM") "(Intercept)" else z_vars
    out <- matrix(0, nrow = 3, ncol = length(start_v))
    colnames(out) <- c(
      unlist(lapply(seq_len(.J), function(j) {
        paste0(c("sigv", "sigu", .bnames), "_class", j)
      })),
      unlist(lapply(seq_len(.J - 1L), function(j) {
        paste0("logit_", .qnames, "_class", j)
      }))
    )
    ## Positivity on the two scales of every block; the betas and the logit
    ## coefficients are unrestricted.
    lower_bob <- c(
      rep(c(.Machine$double.eps, .Machine$double.eps, rep(-Inf, n_x_vars)), .J),
      rep(-Inf, (.J - 1L) * n_q)
    )
  }
  if (model_name %in% c("ZISF")) {
    start_v <- start_v_zisf
    out <- matrix(0, nrow = 3, ncol = length(start_v))
    colnames(out) <- c("gamma", "sigv", "sigu", c(names(plm_lm$coefficients)))
    lower_bob <- c(-Inf, rep(.Machine$double.eps, 2), rep(-Inf, n_x_vars))
  }
  if (model_name %in% c("ZISF_Z")) {
    start_v_zisfz <- if (is.na(beta_0_st)) {
      unname(c(sigma_v, sigma_u, beta_hat, rep(mu, n_z_vars)))
    } else {
      unname(c(sigma_v, sigma_u, beta_0, beta_hat, rep(mu, n_z_vars)))
    }
    start_v <- start_v_zisfz
    out <- matrix(0, nrow = 3, ncol = length(start_v))
    colnames(out) <- c("sigv", "sigu", c(names(plm_lm$coefficients)), z_vars)
    lower_bob <- c(rep(.Machine$double.eps, 2), rep(-Inf, n_x_vars + n_z_vars))
  }
  if (model_name %in% c("NHN")) {
    start_v <- start_v_nhn
    out <- matrix(0, nrow = 3, ncol = length(start_v))
    colnames(out) <- c("lambda", "sigma", c(names(plm_lm$coefficients)))
    lower_bob <- c(rep(.Machine$double.eps, 2), rep(-Inf, n_x_vars))
  }
  if (model_name == "NE") {
    start_v <- start_v_ne
    out <- matrix(0, nrow = 3, ncol = length(start_v))
    colnames(out) <- c("sigv", "sigu", c(names(plm_lm$coefficients)))
    lower_bob <- c(rep(.Machine$double.eps, 2), rep(-Inf, n_x_vars))
  }
  if (model_name == "NR") {
    start_v <- start_v_nr
    out <- matrix(0, nrow = 3, ncol = length(start_v))
    colnames(out) <- c("sigv", "sigu", c(names(plm_lm$coefficients)))
    lower_bob <- c(rep(.Machine$double.eps, 2), rep(-Inf, n_x_vars))
  }
  ## NU / NGE use residual-scaled starting values rather than the flat 0.1 the
  ## older models use.
  if (model_name %in% c("NLN", "NW")) {
    s_eps <- stats::sd(epsilon_hat)
    sv_st <- max(0.5 * s_eps, 1e-3)
    su_st <- max(s_eps, 1e-3)
    th_st <- if (model_name == "NLN") log(max(s_eps, 1e-3)) else 1
    start_v <- if (is.na(beta_0_st)) {
      unname(c(sv_st, su_st, th_st, beta_hat))
    } else {
      unname(c(sv_st, su_st, th_st, beta_0, beta_hat))
    }
    out <- matrix(0, nrow = 3, ncol = length(start_v))
    colnames(out) <- c(
      "sigv", "sigu", if (model_name == "NLN") "mu" else "k",
      c(names(plm_lm$coefficients))
    )
    lower_bob <- c(
      rep(.Machine$double.eps, 2),
      if (model_name == "NLN") -Inf else .Machine$double.eps,
      rep(-Inf, n_x_vars)
    )
  }
  if (model_name %in% c("NU", "NGE")) {
    s_eps <- stats::sd(epsilon_hat)
    sv_st <- max(0.5 * s_eps, 1e-3)
    su_st <- max(if (model_name == "NU") 2 * s_eps else s_eps, 1e-3)
    start_v <- if (is.na(beta_0_st)) {
      unname(c(sv_st, su_st, beta_hat))
    } else {
      unname(c(sv_st, su_st, beta_0, beta_hat))
    }
    out <- matrix(0, nrow = 3, ncol = length(start_v))
    colnames(out) <- c(
      "sigv", if (model_name == "NU") "theta" else "sigu",
      c(names(plm_lm$coefficients))
    )
    lower_bob <- c(rep(.Machine$double.eps, 2), rep(-Inf, n_x_vars))
  }
  if (model_name == "THT") {
    start_v <- start_v_t
    out <- matrix(0, nrow = 3, ncol = length(start_v))
    ## Order is (sigma_u, sigma_v), NOT (sigma_v, sigma_u): start_v_t above is
    ## built as c.
    colnames(out) <- c("sigu", "sigv", "a", c(names(plm_lm$coefficients)))
    ## 2.05 on the df for the same reason as in lower.start(): below 2 the
    ## skew-t has no variance, and below 1 no mean at all.
    lower_bob <- c(rep(.Machine$double.eps, 2), 2.05, rep(-Inf, n_x_vars))
  }
  if (model_name == "tHN") {
    ## Residual-scaled starts, not the flat 0.1 the older models use: both
    ## scales are measured in the units of y.
    s_eps <- stats::sd(epsilon_hat)
    sv_st <- max(0.5 * s_eps, 1e-3)
    su_st <- max(s_eps, 1e-3)
    start_v_thn <- if (is.na(beta_0_st)) {
      unname(c(sv_st, su_st, a_start, beta_hat))
    } else {
      unname(c(sv_st, su_st, a_start, beta_0, beta_hat))
    }
    start_v <- start_v_thn
    out <- matrix(0, nrow = 3, ncol = length(start_v))
    ## Order is (sigma_v, sigma_u, nu) and MUST match sfm.R's tHN likelihood,
    ## which reads sig_v <- x[1]; sig_u <- x[2]; nu <- x[3].
    colnames(out) <- c("sigv", "sigu", "nu", c(names(plm_lm$coefficients)))
    ## nu floor of 2.05: the t has no variance at nu <= 2, and the quadrature in
    ## .log_d_thn() is only validated down to 2.05.
    lower_bob <- c(rep(.Machine$double.eps, 2), 2.05, rep(-Inf, n_x_vars))
  }
  if (model_name == "NG") {
    start_v <- start_v_ng
    out <- matrix(0, nrow = 3, ncol = length(start_v))
    colnames(out) <- c("sigv", "sigu", "mu", c(names(plm_lm$coefficients)))
    lower_bob <- c(rep(.Machine$double.eps, 3), rep(-Inf, n_x_vars))
  }
  if (model_name == "TSL") {
    start_v <- start_v_tsl
    out <- matrix(0, nrow = 3, ncol = length(start_v))
    colnames(out) <- c("sigv", "sigu", "lambda", c(names(plm_lm$coefficients)))
    lower_bob <- c(rep(.Machine$double.eps, 3), rep(-Inf, n_x_vars))
  }
  if (model_name == "NNAK") {
    start_v <- start_v_nnak
    out <- matrix(0, nrow = 3, ncol = length(start_v))
    colnames(out) <- c("sigv", "sigu", "mu", c(names(plm_lm$coefficients)))
    lower_bob <- c(rep(.Machine$double.eps, 3), rep(-Inf, n_x_vars))
  }
  if (model_name == "NTN") {
    start_v <- start_v_ntn
    out <- matrix(0, nrow = 3, ncol = length(start_v))
    colnames(out) <- c("lambda", "sigma", "mu", c(names(plm_lm$coefficients)))
    lower_bob <- c(rep(.Machine$double.eps, 2), rep(-Inf, n_x_vars + 1))
  }

  if (isTRUE(is.numeric(start_val))) {
    start_v <- start_val
  }

  rownames(out) <- c("par", "st_err", "t-val")

  results <- list(
    plm_lm, beta_hat, epsilon_hat, beta_0_st, sigma_u, sigma_v, mu,
    beta_0, lambda, sigma, start_v_ntn, start_v_ng, start_v_nnak,
    start_v_t, start_v_ne, start_v_nr, start_v_nhn, start_v_tsl, start_v, out,
    lower_bob
  )

  names(results) <- c(
    "plm_lm", "beta_hat", "epsilon_hat", "beta_0_st", "sigma_u", "sigma_v", "mu",
    "beta_0", "lambda", "sigma", "start_v_ntn", "start_v_ng", "start_v_nnak",
    "start_v_t", "start_v_ne", "start_v_nr", "start_v_nhn", "start_v_tsl", "start_v", "out",
    "lower_bob"
  )

  return(results)
}

start_panel <- function(formula_x, data, model_name, start_val, intercept, x_vars_vec,
                        individual = NULL, collinear_chk = NULL) {
  sfa_eps <- sfa_alp <- exp_eta <- exp_u <- sigma_v <- sigma_u <-
    sigma_r <- sigma_h <- beta_0 <- lambda <- sigma <- start_v <- out <-
    plm_gtre <- beta_hat <- alpha_hat <- epsilon_hat <- beta_0_st <- beta_se <- NULL
  beta_se_named <- NULL

  plm_tfe <- plm_fd <- NULL

  if (isFALSE(is.numeric(start_val))) {
    if (model_name %in% c("TFE", "TFE_WMLE", "FD")) {
      plm_tfe <- plm(formula_x, data, effect = "individual", model = "within")
      plm_fd <- plm(formula_x, data, effect = "individual", model = "pooling")
    } else {
      ## Guard the random-effects starting-value regression against a
      ## rank-deficient BETWEEN-individual design. The reduction happens at
      ## COLUMN granularity (see .re_start_design): dropping whole terms cannot
      ## express "keep 13 of factor(year)'s 24 dummies", so it left the singular
      ## design untouched and plm::ercomp() then failed in solve().
      formula_start <- formula_x
      data_start <- data
      start_map <- NULL
      if (!is.null(collinear_chk) && length(collinear_chk$between_drop)) {
        if (!isTRUE(collinear_chk$already_warned)) {
          warning(.collinearity_message(collinear_chk, "start_only"), call. = FALSE)
        }
        rd <- .re_start_design(formula_x, data, collinear_chk$between_drop)
        if (!is.null(rd)) {
          formula_start <- rd$formula
          data_start <- rd$data
          start_map <- rd$map
        }
      }

      plm_gtre <- plm(formula_start, data_start, effect = "individual", model = "random")
      beta_hat_raw <- if (isTRUE(intercept == 0)) {
        plm(formula_start, data_start, effect = "individual")$coefficients
      } else {
        plm_gtre$coefficients[-c(1)]
      }
      if (!is.null(start_map)) beta_hat_raw <- .remap_start_names(beta_hat_raw, start_map)
      ## Re-expand to the FULL requested coefficient vector.
      beta_target <- if (isTRUE(intercept == 0)) x_vars_vec else x_vars_vec[x_vars_vec != "(Intercept)"]
      beta_hat <- .expand_start_beta(beta_hat_raw, beta_target, formula_x, data)
      alpha_hat <- ranef(plm_gtre)
      epsilon_hat <- plm_gtre$residuals
      beta_0_st <- if (isTRUE(intercept == 0)) {
        NA
      } else {
        plm_gtre$coefficients[c(1)]
      }
      beta_se <- as.data.frame(summary(plm_gtre)[1])$coefficients.Std..Error
      ## Named copy, so SEQ1/SEQ2 can align SEs to the REQUESTED columns even
      ## when the starting-value regression was reduced.
      beta_se_named <- stats::setNames(
        summary(plm_gtre)$coefficients[, 2], rownames(summary(plm_gtre)$coefficients)
      )
      if (!is.null(start_map)) beta_se_named <- .remap_start_names(beta_se_named, start_map)
    }
  }

  if (isTRUE(is.numeric(start_val))) {
    start_v <- start_val
  }

  if (isFALSE(is.numeric(start_val)) & model_name %in% c("TRE_Z", "GTRE_Z", "TRE", "GTRE", "GTRE_FML", "GTRE_SEQ1", "GTRE_SEQ2")) {
    ## Variance components for the starting values come from the SAME
    ## sequential decomposition that psfm(model_name = "GTRE_SEQ1") reports.
    fit_eps_st <- .fit_nhn_intercept(as.numeric(epsilon_hat))
    fit_alp_st <- .fit_nhn_intercept(as.numeric(alpha_hat))
    sigma_v <- fit_eps_st$par[1]
    sigma_u <- fit_eps_st$par[2]
    sigma_r <- fit_alp_st$par[1]
    sigma_h <- fit_alp_st$par[2]
    exp_u <- fit_eps_st$par[3]
    exp_eta <- fit_alp_st$par[3]
    beta_0 <- beta_0_st + exp_u + exp_eta
    lambda <- sigma_u / sigma_v
    sigma <- sqrt(sigma_u^2 + sigma_v^2)
    if (model_name %in% c("GTRE_Z", "TRE_Z")) {
      beta_0 <- beta_0_st + exp_u
    }

    if (model_name %in% c("GTRE", "TRE")) {
      if (isTRUE(is.numeric(start_val))) {
        start_v <- start_val
      } else {
        start_v <- if (is.na(beta_0_st)) {
          unname(c(lambda, sigma, sigma_r, sigma_h, beta_hat))
        } else {
          unname(c(lambda, sigma, sigma_r, sigma_h, beta_0, beta_hat))
        }
      }

      out <- matrix(0, nrow = 3, ncol = length(start_v))
      rownames(out) <- c("par", "st_err", "t-val")
      if (model_name == "TRE" & isTRUE(is.numeric(start_val))) {
        colnames(out) <- c("lambda", "sig", "sig_r", x_vars_vec)
      } else {
        colnames(out) <- c("lambda", "sig", "sig_r", "sig_h", x_vars_vec)
      }

      if (model_name == "TRE" & isFALSE(is.numeric(start_val))) {
        out <- out[, -c(4)]
        start_v <- start_v[-c(4)]
      }
    }
  }

  results <- list(
    sfa_eps, sfa_alp, exp_eta, exp_u, sigma_v, sigma_u, sigma_r,
    sigma_h, beta_0, lambda, sigma, start_v, out,
    plm_gtre, beta_hat, alpha_hat, epsilon_hat, beta_0_st, beta_se, beta_se_named,
    plm_tfe, plm_fd
  )

  names(results) <- c(
    "sfa_eps", "sfa_alp", "exp_eta", "exp_u", "sigma_v", "sigma_u", "sigma_r",
    "sigma_h", "beta_0", "lambda", "sigma", "start_v", "out",
    "plm_gtre", "beta_hat", "alpha_hat", "epsilon_hat", "beta_0_st", "beta_se", "beta_se_named",
    "plm_tfe", "plm_fd"
  )
  return(results)
}

start.tfe <- function(formula_x, data, model_name, start_val, intercept, x_vars_vec, gamma, individual, N, y_var, n_x_vars) {
  ## `index = individual` is REQUIRED here and must not be dropped.
  plm_tfe <- plm(formula_x, data, effect = "individual", model = "within", index = individual)
  if (isTRUE(is.numeric(start_val))) {
    start_v <- start_val
  } else {
    beta_hat <- plm_tfe$coefficients[c(x_vars_vec)]
    epsilon_hat <- plm_tfe$residuals
    beta_se <- as.data.frame(summary(plm_tfe)[1])$coefficients.Std..Error[-c(1)]
    sfa_eps <- pcs_c(Y = as.numeric(epsilon_hat))[[1]]$par
    exp_u <- sfa_eps[3]
    sigma_v <- sqrt(sfa_eps[2]^2 / (1 + sfa_eps[1]^2))
    sigma_u <- sigma_v * sfa_eps[1]
    lambda <- sigma_u / sigma_v
    sigma <- sqrt(sigma_u^2 + sigma_v^2)
    start_v <- unname(c(lambda, sigma, beta_hat))
    start_v[1] <- if (gamma == TRUE) {
      sigma_u^2 / sigma^2
    } else {
      start_v[1]
    }
  }

  out <- matrix(0, nrow = 3, ncol = length(start_v))
  rownames(out) <- c("par", "st_err", "t-val")
  colnames(out) <- c("lambda", "sig", c(names(plm_tfe$coefficients)))
  colnames(out)[1] <- if (gamma == TRUE) {
    "gamma"
  } else {
    "lambda"
  }

  indiv <- noquote(as.vector(unique(data[, c(individual)])))
  t <- rep(0, N)
  data_i <- Y <- eps <- data_i_vars <- data_i_vars_dm <- one_t <- I_t <- one_t1 <- I_t1 <- as.list(rep(0, N))

  for (ii in seq_len(N)) {
    data_i[[ii]] <- data[which(data[, c(individual)] == indiv[ii]), ]
    t[ii] <- nrow(data_i[[ii]])
    Y[[ii]] <- .demean(as.numeric(data_i[[ii]][, y_var]))
    data_i_vars[[ii]] <- data.frame(data_i[[ii]][, c(x_vars_vec)])
    ## Precomputed once here (same treatment Y already got above) rather than
    ## inside psfm.R's TFE like.tfe()/fn1().
    data_i_vars_dm[[ii]] <- as.matrix(.demean(data_i_vars[[ii]]))
    one_t[[ii]] <- rep(1, t[[ii]])
    I_t[[ii]] <- diag(t[[ii]])
    one_t1[[ii]] <- rep(1, t[[ii]] - 1)
    I_t1[[ii]] <- diag(t[[ii]] - 1)
  }

  if (gamma == TRUE) {
    upper <- c(1, rep(Inf, n_x_vars + 1))
  } else {
    upper <- NA
  }

  results <- list(
    upper, I_t1, one_t1, I_t, one_t, data_i_vars, data_i_vars_dm, Y, t, data_i,
    eps, indiv, out, start_v
  )
  names(results) <- c(
    "upper", "I_t1", "one_t1", "I_t", "one_t", "data_i_vars", "data_i_vars_dm", "Y", "t", "data_i",
    "eps", "indiv", "out", "start_v"
  )
  return(results)
}

lower.start <- function(start_v, model_name, differ) {
  if (model_name == "TRE") {
    lower1 <- c(rep(.0000001, 3), start_v[-c(1:3)] - differ)
  }
  if (model_name == "GTRE") {
    lower1 <- c(rep(.0000001, 4), start_v[-c(1:4)] - differ)
  }
  if (model_name %in% c("NHN", "NE", "NR", "NTN", "NU", "NGE")) {
    lower1 <- c(rep(.0000001, 2), start_v[-c(1:2)] - differ)
  }
  ## THT: the third parameter is degrees of freedom, not a scale.
  if (model_name == "THT") {
    lower1 <- c(rep(.0000001, 2), 2.05, start_v[-c(1:3)] - differ)
  }
  if (model_name == "tHN") {
    lower1 <- c(rep(.0000001, 2), 2.05, start_v[-c(1:3)] - differ)
  }
  if (model_name %in% c("NG", "NNAK", "NW", "TSL")) {
    lower1 <- c(rep(.0000001, 3), start_v[-c(1:3)] - differ)
  }
  ## NLN's third parameter is a meanlog and is genuinely unbounded below.
  if (model_name == "NLN") {
    lower1 <- c(rep(.0000001, 2), start_v[3] - differ, start_v[-c(1:3)] - differ)
  }
  ## LCM's layout is variable-length, so the positivity pattern is derived from
  ## the block structure recorded on start_v rather than hard-coded by index.
  ## .lcm_pos is attached by zsfm() when it builds the starting vector.
  if (model_name %in% c("LCM", "LCM_Z")) {
    pos <- attr(start_v, "lcm_pos")
    if (is.null(pos)) {
      stop("lower.start(): LCM start_v is missing its `lcm_pos` attribute.",
        call. = FALSE)
    }
    lower1 <- ifelse(pos, .0000001, start_v - differ)
  }
  if (model_name %in% c("ZISF")) {
    lower1 <- c(start_v[1] - differ, rep(.0000001, 2), start_v[-c(1:3)] - differ)
  }
  if (model_name %in% c("ZISF_Z")) {
    lower1 <- c(rep(.0000001, 2), start_v[-c(1:2)] - differ)
  }
  if (model_name %in% c("NHN_Z", "NE_Z")) {
    lower1 <- c(rep(.0000001, 1), start_v[-c(1)] - differ)
  }
  ## Two upper-bound vectors, because the optimizer stages need different
  ## things.
  is_pos <- abs(lower1 - .0000001) < 1e-12
  upper1 <- c(start_v + differ)
  upper1_open <- ifelse(is_pos, Inf, start_v + differ)

  results <- list(upper1, lower1, upper1_open)
  names(results) <- c("upper1", "lower1", "upper1_open")
  return(results)
}

start.time <- function() {
  start_time <- Sys.time()
  return(start_time)
}

end.time <- function(start_time) {
  end_time <- Sys.time()
  total_time <- end_time - start_time
  return(total_time)
}
