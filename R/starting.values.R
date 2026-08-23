start_cs <- function(formula_x, data_orig, x_vars_vec, intercept, model_name, n_x_vars, start_val, n_z_vars, z_vars) {
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
  ## THT's third parameter is the degrees of freedom of the skew-t. Starting it
  ## at 1 (Cauchy -- no mean, no variance) puts the optimizer in a region where
  ## the model's moments do not exist and a long way from any plausible fit. Take
  ## a moment start instead: the excess kurtosis of a t_nu is 6/(nu-4), so
  ## nu ~ 4 + 6/kurtosis of the OLS residuals. Clipped to [3, 30] -- above ~30 the
  ## skew-t is numerically indistinguishable from the skew-normal, so there is no
  ## information left to move on, and the likelihood is flat.
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
  start_v_nnak <- if (is.na(beta_0_st)) {
    unname(c(sigma_v, sigma_u, 0.5, beta_hat))
  } else {
    unname(c(sigma_v, sigma_u, 0.5, beta_0, beta_hat))
  }
  start_v_ne <- if (is.na(beta_0_st)) {
    unname(c(sigma_v, sigma_u, beta_hat))
  } else {
    unname(c(sigma_v, sigma_u, beta_0, beta_hat))
  }
  start_v_nhn <- if (is.na(beta_0_st)) {
    unname(c(lambda, sigma, beta_hat))
  } else {
    unname(c(lambda, sigma, beta_0, beta_hat))
  }
  ## NR is started from the Rayleigh moment equations, not the flat 0.1: from the
  ## flat start it reached a WORSE point than the true parameters in 9 of 14
  ## replications at n = 4000. .nr_start() returns NULL on wrongly skewed
  ## residuals, where the moment equations have no admissible solution, and the
  ## flat start is then used as before. See .nr_start() in matrix_utils.R.
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
  ## NU / NGE use residual-scaled starting values rather than the flat 0.1
  ## the older models use. Both have a one-sided parameter measured on the
  ## scale of the data (theta is the SUPPORT WIDTH of u for the uniform, and
  ## sigma_u the exponential mean for the generalized exponential), so a fixed
  ## 0.1 is a poor start whenever y is not O(1) -- residual dispersion gives
  ## the optimizer a start of roughly the right magnitude on any scale.
  ## Simulated-ML models. Third parameter is the lognormal meanlog (may be
  ## negative, hence the -Inf lower bound) or the Weibull shape (strictly
  ## positive, started at 1 = the exponential special case).
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
    ## built as c(sigma_u, sigma_v, ...) and THT's likelihood in sfm.R reads
    ## sig_u <- x[1]; sig_v <- x[2]. These labels used to be the other way round,
    ## so every THT fit reported each scale parameter under the other's name.
    ## Caught by the convergence sweep: the column labelled "sigv" converged to
    ## the true sigma_u and vice versa.
    colnames(out) <- c("sigu", "sigv", "a", c(names(plm_lm$coefficients)))
    ## 2.05 on the df for the same reason as in lower.start(): below 2 the
    ## skew-t has no variance, and below 1 no mean at all.
    lower_bob <- c(rep(.Machine$double.eps, 2), 2.05, rep(-Inf, n_x_vars))
  }
  if (model_name == "tHN") {
    ## Residual-scaled starts, not the flat 0.1 the older models use: both scales
    ## are measured in the units of y, so 0.1 is a poor start whenever y is not
    ## O(1). Same reasoning as NU/NGE/NLN/NW above.
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
    ## which reads sig_v <- x[1]; sig_u <- x[2]; nu <- x[3]. THT's labels were
    ## once transposed relative to its likelihood, so every THT fit reported each
    ## scale under the other name until a convergence sweep caught it; the
    ## ordering here is asserted in tests/testthat/test-thn.R so that cannot
    ## recur silently.
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
    start_v_t, start_v_ne, start_v_nr, start_v_nhn, start_v, out,
    lower_bob
  )

  names(results) <- c(
    "plm_lm", "beta_hat", "epsilon_hat", "beta_0_st", "sigma_u", "sigma_v", "mu",
    "beta_0", "lambda", "sigma", "start_v_ntn", "start_v_ng", "start_v_nnak",
    "start_v_t", "start_v_ne", "start_v_nr", "start_v_nhn", "start_v", "out",
    "lower_bob"
  )

  return(results)
}

start_panel <- function(formula_x, data, model_name, start_val, intercept, x_vars_vec,
                        individual = NULL, collinear_chk = NULL) {
  sfa_eps <- sfa_alp <- exp_eta <- exp_u <- sigma_v <- sigma_u <-
    sigma_r <- sigma_h <- beta_0 <- lambda <- sigma <- start_v <- out <-
    plm_gtre <- beta_hat <- alpha_hat <- epsilon_hat <- beta_0_st <- beta_se <- NULL

  plm_tfe <- plm_fd <- NULL

  if (isFALSE(is.numeric(start_val))) {
    if (model_name %in% c("TFE", "TFE_WMLE", "FD")) {
      plm_tfe <- plm(formula_x, data, effect = "individual", model = "within")
      plm_fd <- plm(formula_x, data, effect = "individual", model = "pooling")
    } else {
      ## Guard the random-effects starting-value regression against a rank-deficient
      ## BETWEEN-individual design (see .check_collinearity() in matrix_utils.R for
      ## why the pooled design can be full rank while this one is not). Without
      ## this, plm fails inside its error-components step with an uninterpretable
      ## solve(crossprod(ZBeta)) LAPACK error, or -- worse -- silently returns a
      ## garbage inverse and therefore meaningless starting values.
      ## `collinear_chk` is supplied by psfm() when the between-individual design
      ## is rank deficient AND the user chose "start_only" (the default): the
      ## requested formula is still what gets estimated, but the offending columns
      ## are removed from THIS starting-value regression only. psfm() has already
      ## dealt with the "error" and "warn_drop" cases before calling us.
      formula_start <- formula_x
      if (!is.null(collinear_chk) && length(collinear_chk$between_drop)) {
        warning(.collinearity_message(collinear_chk, "start_only"), call. = FALSE)
        keep_terms <- setdiff(
          attr(stats::terms(formula_x), "term.labels"),
          .terms_for_columns(formula_x, data, collinear_chk$between_drop)
        )
        formula_start <- stats::reformulate(if (length(keep_terms)) keep_terms else "1",
          response = all.vars(formula_x)[1]
        )
      }

      plm_gtre <- plm(formula_start, data, effect = "individual", model = "random")
      beta_hat_raw <- if (isTRUE(intercept == 0)) {
        plm(formula_start, data, effect = "individual")$coefficients
      } else {
        plm_gtre$coefficients[-c(1)]
      }
      ## Re-expand to the FULL requested coefficient vector, filling any column the
      ## starting-value regression could not identify with its pooled OLS estimate
      ## (a strictly better start than zero, and available at no extra cost).
      ## NOTE the target name set must match what the ORIGINAL code produced:
      ## x_vars_vec carries "(Intercept)" as its first entry, but the non-zero
      ## -intercept branch drops it via [-1]. Expanding to the wrong set makes
      ## start_v one element too long and blows up on colnames(out) downstream.
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
    }
  }

  if (isTRUE(is.numeric(start_val))) {
    start_v <- start_val
  }

  if (isFALSE(is.numeric(start_val)) & model_name %in% c("TRE_Z", "GTRE_Z", "TRE", "GTRE", "GTRE_FML", "GTRE_SEQ1", "GTRE_SEQ2")) {
    ## Variance components for the starting values come from the SAME sequential
    ## decomposition that psfm(model_name = "GTRE_SEQ1") reports: fit an
    ## intercept-only normal/half-normal to the composite residuals (giving
    ## sigma_v, sigma_u) and another to the individual effects (giving sigma_r,
    ## sigma_h).
    ##
    ## This replaces pcs_c(), which did the same decomposition but seeded its own
    ## optimizer with runif() draws. That made every GTRE-family starting value
    ## partly RANDOM: two identical calls on identical data began from different
    ## points, so run time and occasionally the optimum itself varied for no
    ## modelled reason. .fit_nhn_intercept() is deterministic, starts from
    ## sd(y)/sqrt(2) rather than a random draw, and runs the same staged
    ## optimizer the rest of the package uses.
    ##
    ## It also returns (sigma_v, sigma_u) DIRECTLY rather than as a
    ## lambda/sigma pair needing inversion, so the algebra that used to convert
    ## back -- sigma_v = sqrt(sigma^2/(1+lambda^2)), sigma_u = lambda*sigma_v --
    ## is gone along with its opportunity for error.
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
    plm_gtre, beta_hat, alpha_hat, epsilon_hat, beta_0_st, beta_se,
    plm_tfe, plm_fd
  )

  names(results) <- c(
    "sfa_eps", "sfa_alp", "exp_eta", "exp_u", "sigma_v", "sigma_u", "sigma_r",
    "sigma_h", "beta_0", "lambda", "sigma", "start_v", "out",
    "plm_gtre", "beta_hat", "alpha_hat", "epsilon_hat", "beta_0_st", "beta_se",
    "plm_tfe", "plm_fd"
  )
  return(results)
}

start.tfe <- function(formula_x, data, model_name, start_val, intercept, x_vars_vec, gamma, individual, N, y_var, n_x_vars) {
  ## `index = individual` is REQUIRED here and must not be dropped. By this
  ## point `data` has been through data_proc2(), which returns a plain
  ## data.frame -- the pdata.frame class (and with it the panel index) is gone,
  ## even though the individual columns are still pseries. Called without an
  ## explicit index, plm() falls back to treating the FIRST TWO COLUMNS as the
  ## individual and time index. That silently consumed the response and the
  ## first regressor, producing "empty model" (with a cryptic "'-' not
  ## meaningful for factors" warning) for any data whose first two columns are
  ## not the panel index. The bug was invisible for years because
  ## data_gen_p() happens to put `name` and `year` first.
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
    ## inside psfm.R's TFE like.tfe()/fn1(), which used to call demean() on
    ## this same per-individual regressor data on EVERY likelihood evaluation
    ## even though it doesn't depend on the parameter vector being optimized --
    ## a pure waste recomputed thousands of times across bobyqa/psoptim/optim.
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
  ## THT: the third parameter is degrees of freedom, not a scale. A lower bound
  ## of 1e-7 lets the optimizer walk into df < 1, where the skew-t has no mean and
  ## the efficiency predictor E[u|e] does not exist. Bound it at 2.05 so both the
  ## mean and the variance of the composed error exist throughout the search.
  if (model_name == "THT") {
    lower1 <- c(rep(.0000001, 2), 2.05, start_v[-c(1:3)] - differ)
  }
  if (model_name == "tHN") {
    lower1 <- c(rep(.0000001, 2), 2.05, start_v[-c(1:3)] - differ)
  }
  if (model_name %in% c("NG", "NNAK", "NW")) {
    lower1 <- c(rep(.0000001, 3), start_v[-c(1:3)] - differ)
  }
  ## NLN's third parameter is a meanlog and is genuinely unbounded below.
  if (model_name == "NLN") {
    lower1 <- c(rep(.0000001, 2), start_v[3] - differ, start_v[-c(1:3)] - differ)
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
  ## Two upper-bound vectors, because the optimizer stages need different things.
  ##
  ## `upper1` keeps the original trust region of +/- differ around the CURRENT
  ## point for every parameter. psoptim REQUIRES finite bounds -- it errors with
  ## "fixed bounds must be provided" otherwise -- so the particle swarm must use
  ## this one.
  ##
  ## `upper1_open` relaxes the bound to Inf on the strictly-positive scale
  ## parameters, and is for the final optim() stage, which handles Inf fine.
  ## Applying a trust region to a variance parameter there is an outright trap:
  ## the box is anchored to wherever the PREVIOUS stage happened to stop, so a
  ## stage that drove sigma_v toward zero leaves the next one an upper bound of
  ## 0 + differ, which it pins against and cannot escape.
  ##
  ## That is what broke NR. On identical data it reached log-likelihoods up to
  ## 107 points BELOW NHN -- the same model by a different closed form -- and
  ## returned sigma_v = 0.500 exactly (= 0 + differ) on 2 of 5 test seeds, with
  ## its variance-parameter MSE flat in n. Opening the bound recovered 69 and 58
  ## log-likelihood points on those seeds.
  ##
  ## The positions needing this are exactly those given a MIN_POSITIVE lower
  ## bound above, so they are detected from lower1 rather than re-listed per
  ## model, which keeps the two in step as models are added. Parameters with a
  ## genuine two-sided range (NLN's meanlog, ZISF's gamma) keep the finite box.
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
