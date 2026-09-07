## psfm_bootstrap() Parametric bootstrap for panel stochastic frontier models
## fit with psfm() from the sfa package.

## Both GTRE arms below reach this, and they must say the same thing: the
## "sml" route reports lambda/sigma with sigh and sigr alongside, the "fiml"
## route reports the four raw scales, and the warning is about the model rather
## than about either parameterization.
.warn_boot_boundary <- function(sig_h_hat, sigr_hat) {
  warning("psfm_bootstrap(): this fit has a persistent scale on the zero ",
    "boundary (sigh = ", signif(sig_h_hat, 4), ", sigr = ",
    signif(sigr_hat, 4), "). A parametric bootstrap draws from the fitted ",
    "model, so it will resample from a data-generating process with that ",
    "component absent, and the bootstrap is in any case inconsistent on the ",
    "boundary of the parameter space. The resulting intervals will understate ",
    "the uncertainty in the persistent split rather than represent it. See ",
    "?psfm on what a collapsed scale means -- it is frequently the correct ",
    "maximum likelihood estimate, but it is not a point to bootstrap around.",
    call. = FALSE
  )
}

psfm_bootstrap <- function(psfm_object,
                           numCores,
                           BOOT,
                           individual,
                           h_type = c("auto", "none", "scalar", "parametric"),
                           maxit.psoptim = 1000,
                           seed_offset = 0,
                           write_back = TRUE,
                           pkgs = c("sfa", "Formula", "pbapply", "truncnorm"),
                           inefdec,
                           rand.gtre = NULL,
                           rand.psoptim = NULL,
                           maxit.bobyqa = 1,
                           maxit.optim = 1) {
  h_type <- match.arg(h_type)

  ## Restore the caller's random stream on the way out.
  .rng_state <- .rng_snapshot()
  on.exit(.rng_restore(.rng_state), add = TRUE)

  ## ---- 0. Basic validation -------------------------------------------------
  required_fields <- c("out", "data", "formula", "model_name")
  missing_fields <- setdiff(required_fields, names(psfm_object))
  if (length(missing_fields) > 0) {
    stop("psfm_object is missing required component(s): ",
      paste(missing_fields, collapse = ", "),
      call. = FALSE
    )
  }

  ## `inefdec` says whether the frontier is a production frontier (u enters
  ## with a minus) or a cost frontier. It has no default here, and until now
  ## omitting it failed inside parallel::clusterExport() -- the promise was
  ## forced during serialization, so the message arrived from
  ## postNode/sendData/serialize with nothing pointing at the user's call.
  ##
  ## Requiring it was the wrong contract anyway. The fit already knows: a user
  ## who restates it can restate it WRONG, and bootstrapping a cost frontier
  ## against a production fit is a silent error rather than a loud one. So it
  ## is recovered from psfm_object$call and only asked for when that fails.
  if (missing(inefdec)) {
    .cl <- psfm_object$call
    .arg <- if (!is.null(.cl)) as.list(.cl)[["inefdec"]] else NULL
    inefdec <- if (is.null(.arg)) {
      formals(psfm)$inefdec ## the fit did not name it, so it took psfm()'s default
    } else {
      tryCatch(eval(.arg, envir = parent.frame()), error = function(e)
        stop("psfm_bootstrap(): `inefdec` was not supplied, and the value used ",
          "for the original fit could not be recovered from psfm_object$call ",
          "(it reads `", paste(deparse(.arg), collapse = ""), "`, which no ",
          "longer evaluates here). Pass `inefdec` explicitly -- it must match ",
          "the fit, or the bootstrap will resample from a differently oriented ",
          "frontier.", call. = FALSE))
    }
    if (!isTRUE(inefdec) && !isFALSE(inefdec)) {
      stop("psfm_bootstrap(): recovered `inefdec` is not TRUE or FALSE. ",
        "Pass it explicitly.", call. = FALSE)
    }
  }

  for (.nm in c("numCores", "BOOT", "individual")) {
    if (eval(call("missing", as.name(.nm)))) {
      stop("psfm_bootstrap(): `", .nm, "` has no default and was not supplied.",
        call. = FALSE)
    }
  }

  if (!requireNamespace("Formula", quietly = TRUE)) {
    stop("Package 'Formula' is required to parse the multi-part model formula.", call. = FALSE)
  }
  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop("Package 'parallel' is required to run the bootstrap in parallel.", call. = FALSE)
  }

  model_name <- psfm_object$model_name
  supported_models <- c("GTRE_Z", "TRE_Z", "GTRE", "GTRE_FML", "TRE", "TFE", "TFE_WMLE", "FD")
  if (!(model_name %in% supported_models)) {
    stop("psfm_bootstrap() does not support model_name = '", model_name, "'. ",
      "Supported: ", paste(supported_models, collapse = ", "), ". ",
      "(GTRE_SEQ1/GTRE_SEQ2/SSFE are moment-based, not MLE; PL80/BC92 ",
      "don't expose the $U/$H structure this function relies on.)",
      call. = FALSE
    )
  }

  data <- psfm_object$data
  out <- psfm_object$out
  n_par <- nrow(out)

  if (!(individual %in% names(data))) {
    stop("Column '", individual, "' (the `individual` argument) was not found in psfm_object$data.", call. = FALSE)
  }

  ## ---- 1. Panel bookkeeping (shared across every family) --------------------
  ids <- data[[individual]]
  uniq_ids <- unique(ids)
  n_id <- length(uniq_ids)
  timez <- as.integer(table(factor(ids, levels = uniq_ids))) ## preserves uniq_ids order
  n_obs <- nrow(data)

  form <- Formula::as.Formula(psfm_object$formula)
  n_rhs <- length(form)[2]

  ## Response variable name.
  y_name <- all.vars(formula(form, lhs = 1, rhs = 0))[1]

  ## H is only returned by GTRE/GTRE_Z -- see header note.
  H_available <- model_name %in% c("GTRE", "GTRE_FML", "GTRE_Z")
  if (H_available) {
    n_h <- length(psfm_object$H)
    if (n_h != n_id) {
      stop("length(psfm_object$H) (", n_h, ") does not match the number of unique ",
        "individuals implied by the `individual` column (", n_id, "). ",
        "Check that `individual` is correct and that $H is one value per individual.",
        call. = FALSE
      )
    }
  } else {
    n_h <- 0L
  }

  ## Per-observation efficiency score field name differs by family: the
  ## randeff family (GTRE_Z/TRE_Z/GTRE/TRE) returns it as $U.
  U_field <- if (model_name %in% c("TFE", "TFE_WMLE", "FD")) "exp_u_hat" else "U"

  ## Whether to run the psoptim (particle-swarm) stage during each bootstrap
  ## refit.
  PSopt_use <- if (model_name %in% c("TFE", "TFE_WMLE")) FALSE else TRUE

  ## Row index of this family's "total scale" parameter in $out.
  scale_row <- switch(model_name,
    "GTRE_Z" = ,
    "TRE_Z" = 1, ## sigv
    "GTRE" = ,
    "TRE" = 2, ## sigma (total SD, not lambda itself)
    ## GTRE_FML reports raw scales in a different order, and its transient
    ## noise scale sits after the frontier block rather than at a fixed row --
    ## resolved by name below, since the row index depends on how many
    ## regressors the formula has.
    "GTRE_FML" = NA_integer_,
    "TFE" = ,
    "TFE_WMLE" = 2, ## sig
    "FD" = c(1, 2) ## sig_u2, sig_v2
  )
  degenerate_scale_floor <- 1e-6 ## well above .Machine$double.eps (~2.2e-16),
  ## far below any plausible fitted scale

  data_x <- model.matrix(form, data = data, rhs = 1)
  ## TFE and FD drop the intercept from their frontier (x) block entirely.
  if (model_name %in% c("TFE", "TFE_WMLE", "FD") &&
    attr(terms(formula(form, rhs = 1)), "intercept") == 1) {
    data_x <- data_x[, -1, drop = FALSE]
  }
  beta_x_hat <- NULL ## set inside each family block below

  ## ---- 2.

  if (model_name %in% c("GTRE_Z", "TRE_Z", "GTRE", "GTRE_FML", "TRE")) {
    ## ---- 2a. "randeff" family --------------------------------------------
    if (model_name %in% c("GTRE_Z", "TRE_Z")) {
      ## Original (unchanged) extraction: z-covariate-driven u, optional
      ## z-covariate-driven h (GTRE_Z only).
      if (h_type == "auto") {
        h_type <- if (model_name == "GTRE_Z") "parametric" else "none"
      }
      if (h_type == "parametric" && n_rhs < 3) {
        stop("h_type = 'parametric' requires a 3-part formula (y ~ x | z | h), ",
          "but the model formula only has ", n_rhs, " RHS part(s).",
          call. = FALSE
        )
      }

      data_z <- model.matrix(form, data = data, rhs = 2)
      data_h <- if (n_rhs >= 3) model.matrix(form, data = data, rhs = 3) else NULL

      Kx <- ncol(data_x)
      Kz <- ncol(data_z)
      Kh <- if (h_type == "parametric") ncol(data_h) else 0L

      expected_n_par <- 2 + Kx + Kz + Kh
      if (h_type == "scalar") expected_n_par <- expected_n_par + 1L

      if (n_par != expected_n_par) {
        stop("Row count of psfm_object$out (", n_par, ") does not match the expected layout ",
          "(2 + Kx + Kz", if (h_type != "none") " + Kh" else "", " = ", expected_n_par, "). ",
          "Check that `h_type` matches how this model was actually specified.",
          call. = FALSE
        )
      }

      sigv_row <- 1
      sigr_row <- 2
      x_rows <- (2 + 1):(2 + Kx)
      z_rows <- (2 + Kx + 1):(2 + Kx + Kz)
      h_rows <- if (h_type == "parametric") {
        (2 + Kx + Kz + 1):(2 + Kx + Kz + Kh)
      } else if (h_type == "scalar") {
        2 + Kx + Kz + 1
      } else {
        integer(0)
      }

      beta_x_hat <- out[x_rows, 1]
      beta_z_hat <- out[z_rows, 1]
      beta_h_hat <- if (h_type != "none") out[h_rows, 1] else NULL
      sigv_hat <- out[sigv_row, 1]
      sigr_hat <- out[sigr_row, 1]
    } else if (model_name == "GTRE_FML") {
      ## Same four-component data-generating process as GTRE -- the estimator
      ## differs, the model does not -- but the parameter vector is the four
      ## RAW scales in the order (sigr, sigv, sigh, sigu), after the frontier
      ## block, rather than GTRE's lambda/sigma reparameterization.
      h_type <- "scalar"
      Kx <- ncol(data_x)
      expected_n_par <- Kx + 4L
      if (n_par != expected_n_par) {
        stop("Row count of psfm_object$out (", n_par, ") does not match the expected ",
          "layout for model_name = 'GTRE_FML' (", expected_n_par, "). ",
          "This usually means psfm_object was not actually fit with this model_name.",
          call. = FALSE
        )
      }
      x_rows   <- seq_len(Kx)
      sigr_hat <- out[Kx + 1L, 1]
      sigv_hat <- out[Kx + 2L, 1]
      sig_h_hat <- out[Kx + 3L, 1]
      sig_u_hat <- out[Kx + 4L, 1]
      scale_row <- Kx + 2L ## sigv, the transient noise scale

      ## Identical warning to the GTRE arm, and it matters more here: this is
      ## the estimator psfm(model_name = "GTRE") reaches by default.
      if (isTRUE(.panel_scale_at_bound(sig_h_hat, sqrt(sig_u_hat^2 + sigv_hat^2))) ||
        isTRUE(.panel_scale_at_bound(sigr_hat, sqrt(sig_u_hat^2 + sigv_hat^2)))) {
        .warn_boot_boundary(sig_h_hat, sigr_hat)
      }

      beta_x_hat <- out[x_rows, 1]
      data_z     <- matrix(1, nrow = n_obs, ncol = 1)
      beta_z_hat <- 2 * log(sig_u_hat)
      beta_h_hat <- 2 * log(sig_h_hat)
      data_h     <- NULL
    } else {
      ## Bare GTRE/TRE: homoskedastic u (and h, for GTRE), no z pipe at all.
      h_type <- if (model_name == "GTRE") "scalar" else "none"

      Kx <- ncol(data_x)
      lambda_hat <- out[1, 1]
      sigma_hat <- out[2, 1]
      sigr_hat <- out[3, 1]
      sig_u_hat <- (lambda_hat * sigma_hat) / sqrt(1 + lambda_hat^2)
      sigv_hat <- sig_u_hat / lambda_hat

      if (model_name == "GTRE") {
        sig_h_hat <- out[4, 1]
        ## A parametric bootstrap resamples from the FITTED model, so a fit
        ## whose persistent scale collapsed to zero generates data with no
        ## persistent component at all -- internally consistent, and useless.
        ## Worse, the bootstrap is inconsistent on the boundary of the
        ## parameter space: the resampled distribution piles up at the
        ## boundary and the resulting interval understates the uncertainty
        ## rather than reflecting it.
        ##
        ## This is not an edge case. One of GTRE's two persistent scales
        ## collapses in roughly a third of samples -- see
        ## notes/sigh_investigation_2026-08-29.md -- so this warning is
        ## expected to fire in ordinary use, and is the point.
        .ref_scale <- sigma_hat
        if (isTRUE(.panel_scale_at_bound(sig_h_hat, .ref_scale)) ||
          isTRUE(.panel_scale_at_bound(sigr_hat, .ref_scale))) {
          .warn_boot_boundary(sig_h_hat, sigr_hat)
        }
        x_rows <- 5:(4 + Kx)
        expected_n_par <- 4 + Kx
      } else {
        x_rows <- 4:(3 + Kx)
        expected_n_par <- 3 + Kx
      }
      if (n_par != expected_n_par) {
        stop("Row count of psfm_object$out (", n_par, ") does not match the expected ",
          "layout for model_name = '", model_name, "' (", expected_n_par, "). ",
          "This usually means psfm_object was not actually fit with this model_name.",
          call. = FALSE
        )
      }

      beta_x_hat <- out[x_rows, 1]
      data_z <- matrix(1, nrow = n_obs, ncol = 1)
      beta_z_hat <- 2 * log(sig_u_hat)
      beta_h_hat <- if (model_name == "GTRE") 2 * log(sig_h_hat) else NULL
      data_h <- NULL
    }

    simulate_dgp <- function(b) {
      v <- rnorm(n_obs, 0, sigv_hat)

      sigma_u <- sqrt(exp(as.vector(data_z %*% beta_z_hat)))
      u <- abs(rnorm(n_obs, 0, sigma_u))

      r_i <- rnorm(n_id, 0, sigr_hat)
      r <- rep(r_i, times = timez)

      h <- switch(h_type,
        "none" = 0,
        "scalar" = {
          sigma_h <- sqrt(exp(unname(beta_h_hat)))
          h_i <- abs(rnorm(n_id, 0, sigma_h))
          rep(h_i, times = timez)
        },
        "parametric" = {
          first_idx <- match(uniq_ids, ids)
          sigma_h_i <- sqrt(exp(as.vector(data_h[first_idx, , drop = FALSE] %*% beta_h_hat)))
          h_i <- abs(rnorm(n_id, 0, sigma_h_i))
          rep(h_i, times = timez)
        }
      )

      if (inefdec == FALSE) {
        as.vector(data_x %*% beta_x_hat) + v + u + r + h
      } else {
        as.vector(data_x %*% beta_x_hat) + v - u + r - h
      }
    }
  } else if (model_name %in% c("TFE", "TFE_WMLE")) {
    ## ---- 2b.
    if (is.null(psfm_object$r_hat_m)) {
      stop("psfm_object$r_hat_m not found -- required to bootstrap a TFE fit ",
        "(the individual fixed effects are held fixed at their original ",
        "point estimates, not redrawn; see this function's header comment).",
        call. = FALSE
      )
    }
    r_hat_m_orig <- psfm_object$r_hat_m ## one value per individual, uniq_ids order
    if (length(r_hat_m_orig) != n_id) {
      stop("length(psfm_object$r_hat_m) (", length(r_hat_m_orig), ") does not match ",
        "the number of unique individuals (", n_id, ").",
        call. = FALSE
      )
    }
    r_fixed <- rep(r_hat_m_orig, times = timez)

    ## rownames(out)[1] is "gamma" instead of "lambda" if the original fit
    ## used psfm(..., gamma = TRUE).
    gamma_used <- identical(rownames(out)[1], "gamma")

    Kx <- ncol(data_x)
    lambda_hat <- out[1, 1] ## "gamma" value directly if gamma_used -- see below
    sigma_hat <- out[2, 1]
    x_rows <- 3:(2 + Kx)
    expected_n_par <- 2 + Kx
    if (n_par != expected_n_par) {
      stop("Row count of psfm_object$out (", n_par, ") does not match the expected ",
        "TFE layout (2 + Kx = ", expected_n_par, ").",
        call. = FALSE
      )
    }
    beta_x_hat <- out[x_rows, 1]

    if (gamma_used) {
      ## gamma = sigma_u^2 / sigma^2 (see start.tfe()'s start_v[1] assignment
      ## in starting.values.R) -- invert to sigma_u.
      sig_u_hat <- sqrt(lambda_hat) * sigma_hat
      sig_v_hat <- sigma_hat * sqrt(1 - lambda_hat)
    } else {
      sig_u_hat <- (lambda_hat * sigma_hat) / sqrt(1 + lambda_hat^2)
      sig_v_hat <- sig_u_hat / lambda_hat
    }

    simulate_dgp <- function(b) {
      v <- rnorm(n_obs, 0, sig_v_hat)
      u <- abs(rnorm(n_obs, 0, sig_u_hat))
      if (inefdec == FALSE) {
        as.vector(data_x %*% beta_x_hat) + v + u + r_fixed
      } else {
        as.vector(data_x %*% beta_x_hat) + v - u + r_fixed
      }
    }
  } else if (model_name == "FD") {
    ## ---- 2c.
    if (!requireNamespace("truncnorm", quietly = TRUE)) {
      stop("Package 'truncnorm' is required to bootstrap an FD fit ",
        "(u_i is drawn from a truncated normal).",
        call. = FALSE
      )
    }

    data_z <- model.matrix(form, data = data, rhs = 2)
    ## FD's z-block ALSO drops its intercept.
    if (attr(terms(formula(form, rhs = 2)), "intercept") == 1) {
      data_z <- data_z[, -1, drop = FALSE]
    }
    Kx <- ncol(data_x)
    Kz <- ncol(data_z)

    expected_n_par <- 3 + Kx + Kz
    if (n_par != expected_n_par) {
      stop("Row count of psfm_object$out (", n_par, ") does not match the expected ",
        "FD layout (3 + Kx + Kz = ", expected_n_par, ").",
        call. = FALSE
      )
    }

    sig_u2_hat <- out[1, 1]
    sig_v2_hat <- out[2, 1]
    mu_hat <- out[3, 1]
    x_rows <- 4:(3 + Kx)
    z_rows <- (4 + Kx):(3 + Kx + Kz)
    beta_x_hat <- out[x_rows, 1]
    delta_hat <- out[z_rows, 1]

    h_it <- as.vector(exp(data_z %*% delta_hat)) ## deterministic given the data

    simulate_dgp <- function(b) {
      v <- rnorm(n_obs, 0, sqrt(sig_v2_hat))
      u_i <- truncnorm::rtruncnorm(n_id, a = 0, mean = mu_hat, sd = sqrt(sig_u2_hat))
      u <- rep(u_i, times = timez) * h_it

      if (inefdec == FALSE) {
        as.vector(data_x %*% beta_x_hat) + v + u
      } else {
        as.vector(data_x %*% beta_x_hat) + v - u
      }
    }
  }

  ## ---- 3. Set up output containers ------------------------------------------
  boot_par <- matrix(0, nrow = BOOT, ncol = n_par + 2)
  colnames(boot_par) <- c(rownames(out), "loglik", "hours")
  boot_eff <- matrix(0, nrow = BOOT, ncol = n_obs)
  rownames(boot_par) <- rownames(boot_eff) <- seq_len(BOOT)
  if (H_available) {
    boot_eff_h <- matrix(0, nrow = BOOT, ncol = n_h)
    rownames(boot_eff_h) <- seq_len(BOOT)
    colnames(boot_eff_h) <- uniq_ids
  } else {
    boot_eff_h <- NULL
  }

  ## ---- 4. The per-replication worker function (generic across families) ----
  boot_one <- function(b) {
    set.seed(b + seed_offset)

    data_b <- data
    data_b[[y_name]] <- simulate_dgp(b)

    MOD <- tryCatch(
      sfa::psfm(
        formula       = form,
        model_name    = model_name,
        data          = data_b,
        maxit.psoptim = maxit.psoptim,
        maxit.bobyqa  = maxit.bobyqa,
        maxit.optim   = maxit.optim,
        individual    = individual,
        inefdec       = inefdec,
        PSopt         = PSopt_use,
        optHessian    = FALSE,
        rand.gtre     = rand.gtre,
        rand.psoptim  = rand.psoptim,
        start_val     = out[, 1]
      ),
      error = function(e) e
    )

    if (inherits(MOD, "error")) {
      return(list(sim = b, ok = FALSE, msg = conditionMessage(MOD)))
    }

    MOD_U <- MOD[[U_field]]
    if (length(MOD_U) != n_obs) {
      return(list(
        sim = b, ok = FALSE,
        msg = paste0(
          "length(MOD$", U_field, ") = ", length(MOD_U),
          " does not match n_obs = ", n_obs
        )
      ))
    }
    if (H_available && length(MOD$H) != n_h) {
      return(list(
        sim = b, ok = FALSE,
        msg = paste0(
          "length(MOD$H) = ", length(MOD$H),
          " does not match n_h = ", n_h
        )
      ))
    }

    if (any(MOD$out[scale_row, 1] < degenerate_scale_floor)) {
      return(list(
        sim = b, ok = FALSE,
        msg = paste0(
          "refit landed on a degenerate variance-boundary mode ",
          "(", paste(rownames(MOD$out)[scale_row], collapse = "/"),
          " < ", degenerate_scale_floor, "); excluded rather than ",
          "letting it corrupt the bootstrap SE -- see this function's ",
          "scale_row/degenerate_scale_floor comment"
        )
      ))
    }

    ## MOD$opt is the raw object from whichever optimizer stage psfm()'s own
    ## internal fallback logic ended up using (see e.g.
    obj_val <- if (!is.null(MOD$opt$value)) MOD$opt$value else MOD$opt$fval

    list(
      sim      = b,
      ok       = TRUE,
      eff      = MOD_U,
      eff_h    = if (H_available) MOD$H else NULL,
      par_row  = c(MOD$out[, 1], obj_val, as.numeric(MOD$total_time, units = "hours"))
    )
  }

  ## ---- 5. Run in parallel ----------------------------------------------------
  cl <- parallel::makeCluster(numCores)
  on.exit(parallel::stopCluster(cl), add = TRUE)

  ## PSOCK workers start as fresh R sessions with the default library path,
  ## not the parent's.
  parallel::clusterCall(cl, function(paths) .libPaths(paths), .libPaths())

  parallel::clusterExport(
    cl,
    varlist = c(
      "data", "form", "model_name", "out", "simulate_dgp",
      "H_available", "U_field", "scale_row", "degenerate_scale_floor", "PSopt_use",
      "timez", "uniq_ids", "ids", "n_id", "n_obs", "n_h",
      "y_name", "individual", "maxit.psoptim", "seed_offset",
      "maxit.bobyqa", "maxit.optim", "inefdec", "rand.gtre", "rand.psoptim"
    ),
    envir = environment()
  )
  for (p in pkgs) {
    parallel::clusterCall(cl, function(pkg) library(pkg, character.only = TRUE), p)
  }

  if (requireNamespace("pbapply", quietly = TRUE)) {
    results <- pbapply::pblapply(X = seq_len(BOOT), FUN = boot_one, cl = cl)
  } else {
    message(
      "Package 'pbapply' not installed -- running without a progress bar. ",
      "Install it (install.packages(\"pbapply\")) to see live progress."
    )
    results <- parallel::parLapply(cl = cl, X = seq_len(BOOT), fun = boot_one)
  }

  ## ---- 6. Assemble results, flag failures -----------------------------------
  failures <- vapply(results, function(r) !isTRUE(r$ok), logical(1))

  for (b in seq_len(BOOT)) {
    if (!failures[b]) {
      boot_par[b, ] <- results[[b]]$par_row
      boot_eff[b, ] <- results[[b]]$eff
      if (H_available) boot_eff_h[b, ] <- results[[b]]$eff_h
    } else {
      boot_par[b, ] <- NA
      boot_eff[b, ] <- NA
      if (H_available) boot_eff_h[b, ] <- NA
    }
  }

  failed_idx <- which(failures)
  if (length(failed_idx) > 0) {
    failed_msgs <- vapply(results[failed_idx], function(r) r$msg, character(1))
    warning(length(failed_idx), " of ", BOOT, " bootstrap replications failed to ",
      "re-estimate and were set to NA.\n",
      paste0("  [b=", failed_idx, "]: ", failed_msgs, collapse = "\n"),
      call. = FALSE
    )
  }

  ## ---- 7. Bootstrap SEs / t-values for each parameter in $out ---------------
  boot_se <- apply(boot_par[, seq_len(n_par), drop = FALSE], 2, sd, na.rm = TRUE)
  boot_tval <- out[, 1] / boot_se
  names(boot_se) <- names(boot_tval) <- rownames(out)

  result <- list(
    boot_par   = boot_par,
    boot_eff   = boot_eff,
    boot_eff_h = boot_eff_h,
    se         = boot_se,
    tval       = boot_tval,
    failures   = if (length(failed_idx) > 0) failed_idx else NULL
  )

  if (write_back) {
    model_out <- psfm_object
    model_out$out[, 2] <- boot_se
    model_out$out[, 3] <- boot_tval
    result$model <- model_out
  }

  result
}
