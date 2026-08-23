## ============================================================================
## psfm_bootstrap()
## Parametric bootstrap for panel stochastic frontier models fit with psfm()
## from the sfa package.
##
## Authors:  David H. Bernstein & Chris Parmeter
##
## Covers all psfm() models except the moment-based/non-MLE ones
## (GTRE_SEQ1, GTRE_SEQ2, SSFE) and the error-components frontier models
## (PL80, BC92, which don't expose the internal $U/$H structure this
## function relies on).
##
## ----------------------------------------------------------------------------
## SUPPORTED MODELS, grouped by "family" (each family shares one
## data-generating process shape and one block of extraction code below):
##
##   "randeff" family: GTRE_Z, TRE_Z, GTRE, TRE
##     y* = X*beta_hat + v + [-]u + r + [-]h   (sign per `inefdec`)
##       v ~ N(0, sigma_v_hat^2)                          -- always
##       u ~ |N(0, sigma_u_i^2)|, sigma_u_i either z-covariate-driven
##           (GTRE_Z/TRE_Z: sigma_u_i = sqrt(exp(z_i'delta_hat))) or a single
##           scalar for every observation (bare GTRE/TRE, which have no z
##           pipe at all -- sigma_u_i = sigma_u_hat for all i)
##       r ~ N(0, sigma_r_hat^2), one draw per individual, repeated over time
##       h (GTRE/GTRE_Z only): same z-covariate-driven-or-scalar mechanism as
##           u, but for the persistent/time-invariant inefficiency component,
##           one draw per individual repeated over time. TRE/TRE_Z have no h.
##     GTRE_Z/TRE_Z's extraction logic (row-slicing $out, sign convention,
##     re-draw mechanics) is UNCHANGED from the original version of this
##     function -- only generalized so the SAME simulate_dgp()/boot_one()
##     code also covers bare GTRE/TRE by constructing an equivalent
##     intercept-only z/h "design" (data_z = a column of 1s, beta_z_hat =
##     2*log(sigma_u_hat) so that sqrt(exp(1*beta_z_hat)) = sigma_u_hat
##     reproduces the SAME formula GTRE_Z's code already uses) rather than by
##     writing new simulate-u/h code.
##
##   "tfe" family: TFE
##     Fixed-effects model (Chen-Schmidt-Wang 2014-style): the individual
##     effects r_i are NOT random draws from a distribution here -- they are
##     nuisance parameters recovered post-hoc from the data via a moment
##     correction (`r_hat_m` in psfm.R's TFE branch, exposed on the fitted
##     object as `$r_hat_m`). The standard parametric-bootstrap treatment for
##     a fixed-effects model is to hold these estimated effects FIXED at
##     their point estimates across every bootstrap replication (redrawing
##     them would contradict treating them as fixed, not random) and only
##     redraw the genuinely stochastic v/u terms:
##       y* = X*beta_hat + v + [-]u + r_hat_m_original   (r_hat_m held fixed)
##       v ~ N(0, sigma_v_hat^2), u ~ |N(0, sigma_u_hat^2)| -- both derived
##       from the fitted lambda_hat/sigma_hat via the same
##       sigma_u = lambda*sigma/sqrt(1+lambda^2), sigma_v = sigma_u/lambda
##       formula psfm.R's own TFE/GTRE/TRE post-estimation TE-recovery code
##       already uses.
##
##   "fd" family: FD (Wang & Ho 2010 first-difference estimator)
##     Derived directly from like.fd()'s own likelihood terms in psfm.R
##     (l5/l6/l7), not from data_gen.p.R's y_fd column -- that column's own
##     DGP was already flagged as uncertain in DATA_GENERATION_REFERENCE.md
##     ("u_fd_star folded-vs-truncated-normal question") and additionally
##     reuses a single `mu` value for two conceptually different roles
##     (truncation location AND the z-scaling exponent), which does not
##     match the actual estimator's two separate parameters (`mu` vs. the
##     z-block's own delta). Reading like.fd()'s l5/l6/l7 terms directly
##     shows they are exactly the classical Battese-Coelli closed-form
##     likelihood for a TIME-INVARIANT TRUNCATED-NORMAL random effect u_i ~
##     N+(mu, sigma_u2) (truncated below at 0) integrated out analytically
##     against a GLS residual, applied to first-differenced data with a
##     deterministic time-varying scale h_it = exp(z_it'delta_hat):
##       y_it* = X_it*beta_hat + v_it -/+ u_i * h_it
##       u_i   ~ TruncatedNormal(mean = mu_hat, sd = sqrt(sig_u2_hat), lower = 0),
##               one draw per individual (time-invariant), via
##               truncnorm::rtruncnorm() (already a declared Import)
##       v_it  ~ N(0, sig_v2_hat), i.i.d. per observation
##       h_it  = exp(z_it' %*% delta_hat), deterministic given data
##     psfm()'s own internal first-differencing handles the rest identically
##     to the original fit -- the response supplied here is in LEVELS, same
##     as psfm_object$data, not pre-differenced.
##
## NOT supported (structurally incompatible with this function's approach):
##   GTRE_SEQ1, GTRE_SEQ2, SSFE  -- moment-based/LSDV, not MLE; no $opt to
##     bootstrap around in the same sense (see sfareg_methods.R's logLik.sfareg
##     for the same NA-with-warning treatment these get elsewhere).
##   PL80, BC92 -- the error-components frontier models; $U/$H don't
##     follow this package's own row-layout/field conventions and PL80/BC92
##     don't expose an analogous per-replication-refittable formula/data
##     contract the way the other models here do.
##
## NOTE ON ROW LAYOUT OF psfm_object$out (per family, confirmed from
## starting.values.R / the model's own branch in psfm.R):
##   GTRE_Z/TRE_Z: sigv, sigr, [x-block], [z-block], [h-block if GTRE_Z]
##   GTRE:         lambda, sigma, sigr, sigh, [x-block]
##   TRE:          lambda, sigma, sigr, [x-block]
##   TFE:          lambda (or "gamma" if fit with gamma=TRUE), sig, [x-block]
##   FD:           sig_u2, sig_v2, mu, [x-block], [z-block]
##
## $H (persistent/individual-level TE score) is only returned by GTRE and
## GTRE_Z -- TRE, TRE_Z, TFE, and FD only return $U (per-observation TE).
## boot_eff_h is therefore only built/populated when the model_name is one
## of those two (this was a latent bug in the original version of this
## function: it unconditionally required length(psfm_object$H) == n_id,
## which would have errored for TRE_Z -- already in the "supported" set --
## since TRE_Z's fitted object has no $H at all).
## ============================================================================

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

  ## Restore the caller's random stream on the way out. boot_one() seeds each
  ## replicate deliberately (set.seed(b + seed_offset)) so the bootstrap is
  ## reproducible, but that must not leak into the user's session -- the
  ## snapshot is taken once here rather than per replicate, so the per-replicate
  ## seeding is untouched and results are unchanged.
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

  if (!requireNamespace("Formula", quietly = TRUE)) {
    stop("Package 'Formula' is required to parse the multi-part model formula.", call. = FALSE)
  }
  if (!requireNamespace("parallel", quietly = TRUE)) {
    stop("Package 'parallel' is required to run the bootstrap in parallel.", call. = FALSE)
  }

  model_name <- psfm_object$model_name
  supported_models <- c("GTRE_Z", "TRE_Z", "GTRE", "TRE", "TFE", "TFE_WMLE", "FD")
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

  ## Response variable name. NOTE: all.vars() is not S3-generic, so calling it
  ## directly on a `Formula` object does not dispatch correctly and silently
  ## returns nothing usable. We instead ask Formula for just the LHS as a
  ## plain base-R formula (lhs = 1, rhs = 0), which all.vars() *does* handle.
  y_name <- all.vars(formula(form, lhs = 1, rhs = 0))[1]

  ## H is only returned by GTRE/GTRE_Z -- see header note.
  H_available <- model_name %in% c("GTRE", "GTRE_Z")
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
  ## randeff family (GTRE_Z/TRE_Z/GTRE/TRE) returns it as $U, but TFE and FD
  ## both return it as $exp_u_hat instead (confirmed from each model's own
  ## names(results) assignment in psfm.R) -- not a naming choice this
  ## function should paper over silently, so it's made explicit here and
  ## used generically in boot_one() below via MOD[[U_field]].
  U_field <- if (model_name %in% c("TFE", "TFE_WMLE", "FD")) "exp_u_hat" else "U"

  ## Whether to run the psoptim (particle-swarm) stage during each bootstrap
  ## refit. The original version of this function hardcoded PSopt = TRUE
  ## unconditionally; psfm()'s own default is PSopt = FALSE (bobyqa -> optim
  ## only, psoptim skipped). Testing TFE's bootstrap (2026-07) found this
  ## mattered a lot: with PSopt = TRUE, TFE repeatedly collapsed onto a
  ## degenerate boundary mode (sig pinned to exactly .Machine$double.eps) --
  ## and critically, this reproduced even refitting the SAME real (non-
  ## bootstrapped) data with just a different psoptim random seed, proving
  ## it's psoptim's own wide/undirected exploration finding and getting
  ## stuck in that region, not a problem with the resampled data or this
  ## function's DGP. With PSopt = FALSE, 5/5 refits of that same data across
  ## different seeds converged to essentially identical, sensible values.
  ## TFE therefore uses PSopt = FALSE here, matching psfm()'s own default and
  ## the empirically safer behavior; the randeff family keeps PSopt = TRUE
  ## since it was already tested clean that way (GTRE_Z/TRE_Z/GTRE/TRE all
  ## converged with 0 failures across repeated testing).
  PSopt_use <- if (model_name %in% c("TFE", "TFE_WMLE")) FALSE else TRUE

  ## Row index of this family's "total scale" parameter in $out -- used below
  ## as a defense-in-depth check for a refit that lands on a degenerate
  ## variance-boundary mode anyway (a variance component pinned near its
  ## numerical floor -- composed-error/panel-SFA likelihoods can have a
  ## genuine, non-error boundary mode at a variance component -> 0).
  ## PSopt_use above removes the main trigger found for TFE, but this check
  ## stays on for every family as a safety net, since an unchecked
  ## degenerate refit would silently corrupt the bootstrap SE (a huge,
  ## meaningless sd() driven by a few degenerate draws) rather than just
  ## misreport one replication. Flagged as a per-replication FAILURE
  ## (excluded, listed in $failures) rather than attempting to fix the
  ## underlying MLE's numerical behavior itself, since a genuine boundary
  ## mode is a modeling-strategy question, not a bug.
  scale_row <- switch(model_name,
    "GTRE_Z" = ,
    "TRE_Z" = 1, ## sigv
    "GTRE" = ,
    "TRE" = 2, ## sigma (total SD, not lambda itself)
    "TFE" = ,
    "TFE_WMLE" = 2, ## sig
    "FD" = c(1, 2) ## sig_u2, sig_v2
  )
  degenerate_scale_floor <- 1e-6 ## well above .Machine$double.eps (~2.2e-16),
  ## far below any plausible fitted scale

  data_x <- model.matrix(form, data = data, rhs = 1)
  ## TFE and FD drop the intercept from their frontier (x) block entirely --
  ## matching data.processing.R's own rule (`x_vars_vec <- if(model_name %in%
  ## c("TFE","FD","SSFE") & intercept==1) colnames(data_x)[-1] else
  ## colnames(data_x)`), since TFE's within-transformation and FD's
  ## first-differencing both remove the constant, leaving no separate
  ## intercept parameter in $out. Applying the same drop here keeps Kx (and
  ## therefore beta_x_hat / the frontier term in simulate_dgp) aligned with
  ## $out's actual row count for these two models.
  if (model_name %in% c("TFE", "TFE_WMLE", "FD") &&
    attr(terms(formula(form, rhs = 1)), "intercept") == 1) {
    data_x <- data_x[, -1, drop = FALSE]
  }
  beta_x_hat <- NULL ## set inside each family block below

  ## ==========================================================================
  ## ---- 2. Family-specific extraction: build `simulate_dgp(b)`, a closure
  ##      that returns a length-n_obs simulated response vector for
  ##      replication b. Everything family-specific lives in this block;
  ##      boot_one() below is generic across all families.
  ## ==========================================================================

  if (model_name %in% c("GTRE_Z", "TRE_Z", "GTRE", "TRE")) {
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
    } else {
      ## Bare GTRE/TRE: homoskedastic u (and h, for GTRE), no z pipe at all.
      ## Build an intercept-only "z design" so the SAME u/h simulation formulas
      ## GTRE_Z's code already uses (sigma = sqrt(exp(design %*% beta)))
      ## reproduce the model's actual homoskedastic sigma_u/sigma_h exactly:
      ## with a single intercept column of 1s, sqrt(exp(1*beta_z_hat)) =
      ## sigma_u_hat iff beta_z_hat = 2*log(sigma_u_hat).
      h_type <- if (model_name == "GTRE") "scalar" else "none"

      Kx <- ncol(data_x)
      lambda_hat <- out[1, 1]
      sigma_hat <- out[2, 1]
      sigr_hat <- out[3, 1]
      sig_u_hat <- (lambda_hat * sigma_hat) / sqrt(1 + lambda_hat^2)
      sigv_hat <- sig_u_hat / lambda_hat

      if (model_name == "GTRE") {
        sig_h_hat <- out[4, 1]
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
    ## ---- 2b. "tfe" family: fixed effects, r held fixed at r_hat_m --------
    ## Both true-fixed-effects estimators share this block: Greene's TFE
    ## ("TFE") and Chen-Schmidt-Wang's within MLE ("TFE_WMLE") assume the
    ## same data-generating process and report the same (lambda|gamma, sig,
    ## beta) row layout plus $r_hat_m/$exp_u_hat, which is all this block
    ## reads.
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
    ## used psfm(..., gamma = TRUE) -- detect and pass the same reparameterization
    ## through to the refit so the likelihood being bootstrapped matches exactly.
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
      ## in starting.values.R) -- invert to sigma_u, sigma_v directly rather
      ## than via lambda (avoids a division by a possibly-tiny sqrt(1-gamma)).
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
    ## ---- 2c. "fd" family: time-invariant truncated-normal u_i, scaled by
    ##      a deterministic time-varying z-covariate factor h_it. Derived
    ##      directly from like.fd()'s l5/l6/l7 likelihood terms -- see this
    ##      function's header comment for the full derivation.
    if (!requireNamespace("truncnorm", quietly = TRUE)) {
      stop("Package 'truncnorm' is required to bootstrap an FD fit ",
        "(u_i is drawn from a truncated normal).",
        call. = FALSE
      )
    }

    data_z <- model.matrix(form, data = data, rhs = 2)
    ## FD's z-block ALSO drops its intercept (data.processing.R:
    ## `z_vars_vec <- if(model_name %in% c("TFE","FD") & intercept_z==1)
    ## colnames(data_z)[-1] else colnames(data_z)`) -- same rationale as the
    ## x-block drop above, applied to rhs=2 specifically for FD (TFE has no
    ## z-block at all in practice, since it's capped at 1 RHS part).
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
    ## internal fallback logic ended up using (see e.g. the TFE branch's
    ## `if(optHessian==FALSE & PSopt==FALSE){opt <- bob1}`). optim()'s result
    ## has $value; minqa::bobyqa()'s result (used when that fallback fires --
    ## which it does here whenever PSopt_use is FALSE, since optHessian is
    ## always FALSE in this function's own refit call above) has $fval
    ## instead, no $value at all. MOD$opt$value would then silently be NULL,
    ## which c() drops rather than errors on -- shortening par_row by one and
    ## corrupting every column after it. Coalesce explicitly instead.
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

  ## PSOCK workers start as fresh R sessions with the default library path, not
  ## the parent's. If sfa is installed anywhere other than a default library --
  ## a project library, renv/packrat, a user-set R_LIBS_USER, or the temporary
  ## <pkg>.Rcheck tree that R CMD check installs into -- the workers' library()
  ## call below fails with "there is no package called 'sfa'" and the whole
  ## bootstrap dies in checkForRemoteErrors(). On Unix the workers usually
  ## inherit R_LIBS from the parent's environment and the problem stays hidden;
  ## on Windows they do not, which is where this first showed up (the vignette's
  ## psfm_bootstrap() chunk failed on windows-latest while both Linux runners
  ## passed). Pushing the parent's search path to the workers makes the call
  ## resolve the same way in the parent and the workers on every platform.
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
