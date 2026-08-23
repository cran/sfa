## ------------------------------------------------------------------
## Panel data generator.
##
## See DATA_GENERATION_REFERENCE.md (project root) for the full mapping of
## every column here to the psfm() model_name it's meant to test, the true
## parameter vector each model should recover, and an explicit list of
## which models are NOT yet covered (or only best-effort covered) by any
## column here. That file is the canonical reference; comments below are a
## shorter pointer to it, not a duplicate of it.
##
## `eta`: decay-rate parameter used only by the BC92 (Battese-Coelli 1992)
## column below. Defaults to 0.1.
## ------------------------------------------------------------------
data_gen_p <-
  function(t, N, rand, sig_u, sig_v, sig_r, sig_h, cons, tau = 0.5, mu = 0, beta1, beta2, eta = 0.1,
           b_k90 = 0.05, c_k90 = 0.01, d_k90 = 0.05, e_k90 = -0.005,
           rho_mvtn = 0.5) {
    if (!is.null(rand)) {
      .rng_state <- .rng_snapshot()
      on.exit(.rng_restore(.rng_state), add = TRUE)
      set.seed(rand)
    }

    n <- N * t ## Observations
    x1 <- log(runif(n, 0, 10)) ## x1 variable
    x2 <- log(runif(n, 1, 50)) ## x2 variable
    r <- rep(rnorm(N, 0, sig_r), each = t) ## r_i individual effects
    u <- abs(rnorm(n, 0, sig_u)) ## u_it half-normal error term
    v <- rnorm(n, 0, sig_v) ## v_it error term
    h <- abs(rep(rnorm(N, 0, sig_h), each = t)) ## h_i individual effects
    name <- rep(1:N, each = t) ## name for individuals
    year <- rep(seq(1, t, 1), N) ## general time frame
    cons <- cons ## constant
    tau <- tau ## dependence paramter for tfe
    x1_w <- tau * r + sqrt(1 - tau^2) * x1 ## x1 variable with dependence
    x2_w <- tau * r + sqrt(1 - tau^2) * x2 ## x2 variable with dependence


    ## make a first-difference data set from WangHo2010
    r_fd <- runif(N, min = 0, max = 1)
    x_fd <- rep(0, n)

    for (i in seq_len(N)) {
      B <- t * i
      A <- B - (t - 1)
      x_fd[A:B] <- rnorm(t, mean = r_fd[i], sd = 1)
    }

    z_fd <- rnorm(n, 0, 1)
    ## NOTE: u_fd_star is a folded normal (abs() of a general N(mu,sig_u)
    ## draw), not a half-normal (mu=0 case) or a truncated normal
    ## (rtruncnorm, used elsewhere in this package e.g. data_gen_cs.R's
    ## u_tn). When mu != 0 these are three mathematically different
    ## distributions -- verify against Wang & Ho (2010) before treating
    ## FD's Monte Carlo bias/RMSE numbers as conclusive if mu != 0.
    u_fd_star <- abs(rep(rnorm(N, mean = mu, sd = sig_u), each = t))
    r_fd <- rep(r_fd, each = t)
    u_fd <- exp(mu * z_fd) * u_fd_star

    ## Output -  psfm
    y_gtre <- r - h + cons + beta1 * x1 + beta2 * x2 + v - u
    y_tre <- r + cons + beta1 * x1 + beta2 * x2 + v - u
    y_pcs <- cons + beta1 * x1 + beta2 * x2 + v - u
    y_tfe <- r + beta1 * x1_w + beta2 * x2_w + v - u
    y_fd <- r_fd + beta1 * x_fd + v - u_fd

    y_gtre_nc <- r - h + beta1 * x1 + beta2 * x2 + v - u
    y_tre_nc <- r + beta1 * x1 + beta2 * x2 + v - u
    y_pcs_nc <- beta1 * x1 + beta2 * x2 + v - u

    c_gtre <- r + h + cons + beta1 * x1 + beta2 * x2 + v + u
    c_tre <- r + cons + beta1 * x1 + beta2 * x2 + v + u
    c_pcs <- cons + beta1 * x1 + beta2 * x2 + v + u
    c_tfe <- r + beta1 * x1_w + beta2 * x2_w + v + u

    c_gtre_nc <- r + h + beta1 * x1 + beta2 * x2 + v + u
    c_tre_nc <- r + beta1 * x1 + beta2 * x2 + v + u
    c_pcs_nc <- beta1 * x1 + beta2 * x2 + v + u

    ## GTRE-Z on u and with u and h
    u_gtre <- rep(0, n)
    h_gtre <- rep(0, N)
    z_gtre <- runif(n, 1, 2)
    zp_gtreN <- runif(N, 0.5, 1.5)
    zp_gtre <- rep(zp_gtreN, each = t)

    ## u_it and h_i half-normal error term with z var and cons
    u_gtre <- abs(rnorm(n, 0, sqrt(exp(0.33 + 0.61 * z_gtre))))
    h_gtre <- abs(rnorm(N, 0, sqrt(exp(0.25 + 0.21 * zp_gtreN))))

    h_z <- rep(h_gtre, each = t)

    y_gtre_z <- r - h + cons + beta1 * x1 + beta2 * x2 + v - u_gtre
    y_gtre_zz <- r - h_z + cons + beta1 * x1 + beta2 * x2 + v - u_gtre
    y_tre_z <- r + cons + beta1 * x1 + beta2 * x2 + v - u_gtre

    ## ------------------------------------------------------------------
    ## Time-invariant firm effect.
    ##
    ## SSFE (Schmidt & Sickles 1984 LSDV) and PL80 (Pitt & Lee 1980) both
    ## assume inefficiency is CONSTANT within a firm across time (estimated
    ## two different ways -- fixed-effects/LSDV for SSFE, ML random effects
    ## for PL80 -- but the same underlying data-generating assumption).
    ## ------------------------------------------------------------------
    u_inv <- abs(rep(rnorm(N, 0, sig_u), each = t)) ## time-invariant firm-level inefficiency
    y_ssfe <- cons + beta1 * x1 + beta2 * x2 + v - u_inv

    ## ------------------------------------------------------------------
    ## Time-varying decay: BC92 (Battese & Coelli 1992) lets time-invariant
    ## inefficiency decay (or grow) over time, u_it = u_i * B_it, with u_it
    ## = u_i at the last observed period (B_it=1 there). Matches psfm()'s
    ## "BC92" internal convention, B_it = exp(-eta*(t - Tref)), Tref = the
    ## last time period in the whole panel (not per-firm); t below is the
    ## function's own `t` argument (i.e. Tref), and `year` is the
    ## per-observation time index.
    ## ------------------------------------------------------------------
    u_bc92 <- exp(eta * (t - year)) * u_inv
    y_bc92 <- cons + beta1 * x1 + beta2 * x2 + v - u_bc92

    ## ------------------------------------------------------------------
    ## Kumbhakar (1990) time-varying inefficiency -- targets
    ## psfm(model_name = "K1990") and "K1990modified".
    ##
    ## Same error-components structure as PL80/BC92: one time-invariant
    ## half-normal draw per firm, scaled by a deterministic path B_it. Only
    ## B_it differs, exactly as .build_Bit() in matrix_utils.R defines it:
    ##
    ##   K1990          B_it = (1 + exp(b*t + c*t^2))^-1
    ##   K1990modified  B_it = 1 + d*(t - T_i) + e*(t - T_i)^2
    ##
    ## `year` here is the within-firm period 1..t and `t` the panel length, so
    ## (year - t) is the (t - T_i) that K1990modified needs. The defaults are
    ## chosen to keep B_it strictly positive and monotone over the whole panel:
    ## at t = 10, K1990 runs 0.485 -> 0.182 and K1990modified 0.145 -> 1.
    ## ------------------------------------------------------------------
    B_k90 <- 1 / (1 + exp(b_k90 * year + c_k90 * year^2))
    u_k90 <- B_k90 * u_inv
    y_k1990 <- cons + beta1 * x1 + beta2 * x2 + v - u_k90

    B_k90m <- 1 + d_k90 * (year - t) + e_k90 * (year - t)^2
    u_k90m <- B_k90m * u_inv
    y_k1990m <- cons + beta1 * x1 + beta2 * x2 + v - u_k90m

    ## Pitt and Lee (1981) Model III: inefficiency is a T-vector drawn from a
    ## multivariate normal TRUNCATED to the negative orthant, equicorrelated
    ## across periods within a firm. Generated LAST so the RNG stream feeding
    ## every column above is untouched -- the same discipline the NR column
    ## needed in data_gen_cs().
    Sig_mvtn <- sig_u^2 * ((1 - rho_mvtn) * diag(t) + rho_mvtn)
    U_mvtn <- tmvtnorm::rtmvnorm(N,
      mean = rep(0, t), sigma = Sig_mvtn,
      lower = rep(-Inf, t), upper = rep(0, t),
      algorithm = "gibbs", burn.in.samples = 100
    )
    u_mvtn <- as.numeric(t(U_mvtn)) ## firm-major, matching the name/year order
    y_pl_mvtn <- cons + beta1 * x1 + beta2 * x2 + v + u_mvtn ## u <= 0 already

    data_trial <- as.data.frame(cbind(
      name, year, cons, x1, x1_w, x2, x2_w, u, v, r, y_gtre, y_tre, y_tfe, h, y_gtre_nc, y_tre_nc, y_pcs, y_pcs_nc, c_gtre, c_tre, c_tfe, c_gtre_nc, c_tre_nc, c_pcs, c_pcs_nc, y_fd, x_fd, u_fd_star, z_fd, r_fd, u_fd, u_gtre, z_gtre, y_gtre_z, y_tre_z, y_gtre_zz, zp_gtre, h_z,
      u_inv, y_ssfe, u_bc92, y_bc92,
      B_k90, u_k90, y_k1990, B_k90m, u_k90m, y_k1990m,
      u_mvtn, y_pl_mvtn
    ))
    colnames(data_trial) <- c(
      "name", "year", "cons", "x1", "x1_w", "x2", "x2_w", "u", "v", "r", "y_gtre", "y_tre", "y_tfe", "h", "y_gtre_nc", "y_tre_nc", "y_pcs", "y_pcs_nc", "c_gtre", "c_tre", "c_tfe", "c_gtre_nc", "c_tre_nc", "c_pcs", "c_pcs_nc", "y_fd", "x_fd", "u_fd_star", "z_fd", "r_fd", "u_fd", "u_gtre", "z_gtre", "y_gtre_z", "y_tre_z", "y_gtre_zz", "zp_gtre", "h_z",
      "u_inv", "y_ssfe", "u_bc92", "y_bc92",
      "B_k90", "u_k90", "y_k1990", "B_k90m", "u_k90m", "y_k1990m",
      "u_mvtn", "y_pl_mvtn"
    )
    data_rand <- pdata.frame(data_trial, c("name", "year"))

    return(data_rand)
  }
