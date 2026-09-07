## Cross-sectional data generator.
data_gen_cs <- function(N, rand, sig_u, sig_v, cons, beta1, beta2, a, mu, sig_w = sig_u,
                        shape_g = 2, m_nak = 1, mu_ln = -0.5, k_w = 1.5,
                        lam_tsl = 1.5) {
  if (!is.null(rand)) {
    .rng_state <- .rng_snapshot()
    on.exit(.rng_restore(.rng_state), add = TRUE)
    set.seed(rand)
  }

  n <- N ## Observations
  x1 <- log(runif(n, 2, 10)) ## x1 variable
  x2 <- log(runif(n, 3, 50)) ## x2 variable

  ## Normal half normal -- targets sfm(model_name="NHN").
  u <- abs(rnorm(n, 0, sig_u)) ## u_it half-normal error term
  v <- rnorm(n, 0, sig_v)

  ## Normal truncated normal -- targets sfm(model_name="NTN").
  u_tn <- rtruncnorm(n = n, mean = mu, sd = sig_u, a = 0)

  ## NHN with Z -- targets sfm(model_name="NHN_Z", formula = y_pcs_z~x1+x2|z).
  z <- runif(n, 1, 2)
  uz <- abs(rnorm(n, 0, exp(0.9 + 0.6 * z))) ## u_it half-normal error term with z var and cons

  ## Normal exponential -- targets sfm(model_name="NE").
  phi <- 1 / sig_u
  u_e <- rexp(n, rate = phi) ## u_it exponential

  ## Normal exponential with Z -- targets sfm(model_name="NE_Z", formula =
  ## y_pcs_ez~x1+x2|z).
  uz_e <- rexp(n, rate = 1 / exp(0.9 + 0.6 * z))

  ## T half T -- targets sfm(model_name="THT").
  u_t <- sig_u * abs(rt(n, a)) ## u_it half-T error term
  v_t <- sig_v * rt(n, a)

  ## Cauchy half Cauchy -- NOT currently matched to any sfm()/psfm() model (no
  ## "half-Cauchy" model_name exists).
  u_c <- abs(rcauchy(n, 0, sig_u)) ## u_it half-cauchy error term
  v_c <- rcauchy(n, 0, sig_v)

  ## Normal uniform -- targets sfm(model_name="NU") (Li, 1996).
  u_u <- runif(n, min = 0, max = sig_u) ## u_it uniform error term

  name <- 1:N ## name for individuals (cross-sectional: one obs each)
  cons <- cons ## constant

  ## Two-tier components (ttsfm()).

  ## Homoskedastic two-tier, exponential -- targets ttsfm(model_name="TTNE")
  ## and ttsfm(model_name="TTNLS"), formula = y_ttne~x1+x2 (no pipes).
  w_tt <- rexp(n, rate = 1 / sig_w) ## w_it exponential, second one-sided component
  y_ttne <- cons + beta1 * x1 + beta2 * x2 + v + w_tt - u_e

  ## Homoskedastic two-tier, half-normal -- targets
  ## ttsfm(model_name="TTHN"), formula = y_tthn~x1+x2 (no pipes).
  w_tt_hn <- abs(rnorm(n, 0, sig_w)) ## w_it half-normal, second one-sided component
  y_tthn <- cons + beta1 * x1 + beta2 * x2 + v + w_tt_hn - u

  ## Heteroskedastic (z/zp) two-tier, half-normal -- targets
  ## ttsfm(model_name="TTHN", formula = y_tthn_z~x1+x2|z|zp).
  zp <- runif(n, 1, 2)
  wz_hn <- abs(rnorm(n, 0, exp(0.7 + 0.5 * zp)))
  y_tthn_z <- cons + beta1 * x1 + beta2 * x2 + v + wz_hn - uz

  ## Zero-inefficiency (zsfm()).

  ## ZISF (no z) -- targets zsfm(model_name="ZISF", formula=y_zisf~x1+x2).
  gamma0 <- 0.3
  prob0 <- exp(-abs(gamma0))
  eff_ind <- rbinom(n, 1, prob0) ## 1 = efficient (u forced to 0)
  y_zisf <- cons + beta1 * x1 + beta2 * x2 + v - u * (1 - eff_ind)

  ## ZISF_Z -- targets zsfm(model_name="ZISF_Z", formula=y_zisf_z~x1+x2|z).
  prob_z_true <- plogis(0.9 + 0.6 * z)
  eff_ind_z <- rbinom(n, 1, prob_z_true)
  y_zisf_z <- cons + beta1 * x1 + beta2 * x2 + v - u * (1 - eff_ind_z)

  ## Output - sfm
  y_pcs <- cons + beta1 * x1 + beta2 * x2 + v - u
  y_pcs_z <- cons + beta1 * x1 + beta2 * x2 + v - uz
  y_pcs_t <- cons + beta1 * x1 + beta2 * x2 + v_t - u_t
  y_pcs_tn <- cons + beta1 * x1 + beta2 * x2 + v - u_tn
  y_pcs_e <- cons + beta1 * x1 + beta2 * x2 + v - u_e
  y_pcs_ez <- cons + beta1 * x1 + beta2 * x2 + v - uz_e
  y_pcs_c <- cons + beta1 * x1 + beta2 * x2 + v_c - u_c
  y_pcs_u <- cons + beta1 * x1 + beta2 * x2 + v - u_u

  ## Normal + Cauchy - half-normal
  y_pcs_w <- cons + beta1 * x1 + beta2 * x2 + v + v_c - u

  ## The five distributions below are drawn HERE, after every pre-existing
  ## random draw, and not next to their siblings above.

  ## Normal gamma -- targets sfm(model_name="NG").
  u_g <- rgamma(n, shape = shape_g, scale = sig_u / shape_g)

  ## Normal Nakagami -- targets sfm(model_name="NNAK").
  u_nak <- sqrt(rgamma(n, shape = m_nak, scale = sig_u^2 / m_nak))

  ## Normal truncated skew-Laplace -- targets sfm(model_name="TSL").
  u_tsl <- {
    out <- numeric(0)
    while (length(out) < n) {
      cand <- rexp(2 * n, rate = 1 / sig_u)
      keep <- runif(length(cand)) < (1 - exp(-lam_tsl * cand / sig_u) / 2)
      out <- c(out, cand[keep])
    }
    out[seq_len(n)]
  }

  ## Normal generalized-exponential -- targets sfm(model_name="NGE").
  u_ge <- -log(1 - sqrt(runif(n))) * sig_u

  ## Normal lognormal -- targets sfm(model_name="NLN").
  u_ln <- rlnorm(n, meanlog = mu_ln, sdlog = sig_u)

  ## Normal Weibull -- targets sfm(model_name="NW").
  u_w <- rweibull(n, shape = k_w, scale = sig_u)

  ## Outcomes for the five distributions added above. Each pairs with exactly one
  ## model_name -- see DATA_GENERATION_REFERENCE.md for the true-value triples.
  y_pcs_g <- cons + beta1 * x1 + beta2 * x2 + v - u_g ## NG
  y_pcs_nak <- cons + beta1 * x1 + beta2 * x2 + v - u_nak ## NNAK
  y_pcs_tsl <- cons + beta1 * x1 + beta2 * x2 + v - u_tsl ## TSL
  y_pcs_ge <- cons + beta1 * x1 + beta2 * x2 + v - u_ge ## NGE
  y_pcs_ln <- cons + beta1 * x1 + beta2 * x2 + v - u_ln ## NLN
  y_pcs_wb <- cons + beta1 * x1 + beta2 * x2 + v - u_w ## NW

  ## Skew-t -- targets sfm(model_name="THT"). Use THIS column, not y_pcs_t.
  lam_st <- rgamma(n, shape = a / 2, rate = a / 2) ## common mixing, eq (5)
  v_st <- rnorm(n, 0, sig_v) / sqrt(lam_st) ## t noise,    scale sig_v
  u_st <- abs(rnorm(n, 0, sig_u)) / sqrt(lam_st) ## half-t ineff, scale sig_u
  y_pcs_st <- cons + beta1 * x1 + beta2 * x2 + v_st - u_st ## THT

  ## Student-t--half-normal -- targets sfm(model_name="tHN").
  v_thn <- sig_v * rt(n, a)
  u_thn <- abs(rnorm(n, 0, sig_u))
  y_pcs_thn <- cons + beta1 * x1 + beta2 * x2 + v_thn - u_thn ## tHN

  ## Normal-Rayleigh -- targets sfm(model_name="NR").
  u_r <- (sig_u / sqrt(2)) * sqrt(rnorm(n)^2 + rnorm(n)^2)
  y_pcs_r <- cons + beta1 * x1 + beta2 * x2 + v - u_r ## NR

  data_trial <- as.data.frame(cbind(
    name, cons, x1, x2, u, uz, v, u_t, v_t, u_c, v_c, u_e, u_u, u_tn,
    y_pcs, y_pcs_t, y_pcs_e, y_pcs_c, y_pcs_u, y_pcs_z, y_pcs_w, y_pcs_tn, z,
    uz_e, y_pcs_ez,
    w_tt, y_ttne, w_tt_hn, y_tthn, zp, wz_hn, y_tthn_z,
    eff_ind, y_zisf, prob_z_true, eff_ind_z, y_zisf_z,
    u_g, y_pcs_g, u_nak, y_pcs_nak, u_tsl, y_pcs_tsl, u_ge, y_pcs_ge, u_ln, y_pcs_ln, u_w, y_pcs_wb,
    lam_st, u_st, v_st, y_pcs_st, v_thn, u_thn, y_pcs_thn, u_r, y_pcs_r
  ))
  colnames(data_trial) <- c(
    "name", "cons", "x1", "x2", "u", "uz", "v", "u_t", "v_t", "u_c", "v_c", "u_e", "u_u", "u_tn",
    "y_pcs", "y_pcs_t", "y_pcs_e", "y_pcs_c", "y_pcs_u", "y_pcs_z", "y_pcs_w", "y_pcs_tn", "z",
    "uz_e", "y_pcs_ez",
    "w_tt", "y_ttne", "w_tt_hn", "y_tthn", "zp", "wz_hn", "y_tthn_z",
    "eff_ind", "y_zisf", "prob_z_true", "eff_ind_z", "y_zisf_z",
    "u_g", "y_pcs_g", "u_nak", "y_pcs_nak", "u_tsl", "y_pcs_tsl", "u_ge", "y_pcs_ge",
    "u_ln", "y_pcs_ln", "u_w", "y_pcs_wb",
    "lam_st", "u_st", "v_st", "y_pcs_st", "v_thn", "u_thn", "y_pcs_thn",
    "u_r", "y_pcs_r"
  )

  data_rand <- data.frame(data_trial)
  data_rand$con <- 1
  return(data_rand)
}
