## ------------------------------------------------------------------
## Cross-sectional data generator.
##
## See DATA_GENERATION_REFERENCE.md (project root) for the full mapping
## of every column here to the sfm()/zsfm()/ttsfm() model_name it's meant
## to test, the true parameter vector each model should recover, and an
## explicit list of which models are NOT yet covered by any column here.
## That file is the canonical reference; comments below are a shorter
## pointer to it, not a duplicate of it.
##
## `sig_w`: scale of the second one-sided component used only by the
## two-tier (TTNE/TTHN/TTNLS) columns below. Defaults to sig_u.
## ------------------------------------------------------------------
data_gen_cs <- function(N, rand, sig_u, sig_v, cons, beta1, beta2, a, mu, sig_w = sig_u,
                        shape_g = 2, m_nak = 1, mu_ln = -0.5, k_w = 1.5) {
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
  ## sigma_u given z is exp(0.9 + 0.6*z) directly (NOT sqrt(exp(...)) -- see
  ## the reference doc for why this differs from psfm()'s GTRE_Z/TRE_Z
  ## convention). True z-coefficients to recover: intercept 0.9, slope 0.6.
  z <- runif(n, 1, 2)
  uz <- abs(rnorm(n, 0, exp(0.9 + 0.6 * z))) ## u_it half-normal error term with z var and cons

  ## Normal exponential -- targets sfm(model_name="NE").
  phi <- 1 / sig_u
  u_e <- rexp(n, rate = phi) ## u_it exponential

  ## Normal exponential with Z -- targets sfm(model_name="NE_Z",
  ## formula = y_pcs_ez~x1+x2|z). Same exp(0.9+0.6*z) convention as uz above
  ## (NE_Z's sigma_u_fun in sfm.R uses the identical direct-exp link), just
  ## applied to an exponential rate instead of a half-normal sd.
  uz_e <- rexp(n, rate = 1 / exp(0.9 + 0.6 * z))

  ## T half T -- targets sfm(model_name="THT"). Shares the single "a"
  ## degrees-of-freedom parameter across the t-distributed noise and
  ## half-t inefficiency, matching THT's likelihood in sfm.R.
  u_t <- sig_u * abs(rt(n, a)) ## u_it half-T error term
  v_t <- sig_v * rt(n, a)

  ## Cauchy half Cauchy -- NOT currently matched to any sfm()/psfm() model
  ## (no "half-Cauchy" model_name exists). Kept for robustness/contaminated-
  ## noise sensitivity checks outside the main package's testable surface.
  u_c <- abs(rcauchy(n, 0, sig_u)) ## u_it half-cauchy error term
  v_c <- rcauchy(n, 0, sig_v)

  ## Normal uniform -- targets sfm(model_name="NU") (Li, 1996). u is uniform on
  ## [0, sig_u], so the reported `theta` (the upper support point) IS sig_u.
  ## NOTE this column predates NU, which was added in 1.1.2; it was previously
  ## commented as having no matching model.
  u_u <- runif(n, min = 0, max = sig_u) ## u_it uniform error term

  name <- 1:N ## name for individuals (cross-sectional: one obs each)
  cons <- cons ## constant

  ## ------------------------------------------------------------------
  ## Two-tier components (ttsfm()).
  ##
  ## Matches ttsfm.R's actual likelihoods: TTNE's fn() exponentiates
  ## alpha = e/sigu + sigv^2/(2*sigu^2) and
  ## a = sigv^2/(2*sigw^2) - e/sigw where e = Y - X*beta, which is exactly
  ## the closed-form composed-error density for
  ##   e = v + w - u,   v ~ N(0, sigv^2),  u ~ Exponential(mean = sigu),
  ##                                        w ~ Exponential(mean = sigw)
  ## (Polachek & Yoon 1987-style two-tier). TTHN's fn() (theta1/theta2/s/
  ## omega1/omega2 reparameterization) is the same e = v + w - u structure
  ## with u, w ~ Half-Normal(sigu), Half-Normal(sigw) instead. TTNLS
  ## (nonlinear least squares) uses E[u]=sigu, E[w]=sigw directly -- exactly
  ## the exponential case's first moments -- so it targets the same
  ## homoskedastic exponential columns as TTNE.
  ## ------------------------------------------------------------------

  ## Homoskedastic two-tier, exponential -- targets ttsfm(model_name="TTNE")
  ## and ttsfm(model_name="TTNLS"), formula = y_ttne~x1+x2 (no pipes).
  w_tt <- rexp(n, rate = 1 / sig_w) ## w_it exponential, second one-sided component
  y_ttne <- cons + beta1 * x1 + beta2 * x2 + v + w_tt - u_e

  ## Homoskedastic two-tier, half-normal -- targets
  ## ttsfm(model_name="TTHN"), formula = y_tthn~x1+x2 (no pipes).
  w_tt_hn <- abs(rnorm(n, 0, sig_w)) ## w_it half-normal, second one-sided component
  y_tthn <- cons + beta1 * x1 + beta2 * x2 + v + w_tt_hn - u

  ## Heteroskedastic (z/zp) two-tier, half-normal -- targets
  ## ttsfm(model_name="TTHN", formula = y_tthn_z~x1+x2|z|zp). sigma_u given z
  ## reuses the same exp(0.9+0.6*z) convention as uz/NHN_Z above; sigma_w
  ## given zp is exp(0.7+0.5*zp) (same direct-exp convention, both pipes are
  ## parameterized identically in ttsfm.R's TTHN fn). True coefficients to
  ## recover: sigma_u ~ z intercept 0.9 / slope 0.6, sigma_w ~ zp intercept
  ## 0.7 / slope 0.5.
  zp <- runif(n, 1, 2)
  wz_hn <- abs(rnorm(n, 0, exp(0.7 + 0.5 * zp)))
  y_tthn_z <- cons + beta1 * x1 + beta2 * x2 + v + wz_hn - uz

  ## ------------------------------------------------------------------
  ## Zero-inefficiency (zsfm()).
  ##
  ## ZISF's likelihood in zsfm.R is a two-component mixture,
  ## prob*N(0,sigv^2) + (1-prob)*[normal/half-normal composed-error density],
  ## i.e. each observation is drawn from an "efficient firm" regime (u=0,
  ## probability prob) or an "inefficient firm" regime (u ~ Half-Normal(sigu),
  ## probability 1-prob) -- the standard zero-inefficiency SFA formulation.
  ## ------------------------------------------------------------------

  ## ZISF (no z) -- targets zsfm(model_name="ZISF", formula=y_zisf~x1+x2).
  ## ZISF's own prob link is prob = exp(-abs(gamma)); using gamma0 = 0.3 here
  ## means the true recoverable gamma is 0.3 (up to the sign ambiguity
  ## |.| already bakes into that link -- either +-0.3 is a correct recovery).
  gamma0 <- 0.3
  prob0 <- exp(-abs(gamma0))
  eff_ind <- rbinom(n, 1, prob0) ## 1 = efficient (u forced to 0)
  y_zisf <- cons + beta1 * x1 + beta2 * x2 + v - u * (1 - eff_ind)

  ## ZISF_Z -- targets zsfm(model_name="ZISF_Z", formula=y_zisf_z~x1+x2|z).
  ## ZISF_Z's prob link is the logistic prob = plogis(z'gamma); reusing the
  ## same z and 0.9/0.6 coefficients as the other z-parameterized columns
  ## above for consistency (even though here they parameterize a mixing
  ## probability, not a variance).
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

  ## ------------------------------------------------------------------
  ## The five distributions below are drawn HERE, after every pre-existing
  ## random draw, and not next to their siblings above. That is deliberate and
  ## must not be 'tidied': inserting an rgamma()/rweibull() call earlier in the
  ## function advances the RNG stream, so every column generated after it comes
  ## out different for the same seed. Placing them last leaves every column that
  ## existed before bit-identical, so earlier Monte Carlo results, the reference
  ## convergence plots and any saved fits all still reproduce.
  ##
  ## Each is matched to the exact parameterization its likelihood in sfm.R uses,
  ## NOT to the textbook default -- they differ in every case, and a mismatch
  ## generates data the estimator cannot recover, which looks like an estimator
  ## bug.
  ## ------------------------------------------------------------------

  ## Normal gamma -- targets sfm(model_name="NG").
  ## sfm.R's NG likelihood carries `- mu*log(sig_u)` and `sig_v/sig_u` inside the
  ## parabolic cylinder function, which is the standard Greene normal-gamma form
  ## with SHAPE = mu and SCALE = sig_u. shape_g is chosen and the scale set to
  ## sig_u/shape_g so that E[u] = shape*scale = sig_u, keeping u on the same
  ## scale as the half-normal column. True values to recover are therefore
  ## sigu = sig_u/shape_g and mu = shape_g, NOT sig_u and 1.
  u_g <- rgamma(n, shape = shape_g, scale = sig_u / shape_g)

  ## Normal Nakagami -- targets sfm(model_name="NNAK").
  ## sfm.R's NNAK likelihood uses sqrt(2*mu*sig_v^2 + sig_u^2) and
  ## lnDv(-2*mu, ...), which identifies mu as the Nakagami SHAPE m and sig_u^2 as
  ## its SPREAD Omega (so sig_u is the root-mean-square of u, E[u^2] = sig_u^2).
  ## A Nakagami(m, Omega) variate is the square root of a Gamma(m, Omega/m).
  ## m = 1 is Rayleigh; m = 0.5 would collapse NNAK onto the half-normal and so
  ## would not test anything NHN does not already cover.
  u_nak <- sqrt(rgamma(n, shape = m_nak, scale = sig_u^2 / m_nak))

  ## Normal generalized-exponential -- targets sfm(model_name="NGE").
  ## sfm.R's NGE likelihood is log(2*lam) + log(T1 - T2) with lam = 1/sig_u,
  ## i.e. the generalized exponential with its SHAPE FIXED AT 2:
  ## f(u) = 2*lam*exp(-lam*u)*(1 - exp(-lam*u)), F(u) = (1 - exp(-lam*u))^2.
  ## Inverting F gives u = -log(1 - sqrt(U))/lam. E[u] = 1.5*sig_u.
  u_ge <- -log(1 - sqrt(runif(n))) * sig_u

  ## Normal lognormal -- targets sfm(model_name="NLN").
  ## sfm.R draws u as exp(shp + sigu*qnorm(h)), i.e. meanlog = the third reported
  ## parameter and sdlog = sig_u. mu_ln = -0.5 with sdlog = 1 puts E[u] at
  ## exp(-0.5 + 0.5) = 1, again matching the half-normal column's scale.
  u_ln <- rlnorm(n, meanlog = mu_ln, sdlog = sig_u)

  ## Normal Weibull -- targets sfm(model_name="NW").
  ## sfm.R draws u as sigu * (-log(1-h))^(1/shp), the Weibull inverse CDF with
  ## SCALE = sig_u and SHAPE = the third reported parameter.
  u_w <- rweibull(n, shape = k_w, scale = sig_u)

  ## Outcomes for the five distributions added above. Each pairs with exactly one
  ## model_name -- see DATA_GENERATION_REFERENCE.md for the true-value triples.
  y_pcs_g <- cons + beta1 * x1 + beta2 * x2 + v - u_g ## NG
  y_pcs_nak <- cons + beta1 * x1 + beta2 * x2 + v - u_nak ## NNAK
  y_pcs_ge <- cons + beta1 * x1 + beta2 * x2 + v - u_ge ## NGE
  y_pcs_ln <- cons + beta1 * x1 + beta2 * x2 + v - u_ln ## NLN
  y_pcs_wb <- cons + beta1 * x1 + beta2 * x2 + v - u_w ## NW

  ## Skew-t -- targets sfm(model_name="THT"). Use THIS column, not y_pcs_t.
  ##
  ## y_pcs_t above builds the composed error from two INDEPENDENT rt() draws,
  ## v_t = sig_v*rt(n,a) and u_t = sig_u*abs(rt(n,a)). That shares the degrees of
  ## freedom between the two components but NOT the mixing variable, and those
  ## are different things. Tancredi (2002) eq (2) requires (eps_i, v_i) to be
  ## JOINTLY bivariate-t, equivalently eq (5): y = h + (eps - z)/sqrt(lambda)
  ## with a SINGLE lambda ~ Gamma(a/2, a/2) common to both components. Only then
  ## is the composed error skew-t, which is what sfm()'s THT likelihood assumes.
  ## Independent draws miss the skew-t 99th percentile by ~0.1 at sig_u = 1,
  ## sig_v = 0.3, a = 5, so THT could not recover its own parameters from
  ## y_pcs_t however correct the likelihood was.
  ##
  ## y_pcs_t is left in place, and left wrong, because removing or reordering it
  ## would advance the RNG stream and change every column generated after it.
  ## It is retained only for backward compatibility -- see the header note above.
  lam_st <- rgamma(n, shape = a / 2, rate = a / 2) ## common mixing, eq (5)
  v_st <- rnorm(n, 0, sig_v) / sqrt(lam_st) ## t noise,    scale sig_v
  u_st <- abs(rnorm(n, 0, sig_u)) / sqrt(lam_st) ## half-t ineff, scale sig_u
  y_pcs_st <- cons + beta1 * x1 + beta2 * x2 + v_st - u_st ## THT

  ## Student-t--half-normal -- targets sfm(model_name="tHN"). Distinct from THT
  ## above: there the two components share ONE scale mixture, so both are t with
  ## the same df and the composed error is skew-t. Here the noise is t and the
  ## inefficiency is half-normal, drawn INDEPENDENTLY, so the tails differ and
  ## the composed density has no closed form.
  ##   v ~ sig_v * t_a          (heavy-tailed noise)
  ##   u ~ |N(0, sig_u^2)|      (conventional one-sided inefficiency)
  ## True values for sfm(model_name="tHN") are (sigv, sigu, nu) = (sig_v, sig_u, a).
  v_thn <- sig_v * rt(n, a)
  u_thn <- abs(rnorm(n, 0, sig_u))
  y_pcs_thn <- cons + beta1 * x1 + beta2 * x2 + v_thn - u_thn ## tHN

  ## Normal-Rayleigh -- targets sfm(model_name="NR").
  ## PARAMETERIZATION. sfm()'s NR likelihood carries sigma_u such that
  ## E[u^2] = sigma_u^2 -- the same second-raw-moment convention NHN uses, which
  ## keeps sigma_u comparable across the two models. The Rayleigh SCALE is
  ## therefore theta = sig_u/sqrt(2), giving
  ##   E[u]   = sig_u*sqrt(pi)/2  ~= 0.886*sig_u
  ##   Var(u) = (1 - pi/4)*sig_u^2
  ## A Rayleigh variate is the length of a 2-D standard normal vector, which is
  ## written out here rather than hidden in a quantile transform so the scale
  ## convention is visible. True values for sfm(model_name="NR") are
  ## (sigv, sigu) = (sig_v, sig_u), as for the other raw-sigma models.
  u_r <- (sig_u / sqrt(2)) * sqrt(rnorm(n)^2 + rnorm(n)^2)
  y_pcs_r <- cons + beta1 * x1 + beta2 * x2 + v - u_r ## NR

  data_trial <- as.data.frame(cbind(
    name, cons, x1, x2, u, uz, v, u_t, v_t, u_c, v_c, u_e, u_u, u_tn,
    y_pcs, y_pcs_t, y_pcs_e, y_pcs_c, y_pcs_u, y_pcs_z, y_pcs_w, y_pcs_tn, z,
    uz_e, y_pcs_ez,
    w_tt, y_ttne, w_tt_hn, y_tthn, zp, wz_hn, y_tthn_z,
    eff_ind, y_zisf, prob_z_true, eff_ind_z, y_zisf_z,
    u_g, y_pcs_g, u_nak, y_pcs_nak, u_ge, y_pcs_ge, u_ln, y_pcs_ln, u_w, y_pcs_wb,
    lam_st, u_st, v_st, y_pcs_st, v_thn, u_thn, y_pcs_thn, u_r, y_pcs_r
  ))
  colnames(data_trial) <- c(
    "name", "cons", "x1", "x2", "u", "uz", "v", "u_t", "v_t", "u_c", "v_c", "u_e", "u_u", "u_tn",
    "y_pcs", "y_pcs_t", "y_pcs_e", "y_pcs_c", "y_pcs_u", "y_pcs_z", "y_pcs_w", "y_pcs_tn", "z",
    "uz_e", "y_pcs_ez",
    "w_tt", "y_ttne", "w_tt_hn", "y_tthn", "zp", "wz_hn", "y_tthn_z",
    "eff_ind", "y_zisf", "prob_z_true", "eff_ind_z", "y_zisf_z",
    "u_g", "y_pcs_g", "u_nak", "y_pcs_nak", "u_ge", "y_pcs_ge",
    "u_ln", "y_pcs_ln", "u_w", "y_pcs_wb",
    "lam_st", "u_st", "v_st", "y_pcs_st", "v_thn", "u_thn", "y_pcs_thn",
    "u_r", "y_pcs_r"
  )

  data_rand <- data.frame(data_trial)
  data_rand$con <- 1
  return(data_rand)
}
