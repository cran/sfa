# sfa 1.1.5

* **The `intro_to_psfm` vignette now builds in a fraction of the time.** At
  CRAN's request, the models it fits are sized as toy illustrations rather
  than as estimation exercises: the simulated panel is 70 firms over 6
  periods instead of 100 over 10, the simulated-ML fits draw 50 Halton
  points via `halton_num` rather than the default
  `ceiling(sqrt(nrow(data))) + 100`, and the `psfm_bootstrap()` example
  uses 5 replications on a 30-firm panel instead of 10 on 60. Vignette
  rebuild time falls by roughly a factor of five. The reported estimates
  therefore differ from previous versions, and the vignette now says
  plainly that they are not to be read as a serious fit.

* Every model fit in that vignette now sets `rand.gtre` and `rand.psoptim`.
  With `PSopt = TRUE` the particle-swarm stage draws from the session's RNG,
  so the vignette's results previously changed from one build to the next;
  they are now reproducible.

* No change to any R code, to `NAMESPACE`, or to the documented interface.

# sfa 1.1.4

## Breaking change

* **`psfm(model_name = "GTRE")` now defaults to full information maximum
  likelihood.** The four ways of fitting the four-component GTRE model were
  four separate `model_name` values, which made them look like four different
  models rather than four routes to the same one. They are now selected with an
  `estimator` argument, in the same spirit as `sfm()`'s `estimator = c("mle",
  "cols")`:

  - `"fiml"` (the default) -- full information ML through the closed-skew-normal
    representation. Deterministic; requires a balanced panel.
  - `"sml"` -- simulated ML over Halton draws. Handles unbalanced panels. **This
    is what `"GTRE"` meant through 1.1.3.**
  - `"seq1"`, `"seq2"` -- the two-step moment estimators.

  Scripts that pass `model_name = "GTRE"` therefore get a different estimator
  than they did, and are warned once per call. Pass `estimator` explicitly to
  silence it. The names `"GTRE_FML"`, `"GTRE_SEQ1"` and `"GTRE_SEQ2"` are
  unchanged and still select the same three routes directly.

  On an **unbalanced** panel `"fiml"` cannot be fitted. Taking the default
  warns and falls back to `"sml"`, because erroring would make `"GTRE"`
  unusable by default on a whole class of data; asking for `"fiml"` explicitly
  errors instead of silently fitting something else.

## New models

* **`psfm(model_name = "PL80_MVTN")` -- Pitt and Lee's (1981) Model III**, the
  multivariate truncated normal panel likelihood from their Appendix 2.
  Inefficiency varies over time *and* is correlated within a firm:
  `u_i = (u_i1,...,u_iT)'` is drawn from a T-variate normal truncated to the
  negative orthant, where the existing `"PL80"` holds inefficiency fixed over
  time.

  Pitt and Lee derived this likelihood and then set it aside, writing that it
  "is difficult to evaluate since the quantities P0 and P(y_i - x_i beta)
  involve T-dimensional numerical integrals", and estimating Model III by
  Zellner SUR instead. Those quantities are multivariate normal orthant
  probabilities; `mnormt::sadmvn()` evaluates one in about 3 ms at T = 6, so a
  likelihood evaluation costs roughly N+1 of them -- about 0.3 s at N = 100,
  and about 12 s for a whole fit at N = 80, T = 4. What was intractable in
  1981 is merely slow now.

  Sigma is parameterized as **equicorrelated**, `sigma_u^2[(1-rho)I + rho 11']`,
  costing two parameters. Unrestricted Sigma costs `T(T+1)/2` -- 21 at T = 6,
  55 at T = 10 -- every one identified only through orthant probabilities. The
  equicorrelated form captures what the general Sigma was for: dependence of a
  firm's inefficiency across periods, with `rho = 0` giving independence over
  time and `rho -> 1` approaching the time-invariant `"PL80"`. Because the form
  is equicorrelated, every matrix quantity is closed form (Sherman-Morrison,
  verified against brute force to 1e-17), so only the orthant probabilities are
  numerical.

  Requires a **balanced** panel with T >= 2; an unbalanced one errors and points
  at `"PL80"`. `data_gen_p()` gains `y_pl_mvtn` and `u_mvtn` columns and a
  `rho_mvtn` argument to test it, generated last so the RNG stream feeding every
  existing column is untouched. Registered in the convergence framework as
  `PL80_MVTN`. See `?PL80_MVTN`.

  Validated before wiring: profiling the likelihood one parameter at a time
  puts the minimum of the negative log-likelihood at the truth for all five
  parameters, and maximizing it recovers `(sigma_v, sigma_u, rho)` =
  (0.32, 0.78, 0.45) against a truth of (0.30, 0.80, 0.50) at N = 300.

* **`npsfm()`, nonparametric stochastic frontier models.** A fifth entry point,
  alongside `sfm()`, `psfm()`, `zsfm()` and `ttsfm()`, for frontiers whose shape
  is estimated by kernel regression rather than assumed linear. Two estimators:

  - `method = "FLW"` -- Fan, Li and Weersink (1996). Fits `E[y|x]`
    nonparametrically, then recovers the scale parameters from the residuals,
    by maximizing their concentrated likelihood in `lambda` for
    `dist = "hn"` and by inverting central moments for `dist = "exp"`,
    `"gamma"` and `"unif"`.
  - `method = "SVKZ"` -- Simar, Van Keilegom and Zelenyuk (2017). Three
    local-linear regressions give `sigma_u(x)` and `sigma_v(x)` pointwise, so
    both variance components vary with the covariates. No optimizer runs.
    Normal-half normal only.
  - `method = "PSZ"` (alias `"KPST"`) -- Park, Simar and Zelenyuk. Local
    maximum likelihood: the frontier and both log variance components get
    local-linear expansions and the kernel-weighted normal-half normal
    likelihood is maximized in those `3(k+1)` parameters, once per observation.
  - `method = "MY"` -- Martins-Filho and Yao. Iterative local likelihood,
    alternating local frontier fits with a global update of `(lambda, sigma)`.
  - `method = "SZ"` -- Simar and Zelenyuk (2011). Passes an existing smooth
    frontier through an output-oriented DEA to impose monotonicity and
    convexity. Needs the **Benchmarking** package, also in `Suggests`.

  `"PSZ"` and `"MY"` run one numerical optimization per observation (for
  `"MY"`, per observation per iteration), so they are one to two orders of
  magnitude slower than `"FLW"`. Both are seeded from an `"FLW"` fit.

  Ported from Christopher Parmeter's research scripts. Results return as class
  `"npsfareg"` rather than `"sfareg"`: there is no parameter vector with
  standard errors, so `coef()`, `vcov()` and `logLik()` would have nothing to
  return. `fitted()`, `residuals()`, `nobs()`, `print()` and `summary()` are
  provided.

  Kernel regression comes from the **np** package, added to `Suggests` rather
  than `Imports` -- nothing else in `sfa` needs it, and `npsfm()` checks for it
  and stops with an install instruction if it is absent.

  A correction worth recording, because it is easy to repeat: the two
  local-likelihood estimators maximize the *composed-error* likelihood, in
  which the local intercept is the frontier `m(x)` itself. They therefore take
  **no** half-normal mean shift, unlike the least-squares methods, whose first
  stage estimates `E[y|x] = m(x) - E[u]` and does need one. Applying the shift
  to `"PSZ"`/`"MY"` biases the whole frontier up by about `E[u]`; mean absolute
  frontier error at `n = 300` fell from 0.372 to 0.116 (`"PSZ"`) and 0.410 to
  0.070 (`"MY"`) once it was removed.

  Against a simulated nonlinear frontier with `sigma_u = 0.6`, `sigma_v = 0.25`,
  both least-squares estimators converge as `n` grows (6 replications at each size):
  `FLW` recovers `sigma_u` = 0.560, 0.539, 0.595 at `n` = 150, 300, 600, and
  `SVKZ` 0.477, 0.506, 0.559, with mean absolute frontier error falling from
  0.102 to 0.043 and 0.206 to 0.079 respectively. `SVKZ`'s downward bias at
  small `n` is the wrong-skew floor: the share of observations whose local
  third moment has the wrong sign, and whose `sigma_u(x)` is therefore set to
  zero, falls from 21.6% to 0.3% over that range.

## New methods and arguments

* **`sfa_diagnostics()`, `plot()` for `"sfareg"`, and convergence reporting.**
  The numerical hardening was already in place -- staged minimizer, clipping
  constants, analytic gradients where they exist -- but nothing was reported
  back. A fit carried `optim()`'s convergence code, message, evaluation counts
  and Hessian, and `print()`/`summary()` showed none of it, so a fit that
  stopped on the iteration cap printed exactly like a converged one.

  `sfa_diagnostics()` returns the convergence code and what it means, the
  eigenvalue spectrum and condition number of the Hessian, whether it is
  positive definite, which parameters load on its flattest direction, the
  implied parameter correlations, and -- with `keep_objective = TRUE` -- the
  gradient at the reported optimum. `plot()` draws the Hessian spectrum, the
  correlation matrix, a likelihood slice per parameter, and the gradient.
  `print()` and `summary()` now report the convergence code.

  **The code by itself is not diagnostic and is not treated as though it were.**
  Across `NHN`, `NE` and `NTN` at n = 150, 500 and 1500, code 52
  ("ABNORMAL_TERMINATION_IN_LNSRCH") appears routinely alongside a maximum
  relative gradient of ~1e-6 and a positive definite Hessian: the staged
  minimizer had already converged and `L-BFGS-B` could not step away from the
  optimum. The same code on `NTN` at n = 150 came with a relative gradient of
  5e+07 and an indefinite Hessian, a real failure. The verdict therefore
  combines the code with the gradient and the Hessian, and distinguishes
  *benign* from *unverified* (a line-search code with no objective retained, so
  no gradient to settle it) from *failure*. Code 1, the iteration limit, is
  never treated as benign.

  On a single `NNAK` fit the report reproduces what the convergence sweeps
  found only across replications: `mu` and `sigu` correlated at 0.998 -- the
  documented ridge -- with the flattest Hessian direction loading on exactly
  that pair.

  When the Hessian is singular enough that some parameter has no usable
  variance, the correlation report drops those parameters and names them,
  rather than vanishing as a whole -- a diagnostic for ill-conditioning should
  be most informative exactly when conditioning is worst, not least.

* **`sfm(keep_objective = TRUE)`** stores the likelihood on the fitted object so
  the gradient and likelihood slices can be computed after the fact. Off by
  default: a closure carries its enclosing environment, so a fit saved with one
  serializes the estimation data too (about 38 KB to 1.7 MB on a 200-observation
  example). Everything else `sfa_diagnostics()` reports works without it.

## New features

* **Model names are matched without regard to case.** `match.arg()` is case
  sensitive, so `psfm(model_name = "gtre")` used to fail with a list of valid
  names that visibly contained what the user had just typed. All five entry
  points now fold case, for `model_name` and for `npsfm()`'s `method`. No entry
  point has a case collision among its choices -- `sfm()`'s `"THT"` and `"tHN"`
  differ in more than case -- so the canonical spelling is always recoverable.

  Exact matches beat partial ones, which matters because `"GTRE"` is a prefix
  of four other names and must resolve to itself. Genuinely ambiguous partials
  (`"GTRE_S"`, between `GTRE_SEQ1` and `GTRE_SEQ2`) are still rejected rather
  than guessed at, and an unrecognized name now suggests the two closest valid
  choices instead of listing everything.

# sfa 1.1.3

## Breaking change

* **`psfm()`'s optimizer iteration defaults have been raised**, from
  `maxit.bobyqa = 100`, `maxit.psoptim = 10`, `maxit.optim = 10` to
  `5000`, `100` and `1000`. `maxit.nlminb` is now an argument (default `500`);
  it was previously hard-coded at 200 in the `"GTRE_FML"` branch and 500
  elsewhere, and could not be set from the call.

  The old caps were binding rather than merely economical. The
  `"K1990"`/`"K1990modified"` code already carried a note that 100 bobyqa
  evaluations left its seven-parameter fits several log-likelihood units short
  of the optimum purely on the iteration cap, and `"GTRE_FML"` at N = 500,
  T = 10 roughly halves its root-mean-square error against known true values
  once the caps are lifted (0.0080 to 0.0037 and 0.0179 to 0.0073 on two
  draws), for about 1.5 times the run time.

  Existing scripts will get more accurate estimates and slower fits. Pass the
  old values explicitly to restore the previous behaviour.

* **`psfm(model_name = "TFE")` now fits a different estimator.** Through
  version 1.1.2 the name `"TFE"` selected Chen, Schmidt and Wang's (2014)
  *within* maximum-likelihood estimator. It now selects Greene's (2005) *true
  fixed effects* estimator, which is what the name means in the literature.
  The Chen-Schmidt-Wang estimator is unchanged and is now
  `model_name = "TFE_WMLE"`.

  Scripts written against 1.1.2 or earlier that pass `"TFE"` will silently get
  a different estimator, so `psfm()` issues a warning whenever `"TFE"` is used.
  To reproduce earlier results, change the name to `"TFE_WMLE"`.

## New features

* **Corrected ordinary least squares, `sfm(estimator = "cols")`.** The moment
  estimator of Olson, Schmidt and Waldman (1980). OLS is consistent for the
  slopes of a composed-error frontier whatever the one-sided distribution --
  only the intercept is biased, by `E[u]` -- so COLS keeps the OLS slopes,
  inverts the central moments of the OLS residuals for the scale parameters,
  and shifts the intercept up by the implied `E[u]`.

  Implemented for `"NHN"`, `"NE"` and `"NG"`; other models error, because the
  moment inversion is distribution-specific. No optimizer runs and the result
  is deterministic, which makes it a natural robustness check against a
  maximum-likelihood fit that may have settled at a local optimum.

  Wrong-skew samples are reported rather than absorbed: a production frontier
  implies a negative third central moment, and when a sample comes out the
  other way the moment equations have no admissible solution. `sfm()` warns,
  reports `sigu = 0` with the whole residual variance assigned to `sigv`, and
  returns no efficiency predictions -- to be read as *no evidence of
  inefficiency in these data*, not as an estimate of zero.

  Standard errors: the OLS slope standard errors are reported and are valid as
  such. The scale parameters and the corrected intercept carry `NA` by
  default, since neither has a closed-form standard error here and the OLS
  intercept standard error would be wrong (it knows nothing about the sampling
  error of a third-moment estimate). Set `cols_boot` for a nonparametric
  bootstrap covering every parameter, with `rand.cols` to make it
  reproducible.

## New models

* **`sfm(model_name = "tHN")`** -- Student's t--half-normal. Heavy-tailed
  *noise* (`v ~ sigma_v * t_nu`) with a conventional half-normal inefficiency
  term (`u ~ |N(0, sigma_u^2)|`), drawn independently.

  This is **not** `THT`. In `THT` (Tancredi 2002) both components come from one
  shared scale mixture, so both are t with the same degrees of freedom and the
  composed error is a closed-form skew-t. In `tHN` the tails differ, there is no
  closed form, and the density is the convolution
  `f(e) = integral_0^Inf f_v(e+u) f_u(u) du`, evaluated by Gauss-Legendre
  quadrature. `tHN` is therefore the natural parametric comparison for the
  density-power robust estimators (`robust = "mlqe"/"psi"/"mdpd"`), which `THT`
  cannot be, because its inefficiency term is heavy-tailed too. Parameters are
  reported as `(sigv, sigu, nu)` -- the conventional order, not `THT`'s
  inverted one. Returns `exp_u_hat` and `u_hat` by a Bayes rule over the same
  quadrature nodes.

  Two documented properties, both surfaced rather than hidden. **The degrees of
  freedom are weakly identified**: on data simulated from the model at
  `n = 1000` with a true `nu = 5`, the profile log-likelihood moves only about
  0.24 across `nu` from 10 to 100, and peaks near 20. Profile over a grid of
  fixed `nu` rather than reporting one selected value. Because of that flat
  ridge every `tHN` fit runs from several widely separated starts, keeps the
  best, and reports the outcome in `thn_starts`, warning when the starts reach
  different optima. **`sigma_u` can collapse onto zero** on real data, the heavy
  noise tail absorbing the whole one-sided component and leaving mean predicted
  efficiency near one; `sfm()` warns and sets `thn_sigma_u_at_bound` rather than
  bounding `sigma_u` away from zero, because the collapse is a property of the
  model and is the thing a user needs to see.

  The quadrature node count scales with `sigma_u/sigma_v` and is not fixed. A
  fixed 96-node rule is accurate near `lambda = 3` but carries 4% relative
  error at `lambda = 20` and 60% at `lambda = 62` -- and the model does reach
  that region. Note that integrating the density to 1 does not detect this: the
  error redistributes across the support and integrates away, so total mass
  still reads 1.000 while the density is 40% wrong pointwise.

* **`data_gen_cs()` gained `y_pcs_thn`** (with `v_thn` and `u_thn`), the
  matching generator for `tHN`.


* **`psfm(model_name = "TFE")`** -- Greene (2005) true fixed effects. The
  composed-error likelihood with one intercept per individual, estimated as a
  profile likelihood in `(lambda, sigma, beta)` with the firm effects
  concentrated out. Reports the same parameter layout as `"TFE_WMLE"`, plus
  `r_hat_m` (the maximum-likelihood firm effects), `exp_u_hat` and `u_hat`.

  Note that this likelihood always has a supremum on the `sigma_v = 0`
  boundary, because the individual effects are unrestricted. The new argument
  `tfe_lambda_max` (default 100) bounds the search accordingly, and a fit that
  pins at the bound warns. See `?psfm` for the details; this is a property of
  the estimator and is one of the motivations for `"TFE_WMLE"`.

* **`psfm(model_name = "K1990")` and `"K1990modified"`** -- Kumbhakar (1990)
  time-varying inefficiency, `B_it = (1 + exp(bt + ct^2))^-1` and
  `B_it = 1 + d(t - T_i) + e(t - T_i)^2`. These share one likelihood with
  `"PL80"` and `"BC92"`, differing only in `B_it`. `K1990`'s `b` and `c` are
  weakly identified, so the fitted `B_it` path is more interpretable than
  either coefficient on its own.

* **Four new inefficiency distributions in `sfm()`**: `"NU"` (normal-uniform),
  `"NGE"` (normal-generalized exponential), `"NLN"` (normal-lognormal) and
  `"NW"` (normal-Weibull). The last two are estimated by simulated maximum
  likelihood over Halton draws.

## New methods and arguments

* `predict()`, `fitted()` and `residuals()` methods for class `"sfareg"`,
  alongside the existing `coef()`, `vcov()`, `logLik()` and `nobs()`.
  `predict()` accepts `newdata`.

* `psfm()` now accepts an ordinary `data.frame` (or tibble/data.table) as well
  as a `plm::pdata.frame`; the panel index is constructed internally from
  `individual` and the new `time` argument. Previously a plain data frame
  failed with an uninformative `"empty model"` error.

* `psfm(collinear_action =)` controls what happens when the
  *between-individual* design used to build starting values is rank deficient
  -- the situation created by, for example, time dummies, which are estimable
  in a pooled specification but collapse onto the intercept once averaged
  within each unit. `"start_only"` (default) keeps the requested model and
  drops the offending columns from the starting-value regression only;
  `"error"` stops and names them; `"warn_drop"` removes them from the
  estimated model. Previously this surfaced as an opaque
  `solve(crossprod(ZBeta))` LAPACK error inside `plm`.

* `sfm(robust =)` selects a divergence-based robust estimator -- `"mlqe"`,
  `"psi"` or `"mdpd"` -- with sandwich standard errors. Currently implemented
  for `model_name = "NHN"`; other models error rather than silently ignoring
  the argument.

* `sfm()` gained `use.nlminb` and `use.bobyqa` (both `"auto"` by default) and
  `maxit.nlminb`, for control over the optimizer stack described below.

## Estimation and performance

* **`nlminb` now leads the optimizer stack** for the models where it helps,
  ahead of the derivative-free stages, with an analytic gradient supplied for
  `NHN`. Model-by-model defaults are chosen automatically: `NHN`, `NE`, `NTN`
  and `NU` use it, while `NR` and `NGE` (where it degraded the fit) do not.
  Typical cross-sectional fits are several times faster with unchanged
  estimates.

* An `nlminb` stage was added to the `PL80`/`BC92`/`K1990`/`K1990modified`
  branch, worth up to +25 log-likelihood on the seven-parameter models, where
  `psfm()`'s low default `maxit.bobyqa` was stopping short of convergence.
  `BC92` now matches or beats `frontier::sfa()` on parameter, variance and
  efficiency accuracy.

* `PL80` and `BC92` are now estimated by a native maximum-likelihood
  implementation instead of wrapping `frontier::sfa()`. Verified against
  `frontier` (coefficients, log-likelihood and predicted efficiencies) on
  balanced and unbalanced panels and on both production and cost
  specifications before the dependency was removed.

* Efficiency prediction (`exp_u_hat`) is now available for `NE` and `NTN`,
  which previously returned none.

## Bug fixes

* **`psfm(OPG_calc = TRUE)` returned `NA` for every OPG and sandwich standard
  error, and wrote a variable into the global environment.** The OPG "meat"
  matrix was stored with a superassignment, `OPG_meat <<- crossprod(score_mat)`,
  on the assumption that the surrounding `tryCatch({...})` introduced a scope to
  escape from. It does not -- `tryCatch()` evaluates its expression in the
  calling frame -- so `<<-` began its search one frame further out, skipped the
  local `OPG_meat <- NULL` binding entirely, and assigned into `globalenv()`.
  The immediately following `solve(OPG_meat)` therefore still saw `NULL` and
  failed, as did the `MASS::ginv()` fallback, so the OPG standard errors were
  always `NA`; the sandwich errors, which reuse the same matrix, were `NA` with
  them. Both failures were reported through the existing handlers as
  `"OPG matrix singular, using pseudoinverse"` and `"OPG matrix is singular"`,
  which pointed at a rank problem in the data rather than at the scoping bug.
  Changed to a plain `<-`. On a fixed-seed `GTRE_Z` fit the parameter estimates
  and Hessian standard errors are bit-identical before and after, the OPG and
  sandwich errors change from `NA` to finite values, the two warnings stop
  firing, and `OPG_meat` no longer appears in the user's workspace.

* **`psfm_bootstrap()` failed on Windows whenever `sfa` was not installed in a
  default library.** The function distributes work over PSOCK cluster workers,
  which start as fresh R sessions holding the *default* library path rather than
  the parent session's. Where the parent had found `sfa` somewhere else -- a
  project library, `renv`/`packrat`, a user-set `R_LIBS_USER`, or the temporary
  `sfa.Rcheck` tree that `R CMD check` installs into -- the workers' `library()`
  call could not see it, and the bootstrap died in
  `parallel:::checkForRemoteErrors()` with `there is no package called 'sfa'`.
  On Unix the workers generally inherit `R_LIBS` from the parent's environment,
  which hid the bug; on Windows they do not. Where the workers instead found a
  *different*, older copy of `sfa` in a default library, they silently ran the
  bootstrap against that version rather than the one the user had loaded. The
  parent's `.libPaths()` is now pushed to the workers before any package is
  loaded, so both sessions resolve every package identically.

* **`stats::dlnorm` was used without being imported.** The rewritten `"NLN"`
  likelihood and its efficiency block call `dlnorm()`, but `NAMESPACE` did not
  import it, so the call resolved only via the search path rather than the
  package namespace. `R CMD check --as-cran` reported it as an undefined global.
  Now imported explicitly.

* **`zsfm()`'s efficiency predictor used a different mixing probability from
  the likelihood it maximised.** `"ZISF"`'s likelihood sets
  `prob = exp(-abs(gamma))`, which makes it exactly symmetric in `gamma`: `+g`
  and `-g` fit identically and the optimizer may return either. The JLMS block
  then used `exp(-gamma)`, so a negative estimate produced a mixing
  "probability" above 1 and silently invalid `post.prob` and `jlms`. This was
  reachable from an ordinary starting value, not a pathology: seeding at the
  negated estimate returns `gamma = -0.3015` with an identical log-likelihood,
  where the old code computed `prob = 1.3519`. Both places now use
  `exp(-abs(gamma))`.

* **`zsfm()`'s two-component mixture is now formed on the log scale.** It built
  `prob*exp(f1) + (1-prob)*exp(f2)` and took `log(f + 1e-10)`. Both terms
  underflow to zero when an observation is unlikely under either regime, and
  the `1e-10` guard then floors that observation's contribution at `-23.03` --
  which does not merely protect the logarithm, it makes the objective *flat*
  across the whole region beyond the floor, exactly where the optimizer needs a
  gradient. A new internal helper, `.log_add2()`, computes
  `log(exp(a) + exp(b))` without leaving the log scale, and `post.prob` is now
  a ratio taken in logs. Fitted coefficients are unchanged to within optimizer
  path noise (worst discrepancy 1.9e-3 over eight seeds, log-likelihoods
  agreeing to 1e-6).

* **`zsfm(logit = TRUE)` now uses `plogis()`** instead of
  `exp(eta)/(1 + exp(eta))`, which overflows to `Inf/Inf = NaN` once the linear
  predictor passes about 710 -- a value the optimizer can reach while
  searching. Identical wherever the old form was finite.

* **`zsfm()`'s efficiency block now branches on `model_name`** rather than on
  `is.na(n_z_vars)`. The likelihood already branched on the model, and the
  parameter layout is a property of the model; keying the predictor off a
  different condition meant the two could disagree, reading `"ZISF"`'s
  parameters under `"ZISF_Z"`'s layout if `n_z_vars` ever arrived as `0`
  rather than `NA`.

  Note left in place, not changed: the `logit = FALSE` branch computes
  `pnorm(eta)/(1 + pnorm(eta))`, which is bounded above by 0.5 and is not the
  probit link. Correcting it would change results for that option, which is a
  modelling decision rather than a cleanup.

* **`sfm(model_name = "NLN")` was integrating a spike.** Its simulated
  likelihood averaged the normal kernel over lognormal draws, but the kernel is
  only `sigma_v` wide in `u` while the lognormal spreads over decades, so
  nearly every draw landed where the kernel is numerically zero and a handful
  carried the whole integral. Measured against a reference verified two ways
  (adaptive quadrature and a 200,000-point Simpson rule, agreeing to 5e-9), the
  simulated log-likelihood at the true parameters was **226.8 units low** at
  n = 3000 under the default draw count -- which is precisely the unexplained
  228-unit gap recorded against this model in the convergence notes. It was
  simulation error, not a defect in the likelihood.

  `sfm()` now substitutes `u = sigma_v*t - e`, turning the integral into a
  standard-normal expectation of the *smooth* lognormal density truncated to
  `u > 0`. The error at the same draw count is then 0.12, and 0.86 at `Nsim =
  50`. Raising `Nsim` was not an alternative: the old error per observation was
  about -0.076 at both n = 1000 and n = 3000 under the `8*sqrt(n)` rule, so the
  total bias grew linearly in n at the same rate as the log-likelihood itself,
  and closing it needed `Nsim` proportional to n -- quadratic work per
  evaluation. `"NW"`, which uses the same machinery but a far lighter-tailed
  inefficiency term, is unchanged.


* **`sfm(model_name = "NR")` was started from a flat guess it could not
  recover from, and was documented as the wrong model.** Two separate
  problems, one in the estimator and one in everything written about it.

  The estimator: `start_cs()` hard-codes `sigma_u = sigma_v = 0.1` for the
  cross-sectional models, and from there `"NR"` converged to a point with a
  *worse* log-likelihood than the true parameter vector in 9 of 14
  replications at n = 4000. It failed in two modes -- a hard collapse to
  `sigma_u = 1e-7` with `sigma_v` inflated to absorb the spread (a ~50
  log-likelihood deficit), and a partial stall at `sigma_u` around 0.85 with
  the intercept ~0.4 too low. `"NR"` is now started by inverting the Rayleigh
  moment equations instead, which reaches the same optimum a truth-seeded run
  finds in all 14. Because the Rayleigh skewness is a constant, the third
  central moment of the residuals identifies `Var(u)` outright. Wrongly skewed
  residuals leave the moment equations with no admissible solution, and the
  old flat start is then used unchanged.

  The documentation: `"NR"` had been described in this package as "an
  alternative closed-form derivation of the same normal/half-normal composed
  error" as `"NHN"`. That is wrong. `"NR"` is normal-**Rayleigh** -- its coded
  density reproduces a numerical normal-Rayleigh convolution to 1e-8 and
  misses the half-normal one by 86%. The two are different families, not
  reparameterizations: the standardized skewness the inefficiency contributes
  is a different constant in each (-0.631 against -0.995), and no
  transformation of a two-scale family moves a standardized moment. The
  likelihood itself was correct throughout and is unchanged.

  `data_gen_cs()` gains `u_r` and `y_pcs_r` to test `"NR"` against its own
  data-generating process; it had previously been tested against `y_pcs`,
  which it cannot fit. The new columns are appended at the end of the
  function, so every existing column is bit-for-bit unchanged. `sigma_u` is on
  the convention `E[u^2] = sigma_u^2`, matching `"NHN"`, so the Rayleigh scale
  is `sigma_u/sqrt(2)`.

* **`sfm(model_name = "NG")` and `"NNAK"` read their frontier coefficients
  from the wrong position in the parameter vector.** The likelihood closure
  slices the coefficients out by a fixed offset -- `x[3:(n_x+2)]` for models
  with two leading scale parameters, `x[4:(n_x+3)]` for three. `NG` and `NNAK`
  carry three (`sigv`, `sigu`, `mu`) but sat in the two-parameter group, so
  the slice took the right *number* of coefficients starting one slot too
  early.

  The consequences were severe and entirely silent. `mu` was used
  simultaneously as the gamma (or Nakagami) shape *and* as the intercept
  coefficient, so it could not move freely -- which is why the shape appeared
  never to leave its starting value. Every remaining slope was shifted one
  place, and the last coefficient never entered the likelihood at all, so it
  simply kept whatever starting value it was given. The efficiency block used
  the correct offset throughout, so the two halves of the model disagreed
  about which number meant what.

  On the package's own test DGP at n = 4000, `NG` returned `sigma_u = 0` and
  stopped **710 log-likelihood units below the true parameter vector**. Across
  15 fits spanning n = 1000 to 5000, none reached the truth's likelihood.
  After the fix, all 15 do, no fit collapses, and mean RMSE against the truth
  falls from 0.60 to 0.27 and declines with n. `NNAK` improves on the same
  fix -- 6 of 7 fits now beat the truth and standard errors are finite, where
  previously the Hessian was singular in every replication -- but it still
  produces occasional failures and is not yet considered repaired.

  The NG density itself was never wrong: it agrees with numerical convolution
  of the normal and gamma densities to 7e-14 across the whole residual range,
  and that check is now a test.

* **`sfm(model_name = "NG")` also started in the wrong place.** `start_cs()`
  hard-codes `sigma_u = sigma_v = 0.1` and the NG start pinned the shape at 1,
  so the search began from `E[u] = 0.1` against a true 1, with the intercept
  at the raw OLS value -- itself `E[u]` below the frontier. `NG` now builds
  candidate starts from the residual moments and sweeps the shape along the
  `E[u] = mu*sigma_u` ridge, which is the direction the data leave weakly
  determined (Ritter and Simar, 1997), polishing the most promising few before
  choosing. The search is reported in `$ng_starts`.

* **`psfm(model_name = "GTRE_FML")` started its search from the wrong
  intercept.** It seeded `beta_0` at the raw panel-regression intercept, while
  `"GTRE"` and `"TRE"` seed theirs at that intercept plus `E[u] + E[h]`. In a
  composed-error model the regression intercept sits below the frontier by
  `(sigma_u + sigma_h) * sqrt(2/pi)` -- 1.12 at the package's own test DGP --
  so the FIML search began a full unit low, in exactly the direction of the
  `sigma_h = 0` boundary optimum where the model collapses to `"TRE"` and the
  intercept absorbs the missing `E[h]`.

  Fits that fell in were not merely imprecise: on one of six test draws the
  reported solution had a *lower* log-likelihood than the true parameter
  vector (-5267.87 against -5267.80), with `sigma_h = 2e-16`, `beta_0 = 0.155`
  against a true 0.5, and `sigma_r` inflated to 0.307 against a true 0.2 --
  a local optimum, not an estimate.

  `"GTRE_FML"` now also builds a second candidate start from the two-step
  moment estimator (the one `"GTRE_SEQ2"` reports), evaluates the likelihood
  at both and begins from the better, as Colombi (2010) and Colombi, Martini
  and Vittadini (2011) recommend for this likelihood. Across the six test
  draws, boundary collapses go from one to none, mean RMSE against the truth
  falls from 0.049 to 0.023, and all six fits now reach a higher likelihood
  than the truth. The chosen start is reported in `start_search`.

* **`sfm(model_name = "THT")` used the wrong likelihood.** The skew-t density
  of Tancredi (2002, eq. 4) has scale `omega = sqrt(sigma_v^2 + sigma_u^2)`,
  but the implementation evaluated the Student-t factor at the *raw* residual
  and omitted the `1/omega` Jacobian, which pins the scale at 1. Because
  `2*f(e)*G(w(e))` is a valid density for any symmetric `f` and odd `w`
  (Azzalini's lemma), the wrong version still integrated to 1 and still
  produced plausible fits rather than an obvious failure -- but `sigma_u` and
  `sigma_v` were then identified only through the skewing term, and the
  degrees of freedom `a` had to absorb the scale mismatch. Fitted `a` and
  `sigma_v` were consequently inconsistent. Fixed; `THT` now reproduces
  `sn::dst()` to machine precision.

* **`sfm(model_name = "THT")` now reports efficiency.** It previously returned
  no efficiency prediction at all. It now returns `exp_u_hat` = E[exp(-u)|e],
  `u_hat` = E[u|e], and `sd_exp_u_hat`, following Tancredi (2002, section 2.2).
  The conditional density in that paper's eq. (7) is a Student-t truncated to
  the non-negative half-line with `df = a + 1`, location `-e*sigma_u^2/omega^2`
  and scale `sqrt((a + e^2/omega^2)*sigma_v^2*sigma_u^2/(omega^2*(a+1)))`,
  which reduces to the Jondrow et al. (1982) normal predictor as `a` grows.
  Note the behaviour this buys: for a large *positive* residual the
  half-normal model drives predicted efficiency to 1 with near-zero
  uncertainty, whereas the skew-t model reads the point as an outlier, so
  efficiency turns back down and `sd_exp_u_hat` widens.

* **`sfm(model_name = "THT")` degrees of freedom are better started and
  bounded.** The starting value was 1 -- the Cauchy case, which has neither a
  mean nor a variance -- and the lower bound was 1e-7. The start is now a
  moment estimate from the excess kurtosis of the OLS residuals
  (`a ~ 4 + 6/kurtosis`, clipped to `[3, 30]`), and the lower bound is 2.05, so
  the search stays where the composed error has both moments.

* **`data_gen_cs()` gained `y_pcs_st`, which is the column `THT` should be
  tested against.** The existing `y_pcs_t` draws its two error components as
  two *independent* `rt()` variates; that shares the degrees of freedom but not
  the mixing variable, so the composed error is not skew-t and `THT` cannot
  recover its own parameters from it. `y_pcs_st` uses the single common
  `Gamma(a/2, a/2)` mixing of Tancredi eq. (5), with `lam_st`, `u_st` and
  `v_st` also returned. `y_pcs_t` is retained unchanged, because renumbering
  the random draws would alter every column generated after it.

* `nobs()` on an `sfm()`, `zsfm()` or `ttsfm()` fit returned `NA` when called
  from inside a function, because the recorded `data` argument was re-evaluated
  in the caller's frame rather than the one the model was fitted in. This
  propagated silently to `BIC()`, which needs the observation count, while
  `AIC()` kept working.

* `logLik()` on the estimators that are not maximum-likelihood
  (`GTRE_SEQ1`, `GTRE_SEQ2`, `SSFE`) returned an unclassed `NA`, which made
  `AIC()` and `BIC()` return `numeric(0)` -- a missing value that disappeared
  instead of propagating. They now return `NA` as documented.

* `psfm()`'s `TFE` and `SSFE` models silently depended on the first two
  columns of `data` being the panel index, because `plm()` was called without
  an explicit `index`. Any data whose first two columns were something else
  produced `"empty model"`.

* Fixed a `checkSymmetricPositiveDefinite()` failure in the `GTRE` models
  caused by asymmetric `dimnames` on an otherwise symmetric covariance matrix.

* The Halton draws used by the simulated-likelihood models were reshaped in
  column-major order, so each observation drew from a narrow, non-equidistributed
  slice of the sequence.

* A singular Hessian now yields `NA` standard errors rather than aborting the
  whole fit.

## Internal changes

* All 76 `stop()` and `warning()` calls in `R/` now pass `call. = FALSE`; 24 of
  them did not, which made the error output inconsistent between older and newer
  code paths.

* The seven `sapply()` calls in `psfm.R` -- six extracting ridge/method
  diagnostics from the GTRE posterior solver, one building the transient
  efficiency vector -- are now `vapply()` with explicit `numeric(1)` and
  `character(1)` templates, so a change in what the solver returns fails loudly
  instead of silently producing a list column.

## Dependencies

* Removed the dependency on **frontier**, along with eight other packages.
  `Imports` went from 25 packages to 15.

* Lowered the R requirement from `R (>= 4.4.0)` to `R (>= 4.0.0)`.

* Added a `testthat` suite covering the model branches, the numerical helpers,
  the S3 methods and the data generators.


# sfa 1.0.4

* Version released on CRAN, 2026-01-21.
