# sfa 1.2.0

A feature release. In brief:

* **A Wald test for whether the endogeneity correction was needed**,
  `endogeneity_test()`, from Hou, Ramalho and Roseta-Palma (2025). On an
  `ivsfm()` fit it tests `H0: rho = 0` -- that the noise is uncorrelated with
  the reduced-form errors, so a plain `sfm()` fit would have been consistent.
  The null is INTERIOR, since rho lives in the open unit ball, so this is an
  ordinary chi-square and not the chi-bar-square mixture
  `inefficiency_test()` needs. Size 0.046 at a nominal 5% over 1000
  replications; power 1.000 at rho = 0.4 and 0.8.

* **`ivsfm()` now reports standard errors for `rho`.** They were `NA`, on the
  ground that `rho = t/sqrt(1 + t't)` has a non-diagonal Jacobian and a
  delta-method value computed as though it were diagonal would be wrong. The
  first half was right and the conclusion was not: the Jacobian is
  `J = (I - rho rho')/s`, so `Var(rho) = J Var(t) J'` exactly. Checked against
  a numerical Jacobian, and validated end-to-end by the size of the test above
  -- a wrong Jacobian would show up there as a mis-sized test. The full
  covariance travels with the fit as `$vcov_rho`, which the joint test needs.

* **`ttsfm(model_name = "TTHN")` no longer fabricates a log-density it cannot
  compute.** The likelihood needs `D = Phi2(x1, 0; rho1) - Phi2(x2, 0; rho2)`,
  a difference of two bivariate normal CDFs of the same order of magnitude, and
  in parts of the parameter space it cancels completely: at
  `sigma_v = 0.3, sigma_u = 1, sigma_w = 0.2` on a 400-observation draw one
  observation's `D` comes back as exactly 0 while both CDFs are O(0.1).

  `pmax(D, .Machine$double.xmin)` did not repair that -- it invented
  `log(D) = -708.4` for the observation, and that single observation moved the
  summed objective by **715 log-units** in a surface whose real curvature is a
  few units per 0.1 step in log sigma. A cliff that size is what L-BFGS-B
  reports as `ABNORMAL_TERMINATION_IN_LNSRCH`. `D` is now rejected when it
  falls below the rounding error of its own subtraction, which leaves genuine
  small probabilities alone: a `D` of 6e-94 where both CDFs are also ~1e-94 has
  full significance and is kept.

  Three smaller repairs in the same branch: the exponentials for
  `sigma_v`/`sigma_u`/`sigma_w` are clipped at `EXP_CLIP_UPPER` as the `TTNE`
  branch already clipped its own; `rho` is held strictly inside `(-1, 1)`,
  which `sigma_v` small enough does reach; and the substitution of
  `-sqrt(.Machine$double.xmax / n)` -- about `-9.5e152` at n = 200, 140 orders
  of magnitude past the finite penalty it sat beside -- is gone, replaced by
  that penalty. It had also mapped a `+Inf` log-density to a large NEGATIVE
  one, so an unusable draw could be rewarded.

  Measured on eight matched fits: one that previously failed outright now
  converges, and one that previously returned `sigma_w = 0.044` against a truth
  of 1 now returns 1.047 **with a better likelihood**. The other six are
  unchanged. `sigma_v` still collapses toward zero on some draws; that is a
  separate problem and is not fixed.

* **`ttsfm()` no longer throws away a fit because of the optimizer's exit
  code.** It was the only entry point calling `stop()` on ANY non-zero
  convergence code, including L-BFGS-B's code 52,
  `ABNORMAL_TERMINATION_IN_LNSRCH`, which it returns whenever its line search
  meets a discontinuity and usually *having improved* the objective. The guard
  was written when nothing checked `opt$convergence` at all; `opt.optim()` now
  rebuilds at the stage-2 point whenever the value or the Hessian is
  non-finite, so what was left of it only destroyed usable fits. Two of eight
  TTHN test fits returned code 52 with sound estimates where the `stop()` had
  aborted them.

* **What skewness and kurtosis a pair of distributions can actually produce**,
  `moment_range()`, from Papadopoulos and Parmeter (2021). `spec_test()` already
  shipped their formal statistic; this is the range check that comes before it.
  Because the variance share `R` lies in `(0,1)`, the composed error can never
  be more skewed than `u` itself -- Normal-Half-Normal is confined to
  `gamma1 in (-0.995, 0)`, Normal-Exponential to `(-2, 0)` -- so residual
  moments outside that set refute the pair with no test statistic involved.
  That is their Type II failure, alongside the Type I failure of wrong skew.

  Five noise and five inefficiency distributions, including the two whose
  shape is a free parameter: the truncated normal sweeps
  `gamma1(u) in (0, 2)` and `gamma2(u) in (-0.243, 6)` and is the only `u`
  here whose composed error can be platykurtic under normal noise, while the
  gamma is unbounded in both and so can **never** suffer a Type II failure.
  The kurtosis range is not the interval between the two components: with
  leptokurtic noise it dips below both, to `0.674` for Laplace noise with a
  half-normal `u` against the half-normal's own `0.869`.

  **It is a diagnostic, not a test, and the print method says so.** PP2021's
  Tables 1-5 replicate -- their population rows exactly, and Tables 2, 3 and 5
  cell for cell -- and they show the range check over-rejects: with correctly
  specified data at `n = 200` and SNR 3, the sample excess kurtosis leaves its
  theoretical range 26% of the time. Use `spec_test()` for the p-value.

* **A frontier that tolerates the wrong skewness**, `esfm()`, from Hafner,
  Manner and Simar (2018). When the sample skewness of the residuals comes out
  positive -- a small-sample accident, not evidence against the model -- the
  classical MLE collapses: `sigma_u` is exactly zero, every firm is fully
  efficient and the fit says nothing. This model keeps the frontier and widens
  the distribution of `u`, with one parameter `gamma` whose SIZE is the scale
  of inefficiency and whose SIGN chooses the direction of skewness. The
  classical model is nested at `gamma > 0`.

  On a 50-observation wrongly skewed draw, `sfm(model_name = "NHN")` returns
  `sigma_u = 0.0003` and mean efficiency 0.99976; `esfm()` returns
  `gamma = -0.57` with the slopes still near truth.

  Because `gamma = 0` is an **interior** point rather than a boundary, the LR
  test of "no inefficiency" is an ordinary `chi2(1)` -- not the chi-bar-square
  mixture `inefficiency_test()` needs for the classical model. Size over 2000
  replications on the paper's design, at a nominal 5%: 0.178 / 0.082 / 0.054 /
  0.061 / 0.058 for n = 50 / 100 / 200 / 500 / 1000, against the paper's
  0.152 / 0.076 / 0.061 / 0.054 / 0.053. Badly oversized below n = 100, usable
  from n = 200 -- their own conclusion.

  **It is not a free improvement everywhere.** The paper reports smaller bias
  in mean efficiency than the classical model even when the population
  skewness has the correct sign, and that reproduces when inefficiency is not
  small: at `gamma = 0.5` the bias is -0.0001 at n = 100 against the classical
  +0.021. But at `gamma = 0.3`, where 22% of samples are wrongly skewed, this
  model overshoots the other way (-0.040 against +0.036) and is the worse of
  the two at n = 50 (-0.080 against +0.044). A small true `gamma` lets it fit
  sizeable inefficiency of either sign to noise. Use it when the wrong
  skewness is the problem.

* **Kim and Schmidt (2008) two-step test**, `uhet_test()`. The procedure
  applied work actually runs -- fit the frontier, regress the JLMS predictor
  on firm characteristics, test the slopes -- is invalid as usually performed,
  because `u-hat` is a *generated* dependent variable and the first-step
  estimation error enters the second-step variance. The correction is not
  distribution-free, so `"NHN"` and `"NE"` are implemented, the two cases the
  paper works through.

  What it buys: with `z` independent of `x` the corrected and naive tests
  agree, as the theory requires. With `z` strongly correlated with `x`
  (rho = 0.9) the naive standard error is **52% too large** at n = 800 and
  the naive test's size collapses to **0.000** at a nominal 5%, while the
  corrected standard error matches the true sampling spread to within half a
  percent (ratio 0.995) and keeps size 0.032.

  **Caveat, stated plainly:** the corrected test is conservative in these
  simulations. On the paper's own base case it rejects 0.024/0.036/0.028 at
  N = 200/500/1000 against a nominal 0.05, where the paper reports 0.0495 at
  N = 200. The *naive* branch does reproduce the paper's figure for the
  uncorrected test (0.023-0.029 against their 0.0237). The implementation was
  checked against equations (5)-(13), the first-step scores agree with an
  analytic half-normal score to 1e-11, `gamma-hat` is unbiased and Gaussian,
  and an oracle test using the true sampling standard deviation has correct
  size (0.045) -- so the gap is in the variance estimate, not in `gamma-hat`
  or the reference distribution. Two candidate explanations were tested and
  **both refuted**: the homoskedastic second-step variance (the paper's "OLS"
  against its "OLS-H") changes nothing, and neither does using the Hessian
  information rather than the outer product when forming `r_i`. The gap is
  unexplained rather than merely unverified. Read a rejection as meaningful
  and a non-rejection as weak evidence.

  Fits on the wrong-skew boundary are refused rather than answered: every
  `E[u|eps]` is then identical, so there is no estimated inefficiency for `z`
  to explain, and the information matrix is singular there (Waldman 1982).

* **Chen and Wang (2012) centered-residuals moment test**, `cw_test()`. Unlike
  `gof_test()`, it never evaluates the composed-error density or distribution
  function: it compares the empirical characteristic function of the *centred*
  residuals with its theoretical value, which follows from the characteristic
  functions of the two components separately -- in closed form for many pairs
  whose convolution is not, the normal-gamma model being the standard example.

  Centring is what makes it work. A frontier's intercept is not identified
  separately from `E[u]`, so the least-squares intercept estimates
  `alpha - E[u]`; centring cancels it, which is why the test is valid off an
  ordinary regression rather than requiring maximum likelihood. The variance
  corrects for both nuisances -- estimating the scales, and the centring
  itself. Omitting the second understates it and oversizes the test.

  Size at a nominal 5% over 500 replications, half-normal data:

  | test | n | tau=1 | tau=1.5 | tau=2 |
  |---|---|---|---|---|
  | cosine | 500 | 0.070 | 0.054 | 0.050 |
  | cosine | 2000 | 0.046 | 0.044 | 0.050 |
  | sine | 500 | 0.092 | 0.064 | 0.042 |
  | sine | 2000 | 0.068 | 0.062 | 0.064 |

  Power against exponential inefficiency at `tau = 1` is 0.231 (n = 500) and
  0.952 (n = 2000) for the cosine test against 0.108 and 0.706 for the sine
  test. Both of the authors' conclusions reproduce: the cosine test is the
  less sensitive to `tau` and the more powerful, so it is the default.
  **`tau` defaults to a single frequency**: passing several at once gives
  rejection rates of 0.09-0.15 at a nominal 5%, because the cosine moments at
  nearby frequencies are near-collinear and the covariance matrix is then
  near-singular. The authors report the same for their combined `tau`.

  A wrongly signed third moment leaves the moment equations with no admissible
  solution, and that is an error rather than a number.

* **The Gauss-Legendre rule is now cached.** It depends on nothing but the
  node count and was rebuilt by an eigendecomposition on every call --
  hundreds of times per `cw_test()`, since `numDeriv` differences a
  quadrature-based moment function. Caching also keeps `eigen()` out of forked
  workers: macOS's threaded BLAS is not fork-safe and `mclapply()` segfaulted
  inside it before this was added.

* **Goodness-of-fit tests for the assumed inefficiency distribution**,
  `gof_test()`, from Wang, Amsler and Schmidt (2011). The distribution of `u`
  is the assumption in this model least often defended and least often tested,
  and it is testable: hold the normality of `v` fixed and the assumed `u`
  implies a distribution for the composed error, so rejecting that is
  rejecting the assumed `u`. Kolmogorov-Smirnov and Pearson chi-square, both
  on the composed error rather than on `u-hat` -- which is the same test,
  since `u-hat = E[u|eps]` is monotonic in `eps`, but far easier, and it
  avoids the trap the paper warns about: comparing the spread of `u-hat`
  with the assumed density of `u` is a mistake, not a diagnostic, because
  those are different distributions.

  Parameter estimation is not ignored. The default parametric bootstrap copies
  the estimation step exactly -- each replication draws from the fitted model,
  rebuilds the response, **refits**, forms its own residuals and recomputes the
  statistic. The `"asymptotic"` option gives a chi-square p-value that the
  authors note is conservative at the MLE, and **returns `NA` rather than a
  KS p-value**, because with parameters estimated the Kolmogorov distribution
  does not apply and Bai's (2003) transformation is not implemented.

  Size and power over 200 replications (half-normal data for size, exponential
  data fitted as half-normal for power), at nominal 0.10 / 0.05:

  | | n | KS | chi2 |
  |---|---|---|---|
  | size | 200 | 0.065 / 0.025 | 0.085 / 0.065 |
  | size | 800 | 0.070 / 0.040 | 0.065 / 0.010 |
  | power | 200 | 0.505 / 0.395 | 0.235 / 0.145 |
  | power | 800 | 0.985 / 0.965 | 0.810 / 0.655 |
  | power | 2000 | 1.000 / 1.000 | 0.990 / 0.985 |

  Both hold their size, mildly conservatively, and **KS dominates chi-square
  at every sample size** -- the paper's own conclusion, and the reason it is
  the better default.

* **Coelli (1995)'s hypothesis tests, and his COLS standard errors.** New
  `inefficiency_test()` tests H0: no technical inefficiency, reporting the
  one-sided LR test, the naive LR test, the Wald ratio and the third-moment
  test side by side. The point of the paper is that **the two tests most often
  reported are the two with the wrong size**: H0 puts gamma on the boundary, so
  the LR statistic follows the Gourieroux-Holly-Monfort mixture
  `0.5 chi2(0) + 0.5 chi2(1)` (5% critical value 2.71, not 3.84) and the Wald
  ratio is not asymptotically normal at all. Measured here over 400
  replications per cell at a nominal 5%, with sigma_u = 0:

  | n | LR one-sided | LR naive | Wald | M3T |
  |---|---|---|---|---|
  | 100 | 0.060 | 0.030 | 0.230 | 0.055 |
  | 400 | 0.055 | 0.033 | 0.217 | 0.045 |
  | 800 | 0.060 | 0.035 | 0.182 | 0.062 |

  The one-sided LR and third-moment tests hold their size; the naive LR test
  rejects about half as often as it should; the Wald test rejects three to four
  times too often and does not improve with n. Under the alternative the
  one-sided LR test has the better power of the two valid tests at every sample
  size, which is Coelli's recommendation. `NTN` fits get the two-restriction
  mixture `0.25/0.5/0.25` over `chi2(0,1,2)`.

* **`sfm(estimator = "cols")` now reports Coelli's analytic standard errors.**
  It previously returned `NA` for `sigv` and `sigu` unless a bootstrap was
  requested, because the OLS errors are not appropriate for them: those
  parameters are non-linear functions of the residual moments `m2` and `m3`,
  not regression coefficients. Coelli's Appendix 1 (A15-A16) gives the correct
  delta-method variance and it is now filled in for `NHN`; the intercept
  additionally carries the error in its `E[u]` shift. `NE` and `NG` keep `NA`
  and the bootstrap, because the moment inversion -- and therefore the
  variance -- is distribution-specific.

  Checked against the nonparametric bootstrap that was already there, on the
  same data: the analytic/bootstrap ratio for `sigu` is 1.047 at n = 200,
  1.011 at n = 800 and 0.984 at n = 3200, which is an asymptotic
  approximation behaving as one should.

  Coelli's A13 inversion turned out to be the moment inversion the package
  already carried as `.gtre_two_step()` (equal to ten decimal places, not
  merely similar), and his A15 variance is the delta method whose GTRE_SEQ2
  implementation is corrected elsewhere in this release -- **including the
  factor of 2 on the covariance term, which A15 independently confirms**, and
  A16's `mu6 - mu3^2 - 6 mu2 mu4 + 9 mu2^3`, which confirms the other two
  repairs. The central moments now live in one shared helper,
  `.nhn_central_moments()`, rather than being derived twice.

* **A composed-error CDF for every cross-sectional model**, `pcomposed_model()`
  and `composed_cdf()`. `pcomposed()` covers the half-normal only, which is
  enough for copulas but not for goodness-of-fit testing, where the CDF is
  needed for whichever inefficiency distribution the model assumes. All
  thirteen are covered: NHN, NTN, NE, NR, NU, NGE, NLN, NW, NG, NNAK, TSL,
  THT and tHN.

  Three parameterizations do not read the way their names suggest and are
  taken from each likelihood rather than from the label: **`THT` lists `sigu`
  before `sigv`**, the only model that does; **`NG`'s `sigu` is the gamma
  scale** and its `mu` the shape; **`NNAK`'s `sigu` is the Nakagami spread**,
  so `Omega = sigu^2`; and `NLN`'s `mu` is a meanlog.

  **`THT` is not an independent convolution.** Tancredi's composed error is
  skew-t, a scale mixture in which `v` and `u` are divided by the *same*
  `sqrt(V/a)`. Treating it as t noise plus an independent half-normal -- which
  is exactly what `tHN` is -- gets the log density wrong by up to 1.74. The
  first draft did precisely that, and the check below is what caught it.

  Verified three ways: against `pcomposed()` for the half-normal, agreeing to
  3.5e-15 relative including a lower tail at `log F = -209`; against a
  four-million-draw simulation for all thirteen models, maximum absolute error
  4e-4 against a Monte Carlo standard error of 7.5e-4; and by numerically
  differentiating the CDF and comparing with each model's own stored
  log-density, which agrees to about 1e-8 for eleven of them. `NLN` and `NW`
  differ by ~3e-3 because their likelihoods are simulated, not because the
  models disagree.


* **`psfm()` died with "system is exactly singular" when a factor was only
  partly collinear between individuals.** Reported from use: `factor(year)` on
  an unbalanced panel, where the between-individual means of the year dummies
  are rank deficient. The guard added in 1.1.5 detected this correctly and
  printed an accurate warning, then did nothing, because it could only express
  the fix by deleting whole formula **terms** -- and 11 of `factor(year)`'s 24
  dummies being collinear names no term that can be deleted without discarding
  the 13 identified ones too. `plm::ercomp()` was then handed the same singular
  design and failed inside `solve()`. The starting-value design is now reduced
  at **column** granularity (`.re_start_design()`), so the identified dummies
  survive, the dropped ones take their starting values from a pooled OLS fit as
  documented, and the likelihood still uses the formula as written. This was a
  no-op for every partly-collinear factor, which is the common case -- a fully
  collinear term was the only shape the old path could handle.

  Two consequences of that path are fixed with it. `collinear_action =
  "warn_drop"` silently dropped nothing in the same situation, and now says so
  instead. And `"GTRE_SEQ1"`/`"GTRE_SEQ2"` took the frontier standard errors
  from the starting-value fit by **position**, which errored on the length
  mismatch once the design was reduced; they are matched by name, with `NA` for
  columns that fit could not identify rather than a neighbour's value.

* **`"GTRE_SEQ1"` and `"GTRE_SEQ2"` are inconsistent for fixed `T`, and now say
  so.** Both hand `plm`'s random-effects output to a second stage that treats it
  as draws from the latent composite errors. It is not: `alpha_hat` is a
  shrunken BLUP and `eps_hat` a quasi-demeaned residual, which attenuates the
  third central moment by about 0.81 at \eqn{T = 10} on the transient side and
  scales the persistent side by \eqn{c^2} and \eqn{c^3}. No factor tends to 1 as
  \eqn{N} grows at fixed \eqn{T}, so both converge to the wrong constants --
  \eqn{\sigma^2_{hr}} settles near 0.68 against a true 0.80. Feeding the same two
  second stages the latent draws instead recovers the truth at the \eqn{\sqrt{N}}
  rate, which puts the fault in the shared first stage rather than in the
  likelihood or the moment inversion. Each call now warns and points at
  `model_name = "GTRE"`. This is a documentation and diagnosis change: the
  estimators are unchanged, and are the large-T procedures they always were.

* **`"GTRE_SEQ2"`'s standard errors were wrong in four ways.** The variance of
  the third central moment used `mu_3^3` for \eqn{\mu_3^2} and `9 * mu_2` for
  \eqn{9\mu_2^3} (also dimensionally inconsistent); the delta-method
  cross-terms were missing their factor of 2; `beta_0`'s error added a standard
  error to a sum of variances instead of its square; and both
  \eqn{\partial\gamma} derivatives were garbled, squaring neither denominator and
  carrying \eqn{m_3} where \eqn{m_2} belongs. Against the Monte Carlo sampling
  standard deviation the reported error for `gamma_uv` was **16 times too
  large** and `sigmaSq_uv` 1.37 times; corrected, both land within 2%. The
  algebra now lives in `.gtre_two_step_se()` so it is unit-tested, and returns
  `NA` at the wrong-skew boundary where it is not differentiable rather than
  `Inf`. `"GTRE_SEQ2"` also reported a `lambda` computed from a *different*
  (likelihood) decomposition than the sigmas printed beside it; it now comes
  from its own.

* **`ttsfm()` returned an infinite penalty on unusable parameter draws, in
  seven places.** `TTNE`, `TTHN` and `TTNLS` each returned
  `.Machine$double.xmax` when a draw produced non-finite contributions. That is
  the construct `sfm()` was converted away from earlier in this release after
  `"NGE"` died on 3 of 45 fits: `optim()` differences the objective to get its
  gradient, and differencing 1.8e308 overflows to a non-finite value, which does
  not steer the search away -- it aborts the fit. `TTHN` was failing at
  \eqn{n = 400} with `ABNORMAL_TERMINATION_IN_LNSRCH`, which is what L-BFGS-B
  does when handed gradient information it cannot use. All seven now return a
  finite penalty.

  The regression test that exists to catch exactly this construct only ever read
  `sfm.R`, which is how `ttsfm.R` kept all seven for as long as the test
  existed. It now reads every file carrying a likelihood closure.

* **`ttsfm()` no longer calls the bivariate normal CDF once per observation.**
  The `TTHN` density was evaluated with `mapply()` over observations, so one
  `mnormt::pmnorm()` call per observation, twice per likelihood evaluation --
  and with a numerically differenced gradient, \eqn{2n(2p+1)} calls per
  optimiser iteration. The correlation argument varies only through
  \eqn{\sigma_u} and \eqn{\sigma_w}, so under homoskedasticity it is one value
  repeated and the whole vector is a single call. Grouped accordingly, with a
  fallback to the per-observation form when the correlation genuinely varies
  per observation -- there, grouping is *slower*, measured at 0.60x. Results are
  bit-identical to the previous form. The gain is about 1.4x, not more:
  `pmnorm()` still loops internally, so the call count falls but the work does
  not. `TTNLS` separately computed each of its two bivariate terms twice; those
  are now computed once.

* **`selsfm()` gains a second model, `model_name = "kts"`**: Kumbhakar, Tsionas
  and Sipilainen (2009), joint estimation of technology choice and technical
  efficiency. Not a variant of Greene (2010), which this entry point already
  fitted, but a different model. Greene has one frontier observed only for the
  selected, with selection correlated with the *noise*. Here both technologies
  are observed, each with its own frontier and its own two scales, and the
  choice depends on **inefficiency itself** -- `u` lowers output through the
  frontier and shifts the adoption decision through the choice equation.
  Neither two-step order works, which is the paper's own argument for
  single-step maximum likelihood: the choice equation cannot be a probit
  because `u` is unobserved, and fitting the frontiers first ignores the
  endogeneity of the choice. Because `u`'s distribution depends on the regime
  while the regime's probability depends on `u`, the marginal density of `u` is
  a mixture of the two half-normals and the observation's density is the same
  mixture of two one-dimensional integrals. `n_nodes` controls the quadrature.

* **`copsfm()` goes from two copula families to fifteen**: Frank, Clayton,
  Gumbel and Joe join Gaussian and FGM, each of the last three with 90/180/270
  rotations, since Clayton, Gumbel and Joe carry only positive dependence and
  nothing rules out a negative association between noise and inefficiency.
  Every density is verified against the second mixed partial of its own CDF, it
  integrates to 1 over the unit square, and it is exactly 1 at the independence
  parameter.

  **`copsfm()` now warns when you choose a family that does not recover its own
  dependence parameter**, and most of them do not. Fitting each family to 25
  samples generated from itself at n = 400, the share of fits returning the
  independence boundary was 0% for Frank and Clayton but 36% for Gumbel, 40%
  for Joe, 56% for Clayton rotated 270 and 60% for Gumbel rotated 90. On data
  from a Gumbel copula with Spearman 0.685 at n = 2000, *every* family
  including the true one returns the boundary and their log-likelihoods differ
  by less than 0.04: the likelihood is close to flat in the dependence
  parameter. That is a property of the model, not of the implementation. The
  warning quotes the measured rate; the rotations that were not measured say
  so. Prefer `"frank"` or `"clayton"`.

* **`zsfm(model_name = "ZISF")` could be trapped at the boundary where the mass
  point vanishes, returning a much worse fit without saying so.** Since
  `P(fully efficient) = exp(-gamma)`, `gamma -> Inf` is the boundary at which
  the model collapses to an ordinary frontier, and it carries its own local
  optimum. On the package's convergence design the default start fell into it
  on roughly **1 replication in 5000**: one fit returned `gamma = 9.93`
  (`p = 5e-5`) with a log-likelihood **188.7 points below** what a start from
  anywhere in the interior reaches. Nothing in the returned object flagged it.

  `zsfm()` now tries several starting values for `gamma` and keeps the best
  attained objective, so it can only improve a fit. The guard applies to
  `"ZISF"` only: `"ZISF_Z"` has no scalar `gamma` -- its layout is
  `(sigv, sigu, beta, z-block)` with the mixing parameterised through the `z`
  coefficients -- so the same perturbation there would move the noise scale.

  Worth knowing if you have compared this model's convergence behaviour before:
  that single fit inflated its sample size's mean squared error 144-fold and
  flattened the log-MSE-on-log-n slope for `gamma` to zero, while four of the
  five sample sizes were textbook (MSE 0.0047, 0.0019, 0.00095, --, 0.00065).

* **The `sigma_u -> 0` repair listed further down opened a second catastrophic
  cancellation, at `sigma_v -> 0`, which is also fixed here.** To be clear about
  the scope, since it affects whether anyone need re-check old results: this
  defect was introduced by that repair and existed only in this release's own
  development. **No published version of the package has it** -- 1.1.5 and
  earlier have the `sigma_u` defect instead, and `.log_phi_tilt()` does not
  exist in them. It is recorded because the mechanism is worth knowing and
  because the fix changes NE's numerics.

  That repair works through `.log_phi_tilt()`, which returns `log Phi(z) + z^2/2`. For an observation
  with `eps < 0` and `sigma_v` small, `z` goes large *positive*, `log Phi(z)`
  is ~0, and so the returned value is essentially `z^2/2` -- which the NE
  branch then subtracted `eps^2/(2 sigma_v^2)` from. At `sigma_v = 1e-8` both
  are about 5e15, where consecutive doubles are ~1 apart, so the result was
  rounding noise: the summed negative log-likelihood read **-5188 where the
  true value is +4398**, an apparently excellent fit in a narrow band that a
  minimiser walks straight into. The domain guard did not catch it, firing only
  at `sigma_v <= 0` while the band sits at 1e-9 to 1e-7.

  The consequence was not a visible failure. `nlminb` found the band, its
  "only accept an improvement" guard passed the point through because the
  objective had improved, and `bobyqa` and `optim` could not escape the corner;
  `optim`'s bounds, rebuilt around the corrupted point, then pinned the
  frontier coefficients. Such fits came back as ordinary `sfareg` objects with
  `convergence = 0` -- correctly, since a bound-constrained stationary point is
  a legitimate success -- while carrying a gradient in the hundreds.

  The subtraction is analytic: `z^2/2 = eps^2/(2 sv^2) + eps/su + sv^2/(2 su^2)`
  identically, so above `NE_TILT_SWITCH` the expression collapses to the 1.1.5
  tilt with no large intermediate formed. The two forms agree to 1e-13 across
  the switch, and below it the tilt form -- and hence the `sigma_u -> 0` fix --
  is untouched. Over 500 replications at `lambda = 3.33`, n = 1000, the rate of
  fits landing more than 0.3 from the truth falls from **3.8% to 1.0%**, the
  log-likelihood **improves on 18 and worsens on none**.

  The root-n convergence sweep was then re-run in full, 1000 replications at
  each of n = 1000 to 5000, and `"NE"` now **passes on the mean basis, which is
  what it used to fail**: slopes -1.023 (`sigma_v`), -1.044 (`sigma_u`), -1.025,
  -1.033 and -1.009 (the frontier coefficients), every R-squared at least 0.989,
  and no failures in 5000 fits. Before the fix those slopes were -0.488, **+0.262**
  -- positive, i.e. error growing with the sample -- -0.745, -0.443 and -0.165.

  This also resolves a dilemma the convergence notes had recorded and left open.
  They found that reverting the likelihood to its 1.1.5 form fixed the tail but
  would reintroduce the `sigma_u -> 0` defect, and concluded "do not revert, fix
  the interaction instead". Above the switch this fix *is* the 1.1.5 algebra, so
  the re-run reproduces that configuration's slopes digit for digit while
  keeping the `sigma_u` protection. Both boundaries are now covered.

  `"NGE"` shares `.log_phi_tilt()` and was checked for the same defect; it does
  not have it.

* **Robust estimation gains its missing half.** `sfm(robust = )` already fitted
  the Normal--Half-Normal frontier by maximum \eqn{L_q}-likelihood, the power-Psi
  criterion and minimum density-power divergence, but the tuning parameter had
  to be supplied by hand with no guidance and no diagnostic. Four new functions
  close that, from Bernstein, Parmeter and Wright (2026):

  * `hscore_select()` chooses the tuning parameter by minimising the Hyvarinen
    score (Sugasawa and Yonekura, 2021) over a coarse-to-fine grid that includes
    the maximum-likelihood endpoint, so the criterion is free to decline
    robustification. It guards two silent failure modes: candidates lost to
    numerical overflow, and warm-start capture in the degenerate
    `sigma_u -> 0` basin.
  * `hscore()` evaluates the criterion, in logarithms. The natural-scale
    expression raises the fitted density to the power `c - 2`, so one
    observation whose density underflows makes the score non-finite; because the
    surviving candidates are those closest to maximum likelihood, the failure
    biases selection toward too little robustness without any warning.
    `stable = FALSE` reproduces the naive evaluation for comparison.
  * `calibrate_c()` gives the fixed alternative: the tuning value at which a
    residual three scale units out retains a stated share of its clean-data
    influence. The influence ratio is **not** monotone for the
    Fisher-consistency-corrected criteria -- Psi and MDPD cross the target twice
    -- so every root is found and reported in a `"roots"` attribute rather than
    resolved silently.
  * `density_weights()` returns the per-observation weight the estimator
    applies, documented with what it can and cannot detect: it responds to
    response contamination and not to a mis-recorded regressor, which can be
    highly influential and still receive a weight near one.

  The power-Psi and MDPD objectives satisfy `obj_MDPD = (1 + c) * obj_Psi`
  identically, so they share a maximiser and are the same estimator at a common
  `c`. Both remain available as peer methods and `calibrate_c()` returns the
  same value for either.

* **`influence_sfa()`**, an influence-function diagnostic for any `"sfareg"`
  fit made with `keep_objective = TRUE`, after Stead, Wheat and Greene (2023).
  The standard outlier rules do not transfer to this model -- they ignore the
  asymmetry of the composed error, so a genuinely inefficient firm reads as an
  outlier and is discarded. It reports the empirical influence function, a
  Cook's-distance-shaped case influence, and **two** sensitivities. Compare
  models on `sensitivity_std`, the self-standardised one: the raw sup-norm
  depends on the parameterisation, and on a clean sample reads 72.6 for `"NHN"`
  against 1322.5 for `"tHN"`, which says nothing about robustness. On the
  standardised scale, contaminating one response takes `"NHN"` from 13.5 to 32.3
  and `"tHN"` from 13.1 to 17.2.

* **A fifth model-fitting entry point, `lcsfm()`** — the latent class
  stochastic frontier, with `n_class` coexisting technologies.

* **New accessors.** `efficiency()` returns technical efficiency with a choice
  of predictor (`"bc"`, `"jlms"`, `"mode"`) and of scale for the dependent
  variable (`logDepVar`); the mode is the only one that reaches exactly 1, which
  it does wherever the posterior mean of inefficiency is negative.
  `meanefficiency()` returns the model-implied `E[exp(-U)]` for ten
  distributions -- seven in closed form -- together with mean efficiency among
  the most efficient p, which is a property of the fitted model rather than of
  the sample. `simulation_se()` reports how much of a simulated-ML standard
  error is simulation noise rather than sampling noise.

* **`sfm()` gains four arguments.** `weights`/`wscale` for observation weights
  (a weight of 2 gives the same estimates as including the row twice);
  `start_from` to seed a hard model from a simpler fitted one, matched by
  parameter name; `scaling` for the Wang-Schmidt scaling-property model on
  `"NTN"`; and `shapehet`, an **experimental** covariate-dependent shape for
  `"NG"`/`"NNAK"` whose limitations are set out in `?sfm`.

* **`vcov(fit, type = "bhhh")`** returns the outer-product-of-gradients
  covariance, which is defined even where the Hessian is not.

* **Regression tables.** `texreg::screenreg()`, `texreg()` and `htmlreg()` now
  render `"sfareg"` fits directly.

* **Two fixes to things that were shipped but did not work.** The `sandwich`
  methods were defined without being registered, so `sandwich::estfun()`,
  `vcovCL()` and `lmtest::coeftest()` all failed on an installed package --
  every test passed because `devtools::load_all()` exports everything. And
  `sim_seed` silently did nothing under the default `sim_type = "halton"`,
  leaving the cross-sectional simulated-ML models with unrandomized draws;
  they are now shifted like the panel models. Fixing the second is what makes
  `simulation_se()` meaningful.

* **`"NNAK"` no longer fails outright at larger shape values.**
  `gsl::hyperg_U` throws above order ~16 and the vectorized call meant one bad
  element nulled the whole vector. There is now a fallback that evaluates the
  parabolic cylinder function by its integral representation, agreeing with the
  route it replaces to 1e-12 where both work.

* **A seventh entry point, `ivsfm()`** — stochastic frontier estimation when a
  regressor is endogenous, in the sense of being correlated with the noise.
  Three estimators from Amsler, Prokhorov and Schmidt (2016): `"IVLIML"` (full
  maximum likelihood over the frontier and the reduced form jointly), `"IVCF"`
  (the two-step control function of Kutlu 2010) and `"C2SLS"` (corrected 2SLS,
  which is moment-based and so reports no log-likelihood). Supplying `uhet`
  gives the model of Amsler, Prokhorov and Schmidt (2017), in which
  environmental variables entering the inefficiency scale may themselves be
  endogenous. Ignoring endogeneity is not a small matter: on simulated data
  with a correlation of 0.6 between the noise and a regressor's reduced-form
  error, an ordinary `sfm()` fit carries a bias of 0.133 in that regressor's
  coefficient where `ivsfm()` carries 0.002. Note that `"IVCF"`'s conventional
  standard errors are understated unless the correlation is zero, and that the
  case where the regressor is correlated with the *inefficiency* rather than
  the noise is not covered.

* **A sixth entry point, `selsfm()`** — Greene's (2010) stochastic frontier
  with a correction for sample selection, for the case where the units in the
  sample are there because of something correlated with the frontier's own
  noise. The selection and frontier equations are supplied as separate
  `selection` and `frontier` arguments. Estimation is Greene's two-step
  scheme: probit on the full sample, then maximum simulated likelihood on the
  selected subsample, with efficiency scores from the simulated Bayes rule
  over the same draws. Note that the standard errors are conditional on the
  first-stage estimate and do not carry the Murphy–Topel correction; bootstrap
  over both stages when inference on `rho` matters.

* **Six new panel estimators**, all classical rather than maximum likelihood
  and none assuming a distribution for inefficiency: `"CSS"`
  (Cornwell-Schmidt-Sickles 1990), `"LS"` (Lee-Schmidt 1993), `"KSS"`
  (Kneip-Sickles-Song 2012), and `"SSRE"`/`"SSCRE"`, which complete the
  Schmidt-Sickles family that `"SSFE"` already began.

* **Heteroskedasticity in more than one error component**, as named formula
  arguments to `sfm()`: `vhet`, `uhet` and `muhet`. The last, with
  `model_name = "NTN"`, is Battese-Coelli (1995).

* **Three corrections that change numeric output**: `"NE"`/`"NGE"` could
  return a positive log-likelihood with `sigma_u` at zero; `nobs()` counted
  rows supplied rather than rows used, so `BIC()` used the wrong `n`; and
  `"NGE"` aborted outright on about 7% of small samples.

* Plus `z_link` for `sfm()` and `ttsfm()`, `marginal_effects()` for panel
  fits, `efficiency_ci()`, and `sfm(model_name = "TSL")`.

Nothing in this release changes the meaning of an existing argument or
`model_name`. Every addition is a new argument defaulting to the previous
behaviour, or a new `model_name` value.


* **Heteroskedasticity in more than one component: new `vhet`, `uhet` and
  `muhet` arguments to `sfm()`.** `sigma_v`, `sigma_u` and (for `"NTN"`) the
  pre-truncation mean of `u` may now each be driven by their own covariates,
  in any combination:

  ```r
  sfm(y ~ x1 + x2 | z_u, vhet = ~ z_v, model_name = "NHN_Z")
  sfm(y ~ x1 + x2, muhet = ~ z_mu, model_name = "NTN")   # Battese-Coelli 1995
  ```

  Named formulas rather than further pipe segments, because pipe POSITION
  already means different things in different families -- the second segment
  is `sigma_w` in `ttsfm()` but `sigma_h` in `psfm()`'s `GTRE_Z` -- and a
  fourth position would not be readable. The existing `| z` keeps its meaning,
  so nothing that worked before changes. `uhet` is the same specification as
  `| z` and refuses if both are given; it exists so that `"NTN"`, which takes
  no pipe segment, can still have a heteroskedastic `sigma_u`.

  Available for `NHN`, `NHN_Z`, `NE`, `NE_Z` and `NTN`. `muhet` requires
  `NTN`, the only family in which a pre-truncation mean exists. Recovery on
  simulated data with all three blocks non-constant, n = 1500: every one of
  the six parameters within one standard error of truth, and the two z-links
  agree to 1e-6 in log-likelihood with deltas differing by exactly the factor
  of two the reparameterization implies.

  Two reasons this matters beyond completeness. Ignored heteroskedasticity in
  `v` biases the estimated FRONTIER, not just its precision (Caudill, Ford and
  Gropper 1995), because the composed error's mean depends on the scales. And
  `vhet` is one of the few specifications that can absorb apparent wrong
  skewness without forcing `sigma_u` onto the zero boundary.

  These fits carry every positive quantity on the log scale, so the parameter
  vector is unconstrained and no variance can be driven onto a boundary by the
  optimizer. `marginal_effects()` reads them with no extra arguments.
  `estimator = "cols"` and the `robust` divergences are moment-based and
  homoskedastic respectively; both refuse rather than silently ignoring the
  new arguments.

* **`print()` and `summary()` now say when a variance scale is on the zero
  boundary**, for `sfm()`'s `sigma_u` under wrong skewness, for `"tHN"`'s
  `sigma_u` absorbed by heavy-tailed noise, and for `psfm()`'s two persistent
  scales. Each says which of those it is, since they mean different things.

  The condition was already warned about at fit time, but a warning is a
  one-off: it does not survive `saveRDS()`, a fresh session, or a loop that
  suppressed warnings to keep the console readable. What people look at
  afterwards is the printed table, and there a collapsed scale appears as an
  ordinary coefficient of `0.0009` with nothing to say that it is a boundary
  solution, that the split it belongs to is unidentified in this sample, or
  that the intercept moved to absorb the difference.

* **`psfm_bootstrap()` warns when asked to bootstrap a boundary fit.** A
  parametric bootstrap resamples from the *fitted* model, so a fit whose
  persistent scale has collapsed resamples from a data-generating process in
  which that component is absent -- and the bootstrap is inconsistent on the
  boundary of the parameter space in any case. The intervals would understate
  the uncertainty in the persistent split rather than represent it, while
  looking perfectly ordinary.

  Since one of GTRE's two persistent scales collapses in roughly a third of
  samples, this is expected to fire in ordinary use. That is the point.

* **`psfm(model_name = "GTRE")` now reports a persistent scale that has
  collapsed to zero**, in `$sigh_at_bound` and `$sigr_at_bound` and in a
  warning naming the surviving scale. Under **both** estimators: `"fiml"`,
  which is the default and packs its results in a separate branch, and
  `"sml"`.

  `GTRE` has two persistent components and the likelihood cannot always
  separate them in a given sample. When it cannot it merges them: one goes to
  zero, the other absorbs its variation and comes back inflated. On simulated
  data with both genuinely present, **one of the two collapsed in 37% of 87
  replications**, and `sigr -> 0` was three to five times commoner than
  `sigh -> 0`.

  **Whether firms or periods help depends on the estimator, and the two
  answers are opposite.** Crossed (N, T) designs on the same DGP, regressing
  the collapse rate on log N and log T:

  | estimator | fits | log N | log T |
  |---|---|---|---|
  | `"fiml"` (the default) | 600 | **-1.66** (0.18) | -0.93 (0.20) |
  | `"sml"` | 251 | -0.68 (0.26) | **-1.07** (0.27) |

  Under the default, firms buy about 1.8 times what periods do; under
  simulated ML it is the other way round. Restricting the `"fiml"` design to
  N <= 200 -- cell for cell the grid the `"sml"` one used -- leaves it at
  -1.64 (0.22), so this is the estimator and not the design; the difference in
  the N coefficient is significant (z = 2.8) and the difference in the T
  coefficient is not.

  So the practical advice under the default is **add firms**: at T = 4 the
  collapse rate falls 68%, 36%, 18%, 6% across N = 50, 100, 200, 400 and
  reaches 0% at N = 400, T = 16. Under `"sml"` that same row is flat (60%,
  39%, 41%), which is why a short panel is a dead end there however many firms
  it has. `?psfm` now gives both.

  **A user-facing consequence: you can tell in advance.** The quantity that
  governs whether a panel supports the split is not `sigma_h/sigma_r` -- that
  only decides which scale dies -- but the persistent variance in the firm
  means against the transient variance still in them,
  `[sigr^2 + sigh^2 (1-2/pi)] / ([sigu^2 (1-2/pi) + sigv^2] / T)`. Across six
  designs varying all four scales it orders the collapse rate perfectly
  (Spearman -1 under both estimators; logistic coefficient on log(S/N) of
  -1.15, SE 0.20): 54% of fits merge at S/N = 0.67, 31% at 1.30, 6% at 11.7.
  Below about 1, expect half your fits to merge; it takes roughly S/N > 5 to
  get under 10%. `?psfm` and the vignette give the formula and a helper.

  One limit, stated because it is easy to over-read: those elasticities were
  measured at `sigma_r = 0.20` against `sigma_h = 0.40`. On a design where the
  random effect dominates instead (0.50 against 0.40, N = 100, T = 6) the
  *level* ordering reverses -- `"fiml"` collapses on 33% of replications and
  `"sml"` on 20%. Which estimator handles the split better depends on how much
  persistent inefficiency there is to find, not on the estimator alone.

  This is usually the CORRECT maximum likelihood estimate, not a failure: on
  one such replication the boundary solution beat the true parameter vector by
  3.66 log units. It is the panel counterpart of the cross-sectional
  wrong-skewness result (Waldman 1982). So it is reported rather than
  prevented -- bounding either scale away from zero would corrupt exactly the
  samples where the boundary is the answer.

  **Your intercept moved with it**, and the warning now says so. One rule
  covers what a collapse does: it *conserves the persistent variance and
  relocates the persistent mean into the intercept*. Those are two different
  destinations -- `h` is half-normal, carrying both a mean
  `E[h] = sigh*sqrt(2/pi)` and a variance `Var(h) = sigh^2 (1 - 2/pi)`, and
  the model has a separate place for each: the intercept takes the mean, since
  `E[y] = beta0 - E[h] - E[u]`, and the surviving persistent scale takes the
  variance. So the survivor does *not* come back at `sqrt(sigh^2 + sigr^2)`;
  that shortcut predicts 0.447 where the answer is 0.313.

  | | surviving scale | intercept |
  |---|---|---|
  | `sigh -> 0` | `sqrt(sigr^2 + sigh^2(1-2/pi))` = 0.3133 | `-E[h]` = -0.3192 |
  | `sigr -> 0` | `sqrt(sigh^2 + sigr^2/(1-2/pi))` = 0.5197 | +0.0955 |

  Measured over 1000 replications at each of five sample sizes: 0.2973,
  -0.3035, 0.5058, +0.0931 -- within 2.5% to 5.1% with nothing fitted, and all
  four slightly below prediction, which is what selection produces rather than
  a wrong constant. The practical consequence is that a collapsed fit does not
  only mislabel the persistent split; it displaces any level read off `coef()`
  by a known and substantial amount.

  **And it is wrong skewness, measured between firms.** Average the OLS
  residuals within each firm -- that average is the persistent part of the
  error -- and look at their skew. Over 1000 replications the samples whose
  between-firm skew came out the wrong way for a production frontier lost
  `sigma_h` to the boundary **69%** of the time against 3.6% for the rest,
  with an ROC area of **0.955** for that one statistic, and not one of them
  collapsed `sigma_r` instead. So a collapsed `sigma_h` is not an accident to
  retry with different starting values: the sample carries no between-firm
  evidence of persistent inefficiency and the estimate says so.

  Read the sign, not a p-value. The firm mean of residuals contains
  `-ubar_i` as well as `-h_i`, and the mean of T half-normals is skewed the
  same way, so the between-firm skew leans negative even at `sigma_h = 0` --
  82% of the time at T = 4, 63% at T = 10, 53% at T = 30. `skewness_test()`
  is therefore not applicable to between-firm residuals as it stands, and
  `?psfm` says so.

  `psfm()` also gains `keep_objective`, which retains the likelihood for the
  simulated-ML panel models so it can be evaluated away from the optimum. That
  is what distinguishes a weakly identified parameter from a badly estimated
  one, and it is how the above was established. `?psfm` has a worked example.

* **Robust and clustered standard errors, via `sandwich`.** `vcov()` returned
  the inverse Hessian and nothing else, which is valid only if the likelihood
  is correctly specified and the observations are independent. Neither is safe
  to assume in applied frontier work, and clustered errors are routine on
  firm-level panels.

  `bread()` and `estfun()` methods are now registered for `"sfareg"` when
  `sandwich` is installed, which is all that package needs:

  ```r
  fit <- sfm(y ~ x1 + x2, model_name = "NHN", data = d, keep_objective = TRUE)
  lmtest::coeftest(fit, vcov. = sandwich::vcovCL(fit, cluster = d$firm))
  ```

  Two requirements, both of which error clearly rather than returning
  something wrong. The fit needs `keep_objective = TRUE`, since the scores are
  differenced from the retained likelihood; and the `robust` divergence
  estimators are excluded, because they do not maximise a log-likelihood and so
  have no score. `sandwich::vcovHC()` does not work and is not expected to --
  its corrections are built from hat values, which are undefined for a
  nonlinear likelihood whose parameters include variance components.

  The scores are central differences with a per-parameter step. Validated two
  ways at the optimum: against `NHN`'s hand-derived analytic gradient (agreeing
  to 3e-07) and against `numDeriv` (2.9e-07), with the column sums vanishing
  relative to their own scale, which is the identity that makes them scores.

  `sandwich` and `lmtest` are in `Suggests`; neither is required to use the
  package.

* **New function `skewness_test()`: a p-value for wrong skewness.** `sfm()`
  already detected the condition -- `$wrong_skew`, `$sigma_u_at_bound`, and a
  warning -- but could not say whether the skew was wrong by more than sampling
  noise. Wrong skew is the single most common reason an applied SFA fit is
  meaningless, so that gap mattered.

  Two tests: D'Agostino's (1970), which is the one Schmidt and Lin (1984)
  appeal to, and Coelli's (1995) `M3T`. Returns an ordinary `"htest"`, so it
  prints like `t.test()` and carries `m3`, `nobs` and `wrong_skew` alongside
  the usual fields. Takes a fitted `sfm()` object -- using the OLS residuals
  the fit recorded, which is what the tests are defined on, rather than the
  composed residuals -- or a bare numeric vector.

  `"agostino"` is the default because it holds its size where the asymptotic
  form does not. One-sided at a nominal 5%, 4,000 replications of symmetric
  residuals:

  | n | 25 | 50 | 100 | 400 |
  |---|---|---|---|---|
  | `coelli` | 0.031 | 0.037 | 0.043 | 0.050 |
  | `agostino` | 0.045 | 0.046 | 0.046 | 0.051 |

  Coelli's is conservative below about n = 100 and the two agree by n = 400.
  The D'Agostino implementation agrees with `moments::agostino.test` to within
  1e-10 at n = 30/60/200/1000.

  `sfm()` fits now also carry `$ols_residuals`, which is what makes the test
  exact rather than re-derived from the recorded call.

* **`lcsfm()` -- the latent class stochastic frontier** (Greene 2005; Orea and
  Kumbhakar 2004; Caudill 2003), a new exported entry point alongside `sfm()`,
  `psfm()`, `zsfm()` and `ttsfm()`. `J` distinct
  technologies coexist in one sample, each with its own frontier `beta_j` and
  its own `sigma_v,j` and `sigma_u,j`, and which one a firm operates is
  unobserved. Class membership follows a multinomial logit with class `J` as
  the reference, constant under `"LCM"` and parameterized by the second formula
  segment under `"LCM_Z"`. The number of classes is `n_class`, default 2.

  Every firm contributes to every class weighted by its class probability, so
  this is not the same as fitting separate frontiers to known subgroups -- the
  groups are not known, and the uncertainty about them is carried into the
  efficiency predictions. `post.prob` is the `n x J` matrix of posterior class
  probabilities, `jlms_class` the class-conditional JLMS scores, and `jlms`
  their posterior-weighted average.

  It is its own function rather than a `model_name` inside `zsfm()`. The two
  are related -- `zsfm()`'s `"ZISF"` is the restricted two-class case in which
  one class has no inefficiency at all -- but a user looking for a latent class
  model would not think to look inside a function named for zero inefficiency,
  and `n_class` means nothing to `ZISF`, whose two components are not free.

  Recovery on simulated two-class data, `n = 800`, true shares 0.60/0.40:
  `beta_1 = (1.00, 0.48, 0.51)` against a true `(1.0, 0.5, 0.5)`,
  `beta_2 = (2.97, 0.98, 0.21)` against `(3.0, 1.0, 0.2)`, scales
  `0.51`/`0.98` against `0.5`/`1.0`, fitted shares `0.598`/`0.402`, and 92% of
  firms assigned to their true class. On three classes at `n = 1200` the
  figure is 95.5% after matching labels.

  **Class labels are identified only up to permutation.** Starting values come
  from splitting the OLS residuals at their `J`-quantiles and refitting within
  each group, which makes a run reproducible on the same data and avoids the
  saddle point that identical starting components would sit at, but it does
  not pin labels to any external ordering. Match classes before comparing two
  fits. The likelihood-ratio test of `J` against `J+1` is not valid here --
  the null puts a class probability on the boundary -- so use `AIC()`/`BIC()`
  and interpretability instead.

* **`psfm()` gave every firm the SAME simulation draws.** A single `R x 2`
  Halton block was built once and recycled across all `N` firms, in three
  separate places. `sfm()` already does the opposite and says why: Train (2002,
  p. 228) attributes Halton's advantage to its coverage *and to the negative
  correlation it induces across observations*, and the second half only exists
  when observations get different blocks. Sharing one block also makes every
  firm's simulation error the same realization, so it cannot average out as `N`
  grows. Firm `i` now gets its own contiguous block.

  Measured, GTRE by simulated ML, 30 replications at each of N = 50/100/200,
  old and new on identical seeds so the comparison is paired: `sigh` improves
  on **62 of 87 pairs** (Wilcoxon p = 0.015) and `sigr` on 48 of 87
  (p = 0.048), while `beta1` -- the control, since the frontier slope is not
  what the draws integrate over -- is unaffected at 46 of 87 (p = 0.41).
  Collapses of `sigh` to exactly zero halved, 11/87 to 6/87.

  The improvement **grows with `N`**, which is the point: 16/28 at N = 50
  (p = 0.58), 21/29 at N = 100 (p = 0.25), 25/30 at N = 200 (p = 0.012). That
  is the shape the mechanism predicts.

  **It does not fix the GTRE convergence failure.** RMSE(N=50)/RMSE(N=200) for
  `sigh` moves from 0.94 to 1.08 where root-n requires 2.00. Shared draws were
  a contributor, not the cause.

* **`rand.gtre` randomized the draws in a way that destroyed what they were
  for.** The old code drew 9999 random permutations of Halton dimension 1 and
  kept whichever correlated least with dimension 2. A two-dimensional Halton
  sequence is worth having because the *pairs* cover the unit square evenly;
  permuting one column preserves each margin, randomizes the pairing, and
  leaves joint coverage no better than random. Measured at `R = 150` with
  primes 2 and 3 and 1000 discarded: joint discrepancy 0.0217 before the
  shuffle and 0.0550 after, against 0.1400 for purely random pairing -- to
  remove a correlation of 0.008 that was never a problem.

  `rand.gtre` now applies a uniform shift modulo 1 (Tuffin 1996; Train 9.3.4),
  which moves the lattice without disturbing its structure: same measurement,
  0.0250. Results change for anyone who passed `rand.gtre`.

* **`psfm(model_name = "KSS")` -- Kneip, Sickles and Song (2012).** The firm
  effect is a smooth function of time lying in an `L`-dimensional space whose
  basis is estimated from the data:
  `alpha_it = sum_r theta_ir * g_r(t)`. Nothing on CRAN implements this, so it
  is a differentiator rather than catch-up.

  It makes explicit what the other time-varying estimators are special cases
  of: `SSFE` is `L = 1` with a constant basis, `LS` is `L = 1` with the basis
  free, `CSS` is `L = 3` with the basis fixed to `{1, t, t^2}`, and `KSS`
  estimates both. Follows the original's three steps -- cross-sectional
  centering by period, then smoothing of each firm's residual trajectory
  followed by an eigendecomposition of their empirical covariance, then
  loadings by least squares -- iterated with `beta`. Balanced panels only,
  which is what the estimator is defined on; it refuses rather than quietly
  fitting something else.

  Recovers the true rank exactly at `L = 1, 2, 3` on simulated data
  (n = 110, T = 10), with the estimated basis spanning the true factor space
  and `sigma_v` recovered to within 3%.

  Both tuning choices are exposed. `kss_smooth` defaults to GCV; `kss_L`
  defaults to the Bai and Ng (2002) `IC_p2` criterion, which is **not** the
  original paper's threshold rule -- said plainly in `?psfm` because a
  different rule can select a different `L`. `$kss$eigenvalues` is returned so
  the choice can be inspected.

  The automatic search is capped at `floor(T/2)`, and warns if it selects that
  cap. `IC_p2` works because the residual variance flattens out once the real
  factors are in; on a short panel, letting `L` reach `T - 1` lets the factors
  span nearly the whole time dimension, so the variance collapses and no
  penalty competes. Measured at n = 100 over five seeds and two true ranks, an
  uncapped criterion returned its maximum on **all ten** designs at `T = 6` and
  again at `T = 8`. With the cap it recovers the true rank exactly from `T = 8`
  upward. An explicit `kss_L` is the user's own call and is not subject to the
  cap.

  One trap worth knowing: `KSS`'s period centering means its `alpha_hat`
  carries no period mean while `CSS`'s and `LS`'s do, so the two differ by a
  per-period constant *by construction*. Correlating them directly understates
  the agreement badly (0.93 against 0.99 on a design where they agree). Compare
  `u_hat`, which is a within-period contrast and invariant to the shift.

* **`psfm(model_name = "SSRE")` and `"SSCRE"`** complete the Schmidt-Sickles
  (1984) family that `SSFE` starts: the GLS estimator, and its
  correlated-random-effects correction. `SSRE` is more efficient than `SSFE`
  and identifies time-invariant regressors, but assumes the effects are
  uncorrelated with the regressors -- in a frontier model, that a firm's
  inefficiency is unrelated to its input choices.

  `SSCRE` adds the within-firm means of the time-varying regressors (Mundlak
  1978), which models that correlation rather than assuming it away. Its
  slopes then equal `SSFE`'s within slopes **exactly**, which the tests pin at
  1e-6; the coefficients on the added means (reported with a `.mean_` prefix)
  are the Mundlak form of the Hausman test.

* **Two classical time-varying panel estimators: `psfm(model_name = "CSS")`
  and `psfm(model_name = "LS")`.** Cornwell, Schmidt and Sickles (1990) and
  Lee and Schmidt (1993). Both sit between `"SSFE"`, which holds inefficiency
  fixed over time, and the ML panel models, which let it move only by assuming
  a distribution for `u` and a decay path for it: these let the firm effect
  itself vary with time and read inefficiency off it, assuming no distribution
  at all.

  They differ in what they spend on flexibility. `CSS` gives every firm its
  own quadratic in time (`3N` parameters; firms may overtake one another;
  `T_i >= 4` required before a firm contributes to beta). `LS` imposes one
  common temporal pattern scaled per firm, `alpha_it = delta_t * alpha_i`
  (`N + T - 1` parameters; the ranking of firms cannot change, only the
  spread). Cross-over versus no cross-over is the substantive choice between
  them.

  `CSS` residualizes on each firm's own `(1, t, t^2)` rather than building
  `3N` dummy columns, and refuses a regressor spanned by those quadratics (a
  pure time trend, say) instead of reporting it as a small coefficient. `LS`
  is a rank-one factor model and is fitted by alternating least squares, so an
  unbalanced panel needs no special-casing. On simulated data at N = 90,
  T = 8, `LS` recovers the temporal pattern with correlation > 0.99 to truth
  and the firm scales at 0.995; both recover the frontier slopes and the noise
  standard deviation.

  Neither is maximum likelihood, so neither carries an `$opt` component and
  `logLik()`/`AIC()`/`BIC()` return `NA` with a warning, as for `SSFE`.

* **The wrong-skew boundary report fired on a knife-edge.** `sfm()` flags
  `sigma_u_at_bound` when the one-sided scale collapses relative to the
  residual SD, and the threshold was `1e-3`. On the design used to test it the
  collapsed sample sits at `8.6e-4` and the interior ones at `0.47-0.71`, so
  the threshold cleared the case that matters by a factor of **1.16** and the
  others by 470. An optimizer's stopping point moves by more than 16% between
  BLAS implementations: the same fit reported `TRUE` locally and `FALSE` on
  another platform. Now `1e-2`, which sits in the middle of a gap spanning
  nearly three orders of magnitude and leaves every sample at least 11x clear.

  This only ever changes the SECOND of two conditions -- the warning also
  requires wrong skew (`m3 >= 0`), which is a deterministic function of the
  data -- so the practical effect is that the report is stable across
  platforms rather than that it fires more often.

* **`nobs()` counted rows SUPPLIED rather than rows USED, so `BIC()` was
  computed against the wrong n whenever any row was dropped for
  missingness.** `nobs.sfareg()` re-evaluated the `data` argument of the
  recorded call, which of course returns the full frame; `data_proc2()` has
  always dropped incomplete cases before fitting. On a 300-row frame with 7
  missing values `nobs()` returned 300 against 293 actually used. The fitted
  object now records the effective count directly, and `nobs()` falls back to
  the length of a stored per-observation vector before it resorts to
  re-evaluating the call.


* **`sfm(model_name = "NE")` could return a positive log-likelihood with
  `sigma_u = 0` and a divergent `sigma_v`.** The log-density was built as
  `-log(sigma_u) + log Phi(z) + eps/sigma_u + sigma_v^2/(2 sigma_u^2)`. Using
  `pnorm(log.p = TRUE)` for the middle term is correct as far as it goes, but
  that term and the tilt both diverge like `z^2/2` with opposite signs as
  `sigma_u -> 0`, and their sum is a catastrophic cancellation. At
  `sigma_v = 587`, `sigma_u = 1e-7` the two are `-/+1.7252e19`, where
  consecutive doubles are 2048 apart, so the sum came back as rounding noise --
  positive noise, which the optimizer then maximized by running `sigma_u` to
  its lower bound.

  This was not rare and it was silent: scanning 150 samples at
  `lambda = 0.75`, `N = 100`, one returned `sigma_v = 587.4`, `sigma_u = 0`,
  `logLik = +468992`, and another `sigma_v = 3.4e15`, `logLik = +7.9e30`, both
  as ordinary `sfareg` objects with no error and no warning. Across a
  12-cell design at 1,500 replications the rate was 0.74%, reaching 4.4% at
  `lambda = 0.5`, `N = 100`.

  New `.log_phi_tilt()` (`matrix_utils.R`) does the cancellation analytically
  instead. Because
  `z^2/2 = eps^2/(2 sigma_v^2) + eps/sigma_u + sigma_v^2/(2 sigma_u^2)`, the
  tail expansion of `log Phi(z)` cancels the tilt exactly and leaves only
  `-eps^2/(2 sigma_v^2) - log(2 pi)/2 - log(-z) + log1p(-1/z^2 + 3/z^4 - ...)`,
  in which no large intermediate is ever formed. It agrees with the previous
  expression to 8.5e-13 over 75,000 evaluations wherever the previous one was
  trustworthy, integrates to 1, and takes **this** failure rate to 0 in 2,200
  fits.

  To be precise about what that last number does and does not cover: it is the
  rate of the catastrophic-cancellation failure -- a positive log-likelihood
  with `sigma_u` at its bound -- and nothing else. It is not a claim that `NE`
  fits are uniformly well behaved. `convergence/CONVERGENCE_STATUS.md` records
  a separate and still-open defect at the other end of the parameter space: at
  `lambda = 3.33` the root-n sweep shows a tail of 2-4% of replications landing
  far from the truth on some parameter, and unlike the rest it does not shrink
  with `n`. There the scales are near the truth and the frontier SLOPES break,
  with `beta1` and `beta2` trading against each other, so it is a different
  failure from this one and is not addressed by this fix. The two live in
  different regions -- cancellation bites at LOW lambda, the tail at high -- and
  the median fit is unaffected in both.

* **`"NGE"` carried the identical defect and is fixed the same way.** Its
  likelihood is a difference of two exponentially tilted Gaussians, both with
  the same structure as `"NE"`; both now go through `.log_phi_tilt()`, in the
  likelihood and in the post-estimation efficiency block. Agreement with the
  previous expression is 7.1e-10. `"TSL"`, `"TTNE"` and `"NE_Z"` share the
  structure and have **not** been checked.

* **`"NE"` starting values now come from a bias-corrected moment estimator**
  (new `R/ne_start.R`) rather than the flat `sigma_u = sigma_v = 0.1` that
  `start_cs()` hands every cross-sectional model. Minus the mean negative OLS
  residual has an exact asymptotic bias factor `h(lambda)`; `.ne_start()`
  divides it out, shrunk by the share of MSE the bias accounts for. Over an
  88-cell design at 2,000 replications, total MSE of the start against the
  truth: 5.49 for the uncorrected version, 2.86 for the third-moment (COLS)
  inversion, 1.88 for their minimum, 1.82 for this one. Effect on the fitted
  MLE is smaller and confined to `lambda >= 1.5`, where it lowers MSE of
  `sigma_u-hat` by 5-51%; below `lambda = 1` it changes nothing. Derived in
  Bernstein, Parmeter and Wright, "Starting Values for the Normal-Exponential
  Stochastic Frontier Model".

  Note this leaves the flat `0.1` start in place for every *other*
  cross-sectional model, which `PROJECT_STATUS.md` still lists as open.

* **`sfm()` now reports when `sigma_u` sits on the zero boundary under wrong
  skewness, instead of saying nothing.** A one-sided scale at zero when the OLS
  residuals are skewed the wrong way is the *correct* maximum likelihood
  estimate -- the Type I failure of Olson, Schmidt and Waldman (1980) -- not a
  numerical problem. The `estimator = "cols"` path has warned about this for
  some time; the likelihood path did not, so a user got a boundary fit with no
  explanation.

  New `$wrong_skew`, `$sigma_u_at_bound` and `$residual_m3` components, plus a
  warning that says what the boundary means and that the efficiency scores are
  uninformative there. Measured at `lambda = 0.75`, `N = 100` over 600 fits:
  13.2% of samples are wrongly skewed, 17.7% of *those* put `sigma_u` on the
  boundary, and **not one** correctly skewed sample does. That is why the fix
  is to report the boundary rather than bound `sigma_u` away from zero -- a
  bound would corrupt precisely the samples where the boundary is the answer.

* **`"NE"` no longer emits a stream of `NaNs produced` warnings.** Its
  likelihood had no guard against a non-positive scale, so the optimizer
  probing `sigma_u <= 0` evaluated `log(sigma_u)` and warned -- 17 times in a
  single fit on a wrongly skewed sample, burying the warning that mattered. It
  now returns a large finite penalty for out-of-domain parameters, as `"NLN"`
  and `"NW"` already do. **No estimate changes**: this bounds the objective's
  domain, not the estimate.

* **`"NGE"` aborted outright on about 7% of small samples, and `"NU"` could
  have.** Both guarded their scale parameters with
  `return(.Machine$double.xmax)`. `optim()` differences the objective to form a
  gradient, and differencing 1.8e308 overflows to a non-finite value, so the
  fit died with `non-finite finite-difference value` instead of the optimizer
  being steered away. Measured before the change: 3 of 45 `"NGE"` fits at
  `N = 150` failed this way, including at `sigma_u = 1, sigma_v = 0.3`. Both now
  use the same finite penalty the other branches use -- 0 of 45 failures after,
  and `"NU"` estimates are bit-identical across 45 fits, since the penalty only
  ever applies where the parameters are already inadmissible.

* **`ttsfm()` gains `z_link` too**, completing the cross-sectional half of the
  variance-determinant inconsistency. Two-tier models placed `z'delta` on the
  standard deviation, like `sfm()` and unlike `psfm()` and the competing
  packages; `z_link = "var"` aligns them. The default is `"sd"`, and fits under
  it are byte-identical to before -- checked against the previous build with
  the particle swarm seeded, for both `"TTNE"` and `"TTHN"`.

  Both tiers move together: `delta` on `zu` and on `zw` each halve exactly
  between the links, while the maximised log-likelihood and the frontier
  coefficients are unchanged.

* **`marginal_effects()` now works on panel fits**, `psfm(model_name = "TRE_Z")`
  and `psfm(model_name = "GTRE_Z")`. These place `z'delta` on the variance,
  and the link is read from the fit rather than assumed, so their effects are
  directly comparable with an `sfm()` fit -- including one fitted with the new
  `z_link = "var"`.

  `"GTRE_Z"` separates persistent inefficiency from transient, and **both**
  are reported: `component = "u"` for the transient block and
  `component = "h"` for the persistent one. Each block is located by name
  rather than by position, because the `sigma_h` coefficients are the
  *trailing* ones -- a positional rule would have silently reported `sigma_h`
  effects under a `sigma_u` label. Output columns carry the component
  (`dE_u.dz` against `dE_h.dzp`), so a table cannot be misread once separated
  from the call that produced it. The default is `"u"`, so nothing existing
  changes.

* **`sfm()` gains `z_link`, which fixes a real trap in comparing `_Z` fits
  across the package.** `sfm()`'s `"NHN_Z"`/`"NE_Z"` put the
  variance-determinant linear predictor on the standard deviation,
  `sigma_u = exp(z'delta)`, while `psfm()`'s `"TRE_Z"`/`"GTRE_Z"` put it on the
  variance, `sigma_u = sqrt(exp(z'delta))` -- as do the competing packages.
  Since `exp(eta) = exp(eta/2)^2` the two fit the same model and return the same
  `sigma_u`, the same log-likelihood and the same marginal effects, but `delta`
  under the SD link is exactly **half** `delta` under the variance link. Reading
  a `delta` from one family as if it came from the other therefore doubles or
  halves every reported effect.

  `z_link = "var"` puts an `sfm()` fit on the same footing as `psfm()` and
  `sfaR`. **The default is `"sd"`, so no existing result changes.** The
  efficiency predictor and the `z_spec` that `marginal_effects()` reads both
  follow whichever link was used.

* **`sfm()`'s simulated-ML models `"NLN"` and `"NW"` now estimate the composed
  density from two proposals at once, which fixes an accuracy defect that no
  single-proposal scheme can.** The integrand is a product of a normal kernel of
  width `sigma_v` centred at `u = -eps` and the inefficiency density. Drawing
  from the inefficiency quantile -- what both models did through 1.1.5 -- misses
  the kernel whenever it is the narrow factor: at `eps = -4.7` with
  `sigma_v = 0.3` the spike needs a draw beyond the 0.9999 quantile, while 400
  draws reach only about 0.9975, and five observations in 800 carried 72% of
  `"NW"`'s total error. Drawing from the noise instead cures that case and
  breaks the mirror one, wherever the inefficiency density is the narrow factor.

  Both draws are now taken -- half the count each -- and combined by the balance
  heuristic of Veach and Guibas (1995), each draw weighted by which proposal was
  likelier to have produced it. There is no selection rule to get wrong. Total
  error against adaptive quadrature, over an 80-cell parameter grid on a
  300-observation sample, restricted to the cells an optimizer can reach:

  | cells within ... of the optimum | 1.1.5 scheme | now |
  |---|---|---|
  | `"NW"`, 200   | 220, worst cell 172    | **16**, worst cell **5.7** |
  | `"NW"`, 500   | 3633, worst cell 1910  | **40**, worst cell **11** |
  | `"NW"`, 1000  | 17900, worst cell 8553 | **128**, worst cell **64** |
  | `"NLN"`, 200  | 233, worst cell 120    | **11**, worst cell **1.0** |
  | `"NLN"`, 500  | 3564, worst cell 1553  | **49**, worst cell **8.8** |
  | `"NLN"`, 1000 | 15962, worst cell 6639 | **205**, worst cell **63** |

  (Each scheme is measured at its own draw rule, so `"NW"`'s 1.1.5 column uses
  the 400 draws that version would have taken and the new column 200.)

  `"NW"`'s draw rule drops from `max(400, 8*sqrt(n))` to the `max(200,
  3*sqrt(n))` `"NLN"` already used, so at n = 1500 an `"NW"` fit costs about
  what it did in 1.1.5 (roughly 50 seconds) while being far more accurate;
  `"NLN"` costs about 1.5x more, since its draw count did not fall. `Nsim`
  still means the total number of draws.

  The two models now also converge at the same rate, so the old advice that
  `"NLN"` needs far more draws than `"NW"` no longer applies -- that was a
  property of the integration scheme, not the model. At n = 3000, measured
  against quadrature at the true parameters, total simulation error is 3.4
  log-likelihood units for `"NLN"` and 3.0 for `"NW"` at the default, falling
  to 0.07 and 0.06 at `Nsim = 6400`. The comparable figure for the old scheme
  was 226.8.

  The accuracy gain shows up in the estimates. Over 12 replications at
  n = 800, bias and RMSE against the true parameters:

  | | sigma_v | sigma_u | shape | x1 |
  |---|---|---|---|---|
  | `"NW"` 1.1.5 | +0.036 / 0.066 | -0.149 / 0.237 | -0.160 / 0.293 | -0.021 / 0.043 |
  | `"NW"` now   | **-0.006 / 0.044** | **-0.024 / 0.112** | **+0.001 / 0.111** | **+0.001 / 0.045** |
  | `"NLN"` 1.1.5| -0.012 / 0.022 | +0.005 / 0.101 | +0.012 / 0.166 | +0.003 / 0.037 |
  | `"NLN"` now  | -0.004 / 0.035 | **+0.003 / 0.080** | **+0.005 / 0.131** | +0.002 / 0.038 |

  `"NW"` is the clear case: its `sigma_u` and shape bias essentially vanish and
  their RMSE roughly halves. `"NLN"` improves more modestly, and its `sigma_v`
  RMSE is slightly worse -- at 12 replications that difference is within noise
  and should not be read either way.

  The efficiency predictor reuses the weights the likelihood computed, so
  `u_hat` cannot drift away from the density that was maximised.

* **`sfm()` gains user control over its simulated-ML draws**, following the
  practice in Train (2002, ch. 9) and matching the options `sfaR` exposes:
  `sim_type` (`"halton"`, `"sobol"`, `"torus"`, `"uniform"`), `antithetics`,
  `sim_burn`, `sim_scrambling`, `sim_prime` and `sim_seed`. These affect the
  models fitted by simulation, `"NLN"` and `"NW"`.

  `"sobol"` with `sim_scrambling` 1-3 is the scrambled sequence of Bhat
  (2003), which removes cross-dimension correlation while keeping the
  coverage. `"uniform"` exists as a baseline to measure against, not as a
  recommendation -- Bhat (2001), quoted by Train, found 100 Halton draws more
  precise than 1000 pseudorandom ones. `antithetics = TRUE` takes half the
  draws from the sequence and creates the rest as mirror images (Hammersley
  and Morton 1956), which costs nothing.

  **Defaults reproduce the previous behaviour exactly** -- Halton, no
  antithetics, 1000 leading elements discarded -- and this is asserted by a
  test that rebuilds the old construction and compares byte for byte, so no
  existing result moves.

* **The draw construction is now one shared internal function**,
  `.sml_draws()`, rather than being written inline at each site. That matters
  beyond tidiness: the sites had drifted apart, with `sfm()` giving each
  observation its own contiguous block of the sequence while `psfm()` hands
  every firm the same one. Train (ch. 9) attributes the value of a
  low-discrepancy sequence to its coverage *and* to the negative correlation
  it induces across observations, and the second only exists when units get
  different blocks. The shared constructor always blocks by unit and takes a
  `dim` argument for the multi-dimensional panel case, so the other entry
  points can adopt it directly.

* **New model `sfm(model_name = "TSL")`: the normal / truncated skew-Laplace
  frontier (Wang 2012).** Inefficiency has the signed-mixture density
  `f(u) = ((1+lambda)/(sigma_u(2*lambda+1))) [2 exp(-u/sigma_u) -
  exp(-(1+lambda) u/sigma_u)]` on `u >= 0` -- the second exponential enters
  *negatively* -- so it nests the exponential model as `lambda -> 0` while
  allowing a non-monotonic inefficiency density. Reports
  `(sigv, sigu, lambda)` plus the frontier, with both `u_hat` and
  `exp_u_hat`.

  The composed density is evaluated as a difference in log space rather than
  by subtracting the raw terms. Both carry `sigma_v^2/(2 sigma_u^2)`, which
  overflows once `sigma_u` is small relative to `sigma_v` -- on a routine grid
  spanning `sigma_u` in {0.05, 0.1, 0.3} and `sigma_v` in {0.3, 1} the direct
  form returns `NaN` or `Inf` at 9 of 30 points, including `sigma_u = 0.05`,
  `sigma_v = 1`, `eps = 0`, which is not an extreme point. The log form is
  finite throughout.

  Verified two ways: against numerical convolution of the implied `u` density
  with the normal noise (agreeing to 8 decimals across `lambda` in
  {0.5, 1.5, 4}), and against an independent implementation, which it matches
  to an identical log-likelihood and identical coefficients on 8 of 8 samples.
  `lambda` is a shape parameter and is the least sharply identified of the
  three -- read its `t`-value before interpreting it.

* `data_gen_cs()` gains `lam_tsl` and the matching `u_tsl` / `y_pcs_tsl`
  columns. The signed mixture cannot be drawn by picking a component, so it
  is drawn by rejection off its first exponential; the acceptance probability
  never falls below 1/2.

* **New function `efficiency_ci()`: Horrace and Schmidt (1996) confidence
  intervals for individual inefficiency.** `sfa` has always reported point
  predictions of `u_i` -- `u_hat` (Jondrow et al. 1982) and `exp_u_hat`
  (Battese and Coelli 1988) -- and nothing about how sharply either is pinned
  down. Both are posterior means of `u` given the composed residual, and the
  posterior is a truncated normal in closed form, so the interval costs no
  estimation beyond the fit itself. Returns a data frame of
  `u_lower`/`u_hat`/`u_upper` and the corresponding technical-efficiency
  bounds, at any `level`.

  Available for `sfm()`'s `"NHN"`, `"NHN_Z"`, `"NE"` and `"NTN"`, whose
  posterior really is a truncated normal. Other models have a posterior of a
  different shape, and the function says so by name rather than returning a
  number the formulas do not support.

  These intervals condition on the fitted parameters: they do not narrow as
  `n` grows, because they measure the irreducible difficulty of splitting one
  residual into noise and inefficiency, not estimation uncertainty. That is
  worth knowing before reading an efficiency ranking closely -- on a routine
  200-observation half-normal fit the median 90% interval for a single unit
  spans about 0.49 in efficiency, against a total spread of 0.85 across all
  the point predictions.

* **`npsfm(method = "SZ")` now solves its DEA step in the package.** The
  output-oriented envelopment program is computed directly, as one linear
  program per unit, by the new internal `.dea_out()`, covering all four
  returns-to-scale settings. `Suggests` gains `lpSolve`, a general-purpose
  linear-programming solver (pure C, no system libraries).

  Results are unchanged, and were verified against an independent reference
  implementation over all four returns-to-scale settings at n = 30, 80 and
  150, with one to three inputs and one to two outputs: agreement to 3e-12 or
  better throughout.

* **New function `marginal_effects()`: what the `_Z` models' `delta`
  coefficients actually mean for inefficiency.** For a fit whose inefficiency
  scale depends on covariates, it returns the per-observation
  `d E[u]/d z_k` and `d Var[u]/d z_k`, with the average marginal effects
  attached. These are what applied papers report; a `delta` on its own is a
  coefficient in a log link for a scale parameter and is not interpretable in
  the units of either `u` or `z`.

  This also defuses a long-standing trap. `sfm()`'s `"NHN_Z"`/`"NE_Z"` put the
  linear predictor on the standard deviation and `psfm()`'s `"TRE_Z"`/
  `"GTRE_Z"` put it on the variance, so the same `delta` means different
  things in the two families -- the marginal effects differ by a factor of
  two. The effect is on the scale of `u` either way, so reporting it rather
  than the coefficient removes the ambiguity. The returned object records
  which convention was used.

  Currently covers `sfm()`'s `"NHN_Z"` and `"NE_Z"`; `psfm()`'s panel `_Z`
  models are not wired up yet. Standard errors are deliberately not reported --
  see `?marginal_effects` for why.

* `"NHN_Z"` and `"NE_Z"` fits now carry a `$z_spec` component (the
  variance-determinant design, its coefficients, and the link convention),
  which is what `marginal_effects()` reads.

* **`sfm(model_name = "NNAK")` now starts from the method of moments.** The
  normal-Nakagami likelihood has `sigma_u -> 0` as a genuine attractor, and
  the previous hard-coded start of `sigma_u = sigma_v = 0.1` sat next to it.
  Over twelve samples at n = 3000 with a true `sigma_u` of 1, the old start
  drove `sigma_u` to 0.0013 and to 0.0000 on two of them -- inefficiency
  vanishing altogether, with the frontier intercept pushed negative -- for a
  log-likelihood 44 and 45 points worse than the new start reaches. The new
  start was never worse on any of the twelve and was strictly better on seven,
  for a mean gain of 7.9 log-likelihood points.

  The construction follows FronPy (Stead 2024, *Journal of Productivity
  Analysis*, the paper these closed forms come from): invert the half-normal
  moment equations for the two scales and shift the frontier intercept up by
  the implied `E[u]`. The half-normal is the right auxiliary because it is the
  `m = 1/2` member of the Nakagami family. The shape still starts at 0.5, as
  before. When the residuals are skewed the wrong way and the moment equations
  have no admissible solution, the old constants are used as before.

* **`sfm(estimator = "mols")` is accepted as a synonym for
  `estimator = "cols"`.** What the package computes under that name is the
  *modified* OLS moment estimator of Olson, Schmidt and Waldman (1980) -- it
  inverts the second and third central moments of the OLS residuals -- and not
  Winsten's *corrected* OLS, which shifts the intercept by the largest
  residual and estimates no variance parameters. The documentation now names
  it MOLS and draws the distinction. `"cols"` is unchanged and not deprecated.

* Fitted `"sfareg"` objects from the four models above carry a new
  `$u_posterior` component (`mu_star`, `sigma_star`), which is what
  `efficiency_ci()` reads.

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
    convexity.

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

* **The vignette still described `"PL80"` and `"BC92"` as wrappers around
  `frontier::sfa()`.** They have been native closed-form maximum-likelihood
  implementations since 1.1.3, the release that dropped **frontier** as a
  dependency -- so for four releases `vignette("intro_to_psfm")` advertised a
  dependency on a competing package that the package had deliberately removed,
  and that appears in neither `Imports` nor `Suggests`. `?psfm` had it right
  throughout ("natively estimated ... verified against `frontier::sfa()`
  before that dependency was removed"); only the vignette was stale.

  Two related corrections in the same file. The `psfm_bootstrap()` section gave
  the wrong reason for excluding `PL80`/`BC92` -- it is that they do not expose
  the `$U`/`$H` efficiency structure the function reads, not that they wrap
  anything -- and its list of supported models omitted `TFE_WMLE` and
  `GTRE_FML`, the latter being what `psfm(model_name = "GTRE")` returns under
  its default estimator. And the model table, which says "eleven" choices,
  now says twenty-one and points at `?psfm` rather than implying the table is
  exhaustive.

* **`psfm_bootstrap()` refused the estimator `psfm()` reaches by default.**
  Since 1.1.4, `psfm(model_name = "GTRE")` defaults to `estimator = "fiml"`
  and returns an object whose `$model_name` is `"GTRE_FML"`.
  `psfm_bootstrap()` checked that name against a list that did not contain it,
  so the natural two-line sequence -- fit a GTRE model, then bootstrap it --
  stopped with `"psfm_bootstrap() does not support model_name = 'GTRE_FML'"`,
  and the error message's own explanation named only the moment-based and
  wrapped models, none of which applied. Reaching the bootstrap at all
  required knowing to pass `estimator = "sml"` to the fit.

  `"GTRE_FML"` is now supported. The data-generating process is the one the
  `"GTRE"` arm already simulated -- same four-component model, different
  estimator -- so only the parameter layout differs: FIML reports the four raw
  scales `(sigr, sigv, sigh, sigu)` after the frontier block rather than
  `"GTRE"`'s `lambda`/`sigma` reparameterization. `$H` and the boundary
  warning come with it.

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
