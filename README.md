# sfa

<!-- badges: start -->
[![R-CMD-check](https://github.com/davidhbernstein/sfa/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/davidhbernstein/sfa/actions/workflows/R-CMD-check.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/sfa)](https://CRAN.R-project.org/package=sfa)
[![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/sfa)](https://CRAN.R-project.org/package=sfa)
[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://cran.r-project.org/web/licenses/GPL-2)
<!-- badges: end -->

**Stochastic frontier analysis in R.** A single, consistent interface to a wide
range of cross-sectional, panel, latent-class, zero-inefficiency, two-tier,
sample-selection, endogenous-regressor, copula and nonparametric stochastic
frontier models, with a common formula syntax for modelling the variance of each
error component and a common `"sfareg"` result object that works with the
standard R modelling generics.

Beyond fitting, it provides the tools to *choose* among the fifteen
cross-sectional inefficiency distributions rather than assume one, and to check
whether the chosen specification is defensible -- see
[Model selection and diagnostics](#model-selection-and-diagnostics).

Written by David H. Bernstein, Christopher F. Parmeter and Alexander D. Stead.

## Installation

The released version from CRAN:

```r
install.packages("sfa")
```

The development version from GitHub:

```r
# install.packages("remotes")
remotes::install_github("davidhbernstein/sfa")
```

The development version here is ahead of CRAN; see [`NEWS.md`](NEWS.md) for what
has changed, including one deliberate breaking change to `psfm(model_name = "TFE")`.

## Quick start

```r
library(sfa)

## Simulate a cross-section with known true parameters
cs <- data_gen_cs(N = 1000, rand = 1, sig_u = 0.3, sig_v = 0.3,
                  cons = 0.5, beta1 = 0.5, beta2 = 0.5, a = 4, mu = 1)

## Normal-half normal frontier
fit <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = cs)

summary(fit)
coef(fit)              # lambda, sigma, (Intercept), x1, x2 -- see note below
logLik(fit)
head(fit$exp_u_hat)    # technical efficiency, E[exp(-u) | e]
head(fit$med_u_hat)    # median efficiency prediction (NHN only)
```

Which efficiency predictors come back depends on the model: `exp_u_hat`
(Battese and Coelli 1988) is returned by most of `sfm()`'s models, and the
Jondrow et al. (1982) point predictor `u_hat` = `E[u | e]` by `NE`, `NTN`,
`NU`, `NGE`, `NLN`, `NW`, `THT` and `tHN`. See `?sfm` for the full list.

A panel model, with a determinant of the inefficiency variance:

```r
pd <- data_gen_p(t = 10, N = 100, rand = 100, sig_u = 1, sig_v = 0.3,
                 sig_r = 0.2, sig_h = 0.4, cons = 0.5, beta1 = 0.5, beta2 = 0.5)

fit_p <- psfm(y_tre_z ~ x1 + x2 | z_gtre, model_name = "TRE_Z",
              data = pd, individual = "name")
```

## The nine entry points

| Function | Fits | Estimators |
|---|---|---|
| `sfm()` | Cross-sectional frontiers | 15 |
| `psfm()` | Panel frontiers | 21 |
| `lcsfm()` | Latent-class frontiers | 3 |
| `zsfm()` | Zero-inefficiency frontiers | 2 |
| `ttsfm()` | Two-tier frontiers | 3 |
| `selsfm()` | Sample-selection frontiers | 1 |
| `ivsfm()` | Frontiers with endogenous regressors | 3 |
| `copsfm()` | Dependence between the error components | 6 copula families, 15 with rotations |
| `npsfm()` | Nonparametric frontiers | 5 |

All but `npsfm()` return an object of class `"sfareg"`. `npsfm()` returns
`"npsfareg"` instead — a kernel-estimated frontier has no parameter vector with
standard errors, so `coef()`, `vcov()` and `logLik()` would have nothing to
return.

`selsfm()` and `ivsfm()` do not take their equations through the `|` pipes:
`selsfm()` takes `selection` and `frontier` as separate formulas, and `ivsfm()`
takes `formula`, `endogenous` and `instruments`. Both **reject** a `|` segment,
because a pipe already means "variance determinant" everywhere else and reusing
it would give one character two meanings.

### `sfm()` — cross-sectional

| `model_name` | Distribution of `u` |
|---|---|
| `NHN`, `NHN_Z` | half normal (`_Z`: with variance determinants) |
| `NE`, `NE_Z` | exponential (`_Z`: with variance determinants) |
| `NTN` | truncated normal |
| `NR` | Rayleigh |
| `NU` | uniform |
| `NG` | gamma |
| `NNAK` | Nakagami |
| `NGE` | generalized exponential |
| `NLN` | lognormal |
| `NW` | Weibull |
| `tHN` | half normal, with Student-*t* noise |
| `THT` | half *t*, with Student-*t* noise |
| `TSL` | truncated skew-Laplace |

`sfm()` also offers `estimator = "cols"` — corrected OLS (Olson, Schmidt and
Waldman 1980), closed-form and deterministic, for `NHN`, `NE` and `NG` — and
robust divergence-based alternatives to MLE via `robust = "mlqe" | "psi" | "mdpd"`
for `NHN`.

The tuning parameter of those robust criteria does not have to be guessed:
`hscore_select()` chooses it by minimising the Hyvarinen score, `calibrate_c()`
gives a fixed weight-matching alternative, and `density_weights()` shows what
the estimator did to each observation. `influence_sfa()` reports the influence
function of any fit — which observations move it, and whether the specification
lets any single one of them move it without bound.

```r
fit <- sfm(y ~ x1 + x2, data = d, model_name = "NHN")
sel <- hscore_select(fit, method = "mlqe")          # data-driven
calibrate_c(sigma_v = 0.3, sigma_u = 0.6)           # or fixed
sfm(y ~ x1 + x2, data = d, robust = "mlqe", c_mlqe = sel$c)
```

Any of `NHN`, `NE` and `NTN` can additionally take covariates in more than one
error component, as named formulas rather than further pipe segments:

```r
sfm(y ~ x1 + x2 | z_u, vhet = ~ z_v, model_name = "NHN_Z")   # heteroskedastic v
sfm(y ~ x1 + x2, muhet = ~ z_mu, model_name = "NTN")         # Battese–Coelli (1995)
```

`vhet` drives the noise scale, `uhet` the inefficiency scale (the same thing
the `| z` segment does), and `muhet` the pre-truncation mean.

### `psfm()` — panel

| `model_name` | Estimator |
|---|---|
| `TRE`, `TRE_Z` | true random effects (Greene 2005) |
| `GTRE`, `GTRE_Z` | generalized true random effects, four-component |
| `GTRE_FML` | GTRE by full maximum likelihood |
| `GTRE_SEQ1`, `GTRE_SEQ2` | sequential/moment-based GTRE |
| `TFE` | true fixed effects (Greene 2005) |
| `TFE_WMLE` | within MLE (Chen, Schmidt and Wang 2014) |
| `FD` | first differences |
| `SSFE` | Schmidt and Sickles (1984) fixed effects (within) |
| `SSRE`, `SSCRE` | Schmidt–Sickles random effects, and correlated random effects (Mundlak 1978) |
| `CSS` | Cornwell, Schmidt and Sickles (1990), firm-specific quadratic in time |
| `LS` | Lee and Schmidt (1993), one common temporal pattern scaled per firm |
| `KSS` | Kneip, Sickles and Song (2012), data-driven temporal basis |
| `PL80` | Pitt and Lee (1980), time-invariant |
| `BC92` | Battese and Coelli (1992) time decay |
| `K1990`, `K1990modified` | Kumbhakar (1990) time patterns |

The last five are one family: each writes the firm effect as
`alpha_it = sum_r theta_ir * g_r(t)` and reads inefficiency off it as distance
from the best firm, assuming no distribution for inefficiency at all. They
differ only in how much of that structure is assumed rather than estimated —
`SSFE` fixes `L = 1` with a constant basis, `LS` frees the basis, `CSS` fixes
`L = 3` to `{1, t, t^2}`, and `KSS` estimates both. `KSS` needs a balanced
panel; the others do not.

`GTRE_SEQ1`, `GTRE_SEQ2`, `SSFE`, `SSRE`, `SSCRE`, `CSS`, `LS` and `KSS` are
not maximum likelihood, so `logLik()` (and hence `AIC()`/`BIC()`) returns `NA`
for them.

`psfm_bootstrap()` provides a parametric bootstrap for GTRE-family fits,
parallelised over cores.

### `zsfm()` — zero inefficiency

`ZISF` and `ZISF_Z`: a mixture of a fully efficient regime and an inefficient
frontier regime, with the regime probability optionally parameterised by
covariates (`ZISF_Z`).

### `lcsfm()` — latent class

`LCM` and `LCM_Z`: the latent class frontier (Greene 2005; Orea and Kumbhakar
2004) — `n_class` unobserved technologies, each with its own frontier and its
own two scales, mixed by a multinomial logit whose covariates are optional
(`LCM_Z`). Returns posterior class probabilities and posterior-weighted
efficiency alongside the class-conditional predictions.

### `ttsfm()` — two tier

`TTNE` (normal–exponential–exponential), `TTHN` (normal–half normal–half
normal), and `TTNLS` (nonlinear least squares, no distributional assumption
beyond the means of the two one-sided components).

### `selsfm()` — sample selection

Greene's (2010) frontier for the case where the units in the sample are there
for reasons correlated with their inefficiency, so estimating on the selected
sample alone is biased. Estimated in two steps — probit, then simulated maximum
likelihood — and so it takes its two equations as separate arguments rather than
through pipes:

```r
selsfm(selection = participate ~ z1 + z2,
       frontier  = y ~ x1 + x2, data = d)
```

### `ivsfm()` — endogenous regressors

Amsler, Prokhorov and Schmidt (2016): one or more regressors correlated with the
statistical noise. Three estimators of the same model, not three models —
`IVLIML` (full-information maximum likelihood), `IVCF` (two-step control
function) and `C2SLS` (corrected 2SLS):

```r
ivsfm(y ~ x1 + x2, endogenous = ~ x2, instruments = ~ w1 + w2,
      data = d, model_name = "IVLIML")
```

With `uhet` the model becomes that of Amsler, Prokhorov and Schmidt (2017), in
which the environmental variables entering the inefficiency scale may themselves
be endogenous.

`endogeneity_test()` asks whether the correction was needed at all — a Wald test
of `H0: rho = 0`, i.e. that the noise is uncorrelated with the reduced-form
errors and a plain `sfm()` fit would already have been consistent.

```r
fit <- ivsfm(y ~ x1 + x2, endogenous = ~ x2, instruments = ~ w1 + w2,
             data = d, model_name = "IVLIML")
endogeneity_test(fit)
```

### `copsfm()` — dependence between `v` and `u`

Drops the independence assumption between the noise and inefficiency
components, coupling them with a copula and integrating the resulting density by
Gauss–Legendre quadrature (`n_nodes`).

```r
copsfm(y ~ x1 + x2, data = d, copula = "frank")
```

| `copula` | parameter | independence at | dependence it can express |
|---|---|---|---|
| `"gaussian"` | `rho` in (−1, 1) | 0 | both signs, no tail dependence |
| `"fgm"` | `theta` in [−1, 1] | 0 | both signs, but weak — Spearman `rho` = `theta`/3 |
| `"frank"` | `theta` real | 0 | both signs, full range, no tail dependence |
| `"clayton"` | `theta` > 0 | 0 | positive, **lower**-tail dependence |
| `"gumbel"` | `theta` >= 1 | 1 | positive, **upper**-tail dependence |
| `"joe"` | `theta` >= 1 | 1 | positive, heavier upper tail than Gumbel |

Clayton, Gumbel and Joe carry *only positive* dependence. Nothing rules out a
negative association between noise and inefficiency, so each also has `90` and
`270` rotations that reverse the sign (`"clayton270"`), and `180`, the survival
copula, which preserves it.

**Which of these actually work, measured rather than assumed.** Every density
here is verified three ways — it equals the second mixed partial of its own CDF,
it integrates to 1 over the unit square, and it is exactly 1 at the independence
parameter. That establishes the *densities* are right. It does not establish
that the dependence parameter is *recoverable*, and for most of them it is not.
Fitting each family to 25 samples generated from itself at n = 400:

| family | recovers its own `theta`? | collapsed to the independence bound |
|---|---|---|
| `frank` | yes (5.43 against a truth of 5) | **0%** |
| `clayton` | yes (2.24 against 2) | **0%** |
| `gaussian`, `fgm` | yes | 0% |
| `gumbel` | no | **36%** |
| `joe` | no | **40%** |
| `clayton270` | no | **56%** |
| `gumbel90` | no | **60%** |

On data generated from a Gumbel copula with Spearman `rho` = 0.685 at n = 2000,
*every* family — including the true one — returns the independence boundary, and
their log-likelihoods differ by less than 0.04. The likelihood is flat in the
dependence parameter. This is a property of the model, not of the code.

**So: prefer `"frank"` or `"clayton"`, and treat the others as exploratory.**
`copsfm()` warns when you select a family that did not recover, quoting its
measured collapse rate; the rotations that were never measured warn that they
were not, rather than implying either outcome.
More generally, the dependence parameter is estimated imprecisely even when it
is recoverable — on a Gaussian design at n = 600 its sampling standard deviation
is 0.385 against a truth of 0.5 — so a comparison across families is descriptive
rather than evidence for a dependence structure.

### `npsfm()` — nonparametric

Estimates the frontier by kernel regression instead of assuming it linear.

| `method` | Estimator |
|---|---|
| `FLW` | Fan, Li and Weersink (1996). Kernel regression for `E[y\|x]`, then the scale parameters from the residuals. Also supports `dist = "exp"`, `"gamma"`, `"unif"` |
| `SVKZ` | Simar, Van Keilegom and Zelenyuk (2017). Local method of moments; `sigma_u(x)` and `sigma_v(x)` vary with the covariates |
| `PSZ` (alias `KPST`) | Park, Simar and Zelenyuk. Local maximum likelihood |
| `MY` | Martins-Filho and Yao. Iterative local likelihood |
| `SZ` | Simar and Zelenyuk (2011). DEA monotonization of a prior smooth fit |

```r
f <- npsfm(y ~ x1 + x2, data = d, method = "FLW", dist = "hn")
head(fitted(f))       # the estimated frontier
head(f$exp_u_hat)     # technical efficiency
```

`PSZ` and `MY` run one numerical optimization per observation — for `MY`, per
observation per iteration — so expect them to be one to two orders of magnitude
slower than `FLW`. `npsfm()` takes a single-part formula and rejects a `| z`
segment: its heteroskedasticity is nonparametric in the covariates themselves.

Kernel regression comes from [**np**](https://CRAN.R-project.org/package=np),
and `SZ`'s DEA step solves one linear program per unit with
[**lpSolve**](https://CRAN.R-project.org/package=lpSolve). Both are in
`Suggests`, not `Imports`, so they are only required if you actually call
`npsfm()`:

```r
install.packages(c("np", "lpSolve"))
```

## Formula syntax

Variance determinants are supplied in extra pipe-delimited segments:

```
y ~ x1 + x2 | z | zp
```

* the first segment is the frontier,
* the second parameterises the variance of the first one-sided component (`sigma_u`),
* the third parameterises a second one-sided component where the model has one
  (`sigma_w` in two-tier models, `sigma_h` in `GTRE_Z`).

Omitted segments default to `1`, i.e. homoskedastic.

> **The link function differs by model family.** `sfm()`'s `NHN_Z`/`NE_Z` and
> `ttsfm()`'s `TTNE`/`TTHN` use `sigma = exp(z'delta)`, while `psfm()`'s
> `GTRE_Z`/`TRE_Z` use `sigma = sqrt(exp(z'delta))` — that is, `delta`
> parameterises the *variance* rather than the standard deviation. Check which
> convention applies before interpreting a coefficient on `z`.

## Model selection and diagnostics

The package offers fifteen cross-sectional inefficiency distributions. These are
the tools for choosing among them, and for asking whether the choice is
defensible at all.

| Function | Question it answers |
|---|---|
| `TIC()`, `vuong()` | Which of two non-nested specifications fits better, without assuming either is correct |
| `spec_test()`, `spec_test_all()` | Is this *pair* of noise/inefficiency distributions defensible, from OLS residuals alone |
| `moment_range()` | Can this pair *produce* the residuals' skewness and kurtosis at all — the range check that precedes the test |
| `sfma()` | What if the data does not identify one — average over distributions instead of choosing |
| `lcsfm_homogeneity()` | Does a latent-class fit beat a single technology |
| `gof_test()` | Is the assumed *inefficiency distribution* right, holding normal noise fixed |
| `cw_test()` | The same question without ever forming the composed density, from OLS residuals |
| `skewness_test()`, `inefficiency_test()` | Is there evidence of inefficiency at all (the wrong-skew problem) |
| `uhet_test()` | Does inefficiency depend on firm characteristics (validly, from a two-step fit) |
| `endogeneity_test()` | Was the endogeneity correction needed — is `rho` different from zero |
| `esfm()`, `symmetry_test()` | Fit a frontier when the residual skewness has the *wrong* sign |
| `influence_sfa()` | Which observations move the fit, and can any single one move it without bound |
| `hscore_select()`, `calibrate_c()`, `density_weights()` | Choosing and reading the robust-divergence tuning parameter |
| `efficiency()`, `meanefficiency()`, `efficiency_ci()` | Efficiency predictions, model-implied means, and Horrace–Schmidt intervals |
| `marginal_effects()` | Effects of the variance determinants on `E[u]` |
| `simulation_se()` | How much of a simulated-ML standard error is simulation noise |
| `pcomposed()`, `dcomposed()` | Distribution and density of the composed error (half-normal) |
| `pcomposed_model()`, `composed_cdf()` | The same CDF for any of the thirteen cross-sectional models |
| `sfa_diagnostics()` | Convergence and boundary diagnostics for a fit |

```r
## The score-based diagnostics differentiate the fitted likelihood, so the fit
## has to have kept it: pass keep_objective = TRUE. TIC(), vuong() and
## influence_sfa() all need this; the others do not.
fit_hn <- sfm(y ~ x1 + x2, data = d, model_name = "NHN", keep_objective = TRUE)
fit_e  <- sfm(y ~ x1 + x2, data = d, model_name = "NE",  keep_objective = TRUE)
vuong(fit_hn, fit_e)                            # neither assumed correct
influence_sfa(fit_hn)                           # who is driving this fit

spec_test_all(residuals(lm(y ~ x1 + x2, d)))    # before fitting anything
moment_range(residuals(lm(y ~ x1 + x2, d)))    # which pairs are even possible
sfma(y ~ x1 + x2, data = d, models = c("NHN", "NE", "NTN"))

## Is the assumed inefficiency distribution itself defensible? Hold the
## normality of the noise fixed and it implies a distribution for the composed
## error, so testing that is testing the assumption on u.
gof_test(fit_hn, data = d, B = 199)             # KS and Pearson chi-square

## Is there any inefficiency to speak of? The one-sided LR test is the one to
## quote: the Wald ratio and the naive LR test both have the wrong size here,
## because the null sits on a boundary.
inefficiency_test(fit_hn)
```

`spec_test()`, `gof_test()`, `lcsfm_homogeneity()` and `sfma()` default to a
**bootstrap** null rather than the published asymptotic one. Their help pages give the
measured size distortions behind that choice: in two cases the asymptotic limit
is badly mis-sized for the models this package fits, because it is stated for a
restricted specification the package does not impose.

## Working with results

`"sfareg"` objects support the usual generics:

```r
coef(fit); vcov(fit); logLik(fit); nobs(fit); AIC(fit); BIC(fit)
fitted(fit); residuals(fit); predict(fit, newdata = ...)
print(fit); summary(fit)
```

`fit$out` is the source of truth — a **`p x 3`** matrix, one **row per
parameter**, with columns `par`, `st_err` and `t-val`. Index it as
`fit$out[, "par"]`, never `fit$out["par", ]`. Its **row** names vary by model:
several report the `lambda = sigma_u/sigma_v`,
`sigma = sqrt(sigma_u^2 + sigma_v^2)` reparameterisation rather than the raw
scale parameters, so read the names rather than assuming a position.

`npsfm()` fits are the exception. They carry no `out` matrix and no standard
errors, so only `fitted()`, `residuals()`, `nobs()`, `print()` and `summary()`
apply; read the frontier, its gradients and the scale estimates off the returned
object (`$frontier`, `$frontier.grad`, `$sigma.u`, `$sigma.v`).

## Simulating data

`data_gen_cs()` and `data_gen_p()` generate cross-sectional and panel data with
known true parameters. Each returns a data frame with one response column per
model family (`y_pcs`, `y_pcs_z`, `y_pcs_r`, `y_tre_z`, ...), so a given
`model_name` is matched to the column generated under its own assumptions. These
generators are how the package's estimators are checked against known truth,
including `npsfm()`'s — `NPSFM_FLW` and `NPSFM_SVKZ` are registered in the
root-n convergence framework and both pass.

## Included data

| Data set | Description |
|---|---|
| `USUtilities` | Panel of US investor-owned fossil-fuel steam electric utilities, 1986-1999 |
| `FinnishElec` | Cross-section of Finnish electricity distribution firms, averaged over a four-year regulatory period |
| `Indian` | Panel of 14 paddy farmers in Aurepalle, India, 1975-76 to 1984-85 |
| `panel89` | Cross-section of US commercial banks, 1989 (Kumbhakar, Parmeter and Tsionas 2013) |

## Citation

```r
citation("sfa")
```

> Bernstein, D. H., Parmeter, C. F., and Stead, A. D. (2026). *Stochastic
> Frontier Analysis: The sfa Package.* Working Paper.

## License

GPL (>= 2). See [`LICENSE.md`](https://github.com/davidhbernstein/sfa/blob/main/LICENSE.md).
