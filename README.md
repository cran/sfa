# sfa

<!-- badges: start -->
[![R-CMD-check](https://github.com/davidhbernstein/sfa/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/davidhbernstein/sfa/actions/workflows/R-CMD-check.yaml)
[![CRAN status](https://www.r-pkg.org/badges/version/sfa)](https://CRAN.R-project.org/package=sfa)
[![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/sfa)](https://CRAN.R-project.org/package=sfa)
[![License: GPL v2](https://img.shields.io/badge/License-GPL%20v2-blue.svg)](https://cran.r-project.org/web/licenses/GPL-2)
<!-- badges: end -->

**Stochastic frontier analysis in R.** A single, consistent interface to a wide
range of cross-sectional, panel, zero-inefficiency, two-tier and nonparametric
stochastic frontier models, with a common formula syntax for modelling the
variance of each error component and a common `"sfareg"` result object that
works with the standard R modelling generics.

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

## The five entry points

| Function | Fits | Estimators |
|---|---|---|
| `sfm()` | Cross-sectional frontiers | 14 |
| `psfm()` | Panel frontiers | 15 |
| `zsfm()` | Zero-inefficiency (latent-class) frontiers | 2 |
| `ttsfm()` | Two-tier frontiers | 3 |
| `npsfm()` | Nonparametric frontiers | 5 |

The first four return an object of class `"sfareg"`. `npsfm()` returns
`"npsfareg"` instead — a kernel-estimated frontier has no parameter vector with
standard errors, so `coef()`, `vcov()` and `logLik()` would have nothing to
return.

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

`sfm()` also offers `estimator = "cols"` — corrected OLS (Olson, Schmidt and
Waldman 1980), closed-form and deterministic, for `NHN`, `NE` and `NG` — and
robust divergence-based alternatives to MLE via `robust = "mlqe" | "psi" | "mdpd"`
for `NHN`.

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
| `SSFE` | Schmidt and Sickles fixed effects |
| `PL80` | Pitt and Lee (1980), time-invariant |
| `BC92` | Battese and Coelli (1992) time decay |
| `K1990`, `K1990modified` | Kumbhakar (1990) time patterns |

`GTRE_SEQ1`, `GTRE_SEQ2` and `SSFE` are not maximum likelihood, so `logLik()`
(and hence `AIC()`/`BIC()`) returns `NA` for them.

`psfm_bootstrap()` provides a parametric bootstrap for GTRE-family fits,
parallelised over cores.

### `zsfm()` — zero inefficiency

`ZISF` and `ZISF_Z`: a mixture of a fully efficient regime and an inefficient
frontier regime, with the regime probability optionally parameterised by
covariates (`ZISF_Z`).

### `ttsfm()` — two tier

`TTNE` (normal–exponential–exponential), `TTHN` (normal–half normal–half
normal), and `TTNLS` (nonlinear least squares, no distributional assumption
beyond the means of the two one-sided components).

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
and `SZ` additionally needs **Benchmarking**. Both are in `Suggests`, not
`Imports`, so they are only required if you actually call `npsfm()`:

```r
install.packages(c("np", "Benchmarking"))
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

## Working with results

`"sfareg"` objects support the usual generics:

```r
coef(fit); vcov(fit); logLik(fit); nobs(fit); AIC(fit); BIC(fit)
fitted(fit); residuals(fit); predict(fit, newdata = ...)
print(fit); summary(fit)
```

`fit$out` is the source of truth — a `3 x p` matrix of estimates, standard
errors and *t*-values. Its column names vary by model: several models report the
`lambda = sigma_u/sigma_v`, `sigma = sqrt(sigma_u^2 + sigma_v^2)`
reparameterisation rather than the raw scale parameters.

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
