\name{psfm}
\alias{psfm}
\title{psfm}
\description{Function to implement various panel data stochastic frontier estimators}
\usage{psfm(formula, model_name = c("TRE_Z", "GTRE_Z", "TRE",
                    "GTRE", "GTRE_FML", "TFE", "TFE_WMLE", "FD",
                    "GTRE_SEQ1", "GTRE_SEQ2", "SSFE",
                    "PL80", "PL80_MVTN", "BC92", "K1990", "K1990modified",
                    "CSS", "LS", "SSRE", "SSCRE", "KSS"), data,
                    maxit.bobyqa = 5000, maxit.nlminb = 500,
                    maxit.psoptim = 100, maxit.optim = 1000,
                    REPORT = 1, trace = 3, pgtol = 0, individual,
                    halton_num = NULL, start_val = FALSE, gamma = FALSE,
                    PSopt = FALSE, optHessian = TRUE, inefdec = TRUE, Method = "L-BFGS-B",
                    verbose = FALSE, rand.gtre = NULL, rand.psoptim = NULL,
                    OPG_calc = FALSE,
                    estimator = c("fiml", "sml", "seq1", "seq2"),
                    collinear_action = c("start_only", "error", "warn_drop"),
                    time = NULL, tfe_lambda_max = 100,
                    keep_objective = FALSE, kss_L = "auto",
                    kss_smooth = "auto", kss_L_max = 7L, mundlak = NULL)}
\arguments{
  \item{formula}{a symbolic description for the model to be estimated}
  \item{model_name}{Which panel estimator to fit. All are described in \sQuote{Details}; in brief:
\describe{
  \item{\code{"TRE"}, \code{"TRE_Z"}}{Greene's (2005) true random effects, by simulated maximum likelihood. \code{_Z} adds determinants of \eqn{\sigma_u}.}
  \item{\code{"GTRE"}, \code{"GTRE_Z"}}{the four-component model separating persistent from transient inefficiency, by simulated maximum likelihood. \code{_Z} adds determinants.}
  \item{\code{"GTRE_FML"}}{the same four-component model by full information maximum likelihood through its closed-skew-normal representation; deterministic, but requires a balanced panel.}
  \item{\code{"GTRE_SEQ1"}, \code{"GTRE_SEQ2"}}{two-step/moment-based estimators of the four-component model. Not maximum likelihood -- no \code{logLik()}.}
  \item{\code{"TFE"}}{Greene's (2005) true fixed effects. \bold{Changed meaning in 1.1.3} -- see \sQuote{Details}.}
  \item{\code{"TFE_WMLE"}}{Chen, Schmidt and Wang's (2014) within maximum-likelihood estimator.}
  \item{\code{"FD"}}{the first-difference estimator of Wang and Ho (2010).}
  \item{\code{"SSFE"}}{Schmidt and Sickles' (1984) fixed-effects estimator, inefficiency read off the firm effects. Not maximum likelihood -- no \code{logLik()}.}
  \item{\code{"CSS"}}{Cornwell, Schmidt and Sickles (1990): a firm-specific quadratic in time, so inefficiency varies over time and firms may cross over. Not maximum likelihood -- no \code{logLik()}.}
  \item{\code{"LS"}}{Lee and Schmidt (1993): one common temporal pattern scaled per firm, \eqn{\alpha_{it} = \delta_t \alpha_i}. Time-varying but no crossover. Not maximum likelihood -- no \code{logLik()}.}
  \item{\code{"SSRE"}}{Schmidt and Sickles (1984) random effects: the GLS estimator, with inefficiency read off the firm effects' BLUPs. More efficient than \code{"SSFE"} and identifies time-invariant regressors, at the cost of assuming the effects are uncorrelated with the regressors. Not maximum likelihood -- no \code{logLik()}.}
  \item{\code{"SSCRE"}}{Correlated random effects: \code{"SSRE"} augmented with the within-firm means of the time-varying regressors (Mundlak 1978), so that correlation between effects and regressors is modelled rather than assumed away. Not maximum likelihood -- no \code{logLik()}.}
  \item{\code{"KSS"}}{Kneip, Sickles and Song (2012): the firm effect is a smooth function of time lying in an \eqn{L}-dimensional space whose basis is estimated from the data. Nests \code{"SSFE"}, \code{"LS"} and \code{"CSS"}. Balanced panels only. Not maximum likelihood -- no \code{logLik()}.}
  \item{\code{"PL80"}}{Pitt and Lee (1980): time-invariant inefficiency.}
  \item{\code{"BC92"}}{Battese and Coelli (1992): exponential time decay.}
  \item{\code{"PL80_MVTN"}}{Pitt and Lee's (1981) Model III: inefficiency varies over time and is CORRELATED within a firm, through a multivariate normal truncated to the negative orthant. This is the likelihood from their Appendix 2, which they derived but did not use, calling the T-dimensional integrals intractable. Balanced panels only. See \code{\link{PL80_MVTN}}.}
  \item{\code{"K1990"}, \code{"K1990modified"}}{Kumbhakar (1990) time patterns for inefficiency.}
} Matching ignores case, so \code{"nhn"}, \code{"NHN"} and \code{"Nhn"} are the same choice.}
  \item{data}{A panel data set. This may be an ordinary \code{data.frame} (or tibble/data.table) -- the panel index is constructed internally from the \code{individual} and \code{time} arguments -- or a \code{plm::pdata.frame} that already carries an index, which is used as-is.}
  \item{maxit.bobyqa}{Maximum number of iterations for the bobyqa optimization routine}
  \item{maxit.nlminb}{Maximum number of iterations for the nlminb optimization routine.}
  \item{maxit.psoptim}{Maximum number of iterations for the psoptim optimization routine}
  \item{maxit.optim}{Maximum number of iterations for the optim optimization routine}
  \item{REPORT}{reporting parameter}
  \item{trace}{Integer. Tracing level passed through to the optimizer; larger values print more.}
  \item{pgtol}{Numeric. Projected-gradient tolerance passed to \code{optim()}'s \code{"L-BFGS-B"} method. \code{0} uses the optimizer's own default.}
  \item{individual}{individual unit in the regression model}
  \item{halton_num}{number of Halton draws to use in SML models}
  \item{start_val}{starting value (optional)}
  \item{gamma}{Logical. Report the variance parameters in the \eqn{\gamma = \sigma_u^2/(\sigma_u^2 + \sigma_v^2)} parameterization instead of the default. Defaults to \code{FALSE}.}
  \item{PSopt}{use psoptim optimization routine (T or F)}
  \item{optHessian}{Logical. Should a numerically differentiated Hessian matrix be returned while using the optim routine? (for optim routine)}
  \item{inefdec}{Production or cost function}
  \item{Method}{The method to be used for optim. See 'Details' within optim.}
  \item{verbose}{Logical. Print optimization progress messages? Default is \code{FALSE}.}
  \item{rand.psoptim}{Integer. Seed for replication of psoptim.  Default to \code{NULL}.}
  \item{rand.gtre}{Integer. Seed for replication of the gtre model.  Default to \code{NULL}.}
  \item{OPG_calc}{Logical. Should the outer product of gradients be calculated for inference.  Default to \code{FALSE}.}
  \item{estimator}{How to estimate the four-component GTRE model, in the same spirit as \code{\link{sfm}}'s \code{estimator}: one model, several routes to it. Applies to \code{model_name = "GTRE"} only, and is ignored with a warning elsewhere.
\describe{
  \item{\code{"fiml"}}{\bold{the default.} Full information maximum likelihood through the closed-skew-normal representation. Deterministic -- nothing is integrated by Monte Carlo, so repeated fits to the same data agree exactly. Requires a \bold{balanced} panel.}
  \item{\code{"sml"}}{Simulated maximum likelihood over Halton draws. Handles unbalanced panels. This is what \code{model_name = "GTRE"} meant through version 1.1.3.}
  \item{\code{"seq1"}, \code{"seq2"}}{The two-step moment-based estimators. Not maximum likelihood, so they carry no \code{logLik()}.}
}
The older names remain and are unchanged: \code{"GTRE_FML"}, \code{"GTRE_SEQ1"} and \code{"GTRE_SEQ2"} select the same three routes directly.

On an \bold{unbalanced} panel \code{"fiml"} cannot be fitted. If you simply took the default, \code{psfm()} warns and falls back to \code{"sml"} -- erroring would make \code{"GTRE"} unusable by default on a whole class of data. If you asked for \code{"fiml"} explicitly it errors instead, rather than silently fitting a different estimator than the one you named.}
  \item{collinear_action}{One of \code{"start_only"} (default), \code{"error"}, or \code{"warn_drop"}. Controls what happens when the \emph{between-individual} design matrix used to build starting values is rank deficient. This is checked separately from the pooled design because the two can differ: time dummies, for example, are routinely estimable in a pooled specification while collapsing onto the intercept once averaged within each unit, which is exactly the situation that makes the internal random-effects initialization fail. \code{"start_only"} keeps the model you asked for and removes the offending columns from the starting-value regression only, taking their starting values from a pooled OLS fit; \code{"error"} stops before estimation and names the offending columns; \code{"warn_drop"} removes those terms from the estimated model. All three report what was found. Only relevant to the models initialized by a random-effects regression (\code{GTRE}, \code{TRE} and their \code{_Z}/\code{SEQ} variants).}
  \item{time}{Character. Name of a time-period column in \code{data}, used by the error-components frontier models (\code{"PL80"}, \code{"BC92"}, \code{"K1990"}, \code{"K1990modified"}). If \code{NULL} (the default), a sequential within-individual time index is constructed from row order, which assumes \code{data} is already sorted chronologically within each individual.}
  \item{tfe_lambda_max}{Numeric. Upper bound placed on \eqn{\lambda = \sigma_u/\sigma_v} when fitting \code{model_name = "TFE"}; ignored by every other model. Defaults to \code{100}, i.e. \eqn{\sigma_v} is not allowed below about 1\% of \eqn{\sigma}. This is a guard against a genuine degeneracy of Greene's likelihood rather than a convergence tolerance -- see \sQuote{Details}. A fit that pins at the bound warns and should not be read as an estimate of \eqn{\lambda}.}
\item{keep_objective}{Retain the likelihood closure on the fitted object as \code{$objective}, for \code{"GTRE"}, \code{"TRE"}, \code{"GTRE_Z"} and \code{"TRE_Z"} -- the simulated-ML panel models. This lets the log-likelihood be evaluated away from the optimum, which is what profiling a single parameter requires; \code{\link{sfa_diagnostics}} uses it, and it is how the \code{sigh} convergence question was localized (see \sQuote{Details}). \code{FALSE} by default, since the closure captures the data and the draws and so is not small.}
\item{kss_L}{\code{"KSS"} only. The dimension \eqn{L} of the factor space. \code{"auto"} (the default) selects it by the Bai and Ng (2002) \eqn{IC_{p2}} criterion; give an integer to fix it.}
\item{kss_smooth}{\code{"KSS"} only. The smoothing penalty applied to each firm's residual trajectory before the basis is extracted. \code{"auto"} (the default) chooses it by generalized cross-validation; \code{0} turns smoothing off, which reduces the model to an unrestricted factor model.}
\item{kss_L_max}{\code{"KSS"} only. Largest \eqn{L} the automatic criterion may choose. The criterion is capped internally at \eqn{\lfloor T/2 \rfloor} as well, for the reason given in \sQuote{Details}, and warns if it selects that cap. An explicit \code{kss_L} is not subject to the \eqn{\lfloor T/2 \rfloor} cap, only to \eqn{\min(T-1, n-1)}.}
\item{mundlak}{Mundlak adjustment terms: a formula naming the regressors whose \bold{within-firm means} should be added to the frontier, or \code{~ .} for every frontier regressor. \code{NULL} (the default) adds nothing. Only for models that treat the firm effect as random; \code{"TFE"}, \code{"TFE_WMLE"}, \code{"FD"} and \code{"SSFE"} refuse it, because they have already removed the effect. See \sQuote{Details}.}
}
\details{The generalized true random effects model (GTRE, 4-component model) and true random effects models (TRE)  are both estimated by simulated maximum likelihood based on the paper by the Fillipini and Greene (2016, JPA). The TRE_Z and GTRE_Z allow for modeling the u-component of the GTRE and TRE with determinants of inefficiency. The first-difference estimator (FD) of Wang and Ho (2010, JoE) is also available.

\bold{Choosing how GTRE is estimated.} Since version 1.1.4 the four routes to the four-component model are selected with \code{estimator} rather than by picking among four \code{model_name} values that looked like four different models. \code{model_name = "GTRE"} now defaults to \code{estimator = "fiml"}; through 1.1.3 it fitted the simulated-ML estimator, which is \code{estimator = "sml"}. Scripts that pass \code{"GTRE"} therefore get a different estimator than they did and are warned once per call. The explicit names \code{"GTRE_FML"}, \code{"GTRE_SEQ1"} and \code{"GTRE_SEQ2"} are unaffected.

\bold{Two ways of estimating the same four-component model.} \code{"GTRE"} and \code{"GTRE_FML"} both fit
\deqn{y_{it} = \beta_0 + x_{it}'\beta + r_i - h_i + v_{it} - u_{it},}
with \eqn{r_i, v_{it}} two-sided and \eqn{h_i, u_{it}} one-sided (half-normal), separating \emph{persistent} inefficiency \eqn{h_i} from \emph{transient} inefficiency \eqn{u_{it}}. They differ in how the likelihood is obtained.

\code{"GTRE"} uses simulated maximum likelihood over Halton draws. \code{"GTRE_FML"} uses \emph{full information} maximum likelihood through the model's closed-skew-normal representation: with \eqn{A = -[1_T, I_T]}, \eqn{V = \mathrm{diag}(\sigma_h^2, \sigma_u^2 I_T)} and \eqn{\Sigma = \sigma_v^2 I_T + \sigma_r^2 1_T 1_T'}, the firm-level density is available in closed form and nothing is integrated by Monte Carlo. The likelihood is therefore deterministic, and repeated fits to the same data give identical answers.

The \eqn{(T+1)}-dimensional normal CDF this requires would normally be the bottleneck, but its covariance is diagonal-plus-rank-one for this model, so it collapses to a one-dimensional Gauss-Hermite integral -- the same reduction used for \code{"TFE"}.

\code{"GTRE_FML"} reports \bold{raw standard deviations} \code{sigr}, \code{sigv}, \code{sigh}, \code{sigu} (individual effect, noise, persistent inefficiency, transient inefficiency), \emph{not} the \eqn{\lambda}/\eqn{\sigma} reparameterization \code{"GTRE"} uses, and it requires a \bold{balanced} panel -- the closed-skew-normal pieces are built at a single \eqn{T}. It errors rather than silently mis-estimating if the panel is unbalanced; use \code{"GTRE"} in that case.

\bold{Starting values for \code{"GTRE_FML"}.} The closed-skew-normal likelihood carries a boundary optimum at \eqn{\sigma_h = 0}, where the model collapses to \code{"TRE"} and the intercept absorbs the missing \eqn{E[h] = \sigma_h\sqrt{2/\pi}}. \emph{Two different things can put a fit there, and the likelihood tells them apart.} A fit whose boundary solution has a \emph{lower} likelihood than the true parameter vector reached it by getting stuck, and is an optimizer failure; one whose boundary solution has a \emph{higher} likelihood reached it because that is where the maximum is, and is the correct answer for that sample (see \sQuote{When a persistent scale collapses to zero} below). The start search addresses only the first kind. Following Colombi (2010) and Colombi, Martini and Vittadini (2011), who recommend seeding this likelihood from the two-step moment estimates, \code{psfm()} builds two candidate starts -- the random-effects decomposition shared with \code{"GTRE"}, and the two-step estimator reported by \code{"GTRE_SEQ2"} -- evaluates the likelihood at each, and begins from the better one. Which was chosen, and the log-likelihood at both, is returned in \code{start_search}. Supplying \code{start_val} overrides the search.

Note that \eqn{\sigma_r} and \eqn{\sigma_h} enter the density only through \eqn{\sigma_r^2 + \sigma_h^2}; their separation rests entirely on the skewness term. That pair is therefore the least precisely estimated part of the model even when the fit is at the global optimum, which is consistent with the confounding between latent heterogeneity and persistent inefficiency reported in the source papers.

\bold{\code{"GTRE_SEQ1"} and \code{"GTRE_SEQ2"} are large-T estimators.} Both take
the random-effects decomposition from \code{plm} and hand its two output series to a
second stage -- \code{"GTRE_SEQ1"} fits an intercept-only normal--half-normal likelihood
to each, \code{"GTRE_SEQ2"} inverts their second and third moments. The second stage
treats those series as if they were draws from \eqn{v_{it}-u_{it}} and \eqn{r_i-h_i},
and at finite \eqn{T} they are not: \code{alpha_hat} is the shrunken predictor
\eqn{c\,(\alpha_i+\bar e_i)} with \eqn{c = T\sigma_\alpha^2/(T\sigma_\alpha^2+\sigma_\epsilon^2)},
and \code{eps_hat} is the quasi-demeaned residual
\eqn{(1-\theta)\alpha_i + e_{it} - \theta\bar e_i} with
\eqn{\theta = 1-\sigma_\epsilon/\sqrt{T\sigma_\alpha^2+\sigma_\epsilon^2}}.

Relative to the latent draws this attenuates the third central moment on the
\eqn{\epsilon} side by \eqn{(1-\theta/T)^3+(T-1)(-\theta/T)^3} (about 0.81 at
\eqn{T=10}) and scales the \eqn{\alpha} side by \eqn{c^2} and \eqn{c^3}. None of these
factors tends to 1 as \eqn{N} grows with \eqn{T} fixed, so \bold{both estimators
converge to the wrong constants and are inconsistent for fixed \eqn{T}}; each call
warns accordingly. On the package's own convergence design (\eqn{T=10},
\eqn{\sigma_h=0.8}, \eqn{\sigma_r=0.4}) the estimates settle at
\eqn{\gamma_{uv}\approx0.87} against 0.917 and \eqn{\sigma^2_{hr}\approx0.68} against
0.80, and the log MSE-on-log \eqn{N} slopes are near zero rather than \eqn{-1}. The
same routines applied to the latent draws themselves recover the truth at the
\eqn{\sqrt{N}} rate, which locates the problem in the first-stage decomposition rather
than in either second stage. Raising \eqn{T} removes the bias; raising \eqn{N} only
sharpens it. For a consistent alternative use \code{model_name = "GTRE"}, which
estimates all four variance components jointly.

The standard errors these two report are delta-method errors built from the
\emph{iid} moment variances, so they describe sampling variability around each
estimator's probability limit, not around the truth. They are also optimistic on the
\eqn{\epsilon} side, whose \eqn{n} residuals are correlated within firm.

\bold{The two true-fixed-effects estimators, and a change of meaning in \code{"TFE"}.} Both \code{"TFE"} and \code{"TFE_WMLE"} fit the same model,
\deqn{y_{it} = \alpha_i + x_{it}'\beta + v_{it} - u_{it},\quad v \sim N(0,\sigma_v^2),\ u \sim N^+(0,\sigma_u^2),}
with \eqn{\alpha_i} an unrestricted individual effect; they differ in how they estimate it.

\code{"TFE"} is Greene's (2005, JoE/JPA) true fixed effects estimator: the ordinary composed-error likelihood with the \eqn{N} individual intercepts estimated jointly with \eqn{(\lambda, \sigma, \beta)}. This package concentrates the \eqn{\alpha_i} out at every likelihood evaluation (each is the unique root of a strictly decreasing score, solved by safeguarded Newton to machine precision), so the reported parameter vector, its standard errors and the \code{out} layout are the same 2 + K as \code{"TFE_WMLE"}. It is the estimator subject to the \emph{incidental parameters problem}: with \eqn{T} fixed, \eqn{\sigma_u} is biased upward and efficiency scores correspondingly pessimistic.

\code{"TFE_WMLE"} is Chen, Schmidt and Wang's (2014, JoE) within maximum-likelihood estimator, which applies maximum likelihood to the deviations from individual means. The \eqn{\alpha_i} drop out of the likelihood entirely, so the estimator is free of the incidental parameters problem; the firm effects are recovered afterwards from the mean residual with the usual \eqn{\sqrt{2/\pi}\,\sigma_u} adjustment.

The two are \emph{not} nested and their log-likelihoods are not on the same scale -- \code{"TFE"} evaluates the likelihood of the data, \code{"TFE_WMLE"} that of the within-transformed deviations -- so \code{logLik()}, \code{AIC()} and \code{BIC()} must not be used to choose between them.

\bold{Through sfa 1.1.2, \code{model_name = "TFE"} selected the Chen-Schmidt-Wang estimator.} It now selects Greene's, matching standard usage; the Chen-Schmidt-Wang estimator is \code{"TFE_WMLE"}. Existing scripts that pass \code{"TFE"} therefore get a different estimator than they did, and are warned once per call.

\bold{Degeneracy of the Greene likelihood.} Because \eqn{\alpha_i} is unrestricted, Greene's likelihood always has a supremum on the \eqn{\sigma_v \to 0} boundary: setting \eqn{\alpha_i = \max_t (y_{it} - x_{it}'\beta)} makes every composed error weakly negative, \eqn{\Phi(-\lambda e/\sigma) \to 1}, and what remains is the deterministic-frontier likelihood, which is finite and attained. On simulated panels this supremum is reached along a path that first dips below the interior local maximum and then rises past it, so an unconstrained search that escapes the interior basin returns \eqn{\sigma_v \approx 0} together with meaningless standard errors; when \eqn{\lambda} is large and \eqn{T} small, no interior maximum need exist at all. \code{tfe_lambda_max} bounds the search accordingly, and a fit that pins at the bound warns. This is a property of the estimator, not of the implementation, and is one of the motivations for \code{"TFE_WMLE"}.

\bold{Both problems are small-\eqn{T} phenomena.} Simulations in which \eqn{T} grows with \eqn{N} fixed at 100 (truth \eqn{\lambda = 3.333}) recover \eqn{\hat\lambda = } 88.5, 4.8, 3.6, 3.6, 3.5 at \eqn{T = } 10, 20, 30, 40, 50: the fit is pinned at the bound at \eqn{T = 10} and has converged by \eqn{T = 50}. The degeneracy and the incidental-parameters bias dissolve together as \eqn{T} grows, so \code{tfe_lambda_max} matters on short panels and is irrelevant on long ones. \code{"TFE_WMLE"}, having no incidental parameters, is unaffected by \eqn{T} either way.

\code{"SSFE"} is the classical Schmidt and Sickles (1984, JBES) fixed-effects (LSDV) estimator: a standard within/fixed-effects panel regression, with technical inefficiency read directly off the estimated firm effects (each firm's distance from the best-performing firm in the sample) rather than estimated by maximum likelihood. It requires no numerical optimization and returns no log-likelihood -- \code{logLik()}/\code{AIC()}/\code{BIC()} are not defined for it (see \code{\link[sfa]{sfareg-methods}}). Note this is a different model from \code{"TFE"}/\code{"TFE_WMLE"} above, despite the similar names: both of those are maximum-likelihood estimators that allow a separate stochastic noise term on top of the fixed effects, while \code{"SSFE"} attributes all within-firm deviation to inefficiency directly.

\code{"K1990"} and \code{"K1990modified"} (Kumbhakar, 1990) extend the same error-components frontier with richer time paths for inefficiency: \code{"K1990"} uses \eqn{B_{it} = (1 + \exp(b t + c t^2))^{-1}} and \code{"K1990modified"} uses \eqn{B_{it} = 1 + d(t - T_i) + e(t - T_i)^2}, against \code{"BC92"}'s \eqn{\exp(-\eta(t - T_{ref}))} and \code{"PL80"}'s time-invariant \eqn{B_{it} = 1}. All four share one likelihood, differing only in \eqn{B_{it}}. Note that \code{"K1990"}'s \code{b} and \code{c} are weakly identified -- quite different pairs produce near-identical fitted time paths -- so the fitted \eqn{B_{it}} path is more interpretable than either coefficient on its own.

\code{"PL80"} (time-invariant panel inefficiency, Pitt and Lee, 1980, JoE) and \code{"BC92"} (time-varying decay, Battese and Coelli, 1992, JPA) are both natively estimated maximum-likelihood implementations of the same "error components frontier": y_it = x_it'beta + v_it - B_it*u_i (production; sign flips for cost), with v_it iid normal, u_i iid half-normal and constant across time within a firm, and B_it = 1 for \code{"PL80"} or B_it = exp(-eta*(t - Tref)) for \code{"BC92"} (Tref being the last time period in the whole panel). The closed-form log-likelihood is obtained by integrating u_i out of the joint density of (epsilon_i, u_i); it nests the ordinary cross-sectional normal-half-normal likelihood when every firm has a single observation. Verified against \code{frontier::sfa()} (matching coefficients, log-likelihood, and \code{frontier::efficiencies()}-style predicted efficiencies to numerical precision across balanced/unbalanced panels and production/cost specifications) before that dependency was removed. The \code{time} argument above supplies each observation's period for \code{"BC92"}'s decay term; if left \code{NULL}, a sequential within-individual time index is constructed from row order, which assumes \code{data} is already sorted chronologically within each individual.

Note on default iteration counts: \code{psfm()}'s defaults were raised in version 1.1.3, from 100/10/10 for \code{maxit.bobyqa}/\code{maxit.psoptim}/\code{maxit.optim} to 5000/100/1000, and \code{maxit.nlminb} (previously hard-coded) is now an argument. The old values were binding rather than merely economical: the \code{K1990}/\code{K1990modified} branch already carried a note that 100 bobyqa evaluations left its seven-parameter fits several log-likelihood units short of the optimum purely on the iteration cap, and \code{GTRE_FML} at N = 500, T = 10 roughly halves its root-mean-square error against known true values when the caps are lifted, at about 1.5 times the run time. They remain below \code{sfm()}'s defaults because several \code{psfm()} models (\code{GTRE}, \code{GTRE_Z}, \code{TRE}, \code{TRE_Z}) are estimated by simulated maximum likelihood over Halton draws, so each objective evaluation is far more expensive than the closed-form likelihoods used elsewhere. Lower them explicitly if a fast approximate call is wanted.}
\references{Fillipini and Greene (2016, JPA); Colombi, Kumbhakar, Martini and Vittadini (2014, JPA) on the closed-skew-normal four-component likelihood, with Colombi (2010) and Colombi, Martini and Vittadini (2011) on the two-step estimator used to seed it; Greene (2005, JoE and JPA); Wang and Ho (2010, JoE); Chen, Schmidt and Wang (2014, JoE); Kumbhakar (1990, JoE); Pitt and Lee (1980, JoE); Battese and Coelli (1992, JPA); Battese and Coelli (1988, JoE); Jondrow, Lovell, Materov and Schmidt (1982, JoE); Schmidt and Sickles (1984, JBES)}
\author{David H. Bernstein}
\note{Standard errors require optHessian set to TRUE}

\seealso{\code{\link{sfm}} for cross-sectional models, \code{\link{psfm_bootstrap}} for bootstrap standard errors, \code{\link{data_gen_p}} for simulating panel data with known true parameters, and \code{\link{sfareg-methods}} for extracting results.}
\value{
An object of class \code{"sfareg"} containing components that vary by model. All models return:

\item{out}{A matrix with parameter estimates, standard errors, and t-values.}
\item{opt}{A list containing the optimization results from the final optimization procedure (not returned for GTRE_SEQ1 and GTRE_SEQ2).}
\item{total_time}{The total computation time for model estimation.}
\item{start_v}{The starting values used in the optimization (not returned for GTRE_SEQ1 and GTRE_SEQ2).}
\item{model_name}{The name of the panel stochastic frontier model estimated.}
\item{formula}{The formula used in the model specification.}
\item{coefficients}{A vector of estimated parameters.}
\item{std.errors}{A vector of standard errors for the estimated parameters (NA if \code{optHessian = FALSE}).}
\item{t.values}{A vector of t-values for the estimated parameters (NA if \code{optHessian = FALSE}).}
\item{call}{The matched call.}
\item{data}{The data used in estimation.}

\bold{Additional model-specific components:}

\bold{For GTRE_FML model:}
\item{start_search}{A list recording the starting-value search: \code{n_tried}, the log-likelihood at each candidate (\code{loglik}), and which was \code{chosen} (\code{"random-effects"} or \code{"two-step"}). \code{NULL} when \code{start_val} was supplied.}
\item{U}{Predicted transient (time-varying) technical efficiency for each observation.}
\item{H}{Predicted persistent (time-invariant) technical efficiency for each individual, named by individual.}

\bold{For GTRE and GTRE_Z models:}
\item{H}{Predicted time-invariant technical efficiency for each individual.}

\bold{For GTRE_Z model:}
\item{out_opg}{OPG (outer product of gradients) standard errors and t-values along with parameters.}
\item{out_sandwich}{sandwich standard errors and t-values along with parameters.}

\bold{For GTRE, GTRE_Z, TRE and TRE_Z models:}
\item{U}{Predicted time-varying technical efficiency for each observation.}

\bold{For TFE and TFE_WMLE models:}
\item{r_hat_m}{Estimated individual-specific fixed effects, named by individual. For \code{"TFE"} these are the maximum-likelihood estimates of \eqn{\alpha_i} directly; for \code{"TFE_WMLE"} they are recovered after estimation from the mean residual plus \eqn{\sqrt{2/\pi}\,\sigma_u}, since the within likelihood does not identify them.}
\item{exp_u_hat}{Predicted technical efficiency, \eqn{E[\exp(-u)|\epsilon]} (Battese and Coelli, 1988).}

\bold{For TFE model only:}
\item{u_hat}{Predicted inefficiency, \eqn{E[u|\epsilon]} (Jondrow et al., 1982).}

\bold{For FD model:}
\item{u_hat}{Predicted technical efficiency in levels.}
\item{h_hat}{Estimated z heterogeneity function values.}
\item{exp_u_hat}{Predicted technical efficiency.}

\bold{For GTRE_SEQ1 and GTRE_SEQ2 models:}
\item{other_parms}{A matrix of additional parameters (lambda, sigma, beta_0 for SEQ1; sigma_u, sigma_v, sigma_h, sigma_r, lambda, sigma for SEQ2).}

\bold{When a persistent scale collapses to zero.} \code{GTRE} carries two persistent components -- \eqn{\sigma_h}, persistent inefficiency, and \eqn{\sigma_r}, the firm random effect -- and the likelihood cannot always separate them in a given sample. When it cannot, it merges them: one scale goes to zero and the other absorbs its variation, coming back inflated.

This is common rather than exceptional. On simulated data at \eqn{T = 6} with both components genuinely present (\eqn{\sigma_r = 0.20}, \eqn{\sigma_h = 0.40}), one of the two collapsed in \bold{37\% of 87 replications}, and the commoner direction was \eqn{\sigma_r \to 0} rather than \eqn{\sigma_h \to 0}. The share with both interior does improve with \eqn{N} (15, 19 and 21 of about 29 at \eqn{N = 50, 100, 200}) but slowly.

\bold{A collapse is usually the right answer, not a numerical failure.} On one such replication -- whose realized \eqn{\sigma_h} scale was 0.3695, so persistent inefficiency was certainly present -- the boundary solution had a log-likelihood \emph{higher than the true parameter vector} by 3.66, a likelihood ratio of about 39. This is the panel counterpart of the cross-sectional wrong-skewness result (Waldman 1982), where \eqn{\sigma_u = 0} is the correct maximum likelihood estimate rather than an optimizer failure; see \code{\link{skewness_test}}. This is not in tension with the start search described above under \code{"GTRE_FML"}: a boundary fit that beats the truth on likelihood is the estimate, while one that loses to it is a trapped optimizer, and the search exists to remove the second so that what remains is the first.

\bold{What a collapse does to the rest of the fit, exactly.} One rule covers it: \emph{a collapse conserves the persistent variance and relocates the persistent mean into the intercept}. Those are two different destinations. \eqn{h} is half-normal, so it carries both a mean \eqn{E[h] = \sigma_h\sqrt{2/\pi}}{E[h] = sigma_h sqrt(2/pi)} and a variance \eqn{\mathrm{Var}(h) = \sigma_h^2(1 - 2/\pi)}{Var(h) = sigma_h^2 (1 - 2/pi)}, and the model has a separate place for each: the frontier intercept takes the mean, since \eqn{E[y] = \beta_0 - E[h] - E[u]}, and the surviving persistent scale takes the variance. The intuitive shortcut -- that the survivor comes back at \eqn{\sqrt{\sigma_h^2 + \sigma_r^2}}{sqrt(sigma_h^2 + sigma_r^2)} -- ignores that split and is wrong.

\tabular{lll}{
  \tab \bold{surviving scale} \tab \bold{intercept} \cr
  \eqn{\sigma_h \to 0}{sigma_h -> 0} \tab \eqn{\hat\sigma_r = \sqrt{\sigma_r^2 + \sigma_h^2(1-2/\pi)}}{sigma_r_hat = sqrt(sigma_r^2 + sigma_h^2 (1 - 2/pi))} \tab shifts by \eqn{-E[h]}{-E[h]} \cr
  \eqn{\sigma_r \to 0}{sigma_r -> 0} \tab \eqn{\hat\sigma_h = \sqrt{\sigma_h^2 + \sigma_r^2/(1-2/\pi)}}{sigma_h_hat = sqrt(sigma_h^2 + sigma_r^2 / (1 - 2/pi))} \tab shifts by \eqn{+\Delta E[h]}{+Delta E[h]}
}

At \eqn{\sigma_h = 0.40}, \eqn{\sigma_r = 0.20} that predicts 0.3133 and \eqn{-0.3192} in the first row, 0.5197 and \eqn{+0.0955} in the second, with nothing fitted. Over 1000 replications at each of five sample sizes the measurements are 0.2973, \eqn{-0.3035}, 0.5058 and \eqn{+0.0931} -- within 2.5\% to 5.1\%, and all four slightly below prediction, which is what selection produces rather than a wrong constant.

The practical consequence is that a collapsed fit does not merely mislabel the persistent split. It displaces the intercept by a known and substantial amount, so any level read off \code{coef()} is displaced with it.

\bold{What a collapse is telling you about your data.} It is wrong skewness, measured between firms rather than within them. Take the OLS residuals, average them within each firm -- that average is the persistent part of the error -- and look at their skew. Over 1000 replications, the samples whose between-firm skew came out positive (the wrong way for a production frontier) lost \eqn{\sigma_h} to the boundary \bold{69\%} of the time, against 3.6\% for the rest; the area under the ROC curve for that single statistic is \bold{0.955}. Not one wrong-skew sample collapsed \eqn{\sigma_r} instead. So a collapsed \eqn{\sigma_h} is not a numerical accident to be retried with different starting values: this sample carries no between-firm evidence of persistent inefficiency, and the estimate is reporting that.

Note that this is a reading of the sign, not a test. The firm mean of the residuals contains \eqn{-\bar u_i} as well as \eqn{-h_i}, and the mean of \eqn{T} half-normals is itself skewed the same way, so on short panels the between-firm skew leans negative even when \eqn{\sigma_h = 0} -- 82\% of the time at \eqn{T = 4}, 63\% at \eqn{T = 10}, 53\% at \eqn{T = 30}. \code{\link{skewness_test}} is therefore \emph{not} applicable to the between-firm residuals as it stands.

\bold{Whether firms or periods help depends on the estimator, and the two answers are opposite.} Crossed \eqn{(N, T)} designs on the same data-generating process, fitting the collapse rate on \eqn{\log N}{log N} and \eqn{\log T}{log T}:

\tabular{llll}{
  \bold{estimator} \tab \bold{fits} \tab \bold{log N} \tab \bold{log T} \cr
  \code{"fiml"} (the default) \tab 600 \tab \bold{-1.66} (0.18) \tab -0.93 (0.20) \cr
  \code{"sml"} \tab 251 \tab -0.68 (0.26) \tab \bold{-1.07} (0.27)
}

Under the default, firms buy about 1.8 times what periods do; under simulated ML it is the other way round, periods buying about 1.6 times what firms do. Restricting the \code{"fiml"} design to \eqn{N \le 200}{N <= 200}, cell for cell the same grid the \code{"sml"} one used, leaves it at \eqn{-1.64} \eqn{(0.22)}, so this is the estimator and not the design; the difference in the \eqn{N} coefficient is significant \eqn{(z = 2.8)} and the difference in the \eqn{T} coefficient is not.

The practical reading, for the estimator you get by default: \bold{add firms}. At \eqn{T = 4} the collapse rate falls 68\%, 36\%, 18\%, 6\% across \eqn{N = 50, 100, 200, 400}, and reaches 0\% at \eqn{N = 400, T = 16}. Under \code{"sml"} the same row is flat -- 60\%, 39\%, 41\% -- which is why a short panel really is a dead end there however many firms it has. If you must use \code{estimator = "sml"} (an unbalanced panel forces it), then length is what buys the split, and a four-period panel will not give it to you.

\bold{How to tell in advance whether your panel can support the split.} One
quantity governs it, and it is not \eqn{\sigma_h/\sigma_r}. What the likelihood has to do is read the persistent part out of firm means that the transient block contaminates, so write the firm-mean composite as \eqn{r_i - h_i - \bar u_i + \bar v_i} and compare its persistent variance with the transient variance left in it:

\deqn{S/N = \frac{\sigma_r^2 + \sigma_h^2(1 - 2/\pi)}{[\sigma_u^2(1 - 2/\pi) + \sigma_v^2]/T}}{S/N = [sigma_r^2 + sigma_h^2 (1 - 2/pi)] / ([sigma_u^2 (1 - 2/pi) + sigma_v^2] / T)}

Across six designs varying all four scales in different combinations (50 replications each, \eqn{N = 100}, \eqn{T = 6}), this orders the collapse rate \bold{perfectly} -- Spearman \eqn{-1.000} under both estimators, and a logistic coefficient on \eqn{\log(S/N)}{log(S/N)} of \eqn{-1.15} \eqn{(0.20)}:

\tabular{lllll}{
  \bold{S/N} \tab \bold{collapse} \tab \eqn{\sigma_h \to 0}{sigma_h -> 0} \tab \eqn{\sigma_r \to 0}{sigma_r -> 0} \tab \cr
  0.67 \tab 54\% \tab 0.14 \tab 0.40 \tab \cr
  0.90 \tab 50\% \tab 0.10 \tab 0.40 \tab \cr
  0.96 \tab 40\% \tab 0.21 \tab 0.19 \tab \cr
  1.30 \tab 31\% \tab 0.15 \tab 0.17 \tab \cr
  5.20 \tab 8\% \tab 0.06 \tab 0.02 \tab \cr
  11.69 \tab 6\% \tab 0.04 \tab 0.02 \tab
}

As a rule of thumb: below \eqn{S/N \approx 1}{S/N about 1} about half your fits will merge the two components, and it takes roughly \eqn{S/N > 5}{S/N > 5} to get under 10\%.

Two things follow that are easy to get backwards. First, \bold{raising \eqn{\sigma_h/\sigma_r} does not help} -- it only decides which scale dies. Going from a ratio of 2 to 4 by shrinking \eqn{\sigma_r} takes \eqn{\sigma_h \to 0}{sigma_h -> 0} down (0.15 to 0.10) and \eqn{\sigma_r \to 0}{sigma_r -> 0} up (0.17 to 0.40), and the net rate \emph{rises} from 31\% to 50\%, because shrinking \eqn{\sigma_r} removes persistent variance and lowers \eqn{S/N} from 1.30 to 0.90. Holding the ratio at 2 and scaling both components instead takes the rate from 31\% to 8\%. Second, \bold{more transient inefficiency hurts}: \eqn{\sigma_u} enters the denominator, so a noisier transient block contaminates the firm means and makes the persistent split harder.

Note what this quantity does and does not explain. \eqn{T} enters it directly, so it covers the contamination channel by which longer panels help. It says nothing about \eqn{N}, which helps through a different route -- sampling error in what is effectively a deconvolution from \eqn{N} firm-level draws -- which is why the two dimensions have to be measured separately, as above.

One limit on all of that, stated because it is easy to over-read. Those elasticities were measured at \eqn{\sigma_r = 0.20} against \eqn{\sigma_h = 0.40}, where persistent inefficiency is the larger of the two time-invariant components. On a design where the random effect dominates instead (\eqn{\sigma_r = 0.50}, \eqn{\sigma_h = 0.40}, \eqn{N = 100}, \eqn{T = 6}) the \emph{level} ordering reverses: \code{"fiml"} collapses on 33\% of replications and \code{"sml"} on 20\%. Which estimator handles the split better is not a fact about the estimators alone -- it depends on how much persistent inefficiency there is to find.

\code{psfm} therefore reports the condition rather than preventing it, in \code{$sigh_at_bound} and \code{$sigr_at_bound} and in a warning naming the surviving scale. Bounding either parameter away from zero would corrupt precisely the samples where the boundary is the answer. What the warning asks of you is to read the two persistent scales together and to treat the split between them as unidentified in that sample -- their \emph{sum} remains informative even where the division does not.

\bold{Profiling a fitted panel likelihood.} With \code{keep_objective = TRUE} the fitted object carries \code{$objective}, the same negative-summed-log-likelihood the optimizer minimised. Evaluating it away from the optimum is the only way to ask whether a parameter is weakly identified, as distinct from badly estimated:

\preformatted{fit <- psfm(y ~ x1 + x2, "GTRE", d, individual = "id",
            estimator = "sml", keep_objective = TRUE)
ih   <- which(names(coef(fit)) == "sigh")
grid <- seq(0.05, 1, by = 0.05)
ll   <- vapply(grid, function(v) \{
          p <- fit$opt$par; p[ih] <- v; -fit$objective(p) \}, numeric(1))
plot(grid, ll - max(ll), type = "l")}

A flat curve means no estimator could recover the parameter from these data; a sharp one means the difficulty lies elsewhere. For \code{GTRE}'s \code{sigh} the curve spans over 200 log units at \eqn{N = 200} and sharpens with \eqn{N}, which is what ruled out weak identification as the explanation for its convergence behaviour.

\bold{The Schmidt--Sickles family.} \code{"SSFE"}, \code{"SSRE"} and \code{"SSCRE"} are the within, GLS and correlated-random-effects members of one family: all three read inefficiency off a firm effect as distance from the best firm, and none assumes a distribution for \eqn{u}. The choice among them is the standard panel trade-off, with a frontier-specific twist. Random effects is the more efficient estimator and is the only one of the three that identifies coefficients on \emph{time-invariant} regressors, which the within transformation sweeps away -- but it requires the effects to be uncorrelated with the regressors, and in a frontier model that assumption says a firm's inefficiency is unrelated to its input choices, which is close to saying managers do not respond to their own productivity.

\code{"SSCRE"} is the answer to that. Adding the within-firm means of the time-varying regressors (Mundlak 1978) models the correlation instead of assuming it away; the slopes on the original regressors then equal the within estimates \emph{exactly} -- not approximately -- while time-invariant regressors remain identified. The coefficients on the added means are reported with a \code{.mean_} prefix, and testing them jointly against zero is the Mundlak form of the Hausman test: if they are indistinguishable from zero, \code{"SSRE"} was safe. A regressor already constant within firm has a group mean identical to itself, so it contributes no Mundlak term; those are exactly the regressors random effects exists to identify, and they are left alone.

\bold{Kneip, Sickles and Song (2012).} \code{"KSS"} writes the firm effect as \eqn{\alpha_{it} = \sum_{r=1}^{L}\theta_{ir} g_r(t)}, where both the basis functions \eqn{g_r} and the dimension \eqn{L} are estimated from the data rather than assumed. That makes explicit what the other time-varying estimators here are special cases of:

\tabular{lll}{
  \strong{model} \tab \strong{\eqn{L}} \tab \strong{basis} \cr
  \code{"SSFE"} \tab 1 \tab fixed, constant in \eqn{t} \cr
  \code{"LS"}   \tab 1 \tab free \cr
  \code{"CSS"}  \tab 3 \tab fixed to \eqn{\{1, t, t^2\}} \cr
  \code{"KSS"}  \tab estimated \tab estimated
}

Estimation follows the three steps of the original: cross-sectional centering by period (which imposes \eqn{\sum_i u_i(t) = 0} and removes any common function of time), then smoothing of each firm's residual trajectory followed by an eigendecomposition of their empirical covariance to give \eqn{\hat g_r}, then loadings by least squares and the effects rebuilt from them -- iterated with \eqn{\beta} until the residual sum of squares settles. Smoothing before the eigendecomposition is what makes the estimated basis a set of smooth \emph{functions} rather than \eqn{T} arbitrary period effects, and is the whole difference between this and an unrestricted factor model.

Two tuning choices are therefore unavoidable, and both are exposed. \code{kss_smooth} defaults to generalized cross-validation. \code{kss_L} defaults to the Bai and Ng (2002) \eqn{IC_{p2}} criterion -- \bold{not} the threshold rule of the original paper, which is stated here because it matters: \eqn{IC_{p2}} is the standard criterion for this problem and the literature that motivates the estimator treats Bai and Ng (2002) and Bai (2009) as the comparable methods, but a different rule can select a different \eqn{L}. Set \code{kss_L} explicitly to take the choice out of the criterion's hands, and inspect \code{$kss$eigenvalues} to see whether the choice was close.

\bold{Mundlak adjustment terms (\code{mundlak}).} The true random effects model assumes the firm effect is independent of the regressors. In production data it usually is not -- better managers choose different input levels -- and the effect is then split between the slopes and inefficiency, biasing the technology parameters and understating inefficiency. Mundlak's (1978) device models the correlation instead of assuming it away, \eqn{\alpha_i = \pi'\bar{x}_i + \delta_i}. Substituting into the frontier leaves the \emph{same likelihood with a wider design matrix}, which is why this is an argument rather than a new \code{model_name}.

Karagiannis and Kellermann (2019) widen it to \eqn{\alpha_i = \pi'\bar{x}^*_i + \gamma'z_i + \delta_i}, adding environmental variables -- conditions outside the producer's control, such as agro-ecological conditions or the geography an airline serves. Those enter the formula as ordinary time-invariant regressors; \code{mundlak} supplies the piece a user cannot easily build by hand, the within-firm means computed against the same panel index the fit uses. Naming a \emph{subset} of the regressors gives their Hausman-Taylor variant, in which only some regressors are taken to correlate with the effect.

Measured over 40 replications with \eqn{N = 80}, \eqn{T = 6}, on a design where the firm effect is correlated with \eqn{\bar{x}} by construction, and on a control where it is not:

\tabular{lllll}{
  design \tab estimator \tab bias \tab RMSE \cr
  effect correlated \tab plain \tab +0.0479 \tab 0.0570 \cr
  effect correlated \tab \code{mundlak} \tab \bold{-0.0061} \tab \bold{0.0310} \cr
  no correlation \tab plain \tab -0.0033 \tab \bold{0.0253} \cr
  no correlation \tab \code{mundlak} \tab -0.0068 \tab 0.0311
}

So it removes the bias where there is one and nearly halves RMSE -- and where there is no correlation it costs about 23\% in RMSE for nothing. It is not free, and it is not a default.

A variable that is already constant within firm gets no mean column: its own mean is itself, so the two would be exactly collinear. Every \code{mundlak} fit also emits a collinearity warning from the \emph{starting-value} regression, which is expected and harmless -- in a between-individual regression \eqn{x} and \eqn{\bar{x}} are the same column by construction. The starting values fall back to pooled OLS and the likelihood is estimated with the full requested formula.

\bold{The criterion needs a long enough panel, and says so.} \eqn{IC_{p2}} works because the residual variance flattens out once the real factors are in and the penalty then dominates. On a short panel that never happens: if \eqn{L} is allowed up to \eqn{T-1} the factors can span nearly the whole time dimension, the variance collapses toward zero, and no penalty competes. Measured at \eqn{n = 100} over five seeds and two true ranks, an uncapped criterion returned its maximum on \emph{all ten} designs at \eqn{T = 6} and again at \eqn{T = 8}. The automatic search is therefore capped at \eqn{\lfloor T/2 \rfloor}, which recovers the true rank exactly from \eqn{T = 8} upward; at \eqn{T = 6} it is right for a one-factor structure and occasionally one too high for a two-factor one. Whenever the automatic choice lands \emph{on} that cap \code{psfm} warns, because a selection at the boundary is not a free choice -- read it as "at least this many", inspect the eigenvalues, and prefer a longer panel to a larger \code{kss_L_max}.

One consequence of the period centering is worth knowing before comparing models: \code{"KSS"}'s \code{alpha_hat} carries no period mean, while \code{"CSS"}'s and \code{"LS"}'s do. The two therefore differ by a per-period constant by construction and should not be correlated against each other directly. \code{u_hat} is a within-period contrast and is invariant to that shift, so it is the quantity to compare across these models -- and is what the estimators are for.

\bold{Time-varying fixed effects without a distribution for u.} \code{"CSS"} and \code{"LS"} sit between \code{"SSFE"}, which holds inefficiency fixed over time, and the maximum-likelihood panel models, which let it move but only by assuming a distribution for \eqn{u} and a decay path for it. Both instead let the firm effect itself vary with time and read inefficiency off it, so neither assumes any distribution and neither returns a log-likelihood.

They differ in what they are willing to spend on flexibility. Cornwell, Schmidt and Sickles (1990) give every firm its own quadratic in time, \eqn{\alpha_{it} = \theta_{i1} + \theta_{i2}t + \theta_{i3}t^2}: \eqn{3N} parameters, firms may overtake one another, and \eqn{T_i \ge 4} is required before a firm contributes to \eqn{\beta} at all (a firm observed three times is fitted exactly by its own quadratic and its residuals are identically zero -- \code{sfa} warns rather than letting such firms silently dilute the estimate). Lee and Schmidt (1993) impose one common temporal pattern scaled per firm, \eqn{\alpha_{it} = \delta_t\alpha_i}: \eqn{N + T - 1} parameters, and the ranking of firms cannot change, only the spread. Cross-over versus no cross-over is the substantive choice between them, not a technical one.

\code{"CSS"} is fitted by residualizing \eqn{y} and \eqn{X} on each firm's own \eqn{(1, t, t^2)} rather than by building \eqn{3N} dummy columns; a regressor spanned by those quadratics -- a pure time trend, for instance -- is a rank failure and is refused rather than reported as a small coefficient. \code{"LS"} is a rank-one factor model, so it is fitted by alternating least squares, which handles an unbalanced panel with no special-casing; \eqn{\delta_1 = 1} fixes the scale that \eqn{\delta_t\alpha_i} leaves free, with a fallback to a norm constraint when the first period's effect is near zero.

For both, inefficiency is identified only relative to the best firm \emph{in each period} -- \eqn{u_{it} = \max_j \alpha_{jt} - \alpha_{it}} -- which is what makes them time-varying at all; comparing against a single overall best would collapse them back to a time-invariant ranking. Neither absorbs an intercept separately: it goes into the firm effect, as in \code{"SSFE"}, so \code{coef()} returns slopes only.

\bold{For SSFE model:}
\item{alpha_hat}{The estimated firm-specific fixed effects (named by individual).}
\item{u_hat}{Estimated technical inefficiency for each firm, as its distance from the best-performing firm's fixed effect.}
\item{exp_u_hat}{Predicted technical efficiency, \code{exp(-u_hat)}.}

\bold{For the GTRE model:}
\item{sigh_at_bound, sigr_at_bound}{\code{TRUE} when that persistent scale has collapsed to the zero boundary, judged relative to \code{sigma}. Frequently the CORRECT maximum likelihood estimate rather than a failure -- see \sQuote{Details} -- so it is reported rather than suppressed, and \code{psfm} warns naming the other persistent scale, since the two must be read together.}

\bold{For CSS and LS models:}
\item{alpha_hat}{The fitted firm effect for every observation, \eqn{\alpha_{it}} -- time-varying, unlike \code{SSFE}'s one value per firm.}
\item{u_hat, exp_u_hat}{Inefficiency and efficiency per observation, measured against the best firm in the same period.}
\item{frontier_alpha}{The best firm effect in each period, i.e. the benchmark \code{u_hat} is measured from.}
\item{sigma_v, residuals, df.residual}{The residual standard deviation, the residuals net of the firm effect, and the residual degrees of freedom (\eqn{NT - K - 3N} for \code{CSS}, \eqn{NT - K - (N + T - 1)} for \code{LS}).}
\item{theta}{\code{CSS} only. An \eqn{N \times 3} matrix of the firm-specific quadratic coefficients.}
\item{delta, alpha_i, ls_iterations, ls_converged}{\code{LS} only. The common temporal pattern (with \code{delta[1] == 1}), the firm-specific scales, and the state of the alternating least-squares iteration.}

\bold{For SSRE and SSCRE models:}
\item{alpha_hat}{The BLUPs of the firm effects, one per firm. Under \code{SSCRE} these are net of the part the Mundlak means explain.}
\item{theta, ercomp}{The GLS transformation parameter and the full variance-components decomposition from \code{\link[plm]{ercomp}}.}
\item{mundlak_terms}{\code{SSCRE} only. The names of the added group-mean regressors.}
\item{plm_fit}{The underlying \code{\link[plm]{plm}} object, for users who want its own methods.}

\bold{For the KSS model:}
\item{kss}{A list: the selected dimension \code{L}, the smoothing penalty \code{kappa}, the estimated basis (a \eqn{T \times L} matrix of orthonormal columns), the \eqn{n \times L} matrix of \code{loadings}, the full \code{eigenvalues} sequence the dimension was chosen from, and the iteration state.}

\bold{For PL80 and BC92 models:}
\item{exp_u_hat}{Predicted technical efficiency for each observation, E[exp(-B_it*u_i)|epsilon_i] (generalizes Battese and Coelli's (1988, JoE) point estimator to the panel case, pooling all of a firm's observations to predict u_i and scaling by B_it per period).}
}

\note{
The GTRE_SEQ1 and GTRE_SEQ2 models use sequential estimation methods and do not return optimization objects or starting values. The SSFE, CSS, LS, SSRE, SSCRE and KSS models are likewise not maximum-likelihood based and return no optimization object, starting values, or log-likelihood. All panel models require the \code{individual} argument to identify panel units.
}
\examples{
\donttest{
library(sfa)     

## A deliberately small panel with few simulation draws, so the examples
## build quickly. Do not read the estimates as a serious fit. TRE_Z becomes
## unstable below about 60 firms, so this is close to the floor.
data_trial <- data_gen_p(t=6, N=60, rand = 100,
                         sig_u = 1,  sig_v = 0.3,
                         sig_r = .2, sig_h = .4,
                         cons = 0.5, beta1 = 0.5,
                         beta2 = 0.5)

## rand.psoptim is not optional with PSopt = TRUE: the particle swarm draws
## from the session RNG, so results otherwise change between builds.
max_tre_z   <-  psfm(formula    = y_tre_z ~ x1 +x2| z_gtre,
                     model_name = "TRE_Z",
                     data       = data_trial,
                     individual = "name",
                     halton_num = 50,
                     rand.psoptim = 7L,
                     PSopt      = TRUE)

## Time-varying fixed effects, no distribution assumed for u. CSS gives each
## firm its own quadratic in time; LS gives all firms one common pattern.
css_fit     <-  psfm(formula = y_tre ~ x1 + x2, model_name = "CSS",
                     data = data_trial, individual = "name", time = "year")
ls_fit      <-  psfm(formula = y_tre ~ x1 + x2, model_name = "LS",
                     data = data_trial, individual = "name", time = "year")
ls_fit$delta      # the common temporal pattern, normalized to delta[1] = 1

## The Schmidt-Sickles family: within, GLS, and correlated random effects.
## SSCRE's slopes equal SSFE's within slopes exactly (Mundlak 1978); the
## .mean_ coefficients are the Mundlak form of the Hausman test.
ssre_fit    <-  psfm(formula = y_tre ~ x1 + x2, model_name = "SSRE",
                     data = data_trial, individual = "name", time = "year")
sscre_fit   <-  psfm(formula = y_tre ~ x1 + x2, model_name = "SSCRE",
                     data = data_trial, individual = "name", time = "year")

## Kneip, Sickles and Song (2012): the temporal basis and its dimension are
## both estimated. Inspect what the criterion chose before trusting it.
kss_fit     <-  psfm(formula = y_tre ~ x1 + x2, model_name = "KSS",
                     data = data_trial, individual = "name", time = "year")
kss_fit$kss$L                 # selected dimension
kss_fit$kss$eigenvalues       # was the choice close?
}
}
