\name{psfm}
\alias{psfm}
\title{psfm}
\description{Function to implement various panel data stochastic frontier estimators}
\usage{psfm(formula, model_name = c("TRE_Z", "GTRE_Z", "TRE",
                    "GTRE", "GTRE_FML", "TFE", "TFE_WMLE", "FD",
                    "GTRE_SEQ1", "GTRE_SEQ2", "SSFE",
                    "PL80", "PL80_MVTN", "BC92", "K1990", "K1990modified"), data,
                    maxit.bobyqa = 5000, maxit.nlminb = 500,
                    maxit.psoptim = 100, maxit.optim = 1000,
                    REPORT = 1, trace = 3, pgtol = 0, individual,
                    halton_num = NULL, start_val = FALSE, gamma = FALSE,
                    PSopt = FALSE, optHessian = TRUE, inefdec = TRUE, Method = "L-BFGS-B",
                    verbose = FALSE, rand.gtre = NULL, rand.psoptim = NULL,
                    OPG_calc = FALSE,
                    estimator = c("fiml", "sml", "seq1", "seq2"),
                    collinear_action = c("start_only", "error", "warn_drop"),
                    time = NULL, tfe_lambda_max = 100)}
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
  \item{\code{"PL80"}}{Pitt and Lee (1980): time-invariant inefficiency.}
  \item{\code{"BC92"}}{Battese and Coelli (1992): exponential time decay.}
  \item{\code{"PL80_MVTN"}}{Pitt and Lee's (1981) Model III: inefficiency varies over time and is CORRELATED within a firm, through a multivariate normal truncated to the negative orthant. This is the likelihood from their Appendix 2, which they derived but did not use, calling the T-dimensional integrals intractable. Balanced panels only. See \code{\link{PL80_MVTN}}.}
  \item{\code{"K1990"}, \code{"K1990modified"}}{Kumbhakar (1990) time patterns for inefficiency.}
} Matching ignores case, so \code{"nhn"}, \code{"NHN"} and \code{"Nhn"} are the same choice.}
  \item{data}{A panel data set. This may be an ordinary \code{data.frame} (or tibble/data.table) -- the panel index is constructed internally from the \code{individual} and \code{time} arguments -- or a \code{plm::pdata.frame} that already carries an index, which is used as-is.}
  \item{maxit.bobyqa}{Maximum number of iterations for the bobyqa optimization routine}
  \item{maxit.nlminb}{Maximum number of iterations for the nlminb optimization routine. Previously hard-coded (200 in the \code{GTRE_FML} branch, 500 elsewhere) and not settable from the call.}
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
}
\details{The generalized true random effects model (GTRE, 4-component model) and true random effects models (TRE)  are both estimated by simulated maximum likelihood based on the paper by the Fillipini and Greene (2016, JPA). The TRE_Z and GTRE_Z allow for modeling the u-component of the GTRE and TRE with determinants of inefficiency. The first-difference estimator (FD) of Wang and Ho (2010, JoE) is also available.

\bold{Choosing how GTRE is estimated.} Since version 1.1.4 the four routes to the four-component model are selected with \code{estimator} rather than by picking among four \code{model_name} values that looked like four different models. \code{model_name = "GTRE"} now defaults to \code{estimator = "fiml"}; through 1.1.3 it fitted the simulated-ML estimator, which is \code{estimator = "sml"}. Scripts that pass \code{"GTRE"} therefore get a different estimator than they did and are warned once per call. The explicit names \code{"GTRE_FML"}, \code{"GTRE_SEQ1"} and \code{"GTRE_SEQ2"} are unaffected.

\bold{Two ways of estimating the same four-component model.} \code{"GTRE"} and \code{"GTRE_FML"} both fit
\deqn{y_{it} = \beta_0 + x_{it}'\beta + r_i - h_i + v_{it} - u_{it},}
with \eqn{r_i, v_{it}} two-sided and \eqn{h_i, u_{it}} one-sided (half-normal), separating \emph{persistent} inefficiency \eqn{h_i} from \emph{transient} inefficiency \eqn{u_{it}}. They differ in how the likelihood is obtained.

\code{"GTRE"} uses simulated maximum likelihood over Halton draws. \code{"GTRE_FML"} uses \emph{full information} maximum likelihood through the model's closed-skew-normal representation: with \eqn{A = -[1_T, I_T]}, \eqn{V = \mathrm{diag}(\sigma_h^2, \sigma_u^2 I_T)} and \eqn{\Sigma = \sigma_v^2 I_T + \sigma_r^2 1_T 1_T'}, the firm-level density is available in closed form and nothing is integrated by Monte Carlo. The likelihood is therefore deterministic, and repeated fits to the same data give identical answers.

The \eqn{(T+1)}-dimensional normal CDF this requires would normally be the bottleneck, but its covariance is diagonal-plus-rank-one for this model, so it collapses to a one-dimensional Gauss-Hermite integral -- the same reduction used for \code{"TFE"}.

\code{"GTRE_FML"} reports \bold{raw standard deviations} \code{sigr}, \code{sigv}, \code{sigh}, \code{sigu} (individual effect, noise, persistent inefficiency, transient inefficiency), \emph{not} the \eqn{\lambda}/\eqn{\sigma} reparameterization \code{"GTRE"} uses, and it requires a \bold{balanced} panel -- the closed-skew-normal pieces are built at a single \eqn{T}. It errors rather than silently mis-estimating if the panel is unbalanced; use \code{"GTRE"} in that case.

\bold{Starting values for \code{"GTRE_FML"}.} The closed-skew-normal likelihood carries a boundary optimum at \eqn{\sigma_h = 0}, where the model collapses to \code{"TRE"} and the intercept absorbs the missing \eqn{E[h] = \sigma_h\sqrt{2/\pi}}. A fit that falls into it can have a \emph{lower} likelihood than the true parameter vector, so it is an optimizer failure rather than a limit of the data. Following Colombi (2010) and Colombi, Martini and Vittadini (2011), who recommend seeding this likelihood from the two-step moment estimates, \code{psfm()} builds two candidate starts -- the random-effects decomposition shared with \code{"GTRE"}, and the two-step estimator reported by \code{"GTRE_SEQ2"} -- evaluates the likelihood at each, and begins from the better one. Which was chosen, and the log-likelihood at both, is returned in \code{start_search}. Supplying \code{start_val} overrides the search.

Note that \eqn{\sigma_r} and \eqn{\sigma_h} enter the density only through \eqn{\sigma_r^2 + \sigma_h^2}; their separation rests entirely on the skewness term. That pair is therefore the least precisely estimated part of the model even when the fit is at the global optimum, which is consistent with the confounding between latent heterogeneity and persistent inefficiency reported in the source papers.

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

\bold{For SSFE model:}
\item{alpha_hat}{The estimated firm-specific fixed effects (named by individual).}
\item{u_hat}{Estimated technical inefficiency for each firm, as its distance from the best-performing firm's fixed effect.}
\item{exp_u_hat}{Predicted technical efficiency, \code{exp(-u_hat)}.}

\bold{For PL80 and BC92 models:}
\item{exp_u_hat}{Predicted technical efficiency for each observation, E[exp(-B_it*u_i)|epsilon_i] (generalizes Battese and Coelli's (1988, JoE) point estimator to the panel case, pooling all of a firm's observations to predict u_i and scaling by B_it per period).}
}

\note{
The GTRE_SEQ1 and GTRE_SEQ2 models use sequential estimation methods and do not return optimization objects or starting values. The SSFE model is likewise not maximum-likelihood based and returns no optimization object, starting values, or log-likelihood. All panel models require the \code{individual} argument to identify panel units.
}
\examples{
\donttest{
library(sfa)     

data_trial <- data_gen_p(t=10,N=100, rand = 100, 
                         sig_u = 1,  sig_v = 0.3, 
                         sig_r = .2, sig_h = .4, 
                         cons = 0.5, beta1 = 0.5,
                         beta2 = 0.5)

max_tre_z   <-  psfm(formula    = y_tre_z ~ x1 +x2| z_gtre,
                     model_name = "TRE_Z",
                     data       = data_trial,
                     individual = "name",
                     PSopt      = TRUE)
}
}