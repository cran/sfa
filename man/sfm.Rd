\name{sfm}
\alias{sfm}
\title{sfm}
\description{Implementation of the cross-sectional stochastic frontier model across an array of distributional assumptions for both v and u (user specified). For panel models, see the psfm() call.}
\usage{sfm(formula, model_name = c("NHN", "NHN_Z", "NE", "NE_Z", "NR",
    "THT", "NTN", "NG", "NNAK", "NU", "NGE", "NLN", "NW", "tHN",
    "TSL"), data, maxit.bobyqa = 10000, maxit.nlminb = 500,
    maxit.psoptim = 1000, maxit.optim = 1000, REPORT = 1, trace = 2,
    pgtol = 0, start_val = FALSE, PSopt = FALSE, use.nlminb = "auto",
    use.bobyqa = "auto", optHessian = TRUE, inefdec = TRUE, upper =
    NA, Method = "L-BFGS-B", robust = c("mle", "mlqe", "psi",
    "mdpd"), c_mlqe = 0.2, eta = 0.01, alpha = 0.2, verbose = FALSE,
    start_from = NULL, Nsim = "auto", z_link = c("sd", "var"), vhet =
    NULL, uhet = NULL, muhet = NULL, scaling = NULL, shapehet = NULL,
    weights = NULL, wscale = TRUE, sim_type = c("halton", "sobol",
    "torus", "uniform"), antithetics = FALSE, sim_burn = NULL,
    sim_scrambling = 0L, sim_prime = NULL, sim_seed = NULL,
    rand.psoptim = NULL, keep_objective = FALSE, estimator = c("mle",
    "cols", "mols"), cols_boot = 0, rand.cols = NULL)
}
\arguments{
\item{formula}{a symbolic description for the model to be estimated}
\item{model_name}{Which distributional pair to fit. The noise term \eqn{v} is
normal except where stated, and the one-sided inefficiency term \eqn{u} follows:
\describe{
  \item{\code{"NHN"}}{half normal.}
  \item{\code{"NHN_Z"}}{half normal with heteroskedastic \eqn{\sigma_u}, parameterized by the determinants \eqn{z} given in the second part of the formula (\code{y ~ x | z}) as \eqn{\sigma_u = \exp(z'\delta)}.}
  \item{\code{"NE"}}{exponential.}
  \item{\code{"NE_Z"}}{exponential with heteroskedastic \eqn{\sigma_u}, parameterized as for \code{"NHN_Z"}.}
  \item{\code{"NTN"}}{truncated normal.}
  \item{\code{"NR"}}{Rayleigh. Note this is a distinct family from \code{"NHN"}, not a reparameterization of it: the Rayleigh density carries an extra factor of \eqn{u}, and its standardized skewness is fixed at -0.631 against the half normal's -0.996.}
  \item{\code{"NU"}}{uniform.}
  \item{\code{"NG"}}{gamma. See \code{ng_starts} under \sQuote{Value} -- this likelihood determines \eqn{E[u] = \mu\sigma_u} far better than either factor separately.}
  \item{\code{"NNAK"}}{Nakagami, with shape \eqn{m}. Nests \code{"NR"} at \eqn{m = 1}. The shape is weakly identified in practice, so \eqn{\sigma_u} and \eqn{m} are strongly correlated across samples.}
  \item{\code{"NGE"}}{generalized exponential.}
  \item{\code{"NLN"}}{lognormal. Estimated by simulated maximum likelihood; see \code{Nsim} and \sQuote{Details}.}
  \item{\code{"NW"}}{Weibull. Also simulated maximum likelihood.}
  \item{\code{"tHN"}}{half normal, with Student-\emph{t} noise.}
  \item{\code{"TSL"}}{truncated skew-Laplace.}
  \item{\code{"THT"}}{half \emph{t}, with Student-\emph{t} noise.}
} Matching ignores case, so \code{"nhn"}, \code{"NHN"} and \code{"Nhn"} are the same choice.}
\item{data}{A data frame containing the variables named in \code{formula}.}
\item{maxit.bobyqa}{Maximum number of iterations for the bobyqa optimization routine}
\item{maxit.psoptim}{Maximum number of iterations for the psoptim optimization routine}
\item{maxit.optim}{Maximum number of iterations for the optim optimization routine}
\item{REPORT}{reporting parameter}
\item{trace}{Integer. Tracing level passed through to the optimizer; larger values print more.}
\item{pgtol}{Numeric. Projected-gradient tolerance passed to \code{optim()}'s \code{"L-BFGS-B"} method. \code{0} uses the optimizer's own default.}
\item{start_val}{starting value (optional)}
\item{PSopt}{use psoptim optimization routine (T or F)}
\item{maxit.nlminb}{Maximum number of iterations for the nlminb optimization routine.}
\item{use.nlminb}{\code{"auto"} (default), or a logical. Whether to run the \code{stats::nlminb} stage, which precedes bobyqa in the optimizer stack and accepts an analytic gradient where one is available. Under \code{"auto"} it runs only for the models where it has been verified to reach an equal-or-better optimum than the previous bobyqa-based path (\code{NHN}, \code{NE}, \code{NTN}, \code{NU}), where it is roughly 4-8 times faster at identical accuracy.}
\item{use.bobyqa}{\code{"auto"} (default), or a logical. Whether to run the derivative-free bobyqa stage. Under \code{"auto"} it runs for every model \emph{except} those handled by nlminb above, reproducing the previous behaviour exactly. It is kept on by default for the remaining models deliberately: for \code{NR} and \code{NGE} the nlminb path converges to a materially worse optimum, so accuracy, not speed, decides which stage runs.}
\item{optHessian}{Logical. Should a numerically differentiated Hessian matrix be returned while using the optim routine? (for optim routine)}
\item{inefdec}{Production or cost function}
\item{upper}{Vector of upper values for the optim package.}
\item{Method}{The method to be used for optim. See 'Details' within optim.}
\item{robust}{One of \code{"mle"} (default), \code{"mlqe"}, \code{"psi"}, or
\code{"mdpd"}. Selects a robust divergence-based objective in place of
ordinary maximum likelihood -- see Details. Currently only implemented for
\code{model_name = "NHN"}; \code{sfm()} stops with a clear error if a
non-\code{"mle"} value is requested for any other model.}
\item{c_mlqe}{Numeric, default \code{0.20}. The MLqE power parameter
(\code{c} in \code{f^c}; the corresponding \code{q = 1-c}), used only when
\code{robust = "mlqe"}.}
\item{eta}{Numeric, default \code{0.01}. The Psi-likelihood power parameter,
used only when \code{robust = "psi"}.}
\item{alpha}{Numeric, default \code{0.2}. The MDPD (minimum density power
divergence) power parameter, used only when \code{robust = "mdpd"}.}
\item{verbose}{Logical. Print optimization progress messages? Default is \code{FALSE}.}
\item{weights}{Optional observation weights, one per row of the data actually fitted (that is, after missing-data handling; a length mismatch is an error rather than being recycled). They multiply the per-observation log-likelihood contributions. For \bold{frequency} weights this is exact -- a weight of 2 gives the same estimates as including the row twice, which the package tests -- and for \bold{sampling} weights it yields a pseudo-likelihood whose point estimates are consistent but whose Hessian-based standard errors are not. In that case use the \pkg{sandwich} methods this package registers, e.g. \code{sandwich::sandwich(fit)}. Cannot be combined with the robust divergence estimators, which reweight observations themselves.}
\item{wscale}{Rescale \code{weights} to sum to the number of observations. This puts the log-likelihood on the same scale as an unweighted fit, so \code{AIC}/\code{BIC} stay comparable; it is a common factor on the objective and does not move the estimates. Set \code{FALSE} for genuine frequency weights, where the unscaled total is the quantity of interest.}
\item{shapehet}{\bold{Experimental.} Optionally, a one-sided formula parameterizing the \emph{shape} of \eqn{u} rather than its scale, as \eqn{\mu_i = \mu\exp(z_i'\delta)}. Implemented for \code{model_name = "NG"} and \code{"NNAK"}, the two families that have a shape parameter at all. On a design with a true shape effect it improves the log-likelihood by a very large margin (about 10,000 points at \eqn{n = 2000}), so it does capture the variation. \bold{But its parameters do not recover well in testing}: on that design the shape coefficient came back at 0.19 against a true 0.4 while \eqn{\sigma_V} collapsed to its lower bound, which is the signature of a degenerate optimum rather than an estimate. Plain \code{"NG"} recovers cleanly on the same design without a shape effect, so this is specific to the extended model and not to \code{"NG"} itself. Treat the coefficients as exploratory, check \eqn{\sigma_V} for a boundary collapse before reading them, and prefer \code{uhet} where scale heterogeneity is what you actually mean.}
\item{scaling}{Optionally, a one-sided formula of covariates imposing the \bold{scaling property} of Wang and Schmidt (2002) and Alvarez, Amsler, Orea and Schmidt (2006): \eqn{u_i = h(z_i,\delta)\,u^*_i} with \eqn{h = \exp(z_i'\delta)}, so the covariates scale a common draw \eqn{u^*} rather than entering its variance. Both \eqn{\sigma_U} and the pre-truncation mean move by the \emph{same} factor, which holds the shape of the distribution fixed and leaves only its scale free -- the identification argument those papers make. Implemented for \code{model_name = "NTN"} only: for a half-normal \eqn{u}, \eqn{h(z)|N(0,\sigma^2)|} \emph{is} \eqn{|N(0,(h\sigma)^2)|}, so \code{uhet} already fits that model exactly and a separate option would be two names for one thing. Cannot be combined with \code{uhet} or \code{muhet}, since letting those vary separately would undo the restriction. The reported \code{scale.} coefficients are \eqn{\delta}; \code{Zu.(Intercept)} and \code{Zmu.(Intercept)} are the scale and pre-truncation mean of the common draw.}
\item{start_from}{Optionally, a previously fitted \code{"sfareg"} object to take starting values from. Parameters are matched \bold{by name}, not by position: \code{"NE"} reports \code{sigv}/\code{sigu} where \code{"NHN"} reports \code{lambda}/\code{sigma}, so a positional copy would slide a variance ratio into a standard deviation. Anything unmatched keeps the start the package computes for it, and a fit sharing nothing but the frontier coefficients carries only those. Fitting a simpler model first and cascading is the documented remedy for the models with awkward likelihoods: over \bold{forty} samples for \code{"NG"} seeded from \code{"NE"}, the seeded fit was never worse, indistinguishable on thirty-nine, and worth \bold{176.7 log-likelihood points} on one. The rescue is rare and large rather than common and small, which is what makes it worth doing: the sample it rescues is one where the default start converges quietly to a much worse optimum, with nothing in the returned object to say so.}
\item{Nsim}{Total number of simulation draws per observation for the simulated-maximum-likelihood models (\code{"NLN"}, \code{"NW"}); ignored by every closed-form model. \code{"auto"} (the default) uses \code{max(200, ceiling(3*sqrt(n)))}. Half the draws are taken from each of the two proposals described under \sQuote{Details}. Simulated ML is consistent only if the draw count grows with the sample size -- at a fixed \code{Nsim} the simulation bias does not vanish and every parameter except the frontier slopes converges to the wrong value. The previous fixed default of \code{100} did exactly that: at \eqn{n = 3000} it returned \code{(0.362, 0.821, 1.284, 0.328)} against a truth of \code{(0.3, 1.0, 1.5, 0.5)} for \code{"NW"}. A numeric value below the automatic floor warns. Earlier versions warned that \code{"NLN"} needed far more draws than \code{"NW"}; that was a property of the old integration scheme, not of the model. Both now behave alike: at \eqn{n = 3000}, measured against quadrature at the true parameters, total simulation error is 3.4 log-likelihood units for \code{"NLN"} and 3.0 for \code{"NW"} at the default, falling to 0.07 and 0.06 at \code{Nsim = 6400}. The comparable figure for the old scheme was 226.8.}
\item{rand.psoptim}{Integer. seed for replication of psoptim.  Default to \code{NULL}.}
  \item{keep_objective}{Logical. Store the likelihood function on the fitted object, so that \code{\link{sfa_diagnostics}} can compute the gradient at the optimum and draw likelihood slices. Defaults to \code{FALSE}: a closure carries its enclosing environment, so a fit saved with one serializes the estimation data too -- on a 200-observation example the object grows from about 38 KB to 1.7 MB. Everything else \code{sfa_diagnostics()} reports works without it.}
\item{estimator}{\code{"mle"} (default), or \code{"cols"}/\code{"mols"} -- two names for the same thing, see \sQuote{Details} -- which select the modified ordinary least squares moment estimator instead of maximum likelihood: OLS slopes are kept as they are, the second and third central moments of the OLS residuals are inverted for the scale parameters, and the intercept is shifted up by the implied \eqn{E[u]}. Closed form, no optimizer, deterministic. Implemented for \code{"NHN"}, \code{"NE"} and \code{"NG"}; other models error, because the moment inversion is distribution-specific. See \sQuote{Details}.}
  \item{z_link}{Scale on which the variance-determinant linear predictor is placed, for the \code{_Z} models. \code{"sd"} (the default, and \code{sfm}'s historical convention) sets \eqn{\sigma_u = \exp(z'\delta)}; \code{"var"} sets \eqn{\sigma_u = \sqrt{\exp(z'\delta)}}, which is what \code{\link{psfm}} and the competing packages use. The two are reparameterizations of the same model -- they reach the same maximised log-likelihood and the same \code{\link{marginal_effects}} -- but \eqn{\delta} itself differs by a factor of two, so a \eqn{\delta} from \code{sfm} is not comparable with one from \code{psfm} or from \pkg{sfaR} unless the same convention is used. Set \code{z_link = "var"} to put them on one footing. See \sQuote{Details}.}
\item{vhet}{One-sided formula for the determinants of the \emph{noise} scale, e.g. \code{vhet = ~ z1 + z2}, giving \eqn{\sigma_v = \exp(z_v'\delta_v)} under the SD link. \code{NULL} (the default) leaves \eqn{\sigma_v} constant. Available for \code{"NHN"}, \code{"NHN_Z"}, \code{"NE"}, \code{"NE_Z"} and \code{"NTN"}.}
\item{uhet}{One-sided formula for the determinants of \eqn{\sigma_u}. This is the same specification as the \code{| z} segment of \code{formula}, and it is an error to give both; \code{uhet} exists so that \code{"NTN"}, which takes no pipe segment, can still have a heteroskedastic \eqn{\sigma_u}.}
\item{muhet}{One-sided formula for the determinants of the \emph{pre-truncation mean} of \eqn{u}, \eqn{\mu = z_\mu'\delta_\mu}. This is the Battese and Coelli (1995) inefficiency-effects model and requires \code{model_name = "NTN"}, the only family in which \eqn{\mu} exists. Unlike the two scale blocks \eqn{\mu} carries no link, since it may take either sign.}
  \item{sim_type}{Sequence used to build the simulated-ML draws, for the models fitted by simulation (\code{"NLN"}, \code{"NW"}). \code{"halton"} (default) reproduces previous behaviour. \code{"sobol"} accepts \code{sim_scrambling} 1--3, giving the scrambled sequence of Bhat (2003). \code{"torus"} is the Korobov sequence and accepts \code{sim_prime}. \code{"uniform"} is pseudorandom, present as a baseline to measure against rather than as a recommendation. See \sQuote{Details}.}
  \item{antithetics}{Logical. If \code{TRUE}, half the draws are taken from the sequence and the other half created as their mirror images (Hammersley and Morton 1956). Default \code{FALSE}.}
  \item{sim_burn}{Number of leading elements of the sequence discarded. Default \code{NULL}, meaning 1000. The initial elements of different sequences are correlated with one another and discarding them removes it.}
  \item{sim_scrambling}{Integer 0--3, \code{sim_type = "sobol"} only, passed to \code{randtoolbox::sobol}. Non-zero gives a scrambled sequence.}
  \item{sim_prime}{Prime for \code{sim_type = "torus"}. Default \code{NULL} uses the generator's own choice.}
  \item{sim_seed}{Seed for the draws, for reproducibility. The caller's random-number stream is restored afterwards.}
  \item{cols_boot}{Integer. Number of nonparametric bootstrap resamples used to obtain standard errors under \code{estimator = "cols"}. Zero (the default) skips the bootstrap, in which case only the OLS slope standard errors are reported and the scale parameters and corrected intercept carry \code{NA}. Ignored when \code{estimator = "mle"}.}
  \item{rand.cols}{Integer. Seed for the COLS bootstrap, for reproducibility. The caller's random-number stream is restored afterwards. Default \code{NULL}.}
}
\details{The options include the Normal-Half Normal (NHN), Normal-exponential (NE), Student's t-Half t (THT), the Normal-Truncated Normal (NTN), the Normal-Uniform (NU), the Normal-Generalized Exponential (NGE), the Normal-Lognormal (NLN), the Normal-Weibull (NW), the Student's t-half normal (tHN) and the Normal-Truncated Skew-Laplace (TSL).  NHN_Z and NE_Z are extensions for the NHN and NE models that allow for modeling the u-component of those models with determinants of inefficiency.

\bold{Modified ordinary least squares (\code{estimator = "cols"} or \code{"mols"}).} The moment estimator of Olson, Schmidt and Waldman (1980, \emph{Journal of Econometrics}), standardly called MOLS. Note the name: this is \emph{not} Winsten's corrected OLS, which shifts the intercept by the largest residual and estimates no variance parameters at all. The package has shipped this estimator under the spelling \code{"cols"} since it was introduced, so that spelling continues to work and is not deprecated; \code{"mols"} is accepted as an exact synonym so the method can be found under its usual name. OLS is consistent for the \emph{slopes} of a composed-error frontier whatever the one-sided distribution; only the intercept is biased, by \eqn{E[u]}, because the composed error does not have zero mean. COLS therefore keeps the OLS slopes, inverts the central moments of the OLS residuals for the scale parameters, and shifts the intercept up by the implied \eqn{E[u]}. Writing \eqn{m_2}, \eqn{m_3} for the second and third central moments of \eqn{e = v - u}: for \code{"NHN"}, \eqn{m_3 = \sigma_u^3\sqrt{2/\pi}(1-4/\pi)} and \eqn{m_2 = \sigma_v^2 + \sigma_u^2(1-2/\pi)} with \eqn{E[u] = \sigma_u\sqrt{2/\pi}}; for \code{"NE"}, \eqn{m_3 = -2\sigma_u^3} and \eqn{m_2 = \sigma_v^2+\sigma_u^2} with \eqn{E[u] = \sigma_u}; for \code{"NG"}, from the third and fourth cumulants, \eqn{\sigma_u = -k_4/(3k_3)} and \eqn{\mu = -k_3/(2\sigma_u^3)} with \eqn{E[u] = \mu\sigma_u}.

No optimizer runs and the answer is deterministic, which makes COLS a natural robustness check against a maximum-likelihood fit that may have settled at a local optimum. It is less efficient than ML when the distributional assumption holds, and the \code{"NG"} inversion in particular leans on a fourth cumulant and is noticeably noisier than the other two at moderate sample sizes.

\bold{Wrong skew under COLS.} A production frontier implies \eqn{m_3 < 0}. When a sample comes out with \eqn{m_3 \ge 0} the moment equations have no admissible solution: \code{sfm()} warns, reports \code{sigu} as 0 with the whole residual variance assigned to \code{sigv}, and returns no efficiency predictions. Read that as \emph{no evidence of inefficiency in these data}, not as an estimate of zero. This is the Type I failure of Olson, Schmidt and Waldman and is common in small samples; \code{$wrong_skew} and \code{$residual_moments} record it.

\bold{\code{NR} is not \code{NHN}.} \code{"NR"} places a Rayleigh distribution on \eqn{u}, not a half-normal, and the two are separate families rather than reparameterizations of one another: the Rayleigh density carries an extra factor of \eqn{u}, and no rescaling removes it. The practical consequence is that the standardized skewness each contributes to the composed error is a different constant, \eqn{-0.631} against the half-normal's \eqn{-0.995}, so a sample generated under one cannot be fitted by the other. \eqn{\sigma_u} is reported on the second-raw-moment convention \eqn{E[u^2] = \sigma_u^2}, matching \code{"NHN"} and making the two directly comparable; the Rayleigh scale is then \eqn{\sigma_u/\sqrt{2}}, with \eqn{E[u] = \sigma_u\sqrt{\pi}/2} and \eqn{Var(u) = (1-\pi/4)\sigma_u^2}. Starting values come from inverting these moments, since the Rayleigh skewness is a constant and so \eqn{m_3} identifies \eqn{Var(u)} outright. Use \code{data_gen_cs()}'s \code{y_pcs_r} column to test it, not \code{y_pcs}.

\bold{\code{tHN} is not \code{THT}.} In \code{THT} (Tancredi 2002) a single scale mixture is shared by both error components, so the noise and the inefficiency are both Student's t with the same degrees of freedom, and that shared mixture is what makes the composed error a closed-form skew-t. In \code{tHN} the noise is Student's t and the inefficiency is half-normal, drawn \emph{independently}, so the two components have different tail behaviour, there is no closed form, and the density is the convolution \eqn{f(e) = \int_0^\infty f_v(e+u) f_u(u) du} evaluated by Gauss-Legendre quadrature. \code{tHN} is the heavy-tailed-\emph{noise} model, which is what makes it the natural parametric comparison for the density-power robust estimators (\code{robust = "mlqe"}, \code{"psi"}, \code{"mdpd"}); \code{THT} cannot play that role because its inefficiency term is heavy-tailed too. Parameters are reported as \code{(sigv, sigu, nu)}, the conventional order, \emph{not} \code{THT}'s inverted \code{(sigu, sigv, a)}.

Two properties of \code{tHN} are worth knowing before using it. First, the degrees of freedom are weakly identified: on data simulated from the model at \eqn{n = 1000} with a true \eqn{\nu = 5}, the profile log-likelihood varies by only about 0.24 across \eqn{\nu} from 10 to 100 and peaks near 20, so a single reported \eqn{\nu} should not be taken at face value. Profile over a grid of fixed \eqn{\nu} and report the path. Because of that flat ridge \code{sfm()} always fits \code{tHN} from several widely separated starting values, keeps the best, and records the outcome in \code{thn_starts}; it warns when they reach different optima. Second, on real data \eqn{\sigma_u} can collapse onto zero, the heavy noise tail absorbing the entire one-sided component and leaving mean predicted efficiency near one. That is a property of the model rather than a numerical failure, so it is reported rather than suppressed: \code{sfm()} warns and sets \code{thn_sigma_u_at_bound}. Fits are also substantially slower than the closed-form models, since every likelihood evaluation is a quadrature whose node count grows with \eqn{\sigma_u/\sigma_v}.

\code{"NU"} assumes \eqn{u \sim U(0,\theta)} (Li 1996; Nguyen 2010), giving the closed-form composed density \eqn{f(e) = \theta^{-1}[\Phi((e+\theta)/\sigma_v) - \Phi(e/\sigma_v)]}. The reported \code{theta} is the upper bound of the inefficiency support, so predicted inefficiency is bounded by construction. \code{"NGE"} assumes a generalized exponential \eqn{F(u) = (1-e^{-\lambda u})^2}, whose composed density is a difference of two exponentially-tilted normal tail terms; it nests the exponential model and permits a non-monotonic inefficiency density. \code{NGE} reports \code{sigu} on the same scale as \code{NE} (the exponential mean, \eqn{\lambda = 1/\sigma_u}) so the two are directly comparable.

\bold{The variance-determinant scale.} \code{sfm}'s \code{_Z} models place \eqn{z'\delta} on the standard deviation, \eqn{\sigma_u = \exp(z'\delta)}, while \code{\link{psfm}}'s \code{TRE_Z}/\code{GTRE_Z} place it on the variance, \eqn{\sigma_u = \sqrt{\exp(z'\delta)}}, as do the competing packages. Since \eqn{\exp(\eta) = \exp(\eta/2)^2}, fitting the same data under the two conventions returns the same \eqn{\sigma_u}, the same log-likelihood and the same marginal effects, but \eqn{\delta} under the SD link is exactly \bold{half} \eqn{\delta} under the variance link. Reading a \eqn{\delta} from one family as though it came from the other therefore doubles or halves every reported effect. \code{z_link} exists so the two can be aligned; the default preserves \code{sfm}'s own convention so existing results do not move.

\bold{Heteroskedasticity in more than one component.} \code{vhet}, \code{uhet} and \code{muhet} let \eqn{\sigma_v}, \eqn{\sigma_u} and (for \code{"NTN"}) the pre-truncation mean \eqn{\mu} each depend on their own covariates, in any combination. They are named formulas rather than further pipe segments deliberately: pipe \emph{position} already carries different meanings in different model families -- the second segment is \eqn{\sigma_w} in \code{\link{ttsfm}} but \eqn{\sigma_h} in \code{\link{psfm}}'s \code{GTRE_Z} -- and a fourth position would not be readable. The existing \code{| z} segment keeps its meaning, so nothing that worked before changes.

Two reasons to reach for \code{vhet} in particular. Ignored heteroskedasticity in \eqn{v} does not merely cost efficiency: Caudill, Ford and Gropper (1995) show that it biases the estimated frontier, because the composed error's mean depends on the scales. And a \eqn{\sigma_v} that varies with observed characteristics is one of the few specifications that can absorb apparent wrong skewness without forcing \eqn{\sigma_u} onto the zero boundary. \code{muhet} with \code{model_name = "NTN"} is the Battese and Coelli (1995) inefficiency-effects model, in which covariates shift the location of inefficiency rather than its dispersion; that is usually the model wanted when the question is \dQuote{what makes firms inefficient}, since \eqn{\delta_\mu} then has a sign that reads directly.

These fits carry every positive quantity on the log scale, so the parameter vector is unconstrained and no variance can be driven onto a boundary by the optimizer; that is also why this path uses \code{nlminb} followed by an unbounded BFGS rather than the bounded four-stage scaffold the other models use. Coefficients are reported prefixed by block -- \code{Zv.}, \code{Zu.}, \code{Zmu.} -- since the same covariate may appear in more than one. \code{\link{marginal_effects}} reads the \eqn{\sigma_u} block of these fits with no extra arguments. \code{estimator = "cols"} and the \code{robust} divergences are moment-based and homoskedastic respectively, so neither combines with these arguments; both refuse rather than silently ignoring them.

\bold{Simulated-ML draws.} \code{"NLN"} and \code{"NW"} have no closed-form composed density and are fitted by maximum simulated likelihood, averaging over draws built by \code{sim_type}. Two properties are enforced whatever the sequence. Each observation receives its own \emph{contiguous block} of the sequence rather than a strided subsequence: a column-major fill hands observation \eqn{i} the elements \eqn{h_i, h_{i+n}, \ldots} of a van der Corput sequence, which is not equidistributed -- at \eqn{n = 500} and \code{Nsim = 100} the first observation's block spanned only \eqn{[0.50, 0.75]}, integrating over a quarter of the inefficiency distribution. And the draws are clamped away from 0 and 1, since \code{qnorm()} of an exact endpoint is infinite.

Blocking by observation is also what makes the draws worth using: Train (2002, ch. 9) attributes the advantage of a low-discrepancy sequence to its coverage \emph{and} to the negative correlation it induces across observations, and the second only arises when observations receive different blocks. On draw counts, Bhat (2001) found 100 Halton draws more precise than 1000 pseudorandom ones, with 125 Halton draws giving half the simulation error of 1000 random draws; \code{Nsim = "auto"} scales the count with \eqn{n} and warns if a supplied count is too small.

\code{"TSL"} assumes a truncated skew-Laplace inefficiency (Wang 2012, eq. 3.7; the distribution is due to Aryal and Rao 2005, and the model was published as Nguyen, Gupta and Wang 2014), with density \eqn{f(u) = \frac{1+\lambda}{\sigma_u(2\lambda+1)}\left[2e^{-u/\sigma_u} - e^{-(1+\lambda)u/\sigma_u}\right]} for \eqn{u \ge 0}. This is a \emph{signed} mixture of two exponentials -- the second enters negatively -- so it nests the exponential model as \eqn{\lambda \to 0} while allowing a non-monotonic inefficiency density, and the composed density is the corresponding signed combination of two normal-exponential terms. It is evaluated as a difference in log space rather than by subtracting the raw terms: both carry \eqn{\sigma_v^2/(2\sigma_u^2)}, which overflows once \eqn{\sigma_u} becomes small relative to \eqn{\sigma_v}, a region the optimizer routinely visits. \code{lambda} is a shape parameter and is the least sharply identified of the three. \bold{This is a property of the family, not of the fit, and it is worth understanding before reporting \code{lambda} at all.} The shape parameter changes the standardised distribution remarkably little: \eqn{E[U]/\sigma_u} equals 1 in both limits \eqn{\lambda \to 0} and \eqn{\lambda \to \infty}, and its maximum over all \eqn{\lambda} is only \eqn{4-2\sqrt{2} \approx \bold{1.1716}}, attained at \eqn{\lambda = 1/\sqrt{2}}; the standardised skewness moves only from about 1.95 down to 1.79 and back. So the data identify the product \eqn{\sigma_u \times m(\lambda)} far better than they identify \eqn{\sigma_u} and \eqn{\lambda} separately, and the two trade off against each other.

Measured on the package's own convergence design (\eqn{\sigma_u = 1}, \eqn{\lambda = 1.5}, 60 replications at each of \eqn{n = 500, 1000, 2000, 4000}), estimating one parameter with the other held at its true value:

\tabular{lll}{
 \bold{estimated} \tab \bold{log-MSE on log-n slope} \tab \bold{sd at n = 4000} \cr
 \eqn{\sigma_u}, \eqn{\lambda} free \tab -0.65 \tab 0.038 \cr
 \eqn{\sigma_u}, \eqn{\lambda} fixed \tab \bold{-1.12} \tab 0.020 \cr
 \eqn{\lambda}, \eqn{\sigma_u} free \tab -0.38 (median) \tab \bold{26.3} \cr
 \eqn{\lambda}, \eqn{\sigma_u} fixed \tab \bold{-1.07} (median) \tab 0.62
}

Each parameter attains the root-\eqn{n} rate when the other is known and fails it when the other is free, which is the signature of weak scale--shape separation rather than of an optimiser problem. \eqn{\lambda}'s dispersion falls fortyfold once \eqn{\sigma_u} is fixed.

\bold{What to do about it.} Report the functionals rather than the primitives where you can: \eqn{E[U]} converges at the root-\eqn{n} rate even with both parameters free (slope -1.03, sd 0.031 at \eqn{n = 4000}), and \eqn{E[e^{-U}]} nearly so (-0.89). Both are what an efficiency analysis actually uses. If \eqn{\lambda} itself is of interest, prefer profile-based inference over its reported \eqn{t}-value, which is a Wald statistic on a parameter whose sampling distribution is neither symmetric nor tightly concentrated.

\code{"NLN"} (\eqn{u} lognormal, reported as \code{sigu} = sdlog and \code{mu} = meanlog) and \code{"NW"} (\eqn{u} Weibull, reported as \code{sigu} = scale and \code{k} = shape) have no closed-form composed density and are estimated by \emph{simulated} maximum likelihood. The integral \eqn{f(e) = \int_0^\infty f_u(u)\,\phi((e+u)/\sigma_v)/\sigma_v\,du} is a product of two densities, and simulating from either one alone fails wherever the other is the narrow factor: draws from \eqn{f_u} miss the normal kernel for a very inefficient unit, while draws from the kernel miss a sharply peaked \eqn{f_u}. Half the draws are therefore taken from each, and combined by the balance heuristic of Veach and Guibas (1995), which weights every draw by the proposal likelier to have produced it. This needs no rule for choosing between the two and remains accurate in both regimes. The draws are generated once and reused at every optimizer iteration, so the objective is deterministic and reproducible; increase \code{Nsim} to reduce simulation error.

The efficiency predictor reuses the weights the likelihood computed, so \code{u_hat} is the Bayes rule under exactly the density that was maximised.

Outputs include E[exp(-u)|e] given by exp_u_hat, following Battese and Coelli (1988, JoE), where appropriate. As of this version \code{exp_u_hat} is also returned for \code{NE}, \code{NTN}, \code{NU}, \code{NGE}, \code{NLN} and \code{NW}, together with the Jondrow et al. (1982) point predictor \code{u_hat} = E[u|e]; previously the exponential and truncated-normal models returned no efficiency prediction at all.

Note on default iteration counts: \code{sfm()}'s optimization defaults (thousands of iterations) are much higher than \code{psfm()}'s (100/10/10). This is intentional: \code{sfm()}'s likelihoods are closed-form (no simulation), so each objective-function evaluation is cheap and many default iterations cost little; several of \code{psfm()}'s models require simulated maximum likelihood over Halton draws, where each evaluation is far more expensive. See \code{\link{psfm}}'s Details for more.

\strong{Robust divergence-based estimation (\code{robust} argument).} For
\code{model_name = "NHN"}, \code{sfm()} can maximize a robust divergence
objective instead of the ordinary log-likelihood, downweighting the
influence of outlying observations:
\itemize{
  \item \code{"mlqe"}: minimum L\eqn{q} estimation (Lq, Ferrari and Yang
  2010), power parameter \code{c_mlqe} (\eqn{q = 1-c}). Objective
  \eqn{\sum (f(e_i)^c - 1)/c}.
  \item \code{"psi"}: Psi-likelihood (Basu et al.-type divergence family),
  power parameter \code{eta}. Objective
  \eqn{\sum f(e_i)^c/c - n \, I_{1+c}/(1+c)}, where
  \eqn{I_{1+c} = \int f(e;\theta)^{1+c}\,de} is evaluated numerically at the
  current parameter values on every optimizer iteration.
  \item \code{"mdpd"}: minimum density power divergence (Basu et al. 1998),
  power parameter \code{alpha}. Objective
  \eqn{((1+c)/c) \sum f(e_i)^c - n \, I_{1+c}}.
}
All three reduce to \code{"mle"} as their power parameter approaches 0.
\strong{Psi and MDPD produce identical parameter point estimates for the
same power parameter} (\eqn{c = \eta = \alpha}) -- the two objectives are a
positive rescaling of one another, \eqn{\text{MDPD}(c) = (1+c)\cdot
\text{Psi}(c)}, which does not change the argmax; both are exposed for
compatibility with the corresponding paper's three-method terminology, not
because they differ numerically. MLqE is a genuinely distinct estimator.

Because these are M-estimators rather than ordinary MLE, the naive
Hessian-inverse standard error (\code{optHessian = TRUE}'s usual output) is
not statistically valid for them -- it assumes the information-matrix
equality, which these objectives do not satisfy in general. When
\code{robust != "mle"}, \code{sfm()} instead reports a sandwich-form
standard error, \eqn{A^{-1} B A^{-1}}, where \eqn{A} is the Hessian of the
summed robust objective (as already computed for the point estimate) and
\eqn{B = \sum_i g_i g_i'} is built from the per-observation gradient of that
objective (obtained via numerical differentiation,
\code{\link[numDeriv]{jacobian}}). This requires \code{optHessian = TRUE};
\code{std.errors}/\code{t.values} are \code{NA} otherwise, same as for
ordinary MLE.

Formulas verified against an independently audited reference implementation
(cross-checked by direct numerical integration and finite differences),
including the sandwich SE's own verification (the c -> 0 limit against the
ordinary MLE Hessian SE, and a Monte Carlo check against the empirical
spread of repeated fits). Currently only implemented for
\code{model_name = "NHN"} -- other models will error clearly if
\code{robust != "mle"} is requested.}
\author{David H. Bernstein and Alexander Stead}
\note{Standard errors require optHessian set to TRUE. When
\code{robust != "mle"} (\code{model_name = "NHN"} only), a sandwich-form SE
is reported instead of the ordinary Hessian-inverse one -- see Details.}

\seealso{\code{\link{psfm}} for panel models, \code{\link{zsfm}} for zero-inefficiency models, \code{\link{ttsfm}} for two-tier models, \code{\link{data_gen_cs}} for simulating data with known true parameters, and \code{\link{sfareg-methods}} for extracting results.}
\value{
An object of class \code{"sfareg"} containing the following components:

\item{out}{A matrix with parameter estimates, standard errors, and t-values.}
\item{opt}{A list containing the optimization results from the final optimization procedure.}
\item{total_time}{The total computation time for model estimation.}
\item{start_v}{The starting values used in the optimization.}
\item{model_name}{The name of the stochastic frontier model estimated.}
\item{formula}{The formula used in the model specification.}
\item{exp_u_hat}{Predicted technical efficiency, E[exp(-u)|e]. Available for models: NHN, NHN_Z, NR, NG, NNAK, NE, NTN, NU, NGE, NLN, NW, THT and tHN.}
\item{u_hat}{Predicted inefficiency, E[u|e] (Jondrow et al. 1982). Available for NE, NTN, NU, NGE, NLN, NW, THT and tHN.}
\item{wrong_skew, sigma_u_at_bound, residual_m3}{Wrong-skewness diagnostics. \code{residual_m3} is the third central moment of the OLS residuals; \code{wrong_skew} is \code{TRUE} when that is non-negative, which for a production frontier means the moment equations admit no interior solution; \code{sigma_u_at_bound} is \code{TRUE} when the fitted one-sided scale has collapsed to the zero boundary. When both are \code{TRUE} the boundary \emph{is} the maximum likelihood estimate -- the Type I failure of Olson, Schmidt and Waldman (1980), studied by Waldman (1982) -- and not a numerical failure, so \code{sfm} reports it rather than bounding \eqn{\sigma_u} away from zero. Read it as no evidence of inefficiency in these data, and treat the efficiency scores as uninformative. A collapse with correctly skewed residuals would be a different matter and does not occur in testing.}
\item{thn_sigma_u_at_bound}{Logical, \code{tHN} only. \code{TRUE} when the fitted \eqn{\sigma_u} collapsed onto zero, so the model reports essentially no inefficiency. Accompanied by a warning. Treat the efficiency scores as uninformative when this is set.}
\item{thn_starts}{List, \code{tHN} only. Multi-start diagnostic: how many starting values were tried and converged, the log-likelihood reached from each, and \code{n_distinct}, the number of separate optima found. \code{n_distinct > 1} means the likelihood surface is multimodal for these data and the single returned fit is only the best of those found.}
\item{estimator}{\code{"cols"} for a corrected-OLS fit; absent for maximum-likelihood fits.}
\item{wrong_skew}{Logical, \code{estimator = "cols"} only. \code{TRUE} when the OLS residuals were skewed the wrong way and the moment equations had no admissible solution.}
\item{residual_moments}{The second, third and fourth central moments of the OLS residuals, \code{estimator = "cols"} only.}
\item{cols_boot_draws}{The matrix of bootstrap replicates (\code{cols_boot} rows, one column per parameter) when \code{cols_boot > 0}; \code{NULL} otherwise.}
\item{ng_starts}{Returned by \code{NG} only. A record of the starting-value search: how many candidates were built (\code{n_tried}), how many were polished before choosing (\code{n_polished}), the log-likelihood at each candidate (\code{loglik_at_start}), and the best value reached (\code{best}). \code{NULL} when \code{start_val} was supplied. The normal-gamma likelihood determines \eqn{E[u] = \mu\sigma_u} far better than it determines either factor -- the weak identification Ritter and Simar (1997) describe -- so the candidates hold \eqn{E[u]} at a moment estimate and sweep the shape along that ridge rather than searching across it. Inspect this if a normal-gamma fit looks implausible: widely differing log-likelihoods among the candidates indicate the surface is more awkward than usual for that data set.}
\item{sd_exp_u_hat}{Standard deviation of exp(-u) given the residual, following Tancredi (2002, section 2.2). Returned by \code{THT} only. It is the quantity that distinguishes the skew-t frontier from the half-normal one: for a large \emph{positive} residual the half-normal model concentrates exp(-u) at 1 and reports near-zero uncertainty, while the skew-t model treats the point as an outlier, widens this spread and does not push the efficiency estimate to 1.}
\item{med_u_hat}{Predicted technical efficiency (median values). Available only for the NHN model.}
\item{coefficients}{A vector of estimated parameters.}
\item{std.errors}{A vector of standard errors for the estimated parameters (NA if \code{optHessian = FALSE}). Sandwich-form (not the ordinary Hessian-inverse form) when \code{robust != "mle"} -- see Details.}
\item{t.values}{A vector of t-values for the estimated parameters (NA if \code{optHessian = FALSE}).}
\item{call}{The matched call.}
\item{het}{Only present when \code{vhet}, \code{uhet} or \code{muhet} was given: the family fitted, the link used, the sizes of the four parameter blocks, the formulas as supplied, and (for \code{"NTN"}) the fitted per-observation \eqn{\mu}.}
\item{sigma_u, sigma_v}{Only present for a heteroskedastic fit: the fitted per-observation scales.}
\item{v_spec, mu_spec}{Only present for a heteroskedastic fit: the noise and pre-truncation-mean blocks, in the same shape as \code{z_spec}.}
\item{robust}{Only present for \code{model_name = "NHN"}: the \code{robust} argument used to fit this model (\code{"mle"}, \code{"mlqe"}, \code{"psi"}, or \code{"mdpd"}).}
\item{robust_c}{Only present for \code{model_name = "NHN"}: the power parameter (\code{c_mlqe}/\code{eta}/\code{alpha}, whichever applies) used, or \code{NA} when \code{robust = "mle"}.}
}
\examples{
\donttest{
library(sfa)     

cs_data_trial   <- data_gen_cs(N= 1000, rand   = 1,  sig_u  = 0.3, sig_v  = 0.3, 
cons   = 0.5,       beta1  = 0.5,   beta2  = 0.5, a      = 4, mu     = 1)

cs.nhnz     <-  sfm(formula    = y_pcs_z ~ x1 +x2| z,    model_name = "NHN_Z",
                    data       = cs_data_trial,          PSopt      = TRUE)

## Robust divergence-based estimation (NHN only): downweights the influence
## of outlying observations relative to ordinary MLE. std.errors/t.values
## are NA for these methods -- see Details.
cs.nhn.mlqe <-  sfm(formula = y_pcs ~ x1 + x2, model_name = "NHN",
                    data = cs_data_trial, robust = "mlqe", c_mlqe = 0.20)
cs.nhn.psi  <-  sfm(formula = y_pcs ~ x1 + x2, model_name = "NHN",
                    data = cs_data_trial, robust = "psi", eta = 0.20)

## Heteroskedasticity in the noise as well as in the inefficiency: `| z`
## still drives sigma_u, `vhet` drives sigma_v. Coefficients come back
## prefixed by block (Zu., Zv.).
cs.het      <-  sfm(formula = y_pcs_z ~ x1 + x2 | z, model_name = "NHN_Z",
                    data = cs_data_trial, vhet = ~ z)

## Battese and Coelli (1995): covariates shift the LOCATION of inefficiency
## rather than its dispersion. Requires the truncated-normal family.
cs.bc95     <-  sfm(formula = y_pcs ~ x1 + x2, model_name = "NTN",
                    data = cs_data_trial, muhet = ~ z)
}
}
