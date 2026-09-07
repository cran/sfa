\name{selsfm}
\alias{selsfm}
\title{Stochastic Frontier Models with Sample Selection or Technology Choice}
\description{
Two models, chosen by \code{model_name}.

\code{"greene"} (the default) fits Greene's (2010) sample-selection stochastic frontier: a normal--half
normal frontier observed only for the subset of units that select into the
sample, where the unobservables driving selection are correlated with the
\emph{noise} in the frontier. Estimation follows the two-step scheme of the
paper -- a probit for the selection equation on the full sample, then maximum
\emph{simulated} likelihood for the frontier on the selected subsample.

\code{"kts"} fits Kumbhakar, Tsionas and Sipilainen (2009), which is a
\emph{different} model rather than a variant. Both technologies are observed,
each with its own frontier and its own two scales, and the technology choice
depends on \bold{inefficiency itself}: \eqn{u} lowers output through the
frontier and shifts the adoption decision through the choice equation. Fitted
by single-step maximum likelihood over the whole sample.
}
\usage{
selsfm(selection, frontier, data, model_name = c("greene", "kts"),
       n_nodes = 64, Nsim = "auto",
       sim_type = c("halton", "sobol", "torus", "uniform"),
       antithetics = FALSE, seed = NULL, inefdec = TRUE,
       maxit.bobyqa = 10000, maxit.psoptim = 1000, maxit.optim = 1000,
       start_val = FALSE, PSopt = FALSE, optHessian = TRUE,
       Method = "L-BFGS-B", verbose = FALSE, rand.psoptim = NULL)
}
\arguments{
\item{selection}{A two-sided formula for the selection equation, \code{d ~ z1 + z2}.
  The response must be binary (0/1, logical, or a two-level factor) and is
  modelled by probit over \emph{all} rows of \code{data}.}
\item{frontier}{A two-sided formula for the frontier, \code{y ~ x1 + x2}, fitted
  on the selected rows only. The dependent variable need not be observed where
  \code{d == 0}; those rows may be \code{NA}.}
\item{model_name}{\code{"greene"} (default) or \code{"kts"}; see
  \sQuote{Description} and \sQuote{The KTS model} below.}
\item{n_nodes}{Gauss--Legendre nodes for the integral over \eqn{u}, used by
  \code{"kts"} only. The half-normal integral is mapped to \eqn{[0,1]} by its
  own CDF, so no truncation point is needed.}
\item{data}{A \code{data.frame} containing both equations' variables, with the
  unselected rows present. Unlike the rest of the package neither formula takes
  a \code{|} segment -- the two equations are separate arguments -- and a pipe
  is an error rather than being silently ignored.}
\item{Nsim}{Number of simulation draws per selected observation. \code{"auto"}
  (the default) uses \code{max(200, ceiling(3*sqrt(n_selected)))}. Simulated ML
  is consistent only if the draw count grows with the sample size; a numeric
  value below the automatic floor warns.}
\item{sim_type}{Low-discrepancy sequence for the draws, passed to the same
  generator the other simulated-ML models use.}
\item{antithetics}{Mirror the draws to halve the number of independent points.}
\item{seed}{Optional seed for the draws. The RNG state is saved and restored, so
  a fit does not disturb the caller's stream.}
\item{inefdec}{\code{TRUE} (default) for a production frontier, where
  inefficiency is subtracted; \code{FALSE} for a cost frontier.}
\item{maxit.bobyqa, maxit.psoptim, maxit.optim}{Iteration caps for the three
  optimizer stages.}
\item{start_val}{Name the starting-value vector in the returned object.}
\item{PSopt}{Run the particle-swarm stage between BOBYQA and \code{optim}.}
\item{optHessian}{Compute the Hessian, and with it the standard errors.}
\item{Method}{Method passed to \code{optim} for the final stage.}
\item{verbose}{Report optimizer progress.}
\item{rand.psoptim}{Optional seed for the particle-swarm stage.}
}
\details{
The model is Greene's Equation (10):
\deqn{d_i = 1[\alpha'z_i + w_i > 0], \quad w_i \sim N(0,1)}
\deqn{y_i = \beta'x_i + v_i - u_i, \quad \mathrm{observed\ only\ when\ } d_i = 1}
\deqn{u_i = \sigma_U|U_i|, \quad U_i \sim N(0,1), \qquad v_i = \sigma_V V_i, \quad V_i \sim N(0,1)}
\deqn{(w_i, v_i) \sim N_2[(0,0), (1, \rho\sigma_V; \rho\sigma_V, \sigma_V^2)]}

Selection operates through the \emph{noise} \eqn{v}, not through the
inefficiency \eqn{u}. That is the substantive assumption separating this model
from Kumbhakar, Tsionas and Sipilainen (2009), where the selection disturbance
is correlated with \eqn{u} instead; the two answer different questions and are
not nested.

Because \eqn{u} cannot be integrated out in closed form, the log-likelihood is
simulated (Greene's Equation (15)): for each selected observation \eqn{R} draws
of \eqn{|U_{ir}|} are taken and
\deqn{\log L = \sum_{d_i=1} \log \frac{1}{R} \sum_{r=1}^{R}
  \phi\!\left(\frac{v_{ir}}{\sigma_V}\right)\frac{1}{\sigma_V}
  \Phi\!\left(\frac{\rho v_{ir}/\sigma_V + \alpha'z_i}{\sqrt{1-\rho^2}}\right)}
with \eqn{v_{ir} = y_i - \beta'x_i + \sigma_U|U_{ir}|}. Everything is evaluated
in logs and accumulated by log-sum-exp, since the summand underflows in the
tails well before the total does.

\bold{A note on the published formula.} Equations (11) and (13)--(15) of Greene
(2010) write \eqn{\sigma_\varepsilon} in the denominator of the \eqn{\Phi}
argument. It must be \eqn{\sigma_V}: the term being standardized is
\eqn{v_{ir} = y_i - \beta'x_i + \sigma_U|U_{ir}|}, which \emph{is} the noise,
and the conditional distribution of \eqn{w} given \eqn{v} for a bivariate normal
with \eqn{\mathrm{Var}(w) = 1} and correlation \eqn{\rho} is
\eqn{N(\rho v/\sigma_V,\, 1-\rho^2)}. This implementation uses \eqn{\sigma_V}.

\bold{Standard errors are conditional on the first stage.} Following the paper's
two-step scheme, \eqn{\alpha} is taken as given in the frontier likelihood. The
reported standard errors come from the second-stage Hessian and therefore do
\emph{not} include the sampling variability of \eqn{\hat\alpha}; Greene applies
the Murphy--Topel correction at this point, which is not implemented here.
Treat them as a lower bound, and bootstrap over both stages when the precision
matters.
}
\value{
An object of class \code{"sfareg"}. As elsewhere in the package \code{out} is a
\code{p x 3} matrix, one row per parameter, with columns \code{par},
\code{st_err} and \code{t-val}; the parameters are ordered
\code{sigma_u}, \code{sigma_v}, \code{rho}, then the frontier coefficients.
Alongside the usual components:
\item{alpha}{First-stage probit coefficients.}
\item{probit}{The fitted first-stage object.}
\item{jlms}{\eqn{E[u_i \mid \varepsilon_i]} by the simulated Bayes rule of
  Greene's Equation (21), computed over the same draws the likelihood used.}
\item{efficiency}{\code{exp(-jlms)}.}
\item{selection_formula}{The selection equation.}
\item{n_selected, n_total}{Selected and total observation counts.}
\item{Nsim, sim_type}{The draw count and sequence actually used.}
}
\section{The KTS model}{
For technology \eqn{I_i \in \{0,1\}},

\deqn{y_i \mid I_i, u_i \sim N(x_i'\beta_{I_i} - u_i,\ \sigma_{v,I_i}^2), \qquad
P(I_i = 1 \mid u_i) = \Phi(z_i'\gamma + \delta u_i), \qquad
u_i \mid I_i \sim N^+(0, \sigma_{u,I_i}^2).}

\eqn{\delta} is the parameter the model exists for: it is the \emph{dual role}
of inefficiency, and \eqn{\delta > 0} means the less efficient are more likely
to choose technology 1. Neither two-step order works -- the choice equation
cannot be a probit because \eqn{u} is unobserved, and fitting the frontiers
first ignores the endogeneity of the choice -- which is why the paper estimates
in one step.

Because \eqn{u}'s distribution depends on \eqn{I} while \eqn{I}'s probability
depends on \eqn{u}, the marginal density of \eqn{u} is a \emph{mixture} of the
two regimes' half-normals, with weight \eqn{\phi_i} solving the paper's
equation (9). The observation's density is then the same mixture of two
one-dimensional integrals, its equation (14).

\bold{Recovery.} On data simulated from this exact model with
\eqn{\delta = 0.9}, the frontier coefficients, the four scales and
\eqn{\delta} all recover, and improve with the sample:

\tabular{llll}{
 \bold{parameter} \tab \bold{n = 1500} \tab \bold{n = 6000} \tab \bold{truth} \cr
 \eqn{\delta} \tab 1.361 \tab 0.922 \tab 0.900 \cr
 \eqn{\sigma_{u,0}} \tab 0.837 \tab 0.843 \tab 0.800 \cr
 \eqn{\sigma_{u,1}} \tab 0.632 \tab 0.640 \tab 0.600 \cr
 choice intercept \tab -0.511 \tab -0.239 \tab -0.200
}

\eqn{\delta} is the parameter that needs the most data, as one would expect of
a coefficient identified only through the coupling between the two equations,
but it is properly identified: holding the rest at the truth, the mean
log-likelihood falls monotonically away from \eqn{\delta = 0.9} in both
directions, by about 0.069 per observation across \eqn{[0.1, 1.8]}.
}

\references{
Kumbhakar, S.C., Tsionas, E.G. and Sipilainen, T. (2009). Joint estimation of
technology choice and technical efficiency: an application to organic and
conventional dairy farming. \emph{Journal of Productivity Analysis}, 31(3),
151--161.

Greene, W. (2010). A stochastic frontier model with correction for sample
selection. \emph{Journal of Productivity Analysis}, \bold{34}(1), 15--24.

Heckman, J. (1979). Sample selection bias as a specification error.
\emph{Econometrica}, \bold{47}(1), 153--161.
}
\seealso{\code{\link{sfm}} for the frontier without a selection correction.}
\examples{
\donttest{
set.seed(1)
n  <- 800
z1 <- rnorm(n); x1 <- rnorm(n)
w  <- rnorm(n)
d  <- as.numeric(0.2 + 0.8 * z1 + w > 0)
v  <- 0.5 * (0.6 * w + sqrt(1 - 0.6^2) * rnorm(n))
y  <- 0.5 + 0.8 * x1 + v - abs(rnorm(n))
y[d == 0] <- NA
dat <- data.frame(y = y, x1 = x1, z1 = z1, d = d)

fit <- selsfm(selection = d ~ z1, frontier = y ~ x1, data = dat, Nsim = 100)
fit$out
}
}
\keyword{models}
