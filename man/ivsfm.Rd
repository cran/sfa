\name{ivsfm}
\alias{ivsfm}
\title{Stochastic Frontier Models with Endogenous Regressors}
\description{
Fits the stochastic frontier model of Amsler, Prokhorov and Schmidt (2016) in
which one or more regressors is endogenous, in the simultaneous-equation sense:
correlated with the statistical noise \eqn{v}. Three estimators are provided --
full information maximum likelihood, a two-step control function, and corrected
2SLS. With \code{uhet} the model becomes that of Amsler, Prokhorov and Schmidt
(2017), in which environmental variables entering the inefficiency scale may
themselves be endogenous.
}
\usage{
ivsfm(formula, endogenous, instruments, data,
      model_name = c("IVLIML", "IVCF", "C2SLS"), uhet = NULL,
      inefdec = TRUE, maxit.bobyqa = 10000, maxit.psoptim = 1000,
      maxit.optim = 1000, start_val = FALSE, PSopt = FALSE,
      optHessian = TRUE, Method = "L-BFGS-B", verbose = FALSE,
      rand.psoptim = NULL)
}
\arguments{
\item{formula}{The structural frontier, \code{y ~ x1 + x2}, listing all
  regressors including the endogenous ones.}
\item{endogenous}{One-sided formula naming the endogenous variables,
  \code{~ x2}. Each must appear in \code{formula} or in \code{uhet}.}
\item{instruments}{One-sided formula of the \emph{excluded} instruments,
  \code{~ w1 + w2}. The included exogenous regressors are added automatically,
  so only the outside instruments are listed here. There must be at least as
  many as there are endogenous variables.}
\item{data}{A \code{data.frame} holding every variable in the four formulas.}
\item{model_name}{\code{"IVLIML"} (the default), \code{"IVCF"} or
  \code{"C2SLS"}; see \sQuote{Details}.}
\item{uhet}{Optional one-sided formula of environmental variables, giving
  \eqn{\sigma_{U,i} = \sigma_U \exp(q_i'\delta)}. Supplying it selects the
  2017 model; omitting it gives the 2016 one. Variables may appear in both
  \code{uhet} and \code{endogenous}.}
\item{inefdec}{\code{TRUE} (default) for a production frontier, where
  inefficiency is subtracted; \code{FALSE} for a cost frontier.}
\item{maxit.bobyqa, maxit.psoptim, maxit.optim}{Iteration caps for the three
  optimizer stages. Ignored by \code{"C2SLS"}, which is closed-form.}
\item{start_val}{Name the starting-value vector in the returned object.}
\item{PSopt}{Run the particle-swarm stage between BOBYQA and \code{optim}.}
\item{optHessian}{Compute the Hessian, and with it the standard errors.}
\item{Method}{Method passed to \code{optim} for the final stage.}
\item{verbose}{Report optimizer progress.}
\item{rand.psoptim}{Optional seed for the particle-swarm stage.}
}
\details{
Unlike the rest of the package neither \code{formula} nor \code{endogenous} nor
\code{instruments} takes a \code{|} segment. A pipe means \dQuote{variance
determinant} in \code{\link{sfm}} and \code{\link{psfm}}, and reusing it here
for \pkg{ivreg}-style instrument syntax would give the same character two
meanings across entry points, so a pipe is an error rather than being
reinterpreted.

\bold{The model.} Writing \eqn{p_i} for the endogenous block and \eqn{z_i} for
the instruments,
\deqn{y_i = \beta'x_i + v_i - u_i, \qquad p_i = \Pi'z_i + \xi_i}
\deqn{u_i = \sigma_{U,i}|U_i|, \quad U_i \sim N(0,1), \qquad
      (v_i, \xi_i') \sim N(0, \Omega)}
with \eqn{u_i} independent of \eqn{(v_i, \xi_i)}. Endogeneity is the
correlation between \eqn{\xi} and \eqn{v}: the regressors are correlated with
the \emph{noise}, not with the inefficiency.

Conditioning \eqn{v} on \eqn{\xi} turns this into an ordinary normal--half
normal frontier in a shifted residual with a smaller noise variance,
\deqn{\mu_{c,i} = \Sigma_{v\xi}\Sigma_{\xi\xi}^{-1}\xi_i, \qquad
      \sigma_c^2 = \sigma_V^2 - \Sigma_{v\xi}\Sigma_{\xi\xi}^{-1}\Sigma_{\xi v}}
which is what makes the likelihood closed-form and needs no simulation.

\bold{The estimators.}
\describe{
  \item{\code{"IVLIML"}}{APS (2016) Equation (13): maximum likelihood over the
    frontier and the reduced form jointly, so \eqn{\Pi} and
    \eqn{\Sigma_{\xi\xi}} are estimated alongside \eqn{\beta}. Efficient, and
    its standard errors are correct as reported.}
  \item{\code{"IVCF"}}{The two-step control function of Kutlu (2010), described
    in APS (2016) Section 4.4. \eqn{\Pi} and \eqn{\Sigma_{\xi\xi}} are fixed at
    their reduced-form least-squares values and only the frontier block is
    maximized. Cheaper, still consistent, but it discards the information about
    the reduced form carried in the frontier likelihood; its conventional
    standard errors understate the true ones unless \eqn{\Sigma_{v\xi} = 0}.
    Measured over 30 replications at \eqn{n = 1000}: with \eqn{\rho = 0.6}
    the mean reported standard error on the endogenous coefficient is 0.912 of
    the Monte Carlo standard deviation, while \code{"IVLIML"}'s is 1.005; with
    \eqn{\rho = 0} the two are identical, as the theory says they should be.
    Bootstrap, or use \code{"IVLIML"}, when the standard errors matter.}
  \item{\code{"C2SLS"}}{APS (2016) Section 4.1: 2SLS, then \eqn{\sigma_U} and
    \eqn{\sigma_V} from the second and third moments of the residuals, then the
    intercept corrected by \eqn{\sqrt{2/\pi}\,\hat\sigma_U}. The residuals use
    the \emph{actual} endogenous regressors, not their fitted values, which the
    paper flags explicitly. Moment-based rather than likelihood-based, so the
    fit carries no \code{$opt} and \code{\link{logLik}} returns \code{NA} with
    a warning -- the same footing as \code{psfm}'s \code{"GTRE_SEQ1"}. It
    estimates no \eqn{\rho}, since it never models the correlation it corrects
    for.}
}

\bold{Inefficiency.} \code{jlms} implements APS (2016) Equations (14)--(15),
which condition on \eqn{\xi} as well as on \eqn{\varepsilon}. Since \eqn{\xi}
is correlated with \eqn{v} it carries information about \eqn{u} even though
\eqn{u} is independent of \eqn{\xi}, so this is a strictly better predictor
than the plain Jondrow et al. (1982) formula.

\bold{What is assumed, and what is not covered.} \eqn{u} is assumed independent
of \eqn{(v, \xi)}. The case where the reduced-form error is correlated with the
inefficiency as well requires the copula likelihood of APS (2016) Section 4.5
and is \emph{not} implemented; the paper notes there that two intuitive
approaches to that case do not work.
}
\value{
An object of class \code{"sfareg"}. \code{out} is a \code{p x 3} matrix, one
row per parameter, columns \code{par}, \code{st_err} and \code{t-val}.
Parameters are the frontier coefficients, \code{sigma_u}, \code{sigma_v}, any
\code{delta_} terms from \code{uhet}, and one \code{rho_} per endogenous
variable (absent for \code{"C2SLS"}). Alongside the usual components:
\item{Pi}{Reduced-form coefficient matrix.}
\item{Sigma_xi}{Reduced-form error covariance.}
\item{rho}{Correlation of each reduced-form error with \eqn{v}.}
\item{sigma_c}{The conditional noise standard deviation \eqn{\sigma_c}, which
  is smaller than \code{sigma_v} whenever \eqn{\rho \ne 0}.}
\item{b_2sls}{The 2SLS coefficients, which seed the likelihood and are the
  estimate \code{"C2SLS"} corrects.}
\item{wrong_skew}{Whether the 2SLS residuals are skewed the wrong way.}
\item{jlms, efficiency}{\eqn{E[u_i \mid \varepsilon_i, \xi_i]} and
  \code{exp(-jlms)}.}
\item{vcov_rho}{The full \eqn{q \times q} covariance of \eqn{\hat\rho}, which
  \code{\link{endogeneity_test}} needs for the joint Wald statistic.}

\bold{Standard errors for \eqn{\rho}} arrived in 1.2.0. They were previously
reported as \code{NA}, on the ground that \eqn{\rho = t/\sqrt{1+t't}} has a
non-diagonal Jacobian and a delta-method value computed as though it were
diagonal would be wrong. That much was true, but the Jacobian itself is short --
\eqn{J = (I - \rho\rho')/s} with \eqn{s = \sqrt{1+t't}} -- so
\eqn{\mathrm{Var}(\hat\rho) = J\,\mathrm{Var}(\hat t)\,J'} is available
exactly. Checked against a numerical Jacobian in the tests, and validated
end-to-end by the size of \code{\link{endogeneity_test}}, which is 0.046 at a
nominal 5\% over 1000 replications.
}
\references{
Amsler, C., Prokhorov, A. and Schmidt, P. (2016). Endogeneity in stochastic
frontier models. \emph{Journal of Econometrics}, \bold{190}(2), 280--288.

Amsler, C., Prokhorov, A. and Schmidt, P. (2017). Endogenous environmental
variables in stochastic frontier models. \emph{Journal of Econometrics},
\bold{199}(2), 131--140.

Kutlu, L. (2010). Battese-Coelli estimator with endogenous regressors.
\emph{Economics Letters}, \bold{109}(2), 79--81.
}
\seealso{\code{\link{endogeneity_test}} for the Wald test of whether the
correction was needed at all; \code{\link{sfm}} for the frontier without one.}
\examples{
\donttest{
set.seed(1)
n   <- 600
w1  <- rnorm(n); w2 <- rnorm(n); x1 <- rnorm(n)
eta <- rnorm(n)
v   <- 0.5 * (0.6 * eta + sqrt(1 - 0.6^2) * rnorm(n))
x2  <- 0.9 * w1 - 0.7 * w2 + 0.5 * x1 + eta
y   <- 0.5 + 0.8 * x1 - 0.6 * x2 + v - abs(rnorm(n))
dat <- data.frame(y = y, x1 = x1, x2 = x2, w1 = w1, w2 = w2)

fit <- ivsfm(y ~ x1 + x2, endogenous = ~x2, instruments = ~ w1 + w2,
             data = dat)
fit$out

## Ignoring the endogeneity biases the coefficient on x2 toward zero.
sfm(y ~ x1 + x2, model_name = "NHN", data = dat)$out
}
}
\keyword{models}
