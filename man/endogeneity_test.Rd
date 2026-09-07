\name{endogeneity_test}
\alias{endogeneity_test}
\alias{print.sfa_endog_test}
\title{Were the regressors endogenous after all?}
\description{
A Wald test of \eqn{H_0\!: \rho = 0} on an \code{\link{ivsfm}} fit -- the hypothesis that the noise \eqn{v} is uncorrelated with the reduced-form errors of the endogenous regressors, so that the endogeneity correction was not needed and a plain \code{\link{sfm}} fit would have been consistent. This is the Durbin--Wu--Hausman idea carried into the stochastic frontier, and is Equation (25) of Hou, Ramalho and Roseta-Palma (2025).
}
\usage{endogeneity_test(object, level = 0.05)}
\arguments{
\item{object}{A fit from \code{\link{ivsfm}} with \code{model_name = "IVLIML"} or \code{"IVCF"}. \code{"C2SLS"} is a moment estimator and carries no correction parameter, so it is refused.}
\item{level}{Size used for the reported verdict.}
}
\value{
An object of class \code{c("sfa_endog_test", "htest")} carrying \code{statistic} (the Wald \eqn{W}), \code{parameter} (its degrees of freedom, the number of endogenous regressors), \code{p.value}, the estimated \code{rho} with its standard errors and z-ratios, and \code{reject}.
}
\details{
\bold{The statistic.} With \eqn{\hat\rho} the vector of correlations between \eqn{v} and the reduced-form errors,

\deqn{W = \hat\rho' \,[\widehat{\mathrm{Var}}(\hat\rho)]^{-1}\, \hat\rho \;\sim\; \chi^2(q),}

\eqn{q} the number of endogenous regressors. With a single endogenous regressor it is exactly the squared z-ratio already printed for \code{rho_} in the fit.

\bold{Why this is an ordinary chi-square.} \eqn{\rho} lives in the open unit ball, so \eqn{\rho = 0} is an \emph{interior} point of the parameter space. That is the opposite of the situation \code{\link{inefficiency_test}} faces, where \eqn{\sigma_u = 0} sits on a boundary and the null distribution is a chi-bar-square mixture. No mixture is needed here.

\bold{Standard errors for \eqn{\rho}.} \code{ivsfm} parameterises the correlations as \eqn{\rho = t/\sqrt{1 + t't}}, which keeps them inside the unit ball for any real \eqn{t}. Before version 1.2.0 the standard errors of \eqn{\rho} were reported as \code{NA} on the ground that this map has a non-diagonal Jacobian. It does, but the Jacobian is short --

\deqn{\partial \rho_i/\partial t_j = \delta_{ij}/s - t_i t_j / s^3, \qquad s = \sqrt{1 + t't},}

that is \eqn{J = (I - \rho\rho')/s} -- so \eqn{\mathrm{Var}(\hat\rho) = J\,\mathrm{Var}(\hat t)\,J'} by the delta method. At the null \eqn{J = I}, which is what makes the test well calibrated there. \code{ivsfm} now reports these standard errors and carries the full \eqn{\mathrm{Var}(\hat\rho)} as \code{$vcov_rho}, which is what the \emph{joint} test needs.

\bold{Reading a non-rejection.} Failing to reject is not evidence for exogeneity. \eqn{\rho} is identified through the instruments, so weak instruments make \eqn{\hat\rho} imprecise and the test powerless; the print method says so rather than letting the p-value speak alone.
}
\references{
Hou, Z., Ramalho, J. J. S. and Roseta-Palma, C. (2025). Dealing with endogeneity in stochastic frontier models: A comparative assessment of estimators. \emph{Energy Economics} 151, 108922.

Amsler, C., Prokhorov, A. and Schmidt, P. (2016). Endogeneity in stochastic frontier models. \emph{Journal of Econometrics} 190(2), 280-288.

Kutlu, L. (2010). Battese-Coelli estimator with endogenous regressors. \emph{Economics Letters} 109(2), 79-81.
}
\seealso{\code{\link{ivsfm}} for the fit, \code{\link{inefficiency_test}} for the boundary-null case, \code{\link{spec_test}} and \code{\link{moment_range}} for the distributional assumptions.}
\examples{
set.seed(9)
n  <- 600
z  <- rnorm(n); x1 <- rnorm(n); eps <- rnorm(n)
x2 <- z + eps
v  <- 0.6 * eps + sqrt(1 - 0.36) * rnorm(n)
d  <- data.frame(y = 0.5 * x1 + 0.5 * x2 + v - abs(rnorm(n)),
                 x1 = x1, x2 = x2, z = z)

fit <- ivsfm(y ~ x1 + x2, endogenous = ~ x2, instruments = ~ z,
             data = d, model_name = "IVCF")
endogeneity_test(fit)
}
\keyword{htest}
