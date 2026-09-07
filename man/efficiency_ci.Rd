\name{efficiency_ci}
\alias{efficiency_ci}
\title{Confidence Intervals for Individual Inefficiency and Efficiency}
\description{
Horrace and Schmidt (1996) intervals for the inefficiency term \eqn{u_i} of a
fitted stochastic frontier model, and the implied intervals for technical
efficiency \eqn{\exp(-u_i)}.
}
\usage{
efficiency_ci(object, level = 0.95, type = c("both", "u", "te"))
}
\arguments{
  \item{object}{an object of class \code{"sfareg"}, as returned by
    \code{\link{sfm}}.}
  \item{level}{coverage of the interval, a single number strictly between 0
    and 1. Defaults to \code{0.95}.}
  \item{type}{which intervals to return: \code{"both"} (the default),
    \code{"u"} for inefficiency only, or \code{"te"} for technical efficiency
    only.}
}
\details{
Conditional on the fitted parameters, \eqn{u_i \mid e_i} is normal with mean
\code{mu_star} and standard deviation \code{sigma_star}, truncated below at
zero -- the same posterior the Jondrow et al. (1982) and Battese and Coelli
(1988) point predictors average over. Inverting it at \eqn{(1-\alpha)/2} and
\eqn{1-(1-\alpha)/2} gives the bounds in closed form, so the interval costs no
estimation beyond the fit itself. The efficiency bounds follow by monotonicity,
with the endpoints swapped, since \eqn{\exp(-u)} is decreasing in \eqn{u}.

\strong{These intervals condition on the estimated parameters.} They describe
where \eqn{u_i} sits given \eqn{e_i} and a \emph{known} frontier, and make no
allowance for sampling error in the slopes or the variance parameters. Horrace
and Schmidt are explicit about this, and it is why the intervals do not narrow
as the sample grows: they measure the irreducible difficulty of splitting a
single residual into noise and inefficiency, not estimation uncertainty. Read a
wide interval as a warning against reading much into that unit's rank.

Available for the models whose posterior really is a truncated normal:
\code{"NHN"}, \code{"NHN_Z"}, \code{"NE"} and \code{"NTN"}. Every other
\code{model_name} has a posterior of a different shape, for which these
formulas do not hold; the function reports that rather than returning a
misleading number.
}
\value{
A data frame with one row per observation, carrying \code{u_lower},
\code{u_hat}, \code{u_upper} and/or \code{te_lower}, \code{te_hat},
\code{te_upper} according to \code{type}, and a \code{"level"} attribute.
\code{u_hat} is the Jondrow et al. (1982) posterior mean and \code{te_hat} the
Battese and Coelli (1988) score -- exactly the values the fit already reports,
taken from the fitted object rather than recomputed.
}
\references{
Horrace, W. C. and Schmidt, P. (1996). Confidence statements for efficiency
estimates from stochastic frontier models. \emph{Journal of Productivity
Analysis}, \bold{7}, 257--282.

Jondrow, J., Lovell, C. A. K., Materov, I. S. and Schmidt, P. (1982). On the
estimation of technical inefficiency in the stochastic frontier production
function model. \emph{Journal of Econometrics}, \bold{19}, 233--238.

Battese, G. E. and Coelli, T. J. (1988). Prediction of firm-level technical
efficiencies with a generalized frontier production function and panel data.
\emph{Journal of Econometrics}, \bold{38}, 387--399.
}
\seealso{
\code{\link{sfm}}, \code{\link{sfa_diagnostics}}
}
\examples{
dat <- data_gen_cs(N = 200, rand = 1, sig_u = 1, sig_v = 0.3,
                   cons = 0.5, beta1 = 0.5, beta2 = 0.5, a = 1, mu = 0.5)
fit <- sfm(y_pcs ~ x1 + x2, data = dat, model_name = "NHN")

ci <- efficiency_ci(fit, level = 0.90)
head(ci)

## How much of an efficiency ranking is real? Compare the spread of the point
## predictions against the width of a single unit's interval.
diff(range(ci$te_hat))
median(ci$te_upper - ci$te_lower)
}
\keyword{models}
