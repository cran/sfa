\name{meanefficiency}
\alias{meanefficiency}
\title{Model-Implied Mean Efficiency}
\description{
The mean efficiency implied by the fitted distribution of inefficiency,
\eqn{E[\exp(-U)]}, and the mean efficiency among the most efficient \eqn{p} of
that distribution.
}
\usage{
meanefficiency(object,
               p = c(0.01, 0.05, 0.10, 0.25, 0.50, 0.75, 0.90, 0.95, 0.99),
               use_closed_form = TRUE)
}
\arguments{
\item{object}{An \code{"sfareg"} fit.}
\item{p}{Proportions in \eqn{(0, 1]} at which to report supra-percentile mean
  efficiency. \code{p = 1} is the whole distribution and returns the
  unconditional mean.}
\item{use_closed_form}{Use the closed-form expression where one exists. Setting
  it to \code{FALSE} forces numerical integration of the same density, which is
  how the closed forms are checked.}
}
\details{
\bold{A different quantity from the reported average.} \code{print} and
\code{summary} report the average of the per-observation conditional
predictions \eqn{E[\exp(-u_i) \mid \varepsilon_i]}. That is a property of the
sample in hand. \eqn{E[\exp(-U)]} is a property of the fitted \emph{model}, and
it is the number that can be compared across studies.

\bold{Closed forms.} Available for the half-normal
(\eqn{2e^{\sigma_U^2/2}\Phi(-\sigma_U)}), exponential
(\eqn{1/(1+\sigma_U)}), truncated normal, Rayleigh, gamma
(\eqn{(1+\theta)^{-k}}), uniform (\eqn{(1-e^{-\theta})/\theta}) and generalized
exponential (\eqn{2/((\sigma_U+1)(\sigma_U+2))}). The lognormal, Weibull and
Nakagami have none in elementary terms and are integrated numerically;
\code{method} in the returned list says which route was taken.

Every closed form is checked against numerical integration of the same density
in the package's tests, and the density is separately checked against its own
quantile function by simulation. The first catches an algebraic slip; the
second catches a parameter being read off the fit incorrectly, which would
otherwise make both routes agree and both be wrong.

\bold{Supra-percentile means} are always integrated. Efficiency is decreasing
in \eqn{u}, so the most efficient \eqn{p} are the smallest \eqn{p} of \eqn{u}:
the reported value is \eqn{p^{-1}\int_0^{Q(p)} e^{-u}f(u)\,du}. It necessarily
exceeds the unconditional mean and falls as \eqn{p} rises.
}
\value{
A list with \code{model}, \code{distribution}, the \code{parameters} used,
\code{mean_efficiency}, the \code{method} used to obtain it, and \code{supra},
a data frame of \code{p} against mean efficiency.
}
\seealso{\code{\link{efficiency}} for the per-observation predictors.}
\examples{
\donttest{
d <- data_gen_cs(N = 400, rand = 1, sig_u = 1, sig_v = 0.3, cons = 0.5,
                 beta1 = 0.5, beta2 = 0.5, a = 1, mu = 0.5)
fit <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = d)

r <- meanefficiency(fit)
r$mean_efficiency
r$supra
}
}
\keyword{models}
