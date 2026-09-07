\name{calibrate_c}
\alias{calibrate_c}
\title{Fixed robustness parameter by weight matching}
\description{
Chooses the robustness tuning parameter so that a residual a given number of
scale units from the reference point retains a common share of its clean-data
influence, comparably across the robust criteria. This is the transparent
alternative to the data-driven selection of \code{\link{hscore_select}}.
}
\usage{
calibrate_c(sigma_v, sigma_u, method = c("mlqe", "psi", "mdpd", "all"),
            target = 0.10, k = 3, range = c(1e-4, 0.6))
}
\arguments{
\item{sigma_v, sigma_u}{Positive scalars at which to calibrate, typically the
estimates from a preliminary fit or the values of a simulation design.}
\item{method}{Robust criterion, or \code{"all"} for a named vector.}
\item{target}{Share of clean-data influence the reference residual should
retain. Default \code{0.10}.}
\item{k}{Reference residual in scale units. Default \code{3}.}
\item{range}{Search interval for the tuning parameter.}
}
\value{A named numeric vector of calibrated tuning parameters, carrying an
attribute \code{"roots"} that lists every solution found for each criterion.}
\details{
The calibration matches the norm of the observationwise estimating-equation
contribution in the intercept, \eqn{\log\sigma_v} and \eqn{\log\sigma_u}
directions. Matching on the intercept direction alone would be vacuous: there
the criteria share the density-power factor \eqn{f^{c-1}f'}, so it would force
the same tuning value for all of them by construction. They differ through the
Fisher-consistency correction, which enters the scale directions.

Note that \eqn{\sigma=(\sigma_v^2+\sigma_u^2)^{1/2}} is the scale parameter of
the composed density, not the standard deviation of the composed error, which is
\eqn{\{\sigma_v^2+(1-2/\pi)\sigma_u^2\}^{1/2}}. The argument \code{k} is measured
in the former.

\bold{The ratio is not monotone, and there can be two answers.} For the
Fisher-consistency-corrected criteria the influence ratio falls below the target
and is then pulled back up again by the integral correction, so it crosses the
target twice. At \eqn{\sigma_v=0.30}, \eqn{\sigma_u=0.60} the \eqn{\Psi}/MDPD
solutions are \eqn{c=0.2387} and \eqn{c=0.5385}. The \bold{larger} is returned:
it is the calibration reported by Bernstein, Parmeter and Wright, and it lies on
the far side of the dip, so the target down-weighting holds there rather than
being passed through on the way. Both are in the \code{"roots"} attribute; the
non-monotonicity is reported rather than resolved silently. MLqE
has a single solution, \eqn{c=0.2166}.

\code{"psi"} and \code{"mdpd"} always calibrate to the same value: the two
objectives differ by the factor \eqn{(1+c)}, which is constant in the parameters
and cancels from the ratio being matched.
}
\examples{
## the calibration used in Bernstein, Parmeter and Wright
calibrate_c(sigma_v = 0.30, sigma_u = 0.60, method = "all")
}
\references{
Bernstein, D.H., Parmeter, C.F. and Wright, I.A. (2026). On Robust Estimation of the Stochastic Frontier Model. Working paper.
}
\seealso{\code{\link{hscore_select}}, \code{\link{density_weights}}}
