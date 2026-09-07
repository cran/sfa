\name{density_weights}
\alias{density_weights}
\title{Density-power weights from a robust frontier fit}
\description{
Returns the weight each observation receives in the estimating equation of a
maximum \eqn{L_q}-likelihood or minimum density-power divergence fit. The weight
is the fitted density raised to the power \eqn{c}, normalised so the largest is
one.
}
\usage{
density_weights(object, sigma_v = NULL, sigma_u = NULL, c = NULL,
                normalize = TRUE)
}
\arguments{
\item{object}{A fitted frontier model, or a numeric vector of composed
residuals.}
\item{sigma_v, sigma_u}{Scale parameters. Required when \code{object} is a
numeric vector; taken from the fit otherwise.}
\item{c}{Robustness tuning parameter. Taken from the fit when available.}
\item{normalize}{Logical; scale so the maximum weight is one. Default
\code{TRUE}.}
}
\value{A numeric vector of weights in \eqn{(0, 1]}, one per observation.}
\section{What this diagnostic can and cannot find}{
The weight is a function of the fitted density at the observation, so it
responds to observations lying far from the fitted surface. It has no purchase
on an observation that is wrong in a \emph{regressor} and, because the fitted
surface bends toward it, ends up close to that surface. A mis-recorded input can
therefore be highly influential and still receive a weight near one. This is not
a defect of the implementation: the bounded influence of these estimators is
bounded with respect to \emph{response} contamination, conditional on a fixed or
bounded design. Pair this with a leverage diagnostic, and with case deletion,
before concluding that an observation is or is not driving a result.
}
\examples{
set.seed(1)
e <- rnorm(200, 0, 0.3) - abs(rnorm(200, 0, 0.6))
w <- density_weights(e, sigma_v = 0.3, sigma_u = 0.6, c = 0.217)
summary(w)
head(order(w), 10)   # the observations the estimator discounts most
}
\references{
Bernstein, D.H., Parmeter, C.F. and Wright, I.A. (2026). On Robust Estimation of the Stochastic Frontier Model. Working paper.
}
\seealso{\code{\link{hscore_select}}, \code{\link{calibrate_c}}}
