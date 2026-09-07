\name{hscore}
\alias{hscore}
\title{Hyvarinen score for a Normal--Half-Normal frontier fit}
\description{
Evaluates the Hyvarinen score used to choose the robustness tuning parameter of
the maximum \eqn{L_q}-likelihood and minimum density-power divergence
estimators, following Sugasawa and Yonekura (2021). Lower is better.
}
\usage{
hscore(e, sigma_v, sigma_u, c, stable = TRUE)
}
\arguments{
\item{e}{Numeric vector of composed residuals \eqn{y - x'\beta}.}
\item{sigma_v, sigma_u}{Positive scalars, the noise and inefficiency scale
parameters.}
\item{c}{Robustness tuning parameter, \eqn{c = 1-q} for MLqE and
\eqn{c = \alpha} for MDPD. \code{c = 0} is the maximum likelihood endpoint.}
\item{stable}{Logical. \code{TRUE} (default) evaluates in log space;
\code{FALSE} uses the natural-scale expression.}
}
\value{A single numeric value, the mean observationwise score.}
\details{
The score is
\deqn{H(c) = n^{-1}\sum_i \{2 D''(e_i) + D'(e_i)^2\},\qquad D = (f^c-1)/c,}
with derivatives taken with respect to the composed residual.

The score is evaluated in logarithms. The natural-scale expression contains the
fitted density raised to the power \eqn{c-2}; since \eqn{c<1} that exponent is
close to \eqn{-2}, so a single observation whose density underflows makes the
score non-finite. The failure is silent and directional: the tuning candidates
that remain evaluable are those closest to maximum likelihood, so a naive
implementation reports too little robustness. Set \code{stable = FALSE} to
reproduce that behaviour for comparison.
}
\references{
Hyvarinen, A. (2005). Estimation of non-normalized statistical models by score
matching. \emph{Journal of Machine Learning Research} 6, 695--709.

Sugasawa, S. and Yonekura, S. (2021). On selection criteria for the tuning
parameter in robust divergence. \emph{Entropy} 23(9), 1147.

Bernstein, D.H., Parmeter, C.F. and Wright, I.A. (2026). On Robust Estimation of the Stochastic Frontier Model. Working paper.
}
\examples{
set.seed(1)
e <- rnorm(500, 0, 0.3) - abs(rnorm(500, 0, 0.6))
hscore(e, sigma_v = 0.3, sigma_u = 0.6, c = 0.20)

## the two evaluations agree wherever both are computable
hscore(e, 0.3, 0.6, 0.20) - hscore(e, 0.3, 0.6, 0.20, stable = FALSE)

## but only one of them survives an observation in the far tail
e2 <- c(e, -60)
hscore(e2, 0.3, 0.6, 0.15)                  # finite
hscore(e2, 0.3, 0.6, 0.15, stable = FALSE)  # not finite
}
\seealso{\code{\link{hscore_select}}, \code{\link{calibrate_c}},
\code{\link{density_weights}}}
