\name{pcomposed}
\alias{pcomposed}
\alias{dcomposed}
\title{The distribution of the composed error}
\description{
Density and distribution function of the stochastic frontier composed error \eqn{\varepsilon = v - u} (production) or \eqn{v + u} (cost), with \eqn{v \sim N(0, \sigma_v^2)} and \eqn{u \sim N^+(0, \sigma_u^2)} independent. The density is a skew-normal and closed form; the distribution function is not, and is obtained here by quadrature that stays accurate far into the tails.
}
\usage{dcomposed(x, sigma_u, sigma_v, inefdec = TRUE, log = FALSE)

pcomposed(q, sigma_u, sigma_v, inefdec = TRUE, lower.tail = TRUE,
          log.p = FALSE, method = c("quadrature", "simulate"),
          n_nodes = 128, R = 1e6, seed = NULL)}
\arguments{
\item{x, q}{Numeric vector of quantiles.}
\item{sigma_u}{Standard deviation of the pre-truncation normal behind \eqn{u}. Zero is allowed and is the no-inefficiency boundary, where the composed error is just the noise.}
\item{sigma_v}{Standard deviation of the noise, strictly positive.}
\item{inefdec}{\code{TRUE} (the default) for a production frontier, \eqn{\varepsilon = v - u}; \code{FALSE} for a cost frontier, \eqn{\varepsilon = v + u}. The two are mirror images.}
\item{log, log.p}{Return the log of the density or probability. Use \code{log.p = TRUE} whenever the result may be extremely small.}
\item{lower.tail}{If \code{TRUE} (the default) probabilities are \eqn{P(\varepsilon \le q)}, otherwise \eqn{P(\varepsilon > q)}. The upper tail is computed \bold{directly} rather than as \eqn{1 - F}; see \sQuote{Details}.}
\item{method}{\code{"quadrature"} (the default) is deterministic Gauss-Legendre. \code{"simulate"} is the estimator of Amsler, Schmidt and Tsay (2019), kept for comparison.}
\item{n_nodes}{Quadrature nodes, at least 16.}
\item{R, seed}{Draws and optional seed for \code{method = "simulate"}. The seed is restored afterwards, so the caller's random stream is untouched.}
}
\value{A numeric vector the length of \code{x} or \code{q}.}
\details{
\bold{Why this exists.} Every likelihood in \pkg{sfa} evaluates the composed error's density and none needs its distribution function, which is why there was not one. Copula models need it: a Gaussian copula on the composed error contains \eqn{\Phi^{-1}(F(\varepsilon))}, so an \eqn{F} that saturates at 0 or 1 does not merely lose precision, it returns an infinity and takes the copula density with it.

\bold{How it is computed.} Conditioning on \eqn{u} leaves a normal distribution function in closed form, so
\deqn{F(q) = E_u\left[\Phi\left((q + su)/\sigma_v\right)\right],}
with \eqn{s = +1} for a production frontier and \eqn{-1} for a cost frontier. The noise is therefore never drawn or integrated, and the tails inherit the accuracy of \code{\link{pnorm}}, which is reliable twenty standard deviations out. The expectation over \eqn{u} is taken by Gauss-Legendre quadrature on \eqn{u = \sigma_u t/(1-t)}, and accumulated by log-sum-exp so that terms which would underflow individually still contribute.

\bold{The upper tail.} Requesting \code{lower.tail = FALSE} negates the argument of \eqn{\Phi} rather than forming \eqn{1 - F}. At \eqn{q = 30} with \eqn{\sigma_u = \sigma_v = 1} the complement is exactly 1 in double precision while the direct computation returns a finite log-probability below \eqn{-10^{2}}.

\bold{Accuracy.} There is no closed form except at \eqn{\lambda = \sigma_u/\sigma_v} equal to 0 (the normal) or 1, where Amsler, Schmidt and Tsay show \eqn{P(Q) = \Phi(Q/(\sqrt{2}\sigma_u))^2} for the cost case. Against that standard the default quadrature is accurate to about \eqn{10^{-13}} relative at probabilities as small as \eqn{10^{-115}}. It is deliberately used everywhere rather than special-cased at \eqn{\lambda = 1}: the production form of that identity is a complement, which would reintroduce the saturation the function exists to avoid.
}
\references{
Amsler, C., Schmidt, P. and Tsay, W.-J. (2019). Evaluating the cdf of the distribution of the stochastic frontier composed error. \emph{Journal of Productivity Analysis}, 52(1), 29--35.
}
\seealso{\code{\link{copsfm}}, \code{\link{sfm}}}
\examples{
## Density integrates to one.
integrate(function(x) dcomposed(x, 1, 0.5), -Inf, Inf)$value

## Distribution function, and a tail that does not saturate.
pcomposed(c(-2, 0, 2), sigma_u = 1, sigma_v = 0.5)
pcomposed(30, 1, 1, lower.tail = FALSE, log.p = TRUE)

## Checked against the one exact standard there is (cost frontier).
su <- sv <- sqrt(0.5)
c(exact = pnorm(-6 / (sqrt(2) * su))^2,
  computed = exp(pcomposed(-6, su, sv, inefdec = FALSE, log.p = TRUE)))
}
\keyword{distribution}
