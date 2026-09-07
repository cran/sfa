\name{pcomposed_model}
\alias{pcomposed_model}
\alias{composed_cdf}
\title{Composed-error distribution function for any cross-sectional model}
\usage{
pcomposed_model(q, model_name, par, inefdec = TRUE, lower.tail = TRUE,
                log.p = FALSE, n_nodes = 128)

composed_cdf(object, q = NULL, data = NULL, ...)
}
\arguments{
\item{q}{Numeric vector of quantiles. For \code{composed_cdf}, defaults to the
fit's own composed-error residuals.}
\item{model_name}{One of \code{"NHN"}, \code{"NTN"}, \code{"NE"}, \code{"NR"},
\code{"NU"}, \code{"NGE"}, \code{"NLN"}, \code{"NW"}, \code{"NG"},
\code{"NNAK"}, \code{"TSL"}, \code{"THT"}, \code{"tHN"}.}
\item{par}{A \emph{named} parameter vector, as in \code{fit$out[, "par"]}.}
\item{inefdec}{\code{TRUE} for a production frontier, \code{FALSE} for cost.}
\item{lower.tail, log.p}{As elsewhere in R.}
\item{n_nodes}{Quadrature nodes over the inefficiency term.}
\item{object}{An \code{"sfareg"} fit from \code{\link{sfm}}.}
\item{data}{The data the model was fitted to. \code{sfm()} does not retain it,
so it is needed whenever \code{q} is not supplied.}
\item{...}{Passed to \code{pcomposed_model}.}
}
\value{A numeric vector of probabilities (or their logs).}
\description{
\eqn{P(\varepsilon \le q)} for the composed error of any cross-sectional
model, where \code{\link{pcomposed}} covers only the half-normal case.
}
\details{
The probability is computed as an expectation over the inefficiency term,
\deqn{F(q) = E_u\left[F_v\left((q + su)/\sigma_v\right)\right],}
with \eqn{s = +1} for production and \eqn{-1} for cost. Conditioning on
\eqn{u} leaves the noise CDF in closed form, so the noise is never integrated
and the tails inherit \code{pnorm()}'s (or \code{pt()}'s) own accuracy.

Three parameterizations do not read the way their names suggest, and are
taken from each likelihood in \code{sfm()} rather than from the label:
\code{"THT"} lists \code{sigu} \emph{before} \code{sigv}, the only model that
does; \code{"NG"}'s \code{sigu} is the gamma \emph{scale} and its \code{mu}
the \emph{shape}; \code{"NNAK"}'s \code{sigu} is the Nakagami spread, so
\eqn{\Omega = \code{sigu}^2}; and \code{"NLN"}'s \code{mu} is a
\emph{meanlog}.

\bold{\code{"THT"} is not an independent convolution.} Tancredi's (2002)
composed error is skew-\eqn{t}, a scale mixture in which \eqn{v} and \eqn{u}
are divided by the \emph{same} \eqn{\sqrt{V/a}}. Treating it as
\eqn{t} noise plus an independent half-normal -- which is what \code{"tHN"}
actually is -- gets the log density wrong by up to 1.74, so the two models
take different paths here despite looking alike.

Verified two ways: against \code{\link{pcomposed}} for the half-normal case,
agreeing to \eqn{3.5\times10^{-15}} relative including a lower tail at
\eqn{\log F = -209}; and for all thirteen models against a four-million-draw
simulation, with a maximum absolute error of \eqn{4\times10^{-4}} against a
Monte Carlo standard error of \eqn{7.5\times10^{-4}}. Numerically
differentiating this CDF reproduces each model's own stored log-density to
about \eqn{10^{-8}}, except for the simulated-ML models \code{"NLN"} and
\code{"NW"}, where the difference is the simulation error of their
likelihood rather than a disagreement about the model.
}
\seealso{\code{\link{pcomposed}}, \code{\link{dcomposed}}, \code{\link{sfm}}}
\examples{
p <- c(sigv = 0.5, sigu = 1.1)
pcomposed_model(c(-2, -1, 0, 1), "NE", p)

d <- as.data.frame(data_gen_cs(N = 200, rand = 5, sig_u = 1, sig_v = 0.5,
                               cons = 0.5, beta1 = 0.5, beta2 = 0.5,
                               a = 5, mu = 0.1))
f <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = d)
head(composed_cdf(f, data = d))
}
