\name{efficiency}
\alias{efficiency}
\title{Technical Efficiency Scores}
\description{
Extracts observation-specific technical efficiency from a fitted frontier,
choosing among the three usual point predictors and between the two scales the
dependent variable may be on.
}
\usage{
efficiency(object, type = c("bc", "jlms", "mode"), logDepVar = TRUE,
           newdata = NULL)
}
\arguments{
\item{object}{An \code{"sfareg"} fit.}
\item{type}{Which predictor of \eqn{u_i} to use. \code{"bc"} (the default) is
  Battese--Coelli, \eqn{E[\exp(-u_i) \mid \varepsilon_i]}; \code{"jlms"} is
  Jondrow \emph{et al.}, \eqn{\exp(-E[u_i \mid \varepsilon_i])}; \code{"mode"}
  uses the mode of the posterior of \eqn{u_i}.}
\item{logDepVar}{\code{TRUE} (the default) if the dependent variable is on a
  log scale, so that efficiency is \eqn{\exp(-u_i)}; \code{FALSE} if it is on a
  level scale, giving \eqn{1 - u_i / f(x_i)}.}
\item{newdata}{Data to rebuild the fitted frontier from, needed only when
  \code{logDepVar = FALSE}. Supply it when the data used for the fit can no
  longer be recovered from the stored call.}
}
\details{
\bold{Which predictor.} The three differ in what they report about the same
posterior, not in how well they are estimated.
\eqn{E[\exp(-u)] \ge \exp(-E[u])} by Jensen's inequality, so \code{"bc"} is
never below \code{"jlms"}; the gap widens as the posterior spreads out.
\code{"mode"} is the only one that can equal exactly 1, and it does so for every
observation whose posterior mean is negative -- typically much of the efficient
tail. That is a property of the predictor rather than a defect: where the
posterior mean of \eqn{u} is negative the single most likely value really is the
boundary. It is available only for the models whose posterior is a truncated
normal (\code{"NHN"}, \code{"NHN_Z"}, \code{"NE"}, \code{"NTN"}).

\bold{Which scale.} The package elsewhere assumes throughout that the
dependent variable is logged, which is the usual case and makes \eqn{u} a
proportional shortfall. \code{logDepVar = FALSE} instead forms
\eqn{(f(x_i) - u_i)/f(x_i)}. That is a ratio of two estimated quantities and is
less well behaved than the log version: it is undefined where the fitted
frontier is zero and negative where \eqn{u} exceeds it, both of which happen
near zero output. A warning is issued if any score comes back negative or
non-finite, since silence there would be the real failure.
}
\value{
A numeric vector of efficiency scores, one per observation used in the fit.
}
\seealso{
\code{\link{efficiency_ci}} for interval estimates of the same quantity.
}
\examples{
\donttest{
set.seed(2)
n  <- 400
x1 <- rnorm(n); x2 <- rnorm(n)
y  <- 6 + 0.5 * x1 + 0.5 * x2 + rnorm(n, 0, 0.3) - abs(rnorm(n, 0, 0.8))
d  <- data.frame(y = y, x1 = x1, x2 = x2)

fit <- sfm(y ~ x1 + x2, model_name = "NHN", data = d)

summary(efficiency(fit))                      # Battese-Coelli, log scale
summary(efficiency(fit, type = "jlms"))       # never above the above
summary(efficiency(fit, type = "mode"))       # reaches exactly 1
summary(efficiency(fit, logDepVar = FALSE, newdata = d))
}
}
\keyword{models}
