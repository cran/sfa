\name{influence_sfa}
\alias{influence_sfa}
\alias{print.sfa_influence}
\title{Which observations are driving this frontier fit?}
\description{
Empirical influence function and case-influence measures for a fitted
stochastic frontier. Answers \emph{which observations move the estimates}, and
\emph{whether the specification lets any single observation move them without
bound}.
}
\usage{influence_sfa(object, scale = TRUE)}
\arguments{
\item{object}{An \code{"sfareg"} fit. It must retain its likelihood, so fit
with \code{keep_objective = TRUE}; the score matrix cannot otherwise be built.}
\item{scale}{Multiply the influence function by \eqn{n}, so that
\code{influence[i, ]} approximates the effect on the estimate of deleting
observation \eqn{i} and is therefore on the same scale as the coefficients.
Default \code{TRUE}.}
}
\value{
An object of class \code{"sfa_influence"}:
\item{influence}{\eqn{n \times p} matrix, the empirical influence function.}
\item{d}{Length-\eqn{n} case influence \eqn{n\,s_i'Vs_i}, the same quadratic
form as Cook's distance and invariant to reparameterisation.}
\item{sensitivity_std}{\eqn{\sup_i \sqrt{d_i}}, the self-standardised
sensitivity -- \bold{the one to compare across models}.}
\item{sensitivity}{\eqn{\sup_i \|IF_i\|}, the raw sup-norm.}
\item{max_abs_by_parameter}{Largest \eqn{|IF|} in each coordinate.}
}
\details{
The standard outlier rules do not transfer to this model. Least-median-of-squares
and its relatives ignore the \bold{asymmetry} of the composed error, so a
genuinely large \eqn{u_i} -- an inefficient firm, the thing being measured --
reads as an outlier and is discarded, which defeats the analysis.

The influence function is the right object instead. Maximum likelihood is an
M-estimator, so its influence function is a linear transformation of the score,
\eqn{IF_i \propto I(\theta)^{-1}\psi(y_i,\theta)}, and the estimator is
\emph{B-robust} exactly when that is bounded. Both pieces are already available:
\code{\link[=sfareg-methods]{estfun}} supplies the per-observation scores and
\code{\link[=sfareg-methods]{vcov}} the inverse information.

\bold{Which sensitivity to read.} The raw sup-norm depends on how the parameters
happen to be scaled, so it is \emph{not} comparable between models with different
parameter vectors. On a clean sample of 200 it reads 72.6 for \code{"NHN"} and
1322.5 for \code{"tHN"}, which appears to say the robust specification is far
worse; all it reflects is \code{tHN}'s extra \eqn{\nu} on its own scale.
\code{sensitivity_std} measures the same influence function in the information
metric and is invariant to reparameterisation. On that scale, contaminating a
single response takes \code{"NHN"} from 13.5 to 32.3 and \code{"tHN"} from 13.1
to 17.2 -- Stead, Wheat and Greene's finding, and invisible in the raw number.

Neither number has an absolute threshold. The comparison to make is between
specifications on the same data, not either number against a cut-off.

\bold{What to do about a large value.} Nothing, by deletion. A large sensitivity
is a property of the specification: under an unbounded influence function the
next most extreme point simply replaces the one removed. The remedies are a
specification with bounded influence -- Stead, Wheat and Greene's Student's
\eqn{t} noise term, \code{model_name = "tHN"} -- or a robust divergence
criterion, \code{\link{sfm}(robust = )}, whose tuning is chosen by
\code{\link{hscore_select}} or \code{\link{calibrate_c}}.

\bold{What this cannot see.} Like \code{\link{density_weights}}, it is computed
at the fitted surface. An observation mis-recorded in a \emph{regressor} can bend
that surface toward itself and so appear unremarkable here. Pair this with a
leverage check before concluding that an observation is not driving a result.
}
\references{
Stead, A.D., Wheat, P. and Greene, W.H. (2023). Robust maximum likelihood
estimation of stochastic frontier models. \emph{European Journal of Operational
Research}, 309(1), 188--201.

Bernstein, D.H., Parmeter, C.F. and Wright, I.A. (2026). On Robust Estimation of
the Stochastic Frontier Model. Working paper.
}
\seealso{\code{\link{density_weights}}, \code{\link{hscore_select}},
\code{\link{sfm}}, \code{\link{sfareg-methods}}}
\examples{
\donttest{
set.seed(3)
n <- 120
x <- runif(n, 1, 10)
d <- data.frame(y = 1 + 0.5 * log(x) + rnorm(n, 0, 0.3) - abs(rnorm(n, 0, 0.6)),
                x = log(x))

fit <- sfm(y ~ x, data = d, model_name = "NHN", keep_objective = TRUE)
influence_sfa(fit)

## Contaminate one response and compare specifications on the
## self-standardised scale.
d2 <- d; d2$y[1] <- d2$y[1] - 5
a <- sfm(y ~ x, data = d2, model_name = "NHN", keep_objective = TRUE)
b <- sfm(y ~ x, data = d2, model_name = "tHN", keep_objective = TRUE)
c(NHN = influence_sfa(a)$sensitivity_std,
  tHN = influence_sfa(b)$sensitivity_std)
}
}
\keyword{regression}
