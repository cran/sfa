\name{esfm}
\alias{esfm}
\alias{symmetry_test}
\alias{print.esfm}
\title{Stochastic frontier tolerating the wrong skewness}
\usage{
esfm(formula, data, dist = c("halfnormal", "exponential"), inefdec = TRUE,
     start = NULL)

symmetry_test(object)
}
\arguments{
\item{formula, data}{Model and data, as for \code{\link{sfm}}.}
\item{dist}{Base density for the inefficiency term.}
\item{inefdec}{\code{TRUE} for a production frontier, \code{FALSE} for cost.}
\item{start}{Optional numeric start, \code{c(gamma, sigma_v, beta)}.}
\item{object}{An \code{"esfm"} fit.}
}
\value{
An object of class \code{"esfm"} (and \code{"sfareg"}), whose \code{out}
matrix reports \code{gamma}, \code{sigv} and the frontier coefficients, with
\code{Efficiency} holding \eqn{E[\exp(-U) \mid W]} for each firm.
}
\description{
Fits the extended stochastic frontier model of Hafner, Manner and Simar
(2018), in which one parameter \eqn{\gamma} carries both the scale of
inefficiency (its size) and the direction of skewness (its sign), so a
wrongly skewed sample gives a well-defined non-degenerate fit instead of
collapsing.
}
\details{
The classical frontier implies a negatively skewed composed error. When the
sample skewness comes out positive -- common at small \eqn{n} or large
\eqn{\sigma_v/\sigma_u}, and a small-sample accident rather than evidence
against the model -- the MLE collapses: \eqn{\sigma_u} is exactly zero, every
firm is fully efficient, and the fit is uninformative. That is the Type I
failure \code{\link{skewness_test}} reports.

This model keeps the frontier and widens the distribution of \eqn{u}. For
\eqn{\gamma > 0}, \eqn{u} has the classical density with scale
\eqn{|\gamma|}; for \eqn{\gamma < 0} that density is mirrored about zero,
shifted right by \eqn{B = a_0|\gamma|} and truncated to \eqn{[0, B]}, giving
\emph{negative} skewness for \eqn{u} and so positive skewness for the composed
error. The truncation constant \eqn{a_0} is fixed, not estimated: it is chosen
so \eqn{E[u] = k_1|\gamma|} whichever sign \eqn{\gamma} takes, which keeps the
model a one-parameter scale family and avoids the identification difficulties
of bounded-inefficiency models. \eqn{a_0 = 1.3892032925} for the half-normal
and \eqn{1.5936242600} for the exponential.

\bold{The classical model is nested} at \eqn{\gamma > 0}, and \eqn{\gamma = 0}
is an \emph{interior} point rather than a boundary. So the likelihood ratio
test of \eqn{H_0: \gamma = 0} -- no inefficiency, symmetric errors -- is an
ordinary \eqn{\chi^2(1)}, not the chi-bar-square mixture that
\code{\link{inefficiency_test}} needs for the classical model.

\bold{What it does on a wrongly skewed sample.} On a 50-observation draw with
positive residual skewness, \code{sfm(model_name = "NHN")} returns
\eqn{\sigma_u = 0.0003} and mean efficiency 0.99976 -- every firm on the
frontier. \code{esfm()} returns \eqn{\hat\gamma = -0.57} with slopes still
near their true values.

\bold{Size of the LR test}, 2000 replications against the paper's 100,000, on
its design (\eqn{\sigma_v = 1}, \eqn{\gamma = 0}):

\tabular{lrrrrr}{
\eqn{n} \tab 50 \tab 100 \tab 200 \tab 500 \tab 1000 \cr
this package, 5\% \tab 0.178 \tab 0.082 \tab 0.054 \tab 0.061 \tab 0.058 \cr
paper, 5\% \tab 0.152 \tab 0.076 \tab 0.061 \tab 0.054 \tab 0.053
}

The test is badly oversized below \eqn{n = 100} and usable from \eqn{n = 200},
which is the paper's own conclusion.

\bold{Efficiency bias, and a caveat.} The paper reports smaller bias in mean
efficiency than the classical model even when the population skewness has the
correct sign. That reproduces \emph{when inefficiency is not small}: at
\eqn{\gamma = 0.5}, \eqn{\sigma_v = 0.25} the bias in \eqn{E[\exp(-U)]} is
-0.0001 at \eqn{n = 100} against the classical model's +0.021. But at
\eqn{\gamma = 0.3}, where 22\% of samples are wrongly skewed, this model
overshoots the other way: bias -0.040 against the classical +0.036, and it is
the larger of the two at \eqn{n = 50} (-0.080 against +0.044). A small true
\eqn{\gamma} lets the extended model fit sizeable inefficiency of either sign
to noise. \bold{It is the better choice when the wrong skewness is the
problem, not a free improvement everywhere.}
}
\references{
Hafner, C. M., Manner, H. and Simar, L. (2018). The \dQuote{wrong skewness}
problem in stochastic frontier models: A new approach. \emph{Econometric
Reviews} \bold{37}(4), 380--400.

Almanidis, P., Qian, J. and Sickles, R. C. (2014). Stochastic frontier models
with bounded inefficiency. In \emph{Festschrift in Honor of Peter Schmidt}.
}
\seealso{\code{\link{skewness_test}}, \code{\link{inefficiency_test}}, \code{\link{sfm}}}
\examples{
set.seed(3)
n <- 200
lx1 <- rnorm(n, 1.5, 0.3); lx2 <- rnorm(n, 1.8, 0.3)
y <- 0.9 + 0.6 * lx1 + 0.5 * lx2 + rnorm(n, 0, 0.25) - abs(rnorm(n, 0, 0.5))
d <- data.frame(y, lx1, lx2)

f <- esfm(y ~ lx1 + lx2, data = d)
f
symmetry_test(f)
}
