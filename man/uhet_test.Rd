\name{uhet_test}
\alias{uhet_test}
\title{Test whether technical inefficiency depends on firm characteristics}
\usage{
uhet_test(object, z, data = NULL)
}
\arguments{
\item{object}{An \code{"sfareg"} fit from \code{\link{sfm}} with
\code{model_name = "NHN"} or \code{"NE"}, fitted \emph{without} \code{uhet}
and with \code{keep_objective = TRUE}. This is the first step, estimated
under the null.}
\item{z}{A one-sided formula (e.g. \code{~ z1 + z2}) or a numeric matrix of
the characteristics whose effect on inefficiency is being tested.}
\item{data}{The data the model was fitted to.}
}
\value{
A data frame with the corrected and naive statistics, their common degrees of
freedom and p-values. The second-step coefficients and both sets of standard
errors are attached as attributes.
}
\description{
The two-step test of Kim and Schmidt (2008): regress the JLMS predictor on
firm characteristics and test the slopes, with the variance correction that
makes the test valid.
}
\details{
The procedure applied research actually runs is: fit the frontier, take
\eqn{\hat u_i = E[u_i \mid \varepsilon_i]}, regress it on \eqn{z_i}, and test
the slopes. \bold{As usually performed that test is invalid.} \eqn{\hat u} is
a \emph{generated} dependent variable, so the first-step estimation error in
\eqn{\psi} enters the second-step variance; the naive test omits it. Kim and
Schmidt show the omitted term is \eqn{G r_i} with
\eqn{G = E[z_i \nabla_\psi b_i]}, so the naive test is valid only when
\eqn{G = 0} -- which holds when \eqn{z} is independent of \eqn{x}, and
generally fails when they are correlated, the case applied work is usually in.
There is no distribution-free correction, which is why only the two cases the
paper works through are implemented.

\bold{What the correction buys.} With \eqn{z} independent of \eqn{x} the two
tests agree, as the theory says they must. With \eqn{z} strongly correlated
with \eqn{x} (\eqn{\rho = 0.9}) the naive standard error is \bold{52\% too
large} at \eqn{n = 800} and the naive test's size collapses to \bold{0.000} at
a nominal 5\%, while the corrected standard error matches the true sampling
spread to within half a percent (ratio 0.995) and the test retains size 0.032.
A test that never rejects also has no power against local alternatives.

\bold{An honest caveat on size.} The corrected test is \emph{conservative} in
these simulations: on Kim and Schmidt's own base case (scaled exponential,
\eqn{\sigma_v^2 = \lambda = 1}, \eqn{\rho = 0.5}, 1000 replications) it
rejects 0.024, 0.036 and 0.028 of the time at \eqn{N = 200, 500, 1000} against
a nominal 0.05, where the paper reports 0.0495 for its corrected two-step test
at \eqn{N = 200}. The naive branch reproduces the paper's figure for the
\emph{uncorrected} test closely (0.023--0.029 against their 0.0237). The
implementation was checked line by line against equations (5)--(13), the
first-step scores agree with an analytic half-normal score to
\eqn{10^{-11}}, \eqn{\hat\gamma} is unbiased and Gaussian (kurtosis 2.92), and
an oracle test using the true sampling standard deviation has correct size
(0.045) -- so the discrepancy is in the variance estimate, not in
\eqn{\hat\gamma} or the reference distribution: \eqn{cor(|\hat\gamma|, se)}
is 0.375, which compresses the ratio to a standard deviation of 0.93.

Two candidate explanations were tested and \bold{both refuted}. Replacing the
robust second-step variance with the homoskedastic one -- the difference
between the paper's \dQuote{OLS} and \dQuote{OLS-H} variants -- changes
nothing (0.023/0.027/0.036 against 0.022/0.036/0.031). Nor does using the
Hessian information in place of the outer-product one when forming \eqn{r_i}
(0.020/0.026/0.026). The discrepancy is therefore \bold{unexplained}, not
merely unverified. \bold{Treat a rejection as meaningful and a non-rejection
as weak evidence.}

Fits sitting on the wrong-skew boundary are refused rather than answered:
\eqn{\hat u} is then identical for every firm, so there is no estimated
inefficiency for \eqn{z} to explain, and the information matrix is singular
there (Waldman 1982).
}
\references{
Kim, M. and Schmidt, P. (2008). Valid tests of whether technical inefficiency
depends on firm characteristics. \emph{Journal of Econometrics} \bold{144},
409--427.

Wooldridge, J. M. (2002). \emph{Econometric Analysis of Cross Section and
Panel Data}, pp. 139--141.
}
\seealso{\code{\link{sfm}} (its \code{uhet} argument fits the alternative
directly), \code{\link{marginal_effects}}, \code{\link{skewness_test}}}
\examples{
set.seed(2)
n <- 400
x1 <- rnorm(n); x2 <- rnorm(n)
z <- 0.9 * x1 + sqrt(1 - 0.81) * rnorm(n)   # correlated with x: G != 0
y <- 0.5 + 0.5 * x1 + 0.5 * x2 + rnorm(n, 0, 0.5) - abs(rnorm(n, 0, 1))
d <- data.frame(y, x1, x2, z)

f <- sfm(y ~ x1 + x2, model_name = "NHN", data = d, keep_objective = TRUE)
uhet_test(f, ~z, data = d)
}
