\name{cw_test}
\alias{cw_test}
\alias{print.cw_test}
\title{Chen and Wang (2012) centered-residuals moment test}
\usage{
cw_test(x, model_name = NULL, tau = 1, type = c("cosine", "sine"),
        data = NULL, n_nodes = 256L)
}
\arguments{
\item{x}{An \code{"sfareg"} fit, or a numeric vector of regression residuals
-- ordinary least squares residuals are enough.}
\item{model_name}{The distributional pair under the null. Required when
\code{x} is a numeric vector; taken from the fit otherwise. One of
\code{"NHN"}, \code{"NE"}, \code{"NR"}, \code{"NU"}, \code{"NGE"}.}
\item{tau}{Frequency at which the characteristic function is compared. Chen
and Wang suggest a value around 1, and a \emph{single} frequency: combining
several oversizes the test (see Details).}
\item{type}{\code{"cosine"} (recommended) or \code{"sine"}.}
\item{data}{The data the model was fitted to, when \code{x} is a fit.}
\item{n_nodes}{Quadrature nodes for the theoretical moments.}
}
\value{
An object of class \code{"cw_test"}: the statistic \code{D}, its degrees of
freedom and p-value, the moment estimates of \eqn{\sigma_v} and
\eqn{\sigma_u}, and the frequencies used.
}
\description{
Tests the distributional assumption on the composed error using the
characteristic function of the centred residuals, without ever evaluating the
composed-error density or distribution function.
}
\details{
\bold{Why this test rather than \code{\link{gof_test}}.} It never touches the
composed-error PDF or CDF. It compares the empirical characteristic function
of the centred residuals with its theoretical value, and that follows from the
characteristic functions of the two components \emph{separately}, which are in
closed form for many pairs whose convolution is not -- the normal-gamma model
being the standard example of a pair with no closed-form composed density.

\bold{Centring is what makes it work.} A frontier's intercept is not
identified separately from \eqn{E[u]}, so the least-squares intercept
estimates \eqn{\alpha - E[u]} rather than \eqn{\alpha}. Centring the residuals
cancels that inconsistent intercept outright, which is why the test is valid
off an ordinary regression and does not require maximum likelihood. The
variance then has to correct for two nuisances -- estimating the scale
parameters, and the centring itself. Omitting the second understates the
variance and oversizes the test.

\bold{Use one frequency, not several.} Measured here over 500 replications at
a nominal 5\%, half-normal data:

\tabular{lrrrr}{
test \tab \eqn{n} \tab \eqn{\tau=1} \tab \eqn{\tau=1.5} \tab \eqn{\tau=2} \cr
cosine \tab 500 \tab 0.070 \tab 0.054 \tab 0.050 \cr
cosine \tab 2000 \tab 0.046 \tab 0.044 \tab 0.050 \cr
sine \tab 500 \tab 0.092 \tab 0.064 \tab 0.042 \cr
sine \tab 2000 \tab 0.068 \tab 0.062 \tab 0.064
}

Against exponential inefficiency the cosine test's power at \eqn{\tau = 1} is
0.231 at \eqn{n = 500} and 0.952 at \eqn{n = 2000}, against 0.108 and 0.706
for the sine test. Both of Chen and Wang's conclusions reproduce: the cosine
test is the less sensitive to the choice of \eqn{\tau}, and the more powerful.
Passing several frequencies at once gives rejection rates of 0.09 to 0.15 at a
nominal 5\%, because the cosine moments at nearby frequencies are close to
collinear and the covariance matrix is then near-singular; the authors report
the same for their combined \eqn{\tau}.

The moment estimator solves the second- and third-moment equations, which for
\code{"NHN"} and \code{"NE"} is the closed form of the paper's equation (32)
and is the same inversion \code{sfm(estimator = "cols")} uses. A wrongly
signed third moment leaves it with no admissible solution, and that is an
error rather than a number: see \code{\link{skewness_test}}.
}
\references{
Chen, Y.-T. and Wang, H.-J. (2012). Centered-residuals-based moment estimator
and test for stochastic frontier models. \emph{Econometric Reviews}
\bold{31}(6), 625--653.
}
\seealso{\code{\link{gof_test}}, \code{\link{spec_test}}, \code{\link{skewness_test}}}
\examples{
set.seed(4)
n <- 500
X <- matrix(rnorm(n * 3), n, 3)
y <- X \%*\% c(1, 1, 1) + rnorm(n, 0, 1) - abs(rnorm(n, 0, 1))
e <- residuals(lm(y ~ X))

cw_test(e, "NHN")            # true pair
cw_test(e, "NE")             # wrong inefficiency distribution
}
