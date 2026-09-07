\name{inefficiency_test}
\alias{inefficiency_test}
\title{Test whether there is any technical inefficiency}
\usage{
inefficiency_test(object, test = c("lr_1sided", "lr", "wald", "m3t"),
                  level = 0.05)
}
\arguments{
\item{object}{An \code{"sfareg"} fit from \code{\link{sfm}}, with
\code{model_name = "NHN"} or \code{"NTN"}.}
\item{test}{Which statistics to compute; any subset, all by default.}
\item{level}{Significance level used for the \code{reject} column.}
}
\value{
A data frame with one row per test: \code{test}, \code{statistic},
\code{null} (the reference distribution), \code{p.value} and \code{reject}.
The log-likelihoods under both hypotheses are attached as attributes.
}
\description{
Tests \eqn{H_0}: no technical inefficiency, i.e. \eqn{\gamma = 0} in the
Battese-Corra parameterization \eqn{\gamma = \sigma_u^2/(\sigma_u^2 +
\sigma_v^2)}, following Section 3 of Coelli (1995).
}
\details{
\bold{The two tests most often reported are the two with the wrong size.}
\eqn{H_0} puts \eqn{\gamma} on the boundary of the parameter space, so the
usual asymptotics fail. The likelihood ratio statistic is not
\eqn{\chi^2(1)} but the Gourieroux, Holly and Monfort (1982) mixture
\eqn{\tfrac12\chi^2(0) + \tfrac12\chi^2(1)}, whose 5\% critical value is
2.71 rather than 3.84; using the \eqn{\chi^2(1)} value makes the naive LR
test conservative. The Wald ratio \eqn{\hat\gamma/s_{\hat\gamma}} is worse:
it is not asymptotically standard normal at the boundary at all.

Measured in this package on 400 replications per cell, half-normal data with
\eqn{\sigma_u = 0} and a nominal 5\% level:

\tabular{lrrrr}{
\eqn{n} \tab LR one-sided \tab LR naive \tab Wald \tab M3T \cr
100 \tab 0.060 \tab 0.030 \tab 0.230 \tab 0.055 \cr
200 \tab 0.060 \tab 0.033 \tab 0.203 \tab 0.050 \cr
400 \tab 0.055 \tab 0.033 \tab 0.217 \tab 0.045 \cr
800 \tab 0.060 \tab 0.035 \tab 0.182 \tab 0.062
}

The one-sided LR test and the third-moment test hold their size; the naive
LR test rejects about half as often as it should; the Wald test rejects
three to four times too often and does not improve with \eqn{n}. Under the
alternative (\eqn{\sigma_u = \sigma_v = 1}) the one-sided LR test has the
better power of the two valid tests at every sample size, which is Coelli's
recommendation. \bold{The Wald column's apparently high power is not
comparable, because its size is wrong.}

For \code{model_name = "NTN"} the null restricts both \eqn{\gamma} and
\eqn{\mu}, and the mixture becomes
\eqn{\tfrac14\chi^2(0) + \tfrac12\chi^2(1) + \tfrac14\chi^2(2)}.

The \code{"m3t"} statistic is the same one \code{\link{skewness_test}}
reports as \code{test = "coelli"}; it is repeated here so the four tests can
be compared in a single table.
}
\references{
Coelli, T. (1995). Estimators and hypothesis tests for a stochastic frontier
function: A Monte Carlo analysis. \emph{Journal of Productivity Analysis}
\bold{6}, 247--268.

Gourieroux, C., Holly, A. and Monfort, A. (1982). Likelihood ratio test,
Wald test, and Kuhn-Tucker test in linear models with inequality constraints
on the regression parameters. \emph{Econometrica} \bold{50}, 63--80.
}
\seealso{\code{\link{skewness_test}}, \code{\link{spec_test}}, \code{\link{sfm}}}
\examples{
d <- data_gen_cs(N = 300, rand = 7, sig_u = 1, sig_v = 0.5, cons = 0.5,
                 beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.1)
f <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = as.data.frame(d))
inefficiency_test(f)
}
