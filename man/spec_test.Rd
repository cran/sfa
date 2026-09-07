\name{spec_test}
\alias{spec_test}
\alias{spec_test_all}
\alias{print.sfa_spec_test}
\alias{print.sfa_spec_grid}
\title{Is this pair of distributions defensible?}
\description{
Tests an assumed \bold{pair} of distributions for the composed error \eqn{\varepsilon = v \pm u} using only the skewness and excess kurtosis of the OLS residuals. It needs no frontier fit, so it can be run before estimating anything and cannot be contaminated by an optimiser that failed. \code{spec_test_all} runs all twelve supported pairs and ranks them.
}
\usage{spec_test(object, noise = c("normal", "logistic", "laplace", "uniform"),
          inefficiency = c("halfnormal", "exponential", "genexponential"),
          null = c("bootstrap", "asymptotic"),
          B = 999, level = 0.05, seed = NULL)

spec_test_all(object, level = 0.05)}
\arguments{
\item{object}{An \code{"sfareg"} fit, whose \bold{OLS} residuals are used, or a numeric vector of residuals.}
\item{noise}{Assumed distribution of \eqn{v}. Only its excess kurtosis enters, so any symmetric distribution with the same value behaves identically -- testing for Logistic noise is indistinguishable from testing for \eqn{t(9)}.}
\item{inefficiency}{Assumed distribution of \eqn{u}.}
\item{null}{How to get the p-value. \code{"bootstrap"} (the default) simulates from the assumed pair; \code{"asymptotic"} uses the delta-method normal approximation, which is badly over-sized for heavy-tailed \eqn{u} -- see \sQuote{Details}.}
\item{B}{Bootstrap replications. Cheap here: nothing is fitted, so each replicate is \eqn{n} draws and eight moments.}
\item{level}{Size used for the reported verdict.}
\item{seed}{Optional seed for the bootstrap. The RNG state is saved and restored.}
}
\value{
\code{spec_test} returns an object of class \code{c("sfa_spec_test", "htest")}, carrying \code{statistic} (the studentised \eqn{\Gamma}), \code{p.value}, \code{Gamma}, \code{se}, \code{R}, \code{skewness}, \code{excess_kurtosis}, \code{type2_failure} and \code{reject}. \code{spec_test_all} returns a data frame of all twelve pairs sorted by p-value.
}
\details{
\bold{What it tests.} Write \eqn{R = \mathrm{SNR}^2/(1+\mathrm{SNR}^2)} for the share of composed-error variance contributed by \eqn{u}. Papadopoulos and Parmeter (2021, eqs. 5 and 10) give

\deqn{|\gamma_1(\varepsilon)| = \gamma_1(u)\,R^{3/2}, \qquad \gamma_2(\varepsilon) = \gamma_2(v)(1-R)^2 + \gamma_2(u)R^2.}

With both distributions assumed, \eqn{\gamma_1(u)}, \eqn{\gamma_2(u)} and \eqn{\gamma_2(v)} are known constants, so the \emph{skewness alone} identifies \eqn{R} and the kurtosis equation has nothing left to fit. The statistic is the gap,

\deqn{\Gamma = \gamma_2(\varepsilon) - \gamma_2(v)(1-R)^2 - \gamma_2(u)R^2,}

which is zero when the assumed pair is right. Its variance comes from the delta method applied to the central moments \eqn{(m_2, m_3, m_4)}, which needs moments up to the eighth.

\bold{Why the default null is a bootstrap.} The delta-method p-value is not reliable at usual sample sizes, and the failure is systematic rather than marginal. Rejection rates under \emph{correct} specification, 300 replications per cell:

\tabular{lllll}{
 \eqn{u} \tab \eqn{n} \tab SNR \tab asymptotic \tab bootstrap \cr
 half-normal \tab 2000 \tab 1 \tab 0.087 \tab 0.050 \cr
 half-normal \tab 2000 \tab 2 \tab 0.103 \tab 0.057 \cr
 exponential \tab 2000 \tab 1 \tab 0.207 \tab 0.077 \cr
 exponential \tab 500 \tab 2 \tab \bold{0.357} \tab 0.083 \cr
 exponential \tab 2000 \tab 2 \tab 0.307 \tab 0.060
}

Across all 24 cells the asymptotic null rejects \bold{18.7\%} of the time at a nominal 5\%; the bootstrap rejects \bold{5.7\%}. The cause is the eighth moment: for a heavy-tailed \eqn{u} it is estimated far too imprecisely at these sample sizes, and the variance comes out too small. The bootstrap studentises against the same statistic computed under the assumed pair, so it corrects the variance rather than merely relocating the statistic. Power is unaffected -- against wrong pairs at \eqn{n = 2000} it runs 0.93 to 1.00.

\bold{Type II failure.} If \eqn{R > 1} the residuals are more skewed than the assumed \eqn{u} can be -- for instance \eqn{|\gamma_1| > 0.9953} under a half-normal. The pair is then refuted by the skewness alone, and \code{spec_test} says so instead of reporting a p-value that would suggest the question was close.

\bold{Reading \code{spec_test_all}.} These are not sequential tests: no one of them depends on another's outcome, so the procedure implies no multiplicity correction. But picking the largest of twelve p-values \emph{is} a selection, and the ranking should be read as comparative evidence rather than as a significance level.

\bold{Relation to \code{\link{skewness_test}}.} That asks whether the residuals are skewed the way inefficiency makes them -- a test of whether there \emph{is} inefficiency. This asks whether the assumed distributions are the right ones. They answer different questions and are worth running together.

\bold{Provenance.} The statistic here is \emph{reconstructed} from PP2021 equations (5) and (10); the Supplementary Appendix in which PP2023 derive it was not available. The reconstruction is verified two ways in the tests: \eqn{\Gamma} vanishes under correct specification across sixteen designs, and the tabulated constants \eqn{c_{u,1}} and \eqn{c_{u,2}} are reproduced exactly as \eqn{\gamma_1(u)^{-2/3}} and \eqn{1/\gamma_2(u)}. PP2021's own equation (13) is deliberately not used, since PP2023 record that its \eqn{H_k} terms are missing squares.
}
\references{
Papadopoulos, A. and Parmeter, C.F. (2021). Type II failure and specification testing in the Stochastic Frontier Model. \emph{European Journal of Operational Research}, 293(3), 990--1001.

Papadopoulos, A. and Parmeter, C.F. (2023). A specification test for the composed error term in the stochastic frontier model. \emph{Economics Letters}, 233, 111390.
}
\seealso{\code{\link{moment_range}} for the attainable ranges the same two equations imply, \code{\link{skewness_test}}, \code{\link{TIC}}, \code{\link{sfm}}}
\examples{
\donttest{
set.seed(2)
n <- 2000
e <- rnorm(n, 0, 0.6) - abs(rnorm(n, 0, 1))   # normal / half-normal

## The true pair is not rejected.
spec_test(e, noise = "normal", inefficiency = "halfnormal",
          B = 199, seed = 1)

## A wrong one is.
spec_test(e, noise = "laplace", inefficiency = "exponential",
          B = 199, seed = 1)

## All twelve, ranked.
spec_test_all(e)
}
}
\keyword{htest}
