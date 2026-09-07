\name{moment_range}
\alias{moment_range}
\alias{print.sfa_moment_range}
\title{What skewness and excess kurtosis can a pair of distributions produce?}
\description{
Reports, for each pair of noise and inefficiency distributions, the set of skewness and excess kurtosis values the composed error \eqn{\varepsilon = v \pm u} can attain at all, and says which pairs the residual moments rule out. This is the Type II failure check of Papadopoulos and Parmeter (2021): a pair whose attainable set does not contain the sample moments is refuted without any test statistic. It needs no frontier fit.
}
\usage{moment_range(object = NULL,
             noise = c("normal", "logistic", "laplace", "uniform", "t"),
             inefficiency = c("halfnormal", "exponential", "genexponential",
                              "truncatednormal", "gamma"),
             inefdec = NULL, df = 8)}
\arguments{
\item{object}{An \code{"sfareg"} fit, whose \bold{OLS} residuals are used, or a numeric vector of residuals. Omit it to get the attainable ranges alone, as a reference table.}
\item{noise}{Assumed distributions of \eqn{v}; all five by default. Only the excess kurtosis of \eqn{v} enters.}
\item{inefficiency}{Assumed distributions of \eqn{u}; all five by default.}
\item{inefdec}{Orientation. Taken from the fit when there is one, otherwise production. A cost frontier puts the skewness on the positive side.}
\item{df}{Degrees of freedom for \code{noise = "t"}, whose excess kurtosis is \eqn{6/(\nu-4)}. Must exceed 4.}
}
\value{
A data frame of class \code{"sfa_moment_range"}, one row per pair, with the attainable interval for the composed skewness (\code{g1_lo}, \code{g1_hi}) and excess kurtosis (\code{g2_lo}, \code{g2_hi}), logical \code{skew_ok} and \code{kurt_ok}, and a \code{verdict} of \code{"ok"}, \code{"type I (wrong skew)"} or \code{"type II (...)"}. The sample moments, sample size and orientation are attributes.
}
\details{
\bold{The two equations.} With \eqn{R = \mathrm{SNR}^2/(1+\mathrm{SNR}^2)} the share of composed-error variance contributed by \eqn{u}, Papadopoulos and Parmeter (2021, eqs. 5 and 10) give

\deqn{|\gamma_1(\varepsilon)| = \gamma_1(u)\,R^{3/2}, \qquad \gamma_2(\varepsilon) = \gamma_2(v)(1-R)^2 + \gamma_2(u)R^2.}

Because \eqn{R} lies in \eqn{(0,1)}, the composed skewness can never exceed the skewness of \eqn{u} itself, whatever the noise: the Normal-Half-Normal pair is confined to \eqn{\gamma_1(\varepsilon) \in (-0.995, 0)} and the Normal-Exponential pair to \eqn{(-2, 0)}. The kurtosis equation is quadratic in \eqn{R}, so its range is not simply the interval between the two components' values -- with leptokurtic noise it dips below \emph{both} before rising, reaching \eqn{\gamma_2(v)\gamma_2(u)/(\gamma_2(v)+\gamma_2(u))} at \eqn{R = \gamma_2(v)/(\gamma_2(v)+\gamma_2(u))}.

\bold{Which \eqn{u} can fail.} Half-normal, exponential and generalized exponential have fixed skewness and kurtosis, so their attainable sets are the narrowest. The truncated normal is indexed by \eqn{\mu/\sigma_u} and sweeps \eqn{\gamma_1(u) \in (0,2)} and \eqn{\gamma_2(u) \in (-0.243, 6)} -- it is the only one here whose composed error can be platykurtic under normal noise. The gamma has \eqn{\gamma_1(u) = 2/\sqrt{k}} and \eqn{\gamma_2(u) = 6/k}, unbounded as \eqn{k \to 0}, so a gamma inefficiency term can \emph{never} suffer a Type II failure; the price, as PP2021 note, is that \eqn{k < 1} makes the density log-convex.

\bold{This is a diagnostic, not a test.} PP2021 section 2.3 shows the range check over-rejects, and the package's own replication of their Tables 1-5 confirms it: with correctly specified Normal-Half-Normal data at \eqn{n = 200} and \eqn{\mathrm{SNR} = 3}, the sample excess kurtosis falls outside its theoretical range 26\% of the time, and the sample skewness 16\% of the time. Use \code{\link{spec_test}} for the formal statistic, which is built on the same two equations but studentised and bootstrapped.

\bold{Relation to the wrong-skew problem.} A Type I failure -- residual skewness of the wrong sign, Waldman (1982) -- is a property of the data and is reported for every row at once. Type II failures are pair-specific: they are what tells you that the residuals are too skewed, or too heavy-tailed, for the particular \eqn{u} you assumed. Wrong skew is handled by \code{\link{skewness_test}} and, as an estimator, by \code{\link{esfm}}.
}
\references{
Papadopoulos, A. and Parmeter, C. F. (2021). Type II failure and specification testing in the Stochastic Frontier Model. \emph{European Journal of Operational Research} 293(3), 990-1001.

Waldman, D. M. (1982). A stationary point for the stochastic frontier likelihood. \emph{Journal of Econometrics} 18(2), 275-279.
}
\seealso{\code{\link{spec_test}} for the formal test built on the same two equations, \code{\link{skewness_test}} for the Type I failure, \code{\link{esfm}} for an estimator that tolerates it.}
\examples{
## The reference table on its own: no data needed.
moment_range(noise = "normal")

set.seed(4)
d <- data_gen_cs(N = 400, rand = 4, sig_u = 1, sig_v = 0.5, cons = 1,
                 beta1 = 0.5, beta2 = 0.5, a = 1, mu = 0.5)

## Exponential inefficiency, read through a normal-half-normal lens: the
## residuals are more skewed than a half-normal can make them.
fit <- sfm(y_pcs_e ~ x1 + x2, data = d, model_name = "NHN")
moment_range(fit, noise = "normal")
}
\keyword{models}
