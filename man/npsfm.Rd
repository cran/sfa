\name{npsfm}
\alias{npsfm}
\alias{print.npsfareg}
\alias{summary.npsfareg}
\title{Nonparametric Stochastic Frontier Models}
\description{Fits a stochastic frontier whose frontier function is estimated by kernel regression rather than assumed linear. Two estimators are available: the two-step estimator of Fan, Li and Weersink (1996), and the local method-of-moments estimator of Simar, Van Keilegom and Zelenyuk (2017), which additionally lets both variance components vary with the covariates.}
\usage{npsfm(formula, data, method = c("FLW", "SVKZ", "PSZ", "KPST", "MY", "SZ"),
      dist = c("hn", "exp", "gamma", "unif"),
      regtype = c("lc", "ll"), bw.sel = c("cv.ls", "cv.aic"),
      bw = NULL, cost = FALSE, eff = TRUE,
      maxit = 5000, tol = 1e-3, iter = 25,
      rts = c("vrs", "crs", "drs", "irs"),
      prior.fit = NULL, log.form = TRUE, verbose = FALSE)}
\arguments{
\item{formula}{A single-part symbolic description of the frontier, \code{y ~ x1 + x2}. Unlike the package's parametric entry points \code{npsfm()} takes no \code{| z} segment: heteroskedasticity is handled nonparametrically through the covariates themselves under \code{method = "SVKZ"}, and a pipe is an error rather than something silently ignored.}
\item{data}{A data frame containing the variables named in \code{formula}.}
\item{method}{Which estimator to use.
\describe{
  \item{\code{"FLW"}}{Fan, Li and Weersink (1996). Estimates \eqn{E[y|x]} by kernel regression, then recovers the scale parameters from the residuals -- by maximizing their concentrated likelihood in \eqn{\lambda} when \code{dist = "hn"}, and by inverting central moments otherwise. \eqn{\sigma_u} and \eqn{\sigma_v} are constants; only the frontier is smooth.}
  \item{\code{"SVKZ"}}{Simar, Van Keilegom and Zelenyuk (2017). Runs three local-linear regressions -- \eqn{y} on \eqn{x}, then the squared and cubed residuals on \eqn{x} -- and inverts the local moments, so \eqn{\sigma_u(x)} and \eqn{\sigma_v(x)} both vary with the covariates. No optimizer runs. Normal-half normal only.}
  \item{\code{"PSZ"} (or \code{"KPST"})}{Park, Simar and Zelenyuk. Local \emph{maximum likelihood}: at every evaluation point the frontier and both log variance components are given local-linear expansions and the kernel-weighted normal-half normal likelihood is maximized in those \eqn{3(k+1)} parameters. One optimization per observation, so it is much slower than the two above. \eqn{\sigma_u(x)} and \eqn{\sigma_v(x)} vary with the covariates.}
  \item{\code{"MY"}}{Martins-Filho and Yao. Iterative local likelihood: alternates a local-linear fit of the frontier at every evaluation point, holding \eqn{(\lambda,\sigma)} fixed, with a global update of \eqn{(\lambda,\sigma)} from the resulting composed residuals, until the scale parameters stop moving. Cost is (observations \eqn{\times} iterations) optimizations. \eqn{\sigma_u} and \eqn{\sigma_v} are constants; only the frontier is local.}
  \item{\code{"SZ"}}{Simar and Zelenyuk (2011). Not an estimator in its own right: it takes an already-estimated smooth frontier and passes its fitted values through an output-oriented DEA, imposing the monotonicity and (under \code{rts = "vrs"}/\code{"crs"}) convexity that a kernel fit does not guarantee. Supply the prior fit through \code{prior.fit}, or leave it \code{NULL} to fit \code{"SVKZ"} first. Requires the \pkg{Benchmarking} package and is implemented for production frontiers only.}
} Matching ignores case, so \code{"nhn"}, \code{"NHN"} and \code{"Nhn"} are the same choice.}
\item{dist}{Distribution of the one-sided inefficiency term, for \code{method = "FLW"} only: \code{"hn"} (half normal, the default), \code{"exp"} (exponential), \code{"gamma"}, or \code{"unif"} (uniform on \eqn{[0,b]}). Every other method is derived for the half-normal case and errors for anything else.}
\item{regtype}{Kernel regression type passed to \code{np::npregbw}: \code{"lc"} (local constant, the default) or \code{"ll"} (local linear). Applies to \code{"FLW"}; \code{"SVKZ"} is local linear throughout by construction.}
\item{bw.sel}{Bandwidth selection method: \code{"cv.ls"} (least-squares cross-validation, the default) or \code{"cv.aic"}. Ignored when \code{bw} is supplied.}
\item{bw}{Optional numeric vector of bandwidths, one per covariate. When supplied, cross-validation is skipped and these are used directly -- useful for undersmoothing, for sensitivity analysis, or simply to avoid paying for bandwidth selection repeatedly in a simulation.}
\item{cost}{Logical. \code{FALSE} (the default) fits a production frontier, in which inefficiency is subtracted; \code{TRUE} fits a cost frontier.}
\item{eff}{Logical. Compute observation-level inefficiency predictions (\code{u_hat}, \code{exp_u_hat}). Available for \code{dist = "hn"} and \code{dist = "exp"}; the gamma and uniform cases return moment estimates only. Defaults to \code{TRUE}.}
\item{maxit}{Maximum function evaluations for each local optimization under \code{"PSZ"} (passed to \code{minqa::bobyqa}) and maximum iterations for each local fit under \code{"MY"} (passed to \code{optim}). Ignored by the other methods. Defaults to \code{5000}. Raise it if \code{convergence} reports many non-zero codes.}
\item{tol}{Convergence tolerance for \code{"MY"}: the iteration stops when the squared change in \eqn{(\lambda,\sigma)} falls below this. Defaults to \code{1e-3}.}
\item{iter}{Maximum number of outer iterations for \code{"MY"}. Defaults to \code{25}.}
\item{rts}{Returns-to-scale assumption for the DEA step under \code{"SZ"}: \code{"vrs"} (the default), \code{"crs"}, \code{"drs"} or \code{"irs"}. Passed to \code{Benchmarking::dea}.}
\item{prior.fit}{For \code{"SZ"}: a numeric vector of already-estimated frontier values, one per observation, to be monotonized. If \code{NULL} (the default) an \code{"SVKZ"} fit is computed first and used.}
\item{log.form}{For \code{"SZ"}: whether the data are in logs, as is conventional in this literature. When \code{TRUE} (the default) the frontier and covariates are exponentiated before the DEA step and the result is returned to the log scale.}
\item{verbose}{Logical. Report progress through the per-observation loops of \code{"PSZ"} and the outer iterations of \code{"MY"}. Defaults to \code{FALSE}.}
}
\details{
Both estimators relax the parametric frontier of \code{\link{sfm}} while keeping the composed-error structure \eqn{y = m(x) + v - u}. Because a kernel regression of \eqn{y} on \eqn{x} estimates \eqn{E[y|x] = m(x) - E[u]} rather than \eqn{m(x)}, both proceed by fitting that conditional mean and then shifting it back up by an estimate of \eqn{E[u]}.

\bold{Where they differ.} \code{"FLW"} treats \eqn{\sigma_u} and \eqn{\sigma_v} as constants, so its correction is a single number and the fitted gradients are those of the conditional mean. \code{"SVKZ"} estimates \eqn{\sigma_u(x)} from the local third moment, so the correction varies across observations and the frontier gradient picks up an extra term through the chain rule.

\bold{Least squares versus local likelihood.} \code{"FLW"} and \code{"SVKZ"} both begin from a least-squares kernel regression, which estimates \eqn{E[y|x] = m(x) - E[u]} and therefore needs the mean shift described above. \code{"PSZ"} and \code{"MY"} instead maximize the \emph{composed-error} likelihood locally, in which the local intercept is \eqn{m(x)} directly and no shift is applied. They pay for that with one numerical optimization per observation -- for \code{"MY"}, per observation \emph{per iteration} -- so expect them to be one to two orders of magnitude slower than \code{"FLW"}. Both are seeded from an \code{"FLW"} fit.

\bold{Which to use.} On a simulated nonlinear frontier with \eqn{\sigma_u = 0.6}, \eqn{\sigma_v = 0.25} (3 replications), mean absolute frontier error at \eqn{n = 300} was 0.070 for \code{"MY"}, 0.079 for \code{"FLW"}, and 0.116 for both \code{"SVKZ"} and \code{"PSZ"}. \code{"FLW"} is the steadiest for its cost; \code{"MY"} is the most accurate if the run time is acceptable; \code{"SVKZ"} and \code{"PSZ"} earn their keep when \eqn{\sigma_u} genuinely varies with \eqn{x}, which this design does not test.

\bold{Wrong skew.} The identification of \eqn{\sigma_u} rests on the residuals being negatively skewed. \code{"FLW"} inverts a single sample moment and either succeeds or warns. \code{"SVKZ"} inverts a \emph{local} third moment, which is far noisier, and at any point where the estimated skew has the wrong sign the implied \eqn{\sigma_u(x)^3} is negative; following the paper those points are floored at \eqn{\sigma_u(x) = 0} and their contribution to the frontier gradient is set to zero. \code{wrong.skew} records which observations these were. A large share of them means the local third moment is too noisy to be informative and the \code{"SVKZ"} fit should not be trusted -- \code{"FLW"} is much steadier at moderate sample sizes.

\bold{The \code{np} dependency.} Both estimators need kernel regression and bandwidth selection from the \pkg{np} package, which is listed under \code{Suggests} rather than \code{Imports} because nothing else in \pkg{sfa} requires it. \code{npsfm()} checks for it and stops with an install instruction if it is missing. Bandwidth selection by cross-validation is the dominant cost and scales quadratically in the sample size; supply \code{bw} to skip it.
}
\value{
An object of class \code{"npsfareg"}. This is deliberately \emph{not} an \code{"sfareg"} object: there is no parameter vector with standard errors, so \code{coef()}, \code{vcov()} and \code{logLik()} would have nothing meaningful to return. \code{fitted()}, \code{residuals()}, \code{nobs()}, \code{print()} and \code{summary()} are provided. Components:

\item{frontier}{The estimated frontier \eqn{\hat{m}(x)}, i.e. the kernel fit shifted up by the estimated \eqn{E[u]}.}
\item{frontier.grad}{Matrix of estimated frontier gradients, one row per observation and one column per covariate.}
\item{conditional.mean}{The uncorrected kernel fit of \eqn{E[y|x]}, before the mean shift.}
\item{residuals}{Composed residuals measured against the corrected frontier, \eqn{y - \hat{m}(x)}. These are negative up to noise, unlike the centered residuals of the underlying kernel regression.}
\item{mean.correction}{The estimated \eqn{E[u]}: a scalar for \code{"FLW"}, a vector for \code{"SVKZ"}.}
\item{sigma.u, sigma.v}{Estimated scale parameters. Scalars under \code{"FLW"}; vectors of \eqn{\sigma_u(x_i)}, \eqn{\sigma_v(x_i)} under \code{"SVKZ"}.}
\item{lambda, sigma}{Returned by \code{"FLW"} with \code{dist = "hn"} only: \eqn{\lambda = \sigma_u/\sigma_v} and \eqn{\sigma = \sqrt{\sigma_u^2+\sigma_v^2}}.}
\item{theta}{Returned by \code{"FLW"} with \code{dist = "exp"} or \code{"gamma"}: the rate parameter of the one-sided term.}
\item{b}{Returned by \code{"FLW"} with \code{dist = "unif"}: the estimated upper bound of the uniform.}
\item{sigma.u.grad, wrong.skew}{Returned by \code{"SVKZ"}: the gradient of \eqn{\sigma_u(x)}, and a logical vector flagging observations whose local third moment had the wrong sign. \code{"PSZ"} also returns \code{sigma.u.grad} and \code{sigma.v.grad}, the local-linear slopes of the two log variance functions.}
\item{convergence}{Returned by \code{"PSZ"}: the \code{minqa::bobyqa} status code from each observation's local optimization, \code{0} for success. A large share of non-zero codes means \code{maxit} is too low and the local fits should not be trusted.}
\item{iterations, converged, tol.reached}{Returned by \code{"MY"}: how many outer iterations ran, whether the tolerance was met, and the final squared change in \eqn{(\lambda,\sigma)}.}
\item{prior.fit, dea.efficiency, rts}{Returned by \code{"SZ"}: the smooth frontier that was monotonized, the DEA efficiency scores, and the returns-to-scale assumption used.}
\item{u_hat, exp_u_hat}{Jondrow et al. (1982) inefficiency predictions \eqn{E[u|\varepsilon]} and Battese-Coelli (1988) efficiency predictions \eqn{E[\exp(-u)|\varepsilon]}. Returned when \code{eff = TRUE} and \code{dist} is \code{"hn"} or \code{"exp"}.}
\item{bws}{The bandwidth object(s) used: one for \code{"FLW"}, a list of three (\code{r1}, \code{r2}, \code{r3}) for \code{"SVKZ"}.}
\item{method, dist, formula, call, cost, regtype, bw.sel, nobs, total_time, data}{Settings and bookkeeping.}
}
\references{
Fan, Y., Li, Q. and Weersink, A. (1996) 'Semiparametric estimation of stochastic production frontier models', \emph{Journal of Business & Economic Statistics}, 14(4), pp. 460-468.

Simar, L., Van Keilegom, I. and Zelenyuk, V. (2017) 'Nonparametric least squares methods for stochastic frontier models', \emph{Journal of Productivity Analysis}, 47(3), pp. 189-204.

Jondrow, J., Lovell, C.A.K., Materov, I.S. and Schmidt, P. (1982) 'On the estimation of technical inefficiency in the stochastic frontier production function model', \emph{Journal of Econometrics}, 19(2-3), pp. 233-238.

Battese, G.E. and Coelli, T.J. (1988) 'Prediction of firm-level technical efficiencies with a generalized frontier production function and panel data', \emph{Journal of Econometrics}, 38(3), pp. 387-399.
}
\author{Christopher F. Parmeter and David H. Bernstein}
\note{Requires the \pkg{np} package, which is in \code{Suggests}. Bandwidth cross-validation dominates the run time and scales quadratically in the sample size.}
\seealso{\code{\link{sfm}} for parametric cross-sectional frontiers, \code{\link{psfm}} for panel models, and \code{\link{data_gen_cs}} for simulating data with known true parameters.}
\examples{
\donttest{
if (requireNamespace("np", quietly = TRUE)) {
  set.seed(42)
  n  <- 150
  x1 <- runif(n, 1, 4)
  x2 <- runif(n, 1, 4)
  m  <- 1 + 0.6 * log(x1) + 0.4 * sqrt(x2)   # nonlinear frontier
  d  <- data.frame(y = m + rnorm(n, 0, 0.25) - abs(rnorm(n, 0, 0.6)),
                   x1 = x1, x2 = x2)

  ## Fan, Li and Weersink, normal-half normal
  f <- npsfm(y ~ x1 + x2, data = d, method = "FLW", dist = "hn")
  f
  head(fitted(f))
  head(f$exp_u_hat)

  ## Simar, Van Keilegom and Zelenyuk: sigma_u and sigma_v vary with x
  g <- npsfm(y ~ x1 + x2, data = d, method = "SVKZ")
  summary(g$sigma.u)
}
}
}
