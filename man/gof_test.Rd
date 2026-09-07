\name{gof_test}
\alias{gof_test}
\title{Goodness-of-fit test for the assumed inefficiency distribution}
\usage{
gof_test(object, data = NULL, test = c("ks", "chisq"),
         null = c("bootstrap", "asymptotic"), B = 199, cells = 10,
         seed = NULL)
}
\arguments{
\item{object}{An \code{"sfareg"} fit from \code{\link{sfm}}.}
\item{data}{The data the model was fitted to. \code{sfm()} does not retain it,
so it must be supplied or be findable from the fit's \code{call}.}
\item{test}{\code{"ks"}, \code{"chisq"}, or both.}
\item{null}{\code{"bootstrap"} for the parametric bootstrap of Wang, Amsler
and Schmidt, or \code{"asymptotic"}.}
\item{B}{Bootstrap replications. Each one refits the model.}
\item{cells}{Number of equiprobable cells for the chi-square test.}
\item{seed}{Optional seed for the bootstrap.}
}
\value{
A data frame with one row per test: \code{test}, \code{statistic},
\code{null}, \code{p.value} and \code{B} (replications actually used).
}
\description{
Tests the distributional assumption on \eqn{u}, holding the normality of
\eqn{v} fixed, following Wang, Amsler and Schmidt (2011).
}
\details{
The distribution of \eqn{u} is the assumption in this model least often
defended and least often tested, and it is testable: with the noise
distribution maintained, the assumed \eqn{u} implies a distribution for the
composed error, so rejecting that is rejecting the assumed \eqn{u}.

\bold{The test is on \eqn{\varepsilon}, not on \eqn{\hat u}.} Wang, Amsler
and Schmidt point out that \eqn{\hat u = E[u \mid \varepsilon]} is a
monotonic function of \eqn{\varepsilon}, so the KS test is identical either
way and the chi-square test is identical when the cells are defined
conformably -- but \eqn{\varepsilon} is much easier to work with. It is a
mistake, not a diagnostic, to compare the observed spread of \eqn{\hat u}
with the assumed density of \eqn{u}: those are different distributions.

\bold{Parameter estimation cannot be ignored.} Both statistics are evaluated
at \eqn{\hat\theta}, which changes their null distributions. The default
\code{"bootstrap"} handles this by copying the estimation step exactly: each
replication draws a composed error from the fitted model, rebuilds the
response, \emph{refits}, forms its own residuals, and recomputes the
statistic. \code{"asymptotic"} returns a \eqn{\chi^2(k-1-m)} p-value for the
chi-square statistic, which the authors note is \bold{conservative} when
evaluated at the MLE -- that reference distribution belongs to the
minimum-chi-square estimator. There is no asymptotic KS p-value here: with
\eqn{\theta} estimated the Kolmogorov distribution does not apply, and Bai's
(2003) martingale transformation is not implemented, so \code{NA} is returned
rather than a number that would look valid.

\bold{Size and power.} Measured in this package over 200 replications with
\code{B = 49}, half-normal data fitted as half-normal for size and
exponential data fitted as half-normal for power:

\tabular{lrrrrr}{
 \tab \eqn{n} \tab KS .10 \tab KS .05 \tab chi2 .10 \tab chi2 .05 \cr
size \tab 200 \tab 0.065 \tab 0.025 \tab 0.085 \tab 0.065 \cr
size \tab 800 \tab 0.070 \tab 0.040 \tab 0.065 \tab 0.010 \cr
power \tab 200 \tab 0.505 \tab 0.395 \tab 0.235 \tab 0.145 \cr
power \tab 800 \tab 0.985 \tab 0.965 \tab 0.810 \tab 0.655 \cr
power \tab 2000 \tab 1.000 \tab 1.000 \tab 0.990 \tab 0.985
}

Both tests hold their size, mildly conservatively. \bold{The KS test
dominates the chi-square test at every sample size}, which is the paper's own
conclusion and the reason it is the better default: the chi-square statistic
discards the within-cell information that KS uses. Power against the
half-normal/exponential pair -- one of the harder distinctions, since both
are one-parameter and similarly shaped -- is about 0.4 at \eqn{n = 200} and
essentially 1 by \eqn{n = 800}, so a failure to reject in a small sample is
weak evidence.
}
\references{
Wang, W. S., Amsler, C. and Schmidt, P. (2011). Goodness of fit tests in
stochastic frontier models. \emph{Journal of Productivity Analysis}
\bold{35}, 95--118.

Bai, J. (2003). Testing parametric conditional distributions of dynamic
models. \emph{Review of Economics and Statistics} \bold{85}, 531--549.
}
\seealso{\code{\link{spec_test}}, \code{\link{pcomposed_model}}, \code{\link{sfm}}}
\examples{
d <- as.data.frame(data_gen_cs(N = 200, rand = 5, sig_u = 1, sig_v = 0.5,
                               cons = 0.5, beta1 = 0.5, beta2 = 0.5,
                               a = 5, mu = 0.1))
f <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = d)
\donttest{
## Every bootstrap replication REFITS, so keep B small in an example; use
## B = 999 for a real test.
gof_test(f, data = d, B = 39, seed = 1)
}
gof_test(f, data = d, test = "chisq", null = "asymptotic")
}
