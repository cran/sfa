\name{hscore_select}
\alias{hscore_select}
\alias{print.sfa_hscore}
\title{Select the robustness tuning parameter by Hyvarinen score}
\description{
Refits a Normal--Half-Normal stochastic frontier at each candidate value of the
robustness tuning parameter and returns the value minimising the Hyvarinen score
of \code{\link{hscore}}. The search is coarse-to-fine and includes the
maximum-likelihood endpoint \eqn{c=0}, so the criterion is free to decline
robustification altogether.
}
\usage{
hscore_select(object, method = c("mlqe", "psi", "mdpd"),
              range = c(0, 0.60), coarse = 0.01, fine = 0.001,
              half_width = 0.02, multistart = TRUE, verbose = FALSE)
}
\arguments{
\item{object}{A fitted model from \code{\link{sfm}} with
\code{model_name = "NHN"}, used for the data, the formula and the reference
starting values.}
\item{method}{Robust criterion. \code{"psi"} and \code{"mdpd"} are the same
estimator at a common tuning value.}
\item{range}{Numeric length-2 search interval.}
\item{coarse, fine}{Grid steps for the two passes.}
\item{half_width}{Half-width of the refinement window around the coarse
minimum.}
\item{multistart}{Logical; see Details.}
\item{verbose}{Logical, print progress.}
}
\value{
An object of class \code{"sfa_hscore"}: the selected tuning parameter \code{c},
its criterion value \code{hscore}, the whole search \code{path} as a data frame,
the refitted model \code{fit}, and \code{n_finite}, how many candidates were
scorable out of how many tried.
}
\details{
Two failure modes are guarded here, both of which return a converged fit and a
finite criterion value if they are not.

\emph{Silent grid loss.} The criterion is evaluated in log space (see
\code{\link{hscore}}). Evaluated naively it is non-finite wherever any fitted
density underflows, which discards candidates non-randomly -- those nearest
maximum likelihood survive -- and so understates the robustness the data call
for.

\emph{Warm-start capture.} Walking the grid upward and warm-starting each fit
from the previous one is efficient but unguarded: if one fit falls into the
degenerate \eqn{\sigma_u \to 0} basin, every later fit inherits it. With
\code{multistart = TRUE} each candidate is fitted from both the warm start and a
fixed reference start, keeping the better objective.
}
\section{Reading the result}{
A flat region around the minimum is not non-identification in the usual sense; it
means the criterion discriminates weakly among nearby tuning values. The drop
from \eqn{c=0} to the minimum is a descriptive score difference, not a test of
maximum likelihood against a robust alternative. Where the selected value
matters, compare it with the fixed calibration from \code{\link{calibrate_c}} and
check its stability by resampling: the selection can be bimodal across
resamples even when the full-sample criterion clearly prefers a robust solution.
}
\references{
Sugasawa, S. and Yonekura, S. (2021). On selection criteria for the tuning
parameter in robust divergence. \emph{Entropy} 23(9), 1147.

Bernstein, D.H., Parmeter, C.F. and Wright, I.A. (2026). On Robust Estimation of the Stochastic Frontier Model. Working paper.
}
\examples{
\donttest{
set.seed(3)
n <- 200
x <- runif(n, 1, 10)
y <- 1 + 0.5 * log(x) + rnorm(n, 0, 0.3) - abs(rnorm(n, 0, 0.6))
dat <- data.frame(y = y, x = log(x))
fit <- sfm(y ~ x, data = dat, model_name = "NHN")

sel <- hscore_select(fit, method = "mlqe")
sel
plot(sel$path$c, sel$path$hscore, type = "l",
     xlab = "c", ylab = "Hyvarinen score")
abline(v = sel$c, lty = 3)
}
}
\seealso{\code{\link{hscore}}, \code{\link{calibrate_c}},
\code{\link{density_weights}}, \code{\link{sfm}}}
