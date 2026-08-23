\name{sfa_diagnostics}
\alias{sfa_diagnostics}
\alias{plot.sfareg}
\alias{print.sfadiag}
\title{Optimizer Diagnostics for Stochastic Frontier Fits}
\description{Reports how the optimizer actually finished and how well the likelihood identifies the parameters: the convergence code and message, the eigenvalue spectrum and conditioning of the Hessian, the implied parameter correlations, and -- where the objective was retained -- the gradient at the reported optimum and a likelihood slice through each parameter.}
\usage{sfa_diagnostics(object, \dots)

\method{plot}{sfareg}(x, which = 1:4, n_grid = 41, span = 0.25, \dots)

\method{print}{sfadiag}(x, \dots)}
\arguments{
\item{object, x}{An object of class \code{"sfareg"} (for \code{sfa_diagnostics} and \code{plot}), or of class \code{"sfadiag"} (for \code{print}).}
\item{which}{Which panels to draw: \code{1} the Hessian eigenvalue spectrum, \code{2} the parameter correlation matrix, \code{3} a likelihood slice per parameter, \code{4} the gradient. Panels that cannot be drawn are dropped silently; asking for \emph{only} unavailable panels is an error rather than an empty plot.}
\item{n_grid}{Number of points in each likelihood slice. Defaults to \code{41}.}
\item{span}{Half-width of each slice, as a fraction of \code{max(|estimate|, 1)}. Defaults to \code{0.25}.}
\item{\dots}{Currently unused.}
}
\details{
Every stochastic frontier fit in this package runs a staged minimizer and stores \code{optim()}'s convergence code, message, evaluation counts and numerically differentiated Hessian. None of that was reported: a fit that stopped on the iteration cap printed exactly like a converged one. \code{sfa_diagnostics()} is the reporting side of that.

\bold{What the Hessian tells you.} \code{optim()} returns the Hessian of the \emph{negative} log-likelihood, so at a genuine interior optimum it is positive definite. A non-positive eigenvalue means the reported point is not a minimum in that direction. A very small one means the likelihood is nearly flat along it, which is what weak identification looks like numerically -- \code{$hessian$flat_direction} names the parameters that load on that flattest axis, and they are usually the pair the data cannot separate.

\bold{Slices, not profiles.} Panel 3 varies one parameter and holds the others \bold{at their estimates}. A true profile would re-optimize the remaining parameters at every grid point, at the cost of one optimization per point. The slice is the right object for reading curvature at the optimum and for the simplest check there is: the estimate should sit at the peak of its own slice. If it does not, the optimizer stopped early.

\bold{Retaining the objective.} The gradient and the slices need the likelihood itself, which an ordinary fit does not keep. Pass \code{keep_objective = TRUE} to \code{\link{sfm}} to store it. This is off by default because a closure carries its enclosing environment, so a fit saved with one serializes the estimation data alongside the results -- on a 200-observation example the fitted object grows from roughly 38 KB to 1.7 MB.

\bold{The convergence code alone is not diagnostic}, and the flags do not treat it as though it were. Across \code{NHN}, \code{NE} and \code{NTN} at \eqn{n = 150}, \eqn{500} and \eqn{1500}, code \code{52} (\dQuote{ABNORMAL_TERMINATION_IN_LNSRCH}) turns up routinely alongside a maximum relative gradient of about \eqn{10^{-6}} and a positive definite Hessian: the staged minimizer had already found the optimum and the final \code{L-BFGS-B} stage simply could not take a step away from it. The same code on \code{NTN} at \eqn{n = 150} came with a relative gradient of \eqn{5 \times 10^{7}} and an indefinite Hessian, which is a genuine failure. What separates the two is the gradient and the Hessian, so the verdict is built from those and the code is reported as context:

\describe{
  \item{benign}{a line-search code (\code{51}/\code{52}) together with a small relative gradient and a positive definite Hessian. Reported as noise.}
  \item{unverified}{a line-search code with no retained objective, so there is no gradient to settle it. Reported as an open question rather than resolved in the fit's favour -- refit with \code{keep_objective = TRUE}.}
  \item{failure}{anything else, including code \code{1} (the iteration limit), which is never excused: it means the search was still moving when it ran out of budget, whatever the curvature at the stopping point.}
}

\bold{Flags} are raised for a non-zero convergence code, a Hessian that is not positive definite or is worse than \eqn{10^{8}} conditioned, any parameter pair correlated beyond \eqn{|0.95|}, and a large relative gradient. They are phrased as statements about what the fit means, not as generic warnings.
}
\value{
\code{sfa_diagnostics()} returns an object of class \code{"sfadiag"}: a list with \code{model_name}, \code{call}, \code{pnames}, \code{estimates}, \code{convergence} (code, meaning, message, counts, logLik), \code{hessian} (eigenvalues, condition, pos_def, flat_direction), \code{correlation} (matrix, worst_pair, worst_value); \code{convergence} additionally carries \code{benign_nonzero} and \code{unverified_nonzero}, \code{gradient} (gradient, max_abs, relative, max_rel; \code{NULL} without a retained objective), \code{flags}, and \code{has_objective}.

\code{plot.sfareg()} is called for its side effect and returns the \code{"sfadiag"} object invisibly. \code{print.sfadiag()} returns its argument invisibly.
}
\seealso{\code{\link{sfm}}, \code{\link{sfareg-methods}}, \code{\link{summary.sfareg}}}
\author{David H. Bernstein}
\examples{
\donttest{
library(sfa)
d <- data_gen_cs(N = 300, rand = 3, sig_u = 1, sig_v = 0.3,
                 cons = 0.5, beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.5)

## Retaining the objective enables the gradient and the likelihood slices
f <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = d, keep_objective = TRUE)
sfa_diagnostics(f)

## Four-panel diagnostic display
plot(f)

## Works without the objective too, on the Hessian alone
g <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = d)
sfa_diagnostics(g)$hessian$condition
}
}
