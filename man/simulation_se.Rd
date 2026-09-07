\name{simulation_se}
\alias{simulation_se}
\title{Simulation Standard Errors for Simulated-ML Fits}
\description{
Reports how much of a fitted parameter's uncertainty is \emph{simulation} noise
rather than sampling noise, by refitting the same data with independent
randomizations of the draw sequence.
}
\usage{
simulation_se(object, K = 10, seeds = NULL, quiet = FALSE,
              envir = parent.frame())
}
\arguments{
\item{object}{An \code{"sfareg"} fit from a simulated-ML model:
  \code{sfm} with \code{"NLN"} or \code{"NW"}, \code{psfm} with \code{"GTRE"},
  \code{"TRE"}, \code{"TRE_Z"} or \code{"GTRE_Z"}, or \code{selsfm}.}
\item{K}{Number of refits, at least 2.}
\item{seeds}{Optional integer vector of randomization seeds; its length
  overrides \code{K}.}
\item{quiet}{Suppress the per-refit progress messages.}
\item{envir}{Environment in which to re-evaluate the stored call. The default
  is usually right; pass it explicitly if the data cannot otherwise be found.}
}
\details{
Simulated maximum likelihood replaces an integral with an average over draws.
A standard error from the Hessian answers \dQuote{how much would this estimate
move in a new sample}, and says nothing about \dQuote{how much would it move
with the \emph{same} sample and different draws}. Randomizing the sequence
makes the second question answerable: refit \eqn{K} times with independent
shifts and take the standard deviation across refits.

This is possible only because the draws are randomized by a uniform shift. An
unrandomized Halton lattice returns the same answer every time, so the
diagnostic would report exactly zero and mean nothing by it.

\bold{Reading the result.} \code{ratio} is the share of the reported standard
error that comes from having simulated the integral. A large value is an
argument for raising \code{Nsim}, not for collecting more data. On the
package's reference design (\code{"NW"}, \eqn{n = 400}, \code{Nsim = 200}) the
frontier slopes come back at about 0.5\% while the variance and shape
parameters are near 8\% -- which independently reproduces the documented
behaviour that too few draws move everything except the slopes.

Quadrupling \code{Nsim} on that design cuts the simulation standard errors by
roughly a factor of six, faster than the factor of two that independent
pseudorandom draws would give. That is the low-discrepancy sequence earning its
place.
}
\value{
A data frame with one row per parameter: \code{estimate}, \code{sampling_se}
(from the Hessian), \code{simulation_se} (across refits) and their
\code{ratio}. Attributes \code{n_refits} and \code{seed_arg} record how many
refits succeeded and which argument was varied.
}
\seealso{\code{\link{sfa_diagnostics}} for convergence diagnostics on a single fit.}
\examples{
\donttest{
## Kept deliberately small: every refit is a full simulated-ML fit, so the
## cost is K + 1 of them. K = 3 is the smallest that still shows the spread.
d <- data_gen_cs(N = 200, rand = 1, sig_u = 1, sig_v = 0.3, cons = 0.5,
                 beta1 = 0.5, beta2 = 0.5, a = 1, mu = 0.5)
fit <- sfm(y_pcs_wb ~ x1 + x2, model_name = "NW", data = d, Nsim = 200)
simulation_se(fit, K = 3, quiet = TRUE)
}
}
\keyword{models}
