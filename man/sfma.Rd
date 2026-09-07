\name{sfma}
\alias{sfma}
\alias{print.sfma}
\title{Average over inefficiency distributions instead of choosing one}
\description{
Fits several inefficiency distributions to the same data and combines their firm-level inefficiency predictions with data-driven weights. Where \code{\link{TIC}} and \code{\link{vuong}} answer "which distribution?", this answers what to do when the honest reply is "the data does not say" -- selecting one specification and then reporting its efficiencies understates the uncertainty, because the standard errors condition on a choice that was itself estimated.
}
\usage{sfma(formula, data, models = c("NHN", "NE", "NTN", "NR", "NGE"),
     weights = c("sfma", "tic", "aic", "bic", "equal"),
     inefdec = TRUE, quiet = FALSE, ...)}
\arguments{
\item{formula, data, inefdec, ...}{Passed to \code{\link{sfm}} for each candidate.}
\item{models}{Two or more \code{model_name}s to average over. Candidates that fail to fit, or that report no \eqn{E[u \mid \varepsilon]}, are dropped with a message.}
\item{weights}{How to weight the candidates. See \sQuote{Details}.}
\item{quiet}{Suppress the messages about dropped candidates.}
}
\value{
An object of class \code{"sfma"} with components \code{weights}, \code{u_hat} (the averaged \eqn{E[u \mid \varepsilon]}), \code{efficiency}, \code{u_by_model} (the \eqn{n \times S} matrix of per-candidate predictions), \code{fits}, \code{k}, \code{logLik} and, for \code{weights = "sfma"}, \code{criterion}.
}
\details{
\bold{The weighting schemes.}

\code{"sfma"} is the criterion of Parmeter, Wan and Zhang (their eq. 13),
\deqn{C(w) = \|\hat\rho(w) - \hat\rho_{full}\|^2 + n^{1/2}\log(n)\,k'w,}
minimised over the unit simplex. The unobservable target is replaced by the prediction from the most general candidate, exactly as Mallows' \eqn{C_p} replaces an unknown \eqn{\sigma^2} with the full model's estimate; without the penalty all weight would go to the largest candidate by construction. Their Theorem 1 says weights on models that nest the truth go to zero, and they do.

\code{"tic"}, \code{"aic"} and \code{"bic"} are smoothed information-criterion weights in Buckland \emph{et al.}'s (1997) form, \eqn{w_s \propto \exp(-IC_s/2)}. \code{"tic"} is the interesting one, because unlike AIC it does not assume the candidate is correctly specified -- which is the premise of averaging in the first place. \code{"equal"} is the naive average, useful as a baseline.

\bold{\code{"sfma"} needs a unique most general candidate, and refuses without one.} If every candidate has the same number of parameters the penalty \eqn{n^{1/2}\log(n)k'w} is identical across them, so it drops out as a constant and the criterion collapses to "be as close as possible to \eqn{\hat\rho_{full}}" -- minimised by putting weight 1 on whichever model was nominated as full. Measured on exponentially-distributed inefficiency with candidates \code{"NHN"}, \code{"NE"} and \code{"NR"} (all with 5 parameters), that would have placed 0.879 on \code{"NHN"}, the wrong distribution, where \code{"tic"} gave \code{"NE"} 0.934. Rather than break the tie silently, \code{sfma} errors and points at \code{"tic"}. Including a genuinely more general candidate such as \code{"NTN"} makes the criterion well posed.

\bold{Reading the output.} \code{print} reports the per-firm spread of efficiency across candidates. That spread is the reason to average: if the candidates agree, model choice was not costing anything and any single fit would have done; if they disagree, a single-model standard error is not describing the uncertainty a reader cares about.

\bold{Implementation.} The weight problem is a convex quadratic program on the simplex, solved by projected gradient with an exact simplex projection (Duchi \emph{et al.} 2008), so a weight that should be zero is exactly zero. It is checked against brute-force enumeration of the simplex in the tests rather than against itself.
}
\references{
Buckland, S.T., Burnham, K.P. and Augustin, N.H. (1997). Model selection: an integral part of inference. \emph{Biometrics}, 53(2), 603--618.

Duchi, J., Shalev-Shwartz, S., Singer, Y. and Chandra, T. (2008). Efficient projections onto the l1-ball for learning in high dimensions. \emph{ICML}.

Huang, C.J. and Lai, H.-p. (2012). Estimation of stochastic frontier models based on multimodel inference. \emph{Journal of Productivity Analysis}, 38(3), 273--284.

Parmeter, C.F., Wan, A.T.K. and Zhang, X. (2019). Model averaging estimators for the stochastic frontier model. \emph{Journal of Productivity Analysis}, 51(1), 91--103.
}
\seealso{\code{\link{TIC}}, \code{\link{vuong}}, \code{\link{spec_test}}, \code{\link{sfm}}}
\examples{
\donttest{
set.seed(3)
n  <- 300
x1 <- rnorm(n); x2 <- rnorm(n)
## True inefficiency is exponential.
d  <- data.frame(y = 1 + 0.5 * x1 + 0.5 * x2 +
                     rnorm(n, 0, 0.4) - rexp(n, 1),
                 x1 = x1, x2 = x2)

## Information-criterion weights: NE should earn most of the weight.
a <- sfma(y ~ x1 + x2, data = d, models = c("NHN", "NE", "NR"),
          weights = "tic", quiet = TRUE)
a
a$weights

## The Parmeter-Wan-Zhang criterion needs a uniquely largest candidate.
sfma(y ~ x1 + x2, data = d, models = c("NHN", "NE", "NTN"),
     weights = "sfma", quiet = TRUE)
}
}
\keyword{models}
