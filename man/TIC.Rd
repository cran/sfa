\name{TIC}
\alias{TIC}
\alias{vuong}
\alias{print.sfa_vuong}
\title{Choosing between non-nested stochastic frontier specifications}
\description{
Model selection for stochastic frontier models when the candidates are \bold{not nested} and none of them need be correctly specified. \code{TIC} is Takeuchi's information criterion, a version of \code{\link{AIC}} whose penalty is estimated from the data rather than assumed equal to the number of parameters. \code{vuong} is Vuong's pairwise test, which returns a p-value and is allowed to conclude that two specifications are indistinguishable.
}
\usage{TIC(object, detail = FALSE)

vuong(object1, object2, correction = c("none", "aic", "tic"), level = 0.05)}
\arguments{
\item{object, object1, object2}{Objects of class \code{"sfareg"}, fitted with \code{keep_objective = TRUE} (see \sQuote{Details}).}
\item{detail}{If \code{FALSE} (the default) \code{TIC} returns a single number. If \code{TRUE} it returns a list also holding the log-likelihood, the parameter count, the estimated penalty and their ratio.}
\item{correction}{Bias correction applied to the log-likelihood difference before it is standardised. \code{"none"} (the default) is the statistic whose \eqn{N(0,1)} limit Vuong established. \code{"aic"} subtracts the parameter counts and \code{"tic"} the Takeuchi penalties; both shift the mean of the statistic and neither has an established null distribution, so their p-values are flagged as indicative.}
\item{level}{Two-sided size used for the reported critical value and the verdict.}
}
\value{
\code{TIC} returns a single numeric value, or with \code{detail = TRUE} a list with components \code{TIC}, \code{AIC}, \code{logLik}, \code{df}, \code{penalty} and \code{ratio}.

\code{vuong} returns an object of class \code{"sfa_vuong"} with components \code{statistic}, \code{p.value}, \code{n}, \code{lr}, \code{lr_adjusted}, \code{omega}, \code{correction}, \code{penalty}, \code{critical}, \code{level}, \code{models}, \code{favoured} and \code{exact_null}. Lower \code{TIC} is better; \code{favoured} is the name of the winning \code{model_name}, or \code{"neither"}.
}
\details{
\bold{Why not AIC.} \code{sfa} offers many \code{model_name}s for the same data, and \code{\link{AIC}} assumes the fitted model is correctly specified -- exactly the assumption in doubt when choosing between \code{"NHN"}, \code{"NE"}, \code{"NG"} and \code{"NTN"}. Most such pairs are not nested either, so the ordinary likelihood ratio test has no chi-square limit.

\bold{What TIC does.} Takeuchi's criterion replaces AIC's penalty \eqn{d} with \eqn{\mathrm{tr}[H(\hat\theta) I(\hat\theta)^{-1}]}, where \eqn{I} is the sample Fisher information and \eqn{H} the outer product of the per-observation scores. Under correct specification the information matrix equality gives \eqn{H = I}, the trace collapses to \eqn{d}, and TIC equals AIC. The gap between \code{penalty} and \code{df} -- reported as \code{ratio} under \code{detail = TRUE} -- is therefore a readable diagnostic in its own right: a ratio near one says the distributional assumption is not doing visible damage.

\bold{What vuong does.} It tests \eqn{H_0} that the two models are equally close to the truth in the Kullback-Leibler sense. Neither model has to be correct. With \eqn{m_i} the difference of per-observation log-likelihoods, the statistic is \eqn{n^{-1/2}\sum_i m_i / \hat\sigma}, asymptotically standard normal. Because it is a test rather than a criterion it can return \code{"neither"}, which is frequently the honest answer and is not available from AIC.

\bold{Requirement.} Both need per-observation log-likelihoods, so the fits must be made with \code{keep_objective = TRUE}. This is supported by \code{\link{sfm}} (for \code{"NHN"}, \code{"NE"}, \code{"NR"}, \code{"NG"}, \code{"NNAK"}, \code{"THT"}, \code{"NTN"}, \code{"NHN_Z"}, \code{"NE_Z"}, \code{"NU"}, \code{"NGE"}, \code{"NLN"}, \code{"NW"}, \code{"tHN"} and \code{"TSL"}) and by \code{\link{psfm}}. \code{vuong} additionally requires both fits to use the same observations, which it checks by sample size and warns about when the \code{data} arguments differ.

\bold{Estimators that maximise nothing.} \code{psfm}'s \code{"GTRE_SEQ1"}, \code{"GTRE_SEQ2"} and \code{"SSFE"}, and \code{ivsfm}'s \code{"C2SLS"}, are moment-based or FE and have no log-likelihood, so no likelihood-based criterion applies to them and both functions refuse rather than returning a number built on \code{NA}.
}
\references{
Lai, H.-p. and Huang, C.J. (2010). Likelihood ratio tests for model selection of stochastic frontier models. \emph{Journal of Productivity Analysis}, 34(1), 3--13.

Takeuchi, K. (1976). Distribution of informational statistics and a criterion of model fitting. \emph{Suri-Kagaku (Mathematical Sciences)}, 153, 12--18.

Vuong, Q.H. (1989). Likelihood ratio tests for model selection and non-nested hypotheses. \emph{Econometrica}, 57(2), 307--333.
}
\seealso{\code{\link{sfm}}, \code{\link{skewness_test}}, \code{\link{sfa_diagnostics}}, \code{\link{AIC}}}
\examples{
\donttest{
set.seed(42)
n  <- 400
x1 <- rnorm(n); x2 <- rnorm(n)
d  <- data.frame(y = 1 + 0.5 * x1 + 0.5 * x2 +
                     rnorm(n, 0, 0.4) - abs(rnorm(n)),
                 x1 = x1, x2 = x2)

a <- sfm(y ~ x1 + x2, model_name = "NHN", data = d, keep_objective = TRUE)
b <- sfm(y ~ x1 + x2, model_name = "NE",  data = d, keep_objective = TRUE)

## Lower is better, and comparable across non-nested models.
TIC(a)
TIC(b)

## Is the penalty near the parameter count? If so, TIC has little to add.
TIC(a, detail = TRUE)$ratio

## And a test, which may decline to pick a winner.
vuong(a, b)
}
}
\keyword{models}
