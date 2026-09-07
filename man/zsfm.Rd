\name{zsfm}
\alias{zsfm}
\title{Zero-Inefficiency Stochastic Frontier Models}
\description{Fits the zero-inefficiency stochastic frontier model of Kumbhakar, Parmeter and Tsionas (2013), in which the sample is a mixture of fully efficient firms and firms drawn from an ordinary stochastic frontier.}
\usage{zsfm(formula, model_name = c("ZISF", "ZISF_Z"), 
data, maxit.bobyqa = 10000,maxit.psoptim = 1000, maxit.optim = 1000, 
REPORT = 1, trace = 0, pgtol = 0,start_val = FALSE,PSopt = FALSE, 
optHessian = TRUE, inefdec = TRUE, upper = NA, 
Method = "L-BFGS-B",logit = TRUE,verbose=FALSE,rand.psoptim = NULL)}
\arguments{
\item{formula}{a symbolic description for the model to be estimated}
\item{model_name}{Which mixture specification to fit. \code{"ZISF"} holds the probability of being fully efficient constant across the sample; \code{"ZISF_Z"} lets that probability depend on covariates supplied in the second part of the formula (\code{y ~ x | z}). See \sQuote{Details}. Matching ignores case, so \code{"nhn"}, \code{"NHN"} and \code{"Nhn"} are the same choice.}
\item{data}{A data frame containing the variables named in \code{formula}.}
\item{maxit.bobyqa}{Maximum number of iterations for the bobyqa optimization routine}
\item{maxit.psoptim}{Maximum number of iterations for the psoptim optimization routine}
\item{maxit.optim}{Maximum number of iterations for the optim optimization routine}
\item{REPORT}{reporting parameter}
\item{trace}{Integer. Tracing level passed through to the optimizer; larger values print more. \code{0} (the default) is silent.}
\item{pgtol}{Numeric. Projected-gradient tolerance passed to \code{optim()}'s \code{"L-BFGS-B"} method. \code{0} uses the optimizer's own default.}
\item{start_val}{starting value (optional)}
\item{PSopt}{use psoptim optimization routine (T or F)}
\item{optHessian}{Logical. Should a numerically differentiated Hessian matrix be returned while using the optim routine? (for optim routine)}
\item{inefdec}{Production or cost function}
\item{upper}{Vector of upper values for the optim package.}
\item{Method}{The method to be used for optim. See 'Details' within optim.}
\item{logit}{Logical, \code{"ZISF_Z"} only. Link mapping the linear predictor \eqn{z'\gamma} to the probability of full efficiency. \code{TRUE} (the default) uses the logistic link \eqn{\Lambda(z'\gamma)}, which covers \eqn{(0,1)}. \code{FALSE} uses \eqn{\Phi(z'\gamma)/(1 + \Phi(z'\gamma))}, which is \bold{bounded above by 0.5}: under that link no observation can be assigned a probability of full efficiency exceeding one half, whatever the data. This is not the probit analogue of the default and is retained only for backward compatibility; prefer \code{logit = TRUE} unless you specifically want that restriction.}
\item{verbose}{Logical. Print optimization progress messages? Default is \code{FALSE}.}
\item{rand.psoptim}{Integer. seed for replication of psoptim.  Default to \code{NULL}.}
}
\details{The zero-inefficiency stochastic frontier model treats the sample as a mixture of two regimes: with probability \eqn{p} a firm sits exactly on the frontier (\eqn{u = 0}, noise only), and with probability \eqn{1 - p} it is drawn from an ordinary normal-half normal frontier. It is the natural specification when some fraction of the sample is thought to be fully efficient, which an ordinary frontier cannot represent because it forces \eqn{u > 0} almost surely.

Under \code{"ZISF"} the mixing probability is a single constant, parameterized as \eqn{p = \exp(-|\gamma|)}. That parameterization is exactly symmetric in \eqn{\gamma}: \eqn{+\gamma} and \eqn{-\gamma} give identical fits and the optimizer may return either, so the sign of the reported \code{gamma} carries no meaning. Under \code{"ZISF_Z"} the probability varies with covariates \eqn{z} through the link chosen by \code{logit}.

\code{post.prob} reports the posterior probability that each observation belongs to the fully efficient regime, given its residual; \code{jlms} reports the Jondrow et al. (1982) inefficiency prediction.

Note on default iteration counts: like \code{sfm()}, \code{zsfm()}'s likelihood is closed-form, so its optimization defaults (thousands of iterations) are set much higher than \code{psfm()}'s (which uses simulated maximum likelihood for several models and defaults to far fewer iterations for speed). See \code{\link{psfm}}'s Details for more.}
\references{Kumbhakar, S.C., Parmeter, C.F. and Tsionas, E.G. (2013) 'A zero inefficiency stochastic frontier model', \emph{Journal of Econometrics}, 172(1), pp. 66-76. doi:10.1016/j.jeconom.2012.08.021.}
\author{Christopher F. Parmeter and David H. Bernstein}
\note{Standard errors require optHessian set to TRUE}
\seealso{\code{\link{sfm}}, \code{\link{psfm}}, \code{\link{ttsfm}}, \code{\link{sfareg-methods}}, \code{\link{panel89}}}
\value{
An object of class \code{"sfareg"} containing the following components:

\item{out}{A matrix with parameter estimates, standard errors, and t-values.}
\item{opt}{A list containing the optimization results from the final optimization procedure.}
\item{total_time}{The total computation time for model estimation.}
\item{start_v}{The starting values used in the optimization.}
\item{model_name}{The name of the zero-inefficiency stochastic frontier model estimated (\code{"ZISF"} or \code{"ZISF_Z"}).}
\item{formula}{The formula used in the model specification.}
\item{jlms}{Predicted technical efficiency using the Jondrow et al. (1982) conditional mean estimator (JLMS).}
\item{post.prob}{Posterior probabilities of being fully efficient.}
\item{coefficients}{A vector of estimated parameters.}
\item{std.errors}{A vector of standard errors for the estimated parameters (NA if \code{optHessian = FALSE}).}
\item{t.values}{A vector of t-values for the estimated parameters (NA if \code{optHessian = FALSE}).}
\item{call}{The matched call.}
}
\examples{
\donttest{
library(sfa)  

eqz     <- y ~ q1 + q2 + q3 + q4 + q5 + w1 + w2 + w3 + w4 | z

data(panel89)

zsfm(formula    = eqz,
     model_name = "ZISF_Z",
     data       = panel89,
     logit      = TRUE)
}
} 
