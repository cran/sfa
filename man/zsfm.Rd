\name{zsfm}
\alias{zsfm}
\title{Zero-Inflated Stochastic Frontier Model}
\description{Code to use the Zero-Inflated Stochastic Frontier Model}
\usage{zsfm(formula, model_name = c("ZISF", "ZISF_Z"), 
data, maxit.bobyqa = 10000,maxit.psoptim = 1000, maxit.optim = 1000, 
REPORT = 1, trace = 0, pgtol = 0,start_val = FALSE,PSopt = FALSE, 
optHessian, inefdec = TRUE, upper = NA, 
Method = "L-BFGS-B",logit = TRUE,verbose=FALSE,rand.psoptim = NULL)}
\arguments{
\item{formula}{a symbolic description for the model to be estimated}
\item{model_name}{model name for the estimation}
\item{data}{A data set}
\item{maxit.bobyqa}{Maximum number of iterations for the bobyqa optimization routine}
\item{maxit.psoptim}{Maximum number of iterations for the psoptim optimization routine}
\item{maxit.optim}{Maximum number of iterations for the optim optimization routine}
\item{REPORT}{reporting parameter}
\item{trace}{trace}
\item{pgtol}{pgtol}
\item{start_val}{starting value (optional)}
\item{PSopt}{use psoptim optimization routine (T or F)}
\item{optHessian}{Logical. Should a numerically differentiated Hessian matrix be returned while using the optim routine? (for optim routine)}
\item{inefdec}{Production or cost function}
\item{upper}{Vector of upper values for the optim package.}
\item{Method}{The method to be used for optim. See 'Details' within optim.}
\item{logit}{Choice of using logit function}
\item{verbose}{Logical. Print optimization progress messages? Default is \code{FALSE}.}
\item{rand.psoptim}{Integer. seed for replication of psoptim.  Default to \code{NULL}.}
}
\details{Example based on: A zero inefficiency stochastic frontier model, Journal of Econometrics, S. C. Kumbhakar, C. F. Parmeter and E. G. Tsionas, 2013}
\references{S. C. Kumbhakar, C. F. Parmeter and E. G. Tsionas (2013)}
\author{Chris F. Parmeter and David H. Bernstein}
\note{Standard errors require optHessian set to TRUE}
\seealso{panel89}
\value{
An object of class \code{"sfareg"} containing the following components:

\item{out}{A matrix with parameter estimates, standard errors, and t-values.}
\item{opt}{A list containing the optimization results from the final optimization procedure.}
\item{total_time}{The total computation time for model estimation.}
\item{start_v}{The starting values used in the optimization.}
\item{model_name}{The name of the zero-inflated stochastic frontier model estimated (ZISF or ZISF_Z).}
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