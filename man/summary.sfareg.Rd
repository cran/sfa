\name{summary.sfareg}
\alias{summary.sfareg}
\title{Summary method for "sfareg" objects}
\description{Summary of a stochastic frontier fit, as returned by \code{\link{sfm}}, \code{\link{psfm}}, \code{\link{zsfm}} or \code{\link{ttsfm}}.}
\usage{\method{summary}{sfareg}(object, \dots)}
\arguments{
\item{object}{An object of class \code{"sfareg"}, as returned by \code{\link{sfm}}, \code{\link{psfm}}, \code{\link{zsfm}} or \code{\link{ttsfm}}.}
\item{\dots}{Additional arguments passed to other methods}
}
\details{Prints the call, the model name, and the estimated parameters with standard errors and t-values, together with the log-likelihood and the number of observations where the estimator provides them. Which parameters appear depends on the model: several report the \eqn{\lambda = \sigma_u/\sigma_v}, \eqn{\sigma = \sqrt{\sigma_u^2 + \sigma_v^2}} reparameterization rather than the raw scale parameters.}
\value{Invisibly returns its \code{object} argument unchanged, having printed the summary as a side effect.}
\author{David H. Bernstein}
\seealso{\code{\link{print.sfareg}}, \code{\link{sfareg-methods}}}
\examples{
\donttest{
library(sfa)     

cs_data_trial   <- data_gen_cs(N= 1000, rand   = 1,  sig_u  = 0.3, sig_v  = 0.3, 
cons   = 0.5,       beta1  = 0.5,   beta2  = 0.5, a      = 4, mu     = 1)

cs.nhnz     <-  sfm(formula    = y_pcs_z ~ x1 +x2| z,    model_name = "NHN_Z",
                    data       = cs_data_trial,          PSopt      = TRUE)
summary(cs.nhnz)                    
}
}