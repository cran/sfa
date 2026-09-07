\name{print.sfareg}
\alias{print.sfareg}
\title{Print method for "sfareg" objects}
\description{Compact display of a stochastic frontier fit, as returned by \code{\link{sfm}}, \code{\link{psfm}}, \code{\link{zsfm}} or \code{\link{ttsfm}}.}
\usage{\method{print}{sfareg}(x, \dots)}
\arguments{
\item{x}{An object of class \code{"sfareg"}, as returned by \code{\link{sfm}}, \code{\link{psfm}}, \code{\link{zsfm}} or \code{\link{ttsfm}}.}
\item{\dots}{Additional arguments passed to other methods}
}
\details{Prints the model name, the estimated coefficients with their standard errors and t-values (the \code{out} matrix), and the log-likelihood where the estimator has one. For a fuller display, including the frontier and variance-parameter blocks, use \code{\link{summary.sfareg}}.}
\value{No return value, called for side effects}
\author{David H. Bernstein}
\seealso{\code{\link{summary.sfareg}}, \code{\link{sfareg-methods}}}
\examples{
\donttest{
library(sfa)     

cs_data_trial   <- data_gen_cs(N= 1000, rand   = 1,  sig_u  = 0.3, sig_v  = 0.3, 
cons   = 0.5,       beta1  = 0.5,   beta2  = 0.5, a      = 4, mu     = 1)

cs.nhnz     <-  sfm(formula    = y_pcs_z ~ x1 +x2| z,    model_name = "NHN_Z",
                    data       = cs_data_trial,          PSopt      = TRUE)
print(cs.nhnz)
}
}
