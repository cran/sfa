\name{sfa-package}
\alias{sfa-package}
\docType{package}
\title{\packageTitle{sfa}}
\description{\packageDescription{sfa}}
\details{The DESCRIPTION file:
\packageDESCRIPTION{sfa}
\packageIndices{sfa}}
\keyword{package}
\seealso{\code{http://www.davidharrybernstein.com/software}}
\examples{
\donttest{
## Simple application of the generalized true random effects estimator.
library(sfa)

data_trial <- data_gen_p(t=10,N=100,  rand = 100, 
                         sig_u = 1,   sig_v = 0.3, 
                         sig_r = .2,  sig_h = .4, 
                         cons  = 0.5, beta1 = 0.5,
                         beta2 = 0.5)

psfm(formula    = y_gtre ~ x1 + x2,    
     model_name = "GTRE", 
     data       = data_trial,
     individual = "name",
     PSopt      = FALSE)
}               
}                   