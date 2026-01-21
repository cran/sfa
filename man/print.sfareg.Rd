\name{print.sfareg}
\alias{print.sfareg}
\title{sfa Object Summaries}
\description{print function for stochastic frontier models of sfm(), zsfm(), and psfm() calls.}
\usage{\method{print}{sfareg}(x, \dots)}
\arguments{
\item{x}{sfa regression objects of the sfm(), zsfm(), and psfm() calls.}
\item{\dots}{Additional arguments passed to other methods}
}
\details{Allows for the usage of print()}
\value{No return value, called for side effects}
\author{David H. Bernstein}
\examples{
\donttest{
library(sfa)     

cs_data_trial   <- data_gen_cs(N= 1000, rand   = 1,  sig_u  = 0.3, sig_v  = 0.3, 
cons   = 0.5,       beta1  = 0.5,   beta2  = 0.5, a      = 4, mu     = 1)

cs.nhnz     <-  sfm(formula    = y_pcs_z ~ x1 +x2| z,    model_name = "NHN",                  
                    data       = cs_data_trial,          PSopt      = TRUE)
print(cs.nhnz)
}
}