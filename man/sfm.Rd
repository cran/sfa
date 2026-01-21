\name{sfm}
\alias{sfm}
\title{sfm}
\description{Implementation of the cross-sectional stochastic frontier model across an array of distributional assumptions for both v and u (user specified). For panel models, see the psfm() call.}
\usage{sfm(formula, model_name, data,maxit.bobyqa,maxit.psoptim,maxit.optim,REPORT,
trace,pgtol,start_val,PSopt,optHessian,inefdec,upper,Method,eta,alpha,verbose=FALSE,
rand.psoptim=NULL)}
\arguments{
\item{formula}{a symbolic description for the model to be estimated}
\item{model_name}{model name for the estimation includes the: normal-half normal (NHN), 
normal-exponential (NE), student's t-half t (THT),  Normal-Rayleigh (NR), 
and the normal-truncated normal (NTN).}
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
\item{eta}{Parameter used for psi-divergence.}
\item{alpha}{Parameter used for MDPD.}
\item{verbose}{Logical. Print optimization progress messages? Default is \code{FALSE}.}
\item{rand.psoptim}{Integer. seed for replication of psoptim.  Default to \code{NULL}.}
}
\details{The options include the Normal-Half Normal (NHN), Normal-exponential (NE), Student's t-Half t (THT),  and the Normal-Truncated Normal (NTN).  NHN_Z and NE_Z are extensions for the NHN and NE models that allow for modeling the u-component of those models with determinants of inefficiency.

Outputs include E[exp(-u)|e] given by exp_u_hat, following Battese and Coelli (1988, JoE), where appropriate.}
\author{David H. Bernstein and Alexander Stead}
\note{Standard errors require optHessian set to TRUE}

\seealso{see also}
\value{
An object of class \code{"sfareg"} containing the following components:

\item{out}{A matrix with parameter estimates, standard errors, and t-values.}
\item{opt}{A list containing the optimization results from the final optimization procedure.}
\item{total_time}{The total computation time for model estimation.}
\item{start_v}{The starting values used in the optimization.}
\item{model_name}{The name of the stochastic frontier model estimated.}
\item{formula}{The formula used in the model specification.}
\item{exp_u_hat}{Predicted technical efficiency (expected values). Available for models: NHN, NHN_Z, NR, NG, and NNAK.}
\item{med_u_hat}{Predicted technical efficiency (median values). Available only for the NHN model.}
\item{coefficients}{A vector of estimated parameters.}
\item{std.errors}{A vector of standard errors for the estimated parameters (NA if \code{optHessian = FALSE}).}
\item{t.values}{A vector of t-values for the estimated parameters (NA if \code{optHessian = FALSE}).}
\item{call}{The matched call.}
}
\examples{
\donttest{
library(sfa)     

cs_data_trial   <- data_gen_cs(N= 1000, rand   = 1,  sig_u  = 0.3, sig_v  = 0.3, 
cons   = 0.5,       beta1  = 0.5,   beta2  = 0.5, a      = 4, mu     = 1)

cs.nhnz     <-  sfm(formula    = y_pcs_z ~ x1 +x2| z,    model_name = "NHN",                  
                    data       = cs_data_trial,          PSopt      = TRUE)
}
}