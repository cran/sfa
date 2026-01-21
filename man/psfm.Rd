\name{psfm}
\alias{psfm}
\title{psfm}
\description{Function to implement various panel data stochastic frontier estimators}
\usage{psfm(formula, model_name = c("TRE_Z", "GTRE_Z", "TRE",
                    "GTRE", "TFE", "FD", "GTRE_SEQ1", "GTRE_SEQ2"), data,
                    maxit.bobyqa = 100, maxit.psoptim = 10, maxit.optim =
                    10, REPORT = 1, trace = 3, pgtol = 0, individual,
                    halton_num = NULL, start_val = FALSE, gamma = FALSE,
                    PSopt = FALSE, optHessian, inefdec= TRUE, Method = "L-BFGS-B",
                    verbose = FALSE,rand.gtre = NULL, rand.psoptim = NULL)}
\arguments{
  \item{formula}{a symbolic description for the model to be estimated}
  \item{model_name}{model name for the estimation}
  \item{data}{a pdata.frame}
  \item{maxit.bobyqa}{Maximum number of iterations for the bobyqa optimization routine}
  \item{maxit.psoptim}{Maximum number of iterations for the psoptim optimization routine}
  \item{maxit.optim}{Maximum number of iterations for the optim optimization routine}
  \item{REPORT}{reporting parameter}
  \item{trace}{trace}
  \item{pgtol}{pgtol}
  \item{individual}{individual unit in the regression model}
  \item{halton_num}{number of Halton draws to use in SML models}
  \item{start_val}{starting value (optional)}
  \item{gamma}{gamma}
  \item{PSopt}{use psoptim optimization routine (T or F)}
  \item{optHessian}{Logical. Should a numerically differentiated Hessian matrix be returned while using the optim routine? (for optim routine)}
  \item{inefdec}{Production or cost function}
  \item{Method}{The method to be used for optim. See 'Details' within optim.}
  \item{verbose}{Logical. Print optimization progress messages? Default is \code{FALSE}.}
  \item{rand.psoptim}{Integer. Seed for replication of psoptim.  Default to \code{NULL}.}
  \item{rand.gtre}{Integer. Seed for replication of the gtre model.  Default to \code{NULL}.}
}
\details{The generalized true random effects model (GTRE, 4-component model) and true random effects models (TRE)  are both estimated by simulated maximum likelihood based on the paper by the Fillipini and Greene (2016, JPA). The TRE_Z and GTRE_Z allow for modeling the u-component of the GTRE and TRE with determinants of inefficiency. The first-difference estimator (FD) of Wang and Ho (2010, JoE) as well as the True Fixed Effect model estimated by within-maximum likelihood of Chen, Schmidt and Wang (2014, JoE) are also available.}
\references{Fillipini and Greene (2016, JPA); Wang and Ho (2010, JoE); Chen, Schmidt and Wang (2014, JoE)}
\author{David Bernstein}
\note{Standard errors require optHessian set to TRUE}

\seealso{see also}
\value{
An object of class \code{"sfareg"} containing components that vary by model. All models return:

\item{out}{A matrix with parameter estimates, standard errors, and t-values.}
\item{opt}{A list containing the optimization results from the final optimization procedure (not returned for GTRE_SEQ1 and GTRE_SEQ2).}
\item{total_time}{The total computation time for model estimation.}
\item{start_v}{The starting values used in the optimization (not returned for GTRE_SEQ1 and GTRE_SEQ2).}
\item{model_name}{The name of the panel stochastic frontier model estimated.}
\item{formula}{The formula used in the model specification.}
\item{coefficients}{A vector of estimated parameters.}
\item{std.errors}{A vector of standard errors for the estimated parameters (NA if \code{optHessian = FALSE}).}
\item{t.values}{A vector of t-values for the estimated parameters (NA if \code{optHessian = FALSE}).}
\item{call}{The matched call.}
\item{data}{The data used in estimation.}

\bold{Additional model-specific components:}

\bold{For GTRE and GTRE_Z models:}
\item{H}{Predicted time-invariant technical efficiency for each individual.}

\bold{For GTRE, GTRE_Z, TRE and TRE_Z models:}
\item{U}{Predicted time-varying technical efficiency for each observation.}

\bold{For TFE model:}
\item{r_hat_m}{Estimated individual-specific random effects.}
\item{exp_u_hat}{Predicted technical efficiency.}

\bold{For FD model:}
\item{u_hat}{Predicted technical efficiency in levels.}
\item{h_hat}{Estimated z heterogeneity function values.}
\item{exp_u_hat}{Predicted technical efficiency.}

\bold{For GTRE_SEQ1 and GTRE_SEQ2 models:}
\item{other_parms}{A matrix of additional parameters (lambda, sigma, beta_0 for SEQ1; sigma_u, sigma_v, sigma_h, sigma_r, lambda, sigma for SEQ2).}
}

\note{
The GTRE_SEQ1 and GTRE_SEQ2 models use sequential estimation methods and do not return optimization objects or starting values. All panel models require the \code{individual} argument to identify panel units.
}
\examples{
\donttest{
library(sfa)     

data_trial <- data_gen_p(t=10,N=100, rand = 100, 
                         sig_u = 1,  sig_v = 0.3, 
                         sig_r = .2, sig_h = .4, 
                         cons = 0.5, beta1 = 0.5,
                         beta2 = 0.5)

max_tre_z   <-  psfm(formula    = y_tre_z ~ x1 +x2| z_gtre, 
                     model_name = "TRE",                    ## "TRE_Z" also works
                     data       = data_trial,
                     individual = "name",
                     PSopt      = TRUE)
}
}