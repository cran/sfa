\name{lcsfm}
\alias{lcsfm}
\title{Latent Class Stochastic Frontier Models}
\description{Fits the latent class stochastic frontier model of Greene (2005) and Orea and Kumbhakar (2004), in which several unobserved technologies coexist in one sample and each firm contributes to every one of them, weighted by its class probability.}
\usage{lcsfm(formula, model_name = c("LCM", "LCM_Z", "LCM_CN"), 
data, n_class = 2, maxit.bobyqa = 10000, maxit.psoptim = 1000, 
maxit.optim = 1000, REPORT = 1, trace = 0, pgtol = 0, 
start_val = FALSE, PSopt = FALSE, optHessian = TRUE, 
inefdec = TRUE, penalty_c = 0, upper = NA, Method = "L-BFGS-B", 
verbose = FALSE, rand.psoptim = NULL)}
\arguments{
\item{formula}{a symbolic description for the model to be estimated. For \code{"LCM_Z"} the second part gives the covariates driving class membership (\code{y ~ x | z}).}
\item{model_name}{Which specification to fit. \code{"LCM"} holds the class probabilities constant across the sample and lets \bold{every} parameter vary by class; \code{"LCM_Z"} is \code{"LCM"} with the class probabilities depending on covariates supplied in the second part of the formula. \code{"LCM_CN"} is the \bold{contaminated normal} frontier, in which only the \emph{noise scale} varies between components and the frontier, the inefficiency scale and the class probabilities are all common -- see \sQuote{Details}. Matching ignores case, so \code{"lcm"}, \code{"LCM"} and \code{"Lcm"} are the same choice.}
\item{data}{A data frame containing the variables named in \code{formula}.}
\item{n_class}{Integer. The number of latent classes \eqn{J}, at least 2. Defaults to \code{2}. The parameter count grows as \eqn{J(k + 2)} plus \eqn{(J-1)q}, so large \eqn{J} needs a correspondingly large sample; see \sQuote{Details}.}
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
\item{penalty_c}{Tuning constant for Chen \emph{et al.}'s (2001) penalty \eqn{c\log(J^J\prod_j p_j)}, which is added to the log-likelihood so that no class probability can be driven to 0 or 1. The default \code{0} leaves the ordinary likelihood untouched. A positive value is what \code{\link{lcsfm_homogeneity}} needs, and it changes what is maximised: \code{$opt$value} then holds the penalised objective, while \code{$logLik_unpenalised} and \code{$penalty} carry the two pieces separately. Only available for \code{model_name = "LCM"}, since the penalty is defined on a scalar mixing proportion.}
\item{upper}{Vector of upper values for the optim package.}
\item{Method}{The method to be used for optim. See 'Details' within optim.}
\item{verbose}{Logical. Print optimization progress messages? Default is \code{FALSE}.}
\item{rand.psoptim}{Integer. seed for replication of psoptim. Default to \code{NULL}.}
}
\details{
\bold{\code{"LCM_CN"}, the contaminated normal frontier.} A different use of the same machinery: instead of separating technologies, it gives the \emph{noise} heavier tails. The noise density is a scale mixture of normals, \eqn{f_v(v) = \sum_j p_j \phi(v/\sigma_{vj})/\sigma_{vj}}, while the frontier, \eqn{\sigma_u} and the mixing proportions stay common across components. At \eqn{J = 2} this is the scale-contaminated normal, an alternative to a Student-\eqn{t} noise term.

It has a closed form, which is why it is cheap. Since \eqn{f_\epsilon(e) = \int f_v(e + Su) f_u(u)\,du} is linear in \eqn{f_v}, a mixture noise density passes straight through:

\deqn{f_\epsilon(e) = \sum_j p_j\, f_{NHN}(e; \sigma_{vj}, \sigma_u)}

a mixture of ordinary normal/half-normal densities sharing one \eqn{\sigma_u}. No extra integration is required, and the identity is verified against direct numerical integration in the tests.

The parameter vector is \eqn{(\sigma_{v1},\dots,\sigma_{vJ}, \sigma_u, \beta, \mathrm{logit})} rather than one block per class.

\bold{It is also the specification \code{\link{lcsfm_homogeneity}}'''s asymptotic null is stated for}: one scalar parameter differing between components. On \code{"LCM"}, where everything varies, the \eqn{\chi^2_{0:1}} null rejected 63.5\% of the time at a nominal 5\% under a true null; on \code{"LCM_CN"} the same null gives 8.0\%. Still mildly liberal, and the test says so, but usable.

The latent class stochastic frontier supposes that \eqn{J} distinct technologies are present in one sample and that which technology a firm operates is unobserved. Each class \eqn{j} has its own frontier and its own two scales,
\deqn{y_i = x_i'\beta_j + v_{ij} - u_{ij}, \quad v_{ij} \sim N(0, \sigma_{vj}^2), \quad u_{ij} \sim N^+(0, \sigma_{uj}^2),}
and class membership follows a multinomial logit,
\deqn{P(j \mid q_i) = \exp(q_i'\delta_j) / \sum_m \exp(q_i'\delta_m), \quad \delta_J = 0,}
with class \eqn{J} the reference. Every firm contributes to every class, weighted by its class probability, so the log-likelihood mixes the \eqn{J} composed densities rather than assigning firms to groups:
\deqn{\log L_i = \log \sum_j P(j \mid q_i) f_j(\varepsilon_{ij}).}
Under \code{"LCM"} the \eqn{q_i} are a constant, so the class probabilities are the same for every firm; under \code{"LCM_Z"} they are the covariates in the second part of the formula.

This is a different object from fitting separate frontiers to known subgroups: the groups are not known, and the uncertainty about them is carried through into the efficiency predictions. It is also more general than \code{\link{zsfm}}'s zero-inefficiency model, which is the restricted two-class case in which one class has no inefficiency at all.

\code{post.prob} is an \eqn{n \times J} matrix of posterior class probabilities, \code{jlms_class} holds the class-conditional Jondrow et al. (1982) predictions, and \code{jlms} is their posterior-weighted average -- a firm's inefficiency is only defined relative to a frontier, and which frontier it faces is itself uncertain, so reporting the modal class alone would discard that. \code{class} gives the modal class and \code{class_prob} the average prior class shares.

\bold{Label switching.} Class labels in a finite mixture are identified only up to permutation: relabelling the classes and permuting \eqn{\delta} gives the identical likelihood. Starting values are built by splitting the OLS residuals at their \eqn{J}-quantiles and refitting within each group, which makes a run reproducible on the same data -- and avoids the saddle point that identical starting components would sit at, since equal posterior class probabilities give a zero score with respect to the class split -- but it does not pin the labels to any external ordering. Do not compare class 1 across two fits, or against a simulation's true labels, without matching the classes first.

\bold{Choosing \eqn{J}.} The likelihood cannot fall when a class is added, and the usual regularity conditions for a likelihood-ratio test of \eqn{J} against \eqn{J+1} do not hold, since the null puts a class probability on the boundary. Use \code{AIC()}/\code{BIC()}, the separation of the fitted frontiers, and whether the extra class is economically interpretable, rather than a p-value.}
\references{Greene, W. (2005) 'Reconsidering heterogeneity in panel data estimators of the stochastic frontier model', \emph{Journal of Econometrics}, 126(2), pp. 269-303. doi:10.1016/j.jeconom.2004.05.003.

Orea, L. and Kumbhakar, S.C. (2004) 'Efficiency measurement using a latent class stochastic frontier model', \emph{Empirical Economics}, 29(1), pp. 169-183. doi:10.1007/s00181-003-0184-2.

Caudill, S.B. (2003) 'Estimating a mixture of stochastic frontier regression models via the em algorithm: a multiproduct cost function application', \emph{Empirical Economics}, 28(3), pp. 581-598. doi:10.1007/s001810200142.

Jondrow, J., Lovell, C.A.K., Materov, I.S. and Schmidt, P. (1982) 'On the estimation of technical inefficiency in the stochastic frontier production function model', \emph{Journal of Econometrics}, 19(2-3), pp. 233-238. doi:10.1016/0304-4076(82)90004-5.}
\author{Christopher F. Parmeter and David H. Bernstein}
\note{Standard errors require optHessian set to TRUE}
\seealso{\code{\link{lcsfm_homogeneity}}, \code{\link{zsfm}}, \code{\link{sfm}}, \code{\link{psfm}}, \code{\link{ttsfm}}, \code{\link{sfareg-methods}}}
\value{
An object of class \code{"sfareg"} containing the following components:

\item{out}{A matrix with parameter estimates, standard errors, and t-values. One row per parameter, named \code{<par>_class<j>} for the class blocks and \code{logit_<var>_class<j>} for the class-probability coefficients.}
\item{opt}{A list containing the optimization results from the final optimization procedure.}
\item{total_time}{The total computation time for model estimation.}
\item{start_v}{The starting values used in the optimization.}
\item{model_name}{The name of the model estimated (\code{"LCM"} or \code{"LCM_Z"}).}
\item{formula}{The formula used in the model specification.}
\item{jlms}{Posterior-weighted Jondrow et al. (1982) inefficiency predictions.}
\item{post.prob}{An \eqn{n \times J} matrix of posterior class probabilities, whose rows sum to one.}
\item{jlms_class}{An \eqn{n \times J} matrix of class-conditional JLMS predictions.}
\item{class}{Integer vector giving each observation's modal class.}
\item{class_prob}{Average prior class probabilities, one per class.}
\item{n_class}{The number of classes fitted.}
\item{coefficients}{A vector of estimated parameters.}
\item{std.errors}{A vector of standard errors for the estimated parameters (NA if \code{optHessian = FALSE}).}
\item{t.values}{A vector of t-values for the estimated parameters (NA if \code{optHessian = FALSE}).}
\item{call}{The matched call.}
}
\examples{
\donttest{
library(sfa)

## Two coexisting technologies, class probabilities constant.
set.seed(42)
n   <- 400
x1  <- rnorm(n); x2 <- rnorm(n)
cls <- rbinom(n, 1, 0.4) + 1L
b   <- rbind(c(1, 0.5, 0.5), c(3, 1.0, 0.2))
sv  <- c(0.2, 0.2); su <- c(0.5, 1.0)
y   <- b[cls, 1] + b[cls, 2] * x1 + b[cls, 3] * x2 +
       rnorm(n, 0, sv[cls]) - abs(rnorm(n, 0, su[cls]))

fit <- lcsfm(y ~ x1 + x2, model_name = "LCM",
             data = data.frame(y, x1, x2), n_class = 2)
fit$class_prob        # average class shares
head(fit$post.prob)   # posterior class membership, one row per firm
}
}
