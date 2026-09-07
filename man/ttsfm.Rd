\name{ttsfm}
\alias{ttsfm}
\title{Two-Tier Stochastic Frontier Model}
\description{Code to estimate the two-tier stochastic frontier model (2TSF), which decomposes a symmetric error term into two one-sided components (e.g. buyer/seller informational rents) in addition to the usual two-sided noise term.}
\usage{ttsfm(formula, model_name = c("TTNE", "TTHN", "TTNLS"),
data, z_link = c("sd", "var"), maxit.bobyqa = 80000,
maxit.psoptim = 1000, maxit.optim = 1000,
REPORT = 1, trace = 0, pgtol = 0, start_val = FALSE, PSopt = FALSE,
optHessian = TRUE, inefdec = TRUE, upper = NA,
Method = "L-BFGS-B", logit = TRUE, verbose = FALSE, rand.psoptim = NULL)}
\arguments{
\item{formula}{a symbolic description for the model to be estimated. Use \code{y ~ x | zu | zw} to allow determinants \code{zu} of the first one-sided component and \code{zw} of the second.}
\item{model_name}{model name for the estimation: \code{"TTNE"} for the Normal-Exponential-Exponential two-tier model, \code{"TTHN"} for the Normal-Half Normal-Half Normal two-tier model, or \code{"TTNLS"} for a nonlinear-least-squares two-tier estimator with no distributional assumption on \code{u}/\code{w} beyond their means. Matching ignores case, so \code{"nhn"}, \code{"NHN"} and \code{"Nhn"} are the same choice.}
\item{data}{A data frame containing the variables named in \code{formula}.}
\item{z_link}{Scale on which the variance-determinant linear predictors are placed. \code{"sd"} (the default, and \code{ttsfm}'s historical convention) sets \eqn{\sigma_u = \exp(z_u'\delta)} and likewise for \eqn{\sigma_w}; \code{"var"} sets \eqn{\sigma_u = \sqrt{\exp(z_u'\delta)}}, matching \code{\link{psfm}} and the competing packages. The two are reparameterizations of the same model -- identical maximised log-likelihood and identical frontier coefficients -- but \eqn{\delta} under the SD link is exactly \bold{half} \eqn{\delta} under the variance link, so a \eqn{\delta} from \code{ttsfm} is not comparable with one from \code{psfm} or \pkg{sfaR} unless the same convention is used. See \code{\link{sfm}}, which gained the same argument.}
\item{maxit.bobyqa}{Maximum number of iterations for the bobyqa optimization routine}
\item{maxit.psoptim}{Maximum number of iterations for the psoptim optimization routine}
\item{maxit.optim}{Maximum number of iterations for the optim optimization routine}
\item{REPORT}{reporting parameter}
\item{trace}{Integer. Tracing level passed through to the optimizer; larger values print more.}
\item{pgtol}{Numeric. Projected-gradient tolerance passed to \code{optim()}'s \code{"L-BFGS-B"} method. \code{0} uses the optimizer's own default.}
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
\details{The two-tier stochastic frontier model (Polachek and Yoon, 1987) models the composed error as \eqn{v + u - w}, where \code{v} is a symmetric normal noise term and \code{u}, \code{w} are independent one-sided components. \code{"TTNE"} assumes \code{u} and \code{w} are exponential; \code{"TTHN"} assumes they are half-normal; \code{"TTNLS"} makes no distributional assumption and instead estimates \eqn{\sigma_u}/\eqn{\sigma_w} by nonlinear least squares via the scaling property. Determinants of \code{u} and \code{w} can be included via the second and third parts of a two-pipe formula, which parameterize \eqn{\sigma_u} and \eqn{\sigma_w} through a log link.

Like the other composed-error models in this package, all three two-tier likelihoods can have a genuine (non-error) boundary mode where a variance component is estimated near its numerical floor in some samples; spot-check fits before serious use. \code{"TTHN"}'s post-estimation information-deficiency metrics (the \code{metrics} component of the returned object) involve a large number of bivariate-normal-CDF evaluations and are more experimental than the parameter estimates themselves -- \code{NULL} if the underlying calculation fails. For \code{"TTNLS"}, the reported standard error for the \eqn{\sigma_v}-position parameter (present only for layout consistency with \code{"TTNE"}/\code{"TTHN"}, unused in the NLS objective) is always \code{NA} by construction.

Note on default iteration counts: like \code{sfm()}/\code{zsfm()}, \code{ttsfm()}'s likelihoods are closed-form, so its optimization defaults are set much higher than \code{psfm()}'s (which uses simulated maximum likelihood for several models and defaults to far fewer iterations for speed). See \code{\link{psfm}}'s Details for more.}
\references{Polachek, S. and Yoon, B.J. (1987). A Two-Tiered Earnings Frontier Estimation of Employer and Employee Information in the Labor Market. \emph{The Review of Economics and Statistics}.}
\author{David H. Bernstein and Christopher F. Parmeter}
\note{Standard errors require optHessian set to TRUE}
\value{
An object of class \code{"sfareg"} containing the following components:

\item{out}{A matrix with parameter estimates, standard errors, and t-values.}
\item{opt}{A list containing the optimization results from the final optimization procedure.}
\item{total_time}{The total computation time for model estimation.}
\item{start_v}{The starting values used in the optimization.}
\item{model_name}{The name of the two-tier stochastic frontier model estimated (TTNE, TTHN, or TTNLS).}
\item{formula}{The formula used in the model specification.}
\item{coefficients}{A vector of estimated parameters.}
\item{std.errors}{A vector of standard errors for the estimated parameters (NA if \code{optHessian = FALSE}).}
\item{t.values}{A vector of t-values for the estimated parameters (NA if \code{optHessian = FALSE}).}
\item{metrics}{A list of information-deficiency measures (M1, M2, M5, M6, M7, M10; see Polachek and Yoon) evaluated at the final parameter estimates. \code{NULL} if the underlying calculation failed (see Details for the \code{"TTHN"} caveat).}
\item{call}{The matched call.}
}
\examples{
\donttest{
library(sfa)

data <- data_gen_cs(N = 200, rand = 1, sig_u = 0.3, sig_v = 0.2,
                     cons = 2, beta1 = 0.5, beta2 = 0.3, a = 5, mu = 0)

eq <- y_pcs ~ x1 + x2

ttsfm(formula    = eq,
      model_name = "TTNE",
      data       = data)
}
}
