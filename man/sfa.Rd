\name{sfa}
\alias{sfa}
\alias{sfa.fit}
\title{
Fitting stochastic frontier analysis models
}
\description{
\code{sfa} is used to fit stochastic frontier analysis models.
}
\usage{
sfa(formula, data = NULL, intercept = TRUE, fun = "hnormal", pars = NULL, par_mu = NULL, form = "cost", method = "BFGS", ...)
}
\arguments{
  \item{formula}{
an object of class \code{\link{formula}} (or one that can be coerced to that
class): a symbolic description of the model to be fitted.
}
  \item{data}{
a data frame.
}
  \item{intercept}{
logical. If true the model includes intercept.
}
  \item{fun}{
specifies the distribution for the inefficency term u as half-normal ("hnormal"), exponential ("exp"), or truncated-normal ("tnormal").
}
  \item{pars}{
initial values for the parameters to be estimated.
}
  \item{par_mu}{
value for mu in the normal-/truncated-normal case. If mu is known.
}
  \item{form}{
specifies the form of the frontier model as "cost" or "production".
}
  \item{method}{
the method to be used. See \code{\link{optim}} for more details.
}
  \item{\dots}{
ignored.
}
}
\value{
\code{sfa} returns an object of class \code{sfa}:
  \item{y}{response}
  \item{x}{covariables}
  \item{X}{design matrix}
  \item{coef}{coefficients}
  \item{sigmau2}{sigmau2}
  \item{sigmav2}{sigmav2}
  \item{mu}{mu of the truncated-normal distribution (Only if fun = tnormal)}
  \item{par_mu}{NULL if mu is not estimated}
  \item{logLik}{value of the log likelihood function}
  \item{maxlik}{log likelihood function}
  \item{fun}{distribution of the inefficiency term u}
  \item{sc}{specifies the form of the frontier model (-1 = cost, 1 = production)}
  \item{hess}{a symmetric matrix giving an estimate of the Hessian at the solution found (See \code{\link{optim}})}
  \item{ols}{the linear model for the LR-test}
}
\examples{
set.seed(225)
daten <- dgp(n = 100, b = c(1, 2), sc = -1)
test <- sfa(y ~ x, data = daten)
}
\keyword{models}