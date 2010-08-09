\name{methods.sfa}
\alias{methods.sfa}
\alias{coef.sfa}
\alias{print.sfa}
\alias{predict.sfa}
\alias{fitted.sfa}
\alias{logLik.sfa}
\alias{residuals.sfa}
\alias{summary.sfa}
\alias{print.summary.sfa}
\alias{eff.sfa}
\title{
Methods for displaying information about stochastic frontier analysis models
}
\description{
\code{coef.sfa} is used to display the fitted coefficients.
\code{print.sfa} is used to display some information about the fitted SFA.
\code{predict.sfa} is used to predict (new) data with the fitted SFA model.
\code{fitted.sfa} is used to predict the original data with the fitted SFA model.
\code{logLik.sfa} is used to display the value of the log likelihood function.
\code{residuals.sfa} is used to return the residuals of the fitted SFA model.
\code{summary.sfa} is used to calculate the summary result of the SFA.
\code{print.summary.sfa} is used display the summary result of the SFA.
\code{eff.sfa} is used to return the efficiencies of the SFA.
}
\usage{
\method{coef}{sfa}(object, \dots)
\method{print}{sfa}(x, \dots)
\method{predict}{sfa}(object, newdata = NULL, intercept = NULL, \dots)
\method{fitted}{sfa}(object, \dots)
\method{logLik}{sfa}(object, \dots)
\method{residuals}{sfa}(object, \dots)
\method{summary}{sfa}(object, \dots)
\method{print.summary}{sfa} (x, \dots)
\method{eff}{sfa} (object, \dots)
}
\arguments{
  \item{x}{
an object of class sfa
}
  \item{object}{
an object of class sfa
}
  \item{newdata}{
a data frame. If \code{newdata = NULL} then original data will be used.
}
  \item{intercept}{
boolean or \code{NULL}. If \code{intercept = NULL} then the function uses the same
intercept options as specified in sfa.
}
  \item{\dots}{
ignored.
}
}
\examples{
set.seed(225)
daten <- dgp(n = 100, b = c(1, 2), sc = -1)
test <- sfa(y ~ x, data = daten)
coef(test)
print(test)
predict(test)
fitted(test)
logLik(test)
residuals(test)
summary(test)
eff(test)
}
\keyword{methods}