\name{predict.sfareg}
\alias{predict.sfareg}
\alias{fitted.sfareg}
\alias{residuals.sfareg}
\title{Predictions, fitted values and residuals for stochastic frontier fits}
\description{Standard modelling generics for objects of class \code{"sfareg"}, as returned by \code{\link{sfm}}, \code{\link{psfm}}, \code{\link{zsfm}} and \code{\link{ttsfm}}.}
\usage{
\method{predict}{sfareg}(object, newdata = NULL,
        type = c("frontier", "response", "efficiency"), ...)
\method{fitted}{sfareg}(object, ...)
\method{residuals}{sfareg}(object, ...)
}
\arguments{
\item{object}{An object of class \code{"sfareg"}.}
\item{newdata}{Optional data frame containing the frontier regressors at which to predict. If omitted, the data the model was fitted to is used.}
\item{type}{\code{"frontier"} (default) returns the deterministic frontier \eqn{x'\beta}. \code{"response"} shifts it by the predicted inefficiency, giving an estimate of \eqn{E[y \mid x]} rather than of the frontier itself, with the sign following the production/cost convention the model was fitted under. \code{"efficiency"} returns the predicted technical efficiency \eqn{E[\exp(-u) \mid \varepsilon]}.}
\item{...}{Ignored.}
}
\details{Only the first part of a pipe formula (\code{y ~ x | z | zp}) describes the frontier; later parts parameterize variances and take no part in prediction.

Frontier coefficients are identified by matching the fitted coefficient names against the frontier design matrix, so these methods work for every model in the package without depending on where a given model happens to place its auxiliary parameters (\code{lambda}/\code{sigma} for \code{NHN}, \code{sigv}/\code{sigu} for \code{NE}, \code{sigmaSq}/\code{gamma} for \code{PL80}, and so on).

\code{type = "efficiency"} and \code{type = "response"} are available only for the estimation sample, since both require the composed residual and therefore the response; requesting them with \code{newdata} raises an error rather than returning something misleading. \code{"efficiency"} additionally requires a model that returns an efficiency prediction -- see \code{\link{sfm}} for which ones do.}
\value{A numeric vector, one element per observation.}
\seealso{\code{\link{sfm}}, \code{\link{psfm}}, \code{\link{sfareg-methods}}}
\author{David H. Bernstein}
\examples{
\donttest{
library(sfa)
d <- data_gen_cs(N = 200, rand = 1, sig_u = 0.8, sig_v = 0.3,
                 cons = 1, beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.5)
f <- sfm(y_pcs ~ x1 + x2, model_name = "NHN", data = d)

head(fitted(f))
head(residuals(f))
head(predict(f, type = "efficiency"))
head(predict(f, newdata = d[1:5, ]))
}
}
