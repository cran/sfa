\name{LogLik}
\alias{L_hNV}
\alias{L_exp}
\alias{L_trunc}
\alias{L_trunc_mufest}
\title{
The negative log likelihood function of the SFA
}
\description{
The negative log likelihood function is used for estimating the parameters. It varies depending on 
the distribution of the inefficence term u.
\code{L_hNV} is used by halfnormal distribution of u.
\code{L_exp} is used by exponential distribution of u.
\code{L_trunc} is used by truncated normal distribution of u.
\code{L_trunc_mufest} is used by truncated normal distribution of u and constant mu.
}
\usage{
L_hNV(p, y = y, X = X, sc = sc)
L_exp(p, y = y, X = X, sc = sc)
L_trunc(p, y = y, X = X, sc = sc)
L_trunc_mufest(p, mu = mu, y = y, X = X, sc = sc)
}
\arguments{
  \item{p}{
vector with the parameters to estimate
}
  \item{y}{
response
}
  \item{X}{
design matrix of the covariables
}
  \item{sc}{
specifies the form of the frontier model (-1 = cost, 1 = production)
}
  \item{mu}{
if known, the parameter mu
}
}
\value{
returns the value of the log likelihood function
}
\seealso{
\code{\link{sfa}}
}
\keyword{function}