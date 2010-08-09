\name{dgp}
\alias{dgp}
\title{
Sample data generating process
}
\description{
Sample data generating process
}
\usage{
dgp(n, b, intercept = TRUE, sc = -1)
}
\arguments{
  \item{n}{
sample size
}
  \item{b}{
parameter vector
}
  \item{intercept}{
logical, TRUE includes intercept
}
  \item{sc}{
form of the frontier model, -1 for cost frontier model, 1 for production frontier model
}
}
\value{
list
}
\seealso{
\code{\link{sfa}}, \code{\link{rnorm}}, \code{\link{runif}}, \code{\link{abs}}
}
\keyword{data}