\name{Indian}
\alias{Indian}
\docType{data}
\title{Indian}
\description{Panel data on 14 paddy farmers from Aurepalle, India, collected over ten years (1975-76 to 1984-85). Includes farmer characteristics (age, schooling) and production variables (output, land, labor, bullocks, input costs).}
\usage{data("Indian")}
\format{
  A data frame with 273 observations (an unbalanced panel of 34 farmers over 10 years) on the following 10 variables.
  \describe{
    \item{\code{id}}{a numeric vector containing a unique identifier for each farmer}
    \item{\code{yr}}{a numeric vector containing the year of the observation}
    \item{\code{age}}{a numeric vector containing the age of the primary decision maker}
    \item{\code{school}}{a numeric vector containing the number of years of schooling of the primary decision maker}
    \item{\code{yvar}}{a numeric vector containing the natural logarithm of the total value of output (rupees)}
    \item{\code{Lland}}{a numeric vector containing the natural logarithm of the total area of land operated (ha)}
    \item{\code{PIland}}{a numeric vector containing the proportion of land that is irrigated}
    \item{\code{Llabor}}{a numeric vector containing the natural logarithm of the total number of hours of hired and family labour used}
    \item{\code{Lbull}}{a numeric vector containing the natural logarithm of the number of hours of bullock labour used}
    \item{\code{Lcost}}{a numeric vector containing the natural logarithm of the value of inputs including fertilizer, manure, pesticides, machinery, etc.}
  }
}
\source{Battese, G.E. and Coelli, T.J. (1995) 'A model for technical inefficiency effects in a stochastic frontier production function for panel data', Empirical Economics, 20(2), pp. 325-332. doi:10.1007/BF01205442.}
\references{Battese, G.E. and Coelli, T.J. (1992) 'Frontier production functions, technical efficiency and panel data: With application to paddy farmers in India', Journal of Productivity Analysis, 3(1-2), pp. 153-169. doi:10.1007/BF00158774.}
\examples{data(Indian)}
\keyword{datasets}