\name{USUtilities}
\alias{USUtilities}
\docType{data}
\title{USUtilities}
\description{Panel data on U.S. investor-owned fossil fuel-fired steam electric utilities for the period 1986-1999. These data include measures of output, capital, labour and maintenance, and fuel.}
\usage{data("USUtilities")}
\format{
  A data frame with 972 observations (a balanced panel of observations on 81 utilities over 12 years) on the following 7 variables.
  \describe{
    \item{\code{firmID}}{a numeric vector containing a unique firm identifier}
    \item{\code{year}}{a numeric vector containing the year of the observation}
    \item{\code{q}}{a numeric vector containing net steam electric power generation (MWh)}
    \item{\code{K}}{a numeric vector containing capital stock, calculated using a method described by Christensen and Jorgenson (1970)}
    \item{\code{L}}{a numeric vector containing quantity of labor and maintenance, calculated as cost divided by price index}
    \item{\code{F}}{a numeric vector containing quantity of fuel used, calculated as fuel costs divided by fuel price index}
    \item{\code{trend}}{a numeric vector containing an annual time trend (1992=100)}
  }
}
\details{The dataset covers 72 investor-owned utilities after aggregating subsidiaries and excluding plants in states with partial deregulation plans. Data sources include the Energy Information Administration (EIA), Federal Energy Regulatory Commission (FERC), and Bureau of Labor Statistics (BLS). Output is net steam electric generation from fossil fuel-fired boilers.}
\source{Rungsuriyawiboon, S. and Stefanou, S.E. (2007). 'Dynamic Efficiency Estimation: An Application to U.S. Electric Utilities.' Journal of Business & Economic Statistics, 25(2), pp. 226-238. doi:10.1198/073500106000000288}
\references{Christensen, L.R. and Jorgenson, D.W. (1970). 'U.S. Real Product and Real Factor Input, 1928-1967.' Review of Income and Wealth, 16(1), pp. 19-50. doi: 10.1111/j.1475-4991.1970.tb00695.x}
\examples{data(USUtilities)}
\keyword{datasets}