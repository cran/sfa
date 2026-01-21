\name{FinnishElec}
\alias{FinnishElec}
\docType{data}
\title{
FinnishElec
}
\description{
Cross-sectional data on Finnish electricity distribution firms, including annual averages of expenditure and output measures over a four-year regulatory period.
}
\usage{data("FinnishElec")}
\format{
  A data frame with 89 observations on the following 6 variables.
  \describe{
    \item{\code{id}}{a character vector containing a unique identifier for each distribution firm}
    \item{\code{x}}{a numeric vector containing total expenditure (TOTEX*) (1000 Euros)}
    \item{\code{y1}}{a numeric vector containing weighted energy transmitted through the network (GWh of 0.4 kV equivalents)}
    \item{\code{y2}}{a numeric vector containing total length of the network (km)}
    \item{\code{y3}}{a numeric vector containing total number of customers connected to the network}
    \item{\code{z}}{a numeric vector containing the proportion of underground cables in the total network length.}
  }
}
\details{*TOTEX includes capital expenditure (CAPEX), controllable operational expenditure (OPEX), and estimated external cost of interruptions.}
\source{Kuosmanen, T. (2012). 'Stochastic semi-nonparametric frontier estimation of electricity distribution networks: Application of the StoNED method in the Finnish regulatory model.' Energy Economics, 34(6), pp. 2189-2199. doi:10.1016/j.eneco.2012.03.005}
\examples{
data(FinnishElec)
plot(FinnishElec)}
\keyword{datasets}