\name{panel89}
\alias{panel89}
\docType{data}
\title{Panel89}
\description{The dataset is a cross-section of U.S. commercial banks for 1989, extracted from the panel 
dataset used by Kumbhakar, Parmeter and Tsionas (2013) and based on the Federal Reserve Bank of Chicago's Reports of  Condition and Income. It contains detailed cost data with inputs and outputs defined under the  intermediation approach, and input prices constructed as expense-quantity ratios.}
\usage{data("panel89")}
\format{
  A data frame with 4,985 observations on the following 11 variables.
  \describe{
    \item{\code{y}}{a numeric vector containing the natural logarithm of total cost*}
    \item{\code{q1}}{a numeric vector containing the natural logarithm of installment loans}
    \item{\code{q2}}{a numeric vector containing the natural logarithm of real estate loans}
    \item{\code{q3}}{a numeric vector containing the natural logarithm of business loans}
    \item{\code{q4}}{a numeric vector containing the natural logarithm of federal funds sold and securities purchased}
    \item{\code{q5}}{a numeric vector containing the natural logarithm of other assets}
    \item{\code{w1}}{a numeric vector containing the natural logarithm of the price of labour*}
    \item{\code{w2}}{a numeric vector containing the natural logarithm of the price of capital*}
    \item{\code{w3}}{a numeric vector containing the natural logarithm of the price of purchased funds*}
    \item{\code{w4}}{a numeric vector containing the natural logarithm of the price of interest-bearing deposits in total transaction accounts*}
    \item{\code{z}}{a numeric vector containing the natural logarithm of total assets}
  }
}
\details{*The cost and input price variables are normalised by that of a fifth input: the price of interest-bearing deposits in total non-transaction accounts. Total cost is defined as the sum of total expenses for each input. Input prices are derived by dividing the total expense for each input by the corresponding input quantity.}

\source{Kumbhakar, S.C., Parmeter, C.F. and Tsionas, E.G. (2013) 'A zero inefficiency stochastic frontier model', Journal of Econometrics, 172(1), pp. 66-76. doi:10.1016/j.jeconom.2012.08.021.}
\references{Kumbhakar, S.C. and Tsionas, E.G. (2005) 'Measuring technical and allocative inefficiency in the translog cost system: a Bayesian approach', Journal of Econometrics, 126(2), pp. 355-384. doi:10.1016/j.jeconom.2004.05.006.}
\examples{data(panel89)
  plot(panel89)}
\keyword{datasets}