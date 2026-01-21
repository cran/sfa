\name{data_gen_cs}
\alias{data_gen_cs}
\title{Generate Cross-Sectional Data for Stochastic Frontier Analysis}
\description{
  \code{data_gen_cs} generates simulated cross-sectional data based on the stochastic frontier model, allowing for different distributional assumptions for the one-sided technical inefficiency error term (\eqn{u}) and the two-sided idiosyncratic error term (\eqn{v}). The model has the general form:
  \eqn{Y = \beta_0 + \beta_1 x_1 + \beta_2 x_2 + v - u}
  where \eqn{u \geq 0} and represents inefficiency.  All variants are produced so that the user can select those that they want.}
\usage{data_gen_cs(N, rand, sig_u, sig_v, cons, beta1, beta2, a, mu)}
\arguments{
  \item{N}{A single integer specifying the number of observations (cross-sectional units).}
  \item{rand}{A single integer to set the seed for the random number generator, ensuring reproducibility.}
  \item{sig_u}{The standard deviation parameter (\eqn{\sigma_u}) for the base distribution of the one-sided error term \eqn{u}.}
  \item{sig_v}{The standard deviation parameter (\eqn{\sigma_v}) for the base distribution of the two-sided error term \eqn{v}.}
  \item{cons}{The value of the constant term (intercept) in the model.}
  \item{beta1}{The coefficient for the \eqn{x_1} variable.}
  \item{beta2}{The coefficient for the \eqn{x_2} variable.}
  \item{a}{The degrees of freedom parameter for the t half-t distribution (\code{u_t} and \code{v_t}, respectively). Requires the \code{rt} function.}
  \item{mu}{The mean parameter (\eqn{\mu}) for the normal truncated normal distribution (\code{u_tn}). Requires the \code{rtruncnorm} function.}
}
\details{
  The function simulates two explanatory variables, \eqn{x_1} and \eqn{x_2}, as transformations of uniform random variables.
  
  The function generates several different frontier models by combining various distributions for \eqn{u} and \eqn{v}:
  \itemize{
    \item **\eqn{u} Distributions (Inefficiency):** Half-Normal (HN), Truncated Normal (TN), Half-T (HT), Half-Cauchy (HC), Exponential (E), Half-Uniform (HU).
    \item **\eqn{v} Distributions (Idiosyncratic):** Normal (N), t, Cauchy (C).
  }
  
  **Specific Model Outputs (y_pcs variants):**
\itemize{
\item \code{y_pcs}: Normal-Half Normal (N-HN): \eqn{v \sim N(0, \sigma_v^2)}, \eqn{u \sim |N(0, \sigma_u^2)|}.
\item \code{y_pcs_z}: N-HN with Heteroskedastic \eqn{\sigma_u}: \eqn{\sigma_{u,i} = \exp(0.9 + 0.6 Z_i)}, where \eqn{Z} is a uniform variable.
\item \code{y_pcs_t}: T-Half T (T-HT): \eqn{v \sim T(\text{df}=a) \cdot \sigma_v}, \eqn{u \sim |T(\text{df}=a)| \cdot \sigma_u}.
\item \code{y_pcs_tn}: Normal-Truncated Normal (N-TN): \eqn{v \sim N(0, \sigma_v^2)}, \eqn{u \sim TN(\mu, \sigma_u^2)} on \eqn{[0, \infty)}.
\item \code{y_pcs_e}: Normal-Exponential (N-E): \eqn{v \sim N(0, \sigma_v^2)}, \eqn{u \sim Exp(\phi)}, where \eqn{\phi = 1/\sigma_u}.
\item \code{y_pcs_c}: Cauchy-Half Cauchy (C-HC): \eqn{v \sim Cauchy(0, \sigma_v)}, \eqn{u \sim |Cauchy(0, \sigma_u)|}.
\item \code{y_pcs_u}: Normal-Half Uniform (N-HU): \eqn{v \sim N(0, \sigma_v^2)}, \eqn{u \sim U(0, \sigma_u)}.
\item \code{y_pcs_w}: Normal + Cauchy - Half Normal: \eqn{v \sim N(0, \sigma_v^2) + Cauchy(0, \sigma_v)}, \eqn{u \sim |N(0, \sigma_u^2)|}. This introduces a composite \eqn{v} term.
  }
  
  **Note:** The \code{rtruncnorm} function is required for \code{y_pcs_tn} and loads with the package.  In isolation it could be loaded by using \code{library(truncnorm)}.
}
\value{
  A data frame containing \eqn{N} observations with the following columns:
  \item{name}{Individual identifier (simply \eqn{1} to \eqn{N}).}
  \item{cons}{The constant term value.}
  \item{x1}{Simulated explanatory variable \eqn{x_1}.}
  \item{x2}{Simulated explanatory variable \eqn{x_2}.}
  \item{u, uz, u_t, u_c, u_e, u_u, u_tn}{The simulated one-sided error terms under different distributions.}
  \item{v, v_t, v_c}{The simulated two-sided error terms under different distributions.}
  \item{y_pcs, y_pcs_t, y_pcs_e, y_pcs_c, y_pcs_u, y_pcs_z, y_pcs_w, y_pcs_tn}{The dependent variable \eqn{Y} under the corresponding SFA model distributions.}
  \item{z}{The auxiliary variable used for heteroskedasticity in \code{y_pcs_z}.}
  \item{con}{A constant column set to 1, potentially for use in estimation.}
}
\seealso{
  \code{\link[stats]{rnorm}}, \code{\link[stats]{runif}}, \code{\link[stats]{rt}}, \code{\link[stats]{rexp}}, \code{\link[stats]{rcauchy}}, \code{\link[truncnorm]{rtruncnorm}} (if available).
}
\author{David Bernstein}
\examples{

# Generate 100 observations of SFA data
data_sfa <- data_gen_cs(
  N     = 100,
  rand  = 123,
  sig_u = 0.5,
  sig_v = 0.2,
  cons  = 5,
  beta1 = 1.5,
  beta2 = 2.0,
  a     = 5,   # degrees of freedom for T/Half-T
  mu    = 0.1  # mean for Truncated Normal
)

# Display the first few rows of the generated data
head(data_sfa)

# Example of a Normal-Half Normal SFA model data
summary(data_sfa$y_pcs)
plot(density(data_sfa$y_pcs))
}