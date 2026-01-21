\name{data_gen_p}
\alias{data_gen_p}
\title{Generate Panel Data for Stochastic Frontier Analysis}
\description{\code{data_gen_p} generates simulated panel data for estimating various panel stochastic frontier models, including the Generalized True Random Effects (GTRE), True Random Effects (TRE), Pooled Cross-Section (PCS), and True Fixed Effects (TFE) models. The function returns the data as a \code{pdata.frame}.  All variants are produced so that the user can select those that they want.}
\usage{data_gen_p(t, N, rand, sig_u, sig_v, sig_r, sig_h, cons, tau = 0.5, mu = 0, beta1, beta2)}
\arguments{
\item{t}{The number of time periods.}
\item{N}{The number of individuals.}
\item{rand}{A seed for the random number generator to ensure reproducibility.}
\item{sig_u}{The standard deviation (\eqn{\sigma_u}) for the one-sided error component (\eqn{u_{it}}).}
\item{sig_v}{The standard deviation (\eqn{\sigma_v}) for the two-sided error component (\eqn{v_{it}}).}
\item{sig_r}{The standard deviation (\eqn{\sigma_r}) for the two-sided individual effect (\eqn{r_i}).}
\item{sig_h}{The standard deviation (\eqn{\sigma_h}) for the one-sided individual effect (\eqn{h_i}).}
\item{cons}{The constant term (\eqn{\beta_0}) for the frontier models.}
\item{tau}{The dependence parameter (\eqn{\tau}) used for the \code{y_tfe} (TFE) model formulation, default is 0.5.  See Chen, Schmidt, and Wang (2014, Journal of Econometrics).}
\item{mu}{The mean parameter (\eqn{\mu}) used for the Truncated-Normal (TN) component of the \code{y_fd} model with default set to 0. See Wang and Ho (2010, Journal of Econometrics).}
\item{beta1}{The coefficient for the \code{x1} variable (\eqn{\beta_1}).}
\item{beta2}{The coefficient for the \code{x2} variable (\eqn{\beta_2}).}}
\details{A \code{\link[plm:pdata.frame]{pdata.frame}} object with \eqn{N \times t} observations, containing the following columns:
\itemize{
 \item \code{name} Individual identifier.
 \item \code{year} Time period identifier.
 \item \code{cons} The constant term used in the data generation.
 \item \code{x1, x2} Explanatory variables generated from a log-uniform distribution.
 \item \code{x1_w, x2_w} Explanatory variables with dependence parameter \eqn{\tau} and linkage with \eqn{r_i}, used for the TFE model.
 \item \code{u, v, r, h} The generated error and individual effect components.
 \item \code{y_gtre, y_tre, y_pcs, y_tfe} Output variables for the Production Frontier models, including the constant.
 \item \code{y_gtre_nc, y_tre_nc, y_pcs_nc} Output variables for the Production Frontier models, \strong{excluding} the constant.
 \item \code{c_gtre, c_tre, c_pcs, c_tfe} Output variables for the Cost Frontier models, including the constant.
 \item \code{c_gtre_nc, c_tre_nc, c_pcs_nc} Output variables for the Cost Frontier models, \strong{excluding} the constant.
 \item \code{y_fd} Output variable for the first difference model (see Wang and Ho, 2010).
 \item \code{x_fd} Explanatory variable for the \code{y_fd} model.
 \item \code{u_fd_star, z_fd, r_fd, u_fd} Components used to generate \code{y_fd}.
 \item \code{u_gtre, z_gtre, y_gtre_z, y_tre_z} Variables for models with heteroskedastic inefficiency (\eqn{\sigma_{u,i} = \exp(0.9 + 0.6 Z_{i}))}.
 }
  The data is generated based on standard Stochastic Frontier Analysis (SFA) formulations, primarily for a **Production Frontier** where the one-sided error component \eqn{u_{it}} is subtracted:

 \itemize{
   \item \code{y_gtre}: GTRE model: \eqn{y_{it} = \beta_0 + \beta_1 x_{1,it} + \beta_2 x_{2,it} + r_i - h_i + v_{it} - u_{it}}
   \item \code{y_tre}: TRE model: \eqn{y_{it} = \beta_0 + \beta_1 x_{1,it} + \beta_2 x_{2,it} + r_i + v_{it} - u_{it}}
   \item \code{y_pcs}: PCS model: \eqn{y_{it} = \beta_0 + \beta_1 x_{1,it} + \beta_2 x_{2,it} + v_{it} - u_{it}}
   \item \code{y_tfe}: TFE model: \eqn{y_{it} = \beta_1 x_{1,it}^w + \beta_2 x_{2,it}^w + r_i + v_{it} - u_{it}}
   \item \code{y_gtre_z}: GTRE with Heteroskedastic \eqn{u_{it}}: \eqn{\sigma_{u,i} = \exp(0.9 + 0.6 Z_i)}.
 }

 For **Cost Frontier** models, the one-sided error component \eqn{u_{it}} is added (e.g., \code{c_gtre}).

 The error terms are generated as:
 \itemize{
   \item \eqn{r_i \sim N(0, \sigma_r^2)} (individual two-sided effect)
   \item \eqn{h_i \sim |N(0, \sigma_h^2)|} (individual one-sided effect)
   \item \eqn{v_{it} \sim N(0, \sigma_v^2)} (two-sided noise)
   \item \eqn{u_{it} \sim |N(0, \sigma_u^2)|} (one-sided inefficiency)
 }

The First-Difference estimation model (\code{y_fd}) uses a variation where \eqn{r_{i,fd} \sim U(0,1)} and \eqn{u_{it,fd}} is generated using a heteroskedastic truncated-normal structure, reflecting an alternative model type.
 }
\value{A \code{\link[plm:pdata.frame]{pdata.frame}} object containing \eqn{N \times t} observations suitable for Stochastic Frontier Analysis (SFA).}
\references{
Chen, Y., Schmidt, P., & Wang, H. (2014). Consistent estimation of the fixed effects stochastic frontier model. \emph{Journal of Econometrics}, 181(2), 65-76.

Filippini, M., & Greene, W. H. (2016). Persistent and transient productive inefficiency: a maximum simulated likelihood approach. \emph{Journal of Productivity Analysis}, 45, 187-196.

Wang, H., & Ho, C. M. (2010). Estimating fixed-effect panel stochastic frontier models by model transformation. \emph{Journal of Econometrics}, 157(2), 286-296.
}
\author{David Bernstein}
\seealso{data_gen_p, to see all the data generating processes}
\examples{library(sfa) 
# Generate a dataset 
data_trial <- data_gen_p(t=10, N=100, rand = 100, 
                       sig_u = 1,  sig_v = 0.3, 
                       sig_r = .2, sig_h = .4, 
                       cons = 0.5, tau = 0.5,
                       mu= 0.5, beta1 = 0.5,
                       beta2 = 0.5)
 # See the first few rows 
 head(data_trial)
}