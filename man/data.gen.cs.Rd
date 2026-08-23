\name{data_gen_cs}
\alias{data_gen_cs}
\title{Generate Cross-Sectional Data for Stochastic Frontier Analysis}
\description{
  \code{data_gen_cs} generates simulated cross-sectional data based on the stochastic frontier model, allowing for different distributional assumptions for the one-sided technical inefficiency error term (\eqn{u}) and the two-sided idiosyncratic error term (\eqn{v}). The model has the general form:
  \eqn{Y = \beta_0 + \beta_1 x_1 + \beta_2 x_2 + v - u}
  where \eqn{u \geq 0} and represents inefficiency.  All variants are produced so that the user can select those that they want.}
\usage{data_gen_cs(N, rand, sig_u, sig_v, cons, beta1, beta2, a, mu, sig_w = sig_u,
            shape_g = 2, m_nak = 1, mu_ln = -0.5, k_w = 1.5)}
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
  \item{sig_w}{The standard deviation/scale parameter (\eqn{\sigma_w}) for the second one-sided error component \eqn{w}, used to generate two-tier model columns (see \code{\link{ttsfm}}). Defaults to \code{sig_u}.}
  \item{shape_g}{Shape of the gamma inefficiency draw behind \code{y_pcs_g}, which targets \code{sfm(model_name = "NG")}. The scale is set to \code{sig_u/shape_g}, so \eqn{E[u] = \sigma_u} regardless of the shape chosen. \code{sfm()}'s NG likelihood reports the shape as \code{mu} and the scale as \code{sigu}, so the values to recover are \code{shape_g} and \code{sig_u/shape_g}.}
  \item{m_nak}{Nakagami shape \eqn{m} behind \code{y_pcs_nak}, which targets \code{sfm(model_name = "NNAK")}. The spread is \eqn{\Omega = \sigma_u^2}, so \eqn{\sigma_u} is the root-mean-square of \eqn{u}. The default \code{1} is the Rayleigh case; \code{0.5} would collapse the Nakagami onto the half-normal and so duplicate \code{y_pcs}.}
  \item{mu_ln}{Log-scale mean of the lognormal inefficiency draw behind \code{y_pcs_ln}, which targets \code{sfm(model_name = "NLN")}. The log-scale standard deviation is \code{sig_u}. The default \code{-0.5} with \code{sig_u = 1} gives \eqn{E[u] = 1}.}
  \item{k_w}{Weibull shape behind \code{y_pcs_wb}, which targets \code{sfm(model_name = "NW")}. The scale is \code{sig_u}.}
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
\item \code{y_pcs_st}: the column that actually matches \code{\link{sfm}}'s \code{THT}. \code{y_pcs_t} above draws its two components as two \emph{independent} \code{rt()} variates, which shares the degrees of freedom but not the mixing variable, so the composed error is not skew-t and \code{THT} cannot recover it. \code{y_pcs_st} uses the single common \eqn{\lambda \sim \Gamma(a/2, a/2)} of Tancredi (2002), i.e. \eqn{y = h + (\epsilon - z)/\sqrt{\lambda}}. \code{y_pcs_t} is retained unchanged only for backward compatibility.
\item \code{y_pcs_thn}: Student's t--half normal, matching \code{\link{sfm}}'s \code{tHN}: \eqn{v \sim T(\text{df}=a) \cdot \sigma_v} and \eqn{u \sim |N(0,\sigma_u^2)|}, drawn \emph{independently}. Heavy-tailed noise with a conventional one-sided term -- distinct from both \code{y_pcs_t} and \code{y_pcs_st}, in which the inefficiency is heavy-tailed too.
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
  \item{u, uz, u_t, u_c, u_e, u_u, u_tn, u_st, u_thn}{The simulated one-sided error terms under different distributions.}
  \item{u_r}{Rayleigh one-sided error term, for \code{y_pcs_r}. Drawn on the convention \code{sfm}'s \code{"NR"} uses, \eqn{E[u^2] = \sigma_u^2}, so the Rayleigh scale is \eqn{\sigma_u/\sqrt{2}}, \eqn{E[u] = \sigma_u\sqrt{\pi}/2} and \eqn{Var(u) = (1-\pi/4)\sigma_u^2}.}
  \item{v, v_t, v_c, v_st, v_thn}{The simulated two-sided error terms under different distributions.}
\item{lam_st}{The common Gamma(a/2, a/2) scale-mixing variable behind \code{y_pcs_st}; both of that column's components are divided by its square root, which is what makes the composed error skew-t.}
  \item{y_pcs, y_pcs_t, y_pcs_st, y_pcs_thn, y_pcs_e, y_pcs_ez, y_pcs_c, y_pcs_u, y_pcs_z, y_pcs_w, y_pcs_tn}{The dependent variable \eqn{Y} under the corresponding SFA model distributions.}
  \item{y_pcs_r}{Normal-Rayleigh dependent variable, for \code{\link{sfm}}'s \code{"NR"}. Use this column, not \code{y_pcs}: the normal-Rayleigh and normal-half-normal composed errors are different families, not reparameterizations of one another, and \code{"NR"} cannot recover \code{y_pcs}'s parameters. True values are \eqn{(\sigma_v, \sigma_u) = } \code{(sig_v, sig_u)}.}
  \item{z}{The auxiliary variable used for heteroskedasticity in \code{y_pcs_z}, \code{y_pcs_ez}, \code{y_tthn_z}, and \code{y_zisf_z}.}
  \item{uz_e}{Exponential one-sided error term with heteroskedastic scale, for \code{y_pcs_ez}.}
  \item{w_tt, w_tt_hn, wz_hn}{The second one-sided error component (\eqn{w}) for the two-tier columns below, under exponential, homoskedastic half-normal, and heteroskedastic half-normal distributions respectively.}
  \item{zp}{A second auxiliary variable used for heteroskedasticity of \eqn{\sigma_w} in \code{y_tthn_z}.}
  \item{y_ttne}{Homoskedastic two-tier dependent variable (\eqn{v + w - u}, both one-sided components exponential), for \code{\link{ttsfm}}'s \code{"TTNE"}/\code{"TTNLS"}.}
  \item{y_tthn}{Homoskedastic two-tier dependent variable (both one-sided components half-normal), for \code{\link{ttsfm}}'s \code{"TTHN"} with no pipes.}
  \item{y_tthn_z}{Heteroskedastic two-tier dependent variable (\eqn{\sigma_u} a function of \code{z}, \eqn{\sigma_w} a function of \code{zp}), for \code{\link{ttsfm}}'s \code{"TTHN"} with formula \code{~x1+x2|z|zp}.}
  \item{eff_ind, eff_ind_z}{Indicator (1 = efficient, \eqn{u} forced to 0) used to build the zero-inefficiency columns below.}
  \item{prob_z_true}{The true heteroskedastic "efficient regime" probability used to draw \code{eff_ind_z}.}
  \item{y_zisf}{Zero-inefficiency dependent variable with constant efficient-regime probability, for \code{\link{zsfm}}'s \code{"ZISF"}.}
  \item{y_zisf_z}{Zero-inefficiency dependent variable with \code{z}-dependent efficient-regime probability, for \code{\link{zsfm}}'s \code{"ZISF_Z"}.}
  \item{u_g, y_pcs_g}{Gamma one-sided error term and the corresponding dependent variable, for \code{\link{sfm}}'s \code{"NG"}. Shape \code{shape_g}, scale \code{sig_u/shape_g}, so \eqn{E[u] = \sigma_u} whatever shape is chosen.}
  \item{u_nak, y_pcs_nak}{Nakagami one-sided error term and dependent variable, for \code{\link{sfm}}'s \code{"NNAK"}. Shape \code{m_nak}, spread \eqn{\Omega = \sigma_u^2}. The default \code{m_nak = 1} is the Rayleigh case, so this column coincides in distribution with \code{y_pcs_r} at the default.}
  \item{u_ge, y_pcs_ge}{Generalized-exponential one-sided error term and dependent variable, for \code{\link{sfm}}'s \code{"NGE"}, drawn as \eqn{u = -\sigma_u\log(1 - \sqrt{U})} with \eqn{U \sim \mathrm{Unif}(0,1)}.}
  \item{u_ln, y_pcs_ln}{Lognormal one-sided error term and dependent variable, for \code{\link{sfm}}'s \code{"NLN"}. Log-scale mean \code{mu_ln} and log-scale standard deviation \code{sig_u}.}
  \item{u_w, y_pcs_wb}{Weibull one-sided error term and dependent variable, for \code{\link{sfm}}'s \code{"NW"}. Shape \code{k_w}, scale \code{sig_u}.}
  \item{con}{A constant column set to 1, potentially for use in estimation.}
}
\seealso{
  \code{\link{data_gen_p}} for the panel generator, and \code{\link{sfm}}, \code{\link{zsfm}} and \code{\link{ttsfm}} for the estimators these columns are built for.

  \code{\link[stats]{rnorm}}, \code{\link[stats]{runif}}, \code{\link[stats]{rt}}, \code{\link[stats]{rexp}}, \code{\link[stats]{rcauchy}}, \code{\link[truncnorm]{rtruncnorm}} (if available).
}
\author{David H. Bernstein}
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