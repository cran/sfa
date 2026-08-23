\name{PL80_MVTN}
\alias{PL80_MVTN}
\title{Pitt and Lee (1981) Model III: the multivariate truncated normal panel likelihood}
\description{Notes on \code{psfm(model_name = "PL80_MVTN")}, which fits the likelihood derived in Appendix 2 of Pitt and Lee (1981) but never used there.}
\details{
The model is
\deqn{y_{it} = x_{it}'\beta + u_{it} + v_{it}, \qquad u_{it} \le 0,}
with the firm's inefficiency vector \eqn{u_i = (u_{i1},\dots,u_{iT})'} drawn from a \eqn{T}-variate normal \eqn{N(0,\Sigma)} \bold{truncated to the negative orthant}, and \eqn{v_{it}} iid \eqn{N(0,\sigma_v^2)} independent of \eqn{u}. Unlike \code{model_name = "PL80"}, which holds inefficiency fixed over time, here it varies across periods \emph{and} is correlated within a firm -- \eqn{\Sigma} carries that dependence.

\bold{Why it was never used.} Pitt and Lee derived this likelihood and then set it aside, writing that it \dQuote{is difficult to evaluate since the quantities \eqn{P_0} and \eqn{P(y_i - x_i\beta)} involve \eqn{T}-dimensional numerical integrals}, and estimating Model III by Zellner seemingly-unrelated regression instead. Those quantities are orthant probabilities of a multivariate normal. \code{mnormt::sadmvn()} evaluates one in about 3 milliseconds at \eqn{T = 6}, so a likelihood evaluation costs roughly \eqn{N+1} of them -- about 0.3 s at \eqn{N = 100}, and about 12 s for a whole fit at \eqn{N = 80}, \eqn{T = 4}. What was intractable in 1981 is merely slow now.

Writing \eqn{Q^{-1} = \Sigma^{-1} + I/\sigma_v^2} and \eqn{\mu_i = Q\varepsilon_i/\sigma_v^2}, the per-firm log density is
\deqn{\log f(\varepsilon_i) = -\tfrac{T}{2}\log 2\pi - T\log\sigma_v - \tfrac{1}{2}\log|\Sigma| + \tfrac{1}{2}\log|Q| - \log P_0 - \tfrac{1}{2}\left(\tfrac{\varepsilon_i'\varepsilon_i}{\sigma_v^2} - \tfrac{\varepsilon_i'Q\varepsilon_i}{\sigma_v^4}\right) + \log P_i,}
where \eqn{P_0 = \Pr(w \le 0)} for \eqn{w \sim N(0,\Sigma)} is the truncation constant and \eqn{P_i = \Pr(w \le 0)} for \eqn{w \sim N(\mu_i, Q)}.

\bold{\eqn{\Sigma} is equicorrelated, not unrestricted.} It is parameterized as \eqn{\Sigma = \sigma_u^2[(1-\rho)I + \rho \mathbf{1}\mathbf{1}']}, costing two parameters. An unrestricted \eqn{\Sigma} costs \eqn{T(T+1)/2} -- 21 at \eqn{T=6}, 55 at \eqn{T=10} -- on top of \eqn{\beta} and \eqn{\sigma_v}, every one identified only through orthant probabilities. The equicorrelated form captures what the general \eqn{\Sigma} was introduced for: dependence of a firm's inefficiency across periods. \eqn{\rho = 0} gives inefficiency independent over time; \eqn{\rho \to 1} approaches the time-invariant \code{"PL80"} case. Because the form is equicorrelated every matrix quantity above is closed form (Sherman-Morrison), so only the orthant probabilities are numerical.

\bold{Requirements and limits.} A \bold{balanced} panel with \eqn{T \ge 2}: \eqn{\Sigma} is a single \eqn{T \times T} matrix shared by every firm, and with \eqn{T = 1} there is no cross-period dependence for it to describe. An unbalanced panel is an error pointing at \code{"PL80"}. Reported parameters are \code{sigv}, \code{sigu}, \code{rho} and the frontier coefficients.

\code{u_hat} is the posterior mean of \eqn{u} from the Gaussian part, floored at zero, \bold{not} \eqn{E[u_{it}\mid\varepsilon_i]}: the exact conditional mean of a truncated multivariate normal is another \eqn{T}-dimensional integral. It is exact where the truncation does not bind and an approximation where it does.
}
\references{
Pitt, M.M. and Lee, L.-F. (1981) 'The measurement and sources of technical inefficiency in the Indonesian weaving industry', \emph{Journal of Development Economics}, 9(1), pp. 43-64. See Appendix 2 for the likelihood.
}
\seealso{\code{\link{psfm}}, \code{\link{data_gen_p}} for the \code{y_pl_mvtn} column that tests it.}
\author{David H. Bernstein}
\examples{
\donttest{
library(sfa)
d <- data_gen_p(t = 4, N = 60, rand = 5, sig_u = 1, sig_v = 0.3, sig_r = 0.2,
                sig_h = 0.4, cons = 0.5, beta1 = 0.5, beta2 = 0.5)
f <- psfm(y_pl_mvtn ~ x1 + x2, model_name = "PL80_MVTN",
          data = as.data.frame(d), individual = "name")
f$out
}
}
