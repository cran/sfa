\name{skewness_test}
\alias{skewness_test}
\title{Formal tests for wrong skewness in the OLS residuals}
\description{
Tests whether the OLS residuals of a stochastic frontier model are skewed in the direction the model implies. A production frontier has \eqn{e = v - u} with \eqn{u \ge 0}, so the composed error is \bold{negatively} skewed; residuals skewed the other way mean the maximum likelihood estimate of \eqn{\sigma_u} is exactly zero and every efficiency score is one.
}
\usage{skewness_test(object, test = c("agostino", "coelli"),
              alternative = c("auto", "less", "greater", "two.sided"))}
\arguments{
\item{object}{An object of class \code{"sfareg"} from \code{\link{sfm}}, or a numeric vector of residuals. From a fit, the \bold{OLS} residuals recorded at fit time are used -- not the composed residuals of the fitted model, which are not what either test is defined on.}
\item{test}{\code{"agostino"} (the default) is D'Agostino's (1970) skewness test, which is the test Schmidt and Lin (1984) appeal to. \code{"coelli"} is Coelli's (1995) \eqn{M3T}. See \sQuote{Details} for why the default is the former.}
\item{alternative}{\code{"auto"} (the default) runs the one-sided test the model implies: a production frontier predicts negative skew, so the null of zero skew is tested against \eqn{m_3 < 0}. The other three force a direction.}
}
\value{
An object of class \code{"htest"}, so it prints like any other R test. Alongside the usual \code{statistic}, \code{p.value} and \code{estimate} (the sample skewness \eqn{\sqrt{b_1}}) it carries \code{m3}, \code{nobs}, and \code{wrong_skew} -- the last being \code{TRUE} when the third central moment is non-negative, matching \code{sfm}'s own \code{$wrong_skew}.
}
\details{
\bold{What the tests are.} Both start from the third central moment of the OLS residuals. Coelli's \eqn{M3T = m_3/\sqrt{6m_2^3/n}} is asymptotically standard normal under the null of zero skewness. D'Agostino's transforms \eqn{\sqrt{b_1} = m_3/m_2^{3/2}} to approximate normality in finite samples rather than leaning on the asymptotic variance \eqn{6/n}.

\bold{Why \code{"agostino"} is the default.} It holds its size where the asymptotic form does not. Measured here, one-sided at a nominal 5\% over 4,000 replications of symmetric residuals:

\tabular{lllll}{
  \eqn{n} \tab 25 \tab 50 \tab 100 \tab 400 \cr
  \code{coelli} \tab 0.031 \tab 0.037 \tab 0.043 \tab 0.050 \cr
  \code{agostino} \tab 0.045 \tab 0.046 \tab 0.046 \tab 0.051
}

Coelli's is conservative below about \eqn{n = 100}; the two agree by \eqn{n = 400}. The default also matches \pkg{sfaR}'s. The implementation agrees with \code{moments::agostino.test} to within \eqn{10^{-10}}.

\bold{How to read the result.} A small p-value under the default alternative is evidence \emph{for} the frontier specification: the residuals are skewed the way inefficiency makes them. A large p-value is not evidence of a problem by itself -- it may simply mean there is little inefficiency to detect, or too little data to detect it. What matters is the combination with \code{wrong_skew}: residuals skewed the wrong way \emph{and} a fit reporting \code{$sigma_u_at_bound} is the Type I failure of Olson, Schmidt and Waldman (1980), in which case the boundary really is the maximum likelihood estimate (Waldman 1982) and the efficiency scores carry no information. See \code{\link{sfm}}.
}
\references{
Coelli, T. (1995). Estimators and hypothesis tests for a stochastic frontier function: a Monte Carlo analysis. \emph{Journal of Productivity Analysis}, 6(3), 247--268.

D'Agostino, R.B. (1970). Transformation to normality of the null distribution of \eqn{g_1}. \emph{Biometrika}, 57(3), 679--681.

Olson, J.A., Schmidt, P. and Waldman, D.M. (1980). A Monte Carlo study of estimators of stochastic frontier production functions. \emph{Journal of Econometrics}, 13(1), 67--82.

Schmidt, P. and Lin, T.-F. (1984). Simple tests of alternative specifications in stochastic frontier models. \emph{Journal of Econometrics}, 24(3), 349--361.

Waldman, D.M. (1982). A stationary point for the stochastic frontier likelihood. \emph{Journal of Econometrics}, 18(2), 275--279.
}
\seealso{\code{\link{moment_range}} and \code{\link{spec_test}} for the Type II counterpart, \code{\link{sfm}}, \code{\link{sfa_diagnostics}}}
\examples{
\donttest{
cs <- data_gen_cs(N = 200, rand = 1, sig_u = 0.3, sig_v = 0.4, cons = 0.5,
                  beta1 = 0.5, beta2 = 0.5, a = 5, mu = 0.1)
fit <- sfm(y_pcs_e ~ x1 + x2, model_name = "NE", data = cs)

skewness_test(fit)                    # D'Agostino, one-sided
skewness_test(fit, test = "coelli")   # Coelli's M3T

## Also works on a bare residual vector.
skewness_test(residuals(lm(y_pcs_e ~ x1 + x2, data = cs)))
}
}
