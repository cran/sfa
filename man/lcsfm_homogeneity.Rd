\name{lcsfm_homogeneity}
\alias{lcsfm_homogeneity}
\alias{print.sfa_mlrt}
\title{Test a latent class frontier against homogeneity}
\description{
Tests whether a \code{\link{lcsfm}} fit's latent classes are needed at all, against the null that one class describes the data. This is the first question a latent class fit raises and the ordinary likelihood ratio test cannot answer it, because under the null the second class's parameters are unidentified, the mixing proportion sits on a boundary, and the classes are interchangeable.
}
\usage{lcsfm_homogeneity(object, null = c("bootstrap", "chisq01"),
                  B = 199, c = 1, level = 0.05,
                  seed = NULL, quiet = FALSE, envir = parent.frame())}
\arguments{
\item{object}{An \code{"sfareg"} fit from \code{\link{lcsfm}} with \code{model_name = "LCM"} and \code{n_class >= 2}.}
\item{null}{How to obtain the null distribution. \code{"bootstrap"} (the default) is a parametric bootstrap under the fitted one-class model. \code{"chisq01"} uses the published \eqn{\chi^2_{0:1}} limit and \bold{warns}: see \sQuote{Details} for why it does not apply to \code{"LCM"} as implemented.}
\item{B}{Bootstrap replications. Each is two refits, so this dominates the running time. The smallest attainable p-value is \eqn{1/(B+1)}.}
\item{c}{Tuning constant in the penalty \eqn{c\log(J^J\prod_j p_j)}. Chen \emph{et al.} (2001) report their simulations were insensitive to it; Stead, Wheat and Greene use 1 and 5.}
\item{level}{Size used for the reported verdict.}
\item{seed}{Optional integer seed for the bootstrap. The RNG state is saved and restored, so a seeded call does not disturb the caller's stream.}
\item{quiet}{Set \code{TRUE} to suppress the bootstrap progress messages.}
\item{envir}{Environment in which to re-evaluate the fit's \code{data} argument.}
}
\value{
An object of class \code{c("sfa_mlrt", "htest")}, so it prints like any other R test. Alongside \code{statistic} and \code{p.value} it carries \code{null}, \code{B}, \code{boot} (the bootstrap null draws), \code{boot_ok}, \code{penalty_c}, \code{penalty_at_max}, the penalised and unpenalised \eqn{J}-class log-likelihoods, \code{logLik_1}, \code{class_prob}, and \code{reject}.
}
\details{
\bold{The statistic.} \code{\link{lcsfm}} is refitted maximising the \emph{modified} likelihood -- the ordinary log-likelihood plus Chen \emph{et al.}'s penalty \eqn{c\log(J^J\prod_j p_j)}, which diverges to \eqn{-\infty} as any class probability approaches 0 or 1 and so keeps every class occupied. The statistic is twice the difference between that maximum and the maximised log-likelihood of the one-class model, which is \code{sfm(model_name = "NHN")} on the same formula, since a single latent class \emph{is} the ordinary frontier. The penalty is exactly zero at equal class probabilities, so the null model is recovered without a handicap.

\bold{Why the default null is a bootstrap.} Stead, Wheat and Greene report the modified statistic is asymptotically \eqn{\chi^2_{0:1}}, a 50:50 mixture of a point mass at zero and \eqn{\chi^2_1}. That result is stated for the case in which a \emph{scalar} parameter differs between classes -- Zhu and Zhang (2004) require \eqn{\theta_1} and \eqn{\theta_2} scalar, and the paper's own application holds every parameter except the noise scale common across classes.

\code{lcsfm}'s \code{"LCM"} does not do that: it lets \eqn{\sigma_v}, \eqn{\sigma_u} and the whole frontier vector vary by class. The limit therefore does not carry over, and this was measured rather than assumed. Over 200 replications of a one-class data generating process at \eqn{n = 400}:

\tabular{lll}{
  nominal size \tab 10\% \tab 5\% \cr
  actual rejection under \code{"chisq01"} \tab 82.0\% \tab 63.5\%
}

with the sampling distribution tracking \eqn{\chi^2_5} -- exactly the parameter-count difference -- rather than \eqn{\chi^2_{0:1}}, whose median is 0 against an observed 4.21.

A parametric bootstrap makes no appeal to an asymptotic distribution: it simulates from the fitted null and re-runs the whole procedure, so it is valid whatever the parameter structure. Over 60 replications of the same one-class process with \code{B = 99}:

\tabular{lll}{
  nominal size \tab 10\% \tab 5\% \cr
  \code{"chisq01"} \tab 82.0\% \tab 63.5\% \cr
  \code{"bootstrap"} \tab 10.0\% \tab 6.7\%
}

both within a binomial standard error of nominal, with the p-values indistinguishable from uniform (KS against Uniform(0,1), \eqn{p = 0.80}). On a sample that \code{"chisq01"} rejected at \eqn{p = 0.005}, the bootstrap returns \eqn{p = 0.32}, and on genuinely two-class data the statistic is 17 times the largest of 99 null draws.

\bold{Cost.} The bootstrap is \eqn{B+1} pairs of refits. Budget accordingly, and use \code{seed} for reproducibility.

\bold{\code{"LCM_Z"} is refused}, by both this function and \code{lcsfm(penalty_c = )}: it makes the class probabilities depend on covariates, and the penalty is defined on a scalar mixing proportion.
}
\references{
Chen, H., Chen, J. and Kalbfleisch, J.D. (2001). A modified likelihood ratio test for homogeneity in finite mixture models. \emph{Journal of the Royal Statistical Society B}, 63(1), 19--29.

Stead, A.D., Wheat, P. and Greene, W.H. (2023). On hypothesis testing in latent class and finite mixture stochastic frontier models, with application to a contaminated normal-half normal model. \emph{Journal of Productivity Analysis}, 60(1), 37--48.

Zhu, H.-T. and Zhang, H. (2004). Hypothesis testing in mixture regression models. \emph{Journal of the Royal Statistical Society B}, 66(1), 3--16.
}
\seealso{\code{\link{lcsfm}}, \code{\link{TIC}}, \code{\link{skewness_test}}}
\examples{
\donttest{
set.seed(21)
n   <- 300
x1  <- rnorm(n)
cls <- rbinom(n, 1, 0.5)
d   <- data.frame(y = ifelse(cls == 1, 4, 1) + x1 +
                      rnorm(n, 0, 0.4) - abs(rnorm(n, 0, 0.6)),
                  x1 = x1)

f <- lcsfm(y ~ x1, model_name = "LCM", data = d, n_class = 2)

## Small B for the example only; use the default in real work.
lcsfm_homogeneity(f, B = 19, seed = 1, quiet = TRUE)
}
}
\keyword{htest}
