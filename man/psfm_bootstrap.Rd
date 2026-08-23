\name{psfm_bootstrap}
\alias{psfm_bootstrap}
\title{Parametric Bootstrap for \code{psfm()} Panel Stochastic Frontier Models}
\description{
Performs a parametric bootstrap for panel stochastic frontier models fit with
\code{\link[sfa]{psfm}}. Supports \code{model_name \%in\% c("GTRE_Z", "TRE_Z",
"GTRE", "TRE", "TFE", "FD")} -- every \code{psfm()} model except the
moment-based/LSDV estimators (\code{"GTRE_SEQ1"}, \code{"GTRE_SEQ2"},
\code{"SSFE"}, which are not maximum-likelihood and have no \code{$opt} to
bootstrap around in the same sense) and the \code{\link[frontier]{sfa}}-wrapped
models (\code{"PL80"}, \code{"BC92"}, whose \code{$out}/efficiency-score
fields don't follow this package's own conventions). Each supported model's
frontier/variance-component parameter blocks are read directly from
\code{psfm_object$out} and \code{psfm_object$formula}; for each of
\code{BOOT} replications, a new response is simulated from the fitted
parameters using a data-generating process appropriate to that specific
model (see Details), the model is re-estimated in parallel, and bootstrap
standard errors / t-values are returned for every parameter in
\code{psfm_object$out}.
}
\usage{
psfm_bootstrap(psfm_object,
               numCores,
               BOOT,
               individual,
               h_type        = c("auto", "none", "scalar", "parametric"),
               maxit.psoptim = 1000,
               seed_offset   = 0,
               write_back    = TRUE,
               pkgs          = c("sfa", "Formula", "pbapply", "truncnorm"),
               inefdec,
               rand.gtre     = NULL,
               rand.psoptim  = NULL,
               maxit.bobyqa  = 1,
               maxit.optim   = 1)
}
\arguments{
  \item{psfm_object}{A fitted model object of class \code{"sfareg"} returned by
    \code{\link[sfa]{psfm}}, with \code{model_name \%in\% c("GTRE_Z", "TRE_Z",
    "GTRE", "TRE", "TFE", "FD")}. Must contain \code{$out}, \code{$data},
    \code{$formula}, and \code{$model_name}; \code{$H} is additionally
    required for \code{"GTRE"}/\code{"GTRE_Z"} (the only two models that
    return it) and \code{$r_hat_m} is additionally required for
    \code{"TFE"}.}
  \item{numCores}{Integer. Number of cores for the parallel cluster created
    internally via \code{\link[parallel]{makeCluster}}.}
  \item{BOOT}{Integer. Number of bootstrap replications (e.g. 199, 399, 999).}
  \item{individual}{Character. Name of the individual/panel-ID column in
    \code{psfm_object$data} (e.g. \code{"firm.id"}).}
  \item{h_type}{One of \code{"auto"}, \code{"none"}, \code{"scalar"}, or
    \code{"parametric"}. Only meaningful for \code{model_name \%in\%
    c("GTRE_Z", "TRE_Z", "GTRE", "TRE")} (ignored for \code{"TFE"}/\code{"FD"},
    which have no random persistent-inefficiency \code{h} term); controls how
    the individual-specific, time-invariant inefficiency term \code{h} is
    simulated:
    \describe{
      \item{\code{"none"}}{The model has no \code{h} term at all (\code{TRE},
        \code{TRE_Z}).}
      \item{\code{"scalar"}}{\code{h_i ~ abs(N(0, sigma_h))}, a single scalar
        standard deviation -- used automatically for bare \code{"GTRE"} (which
        has no \code{h}-covariates at all, only a single fitted
        \code{sigh}).}
      \item{\code{"parametric"}}{\code{h_i ~ abs(N(0, exp(data_h \%*\% h_vec)))},
        using the third right-hand-side part of the model formula as the
        \code{h}-covariates -- used automatically for \code{"GTRE_Z"}.}
      \item{\code{"auto"}}{(Default) Resolves per \code{model_name}:
        \code{"parametric"} for \code{"GTRE_Z"}, \code{"scalar"} for bare
        \code{"GTRE"}, \code{"none"} for \code{"TRE_Z"}/\code{"TRE"}.}
    }
  }
  \item{maxit.psoptim}{Integer. Passed through to \code{psfm()} as
    \code{maxit.psoptim} on each bootstrap re-estimation. Default \code{1000}.}
  \item{seed_offset}{Integer added to the bootstrap index \code{b} when calling
    \code{set.seed()}, allowing reproducible-but-distinct seeds across multiple
    bootstrap runs. Default \code{0}.}
  \item{write_back}{Logical. If \code{TRUE} (default), the returned list
    includes \code{model}, a copy of \code{psfm_object} with the bootstrap
    standard errors and t-values written into \code{$out[,2]} and
    \code{$out[,3]}. If \code{FALSE}, only the raw bootstrap results are
    returned.}
  \item{pkgs}{Character vector of package names to load on each cluster
    worker via \code{\link[parallel]{clusterCall}}. Default
    \code{c("sfa", "Formula", "pbapply", "truncnorm")} -- \pkg{truncnorm} is
    only actually used when bootstrapping an \code{"FD"} fit, but is included
    in the default so a single default works for every supported model.}
  \item{inefdec}{Logical, with no default; must be supplied. Passed through to
    \code{psfm()} on each bootstrap re-estimation. Should match the
    \code{inefdec} setting used to originally fit \code{psfm_object}.}
  \item{rand.gtre}{Passed through to \code{psfm()} as \code{rand.gtre} on each
    bootstrap re-estimation. Default \code{NULL}.}
  \item{rand.psoptim}{Passed through to \code{psfm()} as \code{rand.psoptim} on
    each bootstrap re-estimation. Default \code{NULL}.}
  \item{maxit.bobyqa}{Integer. Passed through to \code{psfm()} as
    \code{maxit.bobyqa} on each bootstrap re-estimation. Default \code{1}.}
  \item{maxit.optim}{Integer. Passed through to \code{psfm()} as
    \code{maxit.optim} on each bootstrap re-estimation. Default \code{1}.}
}
\details{
Supported models fall into three families, each with its own data-generating
process and \code{$out} row layout (see \code{R/psfm_bootstrap.R}'s header
comment for the full derivation, including how each was cross-checked
against that model's own likelihood in \code{psfm.R}):

\describe{
  \item{\code{"randeff"} family (\code{GTRE_Z}, \code{TRE_Z}, \code{GTRE},
    \code{TRE})}{\code{y* = X \%*\% beta_hat + v + u + r + h} (sign per
    \code{inefdec}), with \code{v ~ N(0, sigma_v^2)}, \code{u} half-normal
    (covariate-driven sigma for \code{GTRE_Z}/\code{TRE_Z}, a single
    homoskedastic scalar for bare \code{GTRE}/\code{TRE}, which have no
    \code{z} pipe at all), \code{r ~ N(0, sigma_r^2)} (one draw per
    individual, repeated over time), and \code{h} per \code{h_type}
    (\code{GTRE}/\code{GTRE_Z} only). \code{$out} row layout: \code{sigv,
    sigr, [x-block], [z-block], [h-block]} for \code{GTRE_Z}/\code{TRE_Z};
    \code{lambda, sigma, sigr, [sigh], [x-block]} for bare
    \code{GTRE}/\code{TRE} (note the different reparameterization -- see
    \code{\link[sfa]{psfm}}'s own model-family sigma-link inconsistency
    note).}
  \item{\code{"tfe"} family (\code{TFE})}{Fixed-effects model: the individual
    effects are held FIXED at \code{psfm_object$r_hat_m}'s original point
    estimates across every replication (not redrawn -- redrawing them would
    contradict treating them as fixed rather than random), and only
    \code{v ~ N(0, sigma_v^2)} and \code{u ~ |N(0, sigma_u^2)|} (both derived
    from the fitted \code{lambda}/\code{sig}, or \code{gamma}/\code{sig} if
    the original fit used \code{psfm(..., gamma = TRUE)}) are redrawn.
    \code{$out} row layout: \code{lambda (or gamma), sig, [x-block]} (no
    intercept row -- the within-transformation removes it).}
  \item{\code{"fd"} family (\code{FD})}{Wang and Ho (2010) first-difference
    estimator: a single time-invariant \code{u_i ~ TruncatedNormal(mean =
    mu, sd = sqrt(sig_u2), lower = 0)} draw per individual (via
    \code{\link[truncnorm]{rtruncnorm}}), scaled by a deterministic
    time-varying factor \code{h_it = exp(z_it \%*\% delta)}, plus
    \code{v_it ~ N(0, sig_v2)}. \code{$out} row layout: \code{sig_u2, sig_v2,
    mu, [x-block], [z-block]} (no intercept row in either the x- or z-block).}
}

\code{psfm_bootstrap} validates that \code{nrow(psfm_object$out)} matches the
expected layout for the resolved model/\code{h_type} before running any
bootstrap replications, and stops with an informative error if it does not.
For \code{GTRE}/\code{GTRE_Z} it similarly validates up front that
\code{length(psfm_object$H)} equals the number of unique individuals implied
by the \code{individual} column (\code{$H} is one value per individual, not
one per observation) -- \code{TRE}/\code{TRE_Z}/\code{TFE}/\code{FD} do not
return \code{$H} at all, so \code{boot_eff_h} is simply \code{NULL} for those
four models (not an error).

The model is re-estimated on each simulated dataset via \code{sfa::psfm()},
executed in parallel across replications using
\code{\link[parallel]{parLapply}} (or \code{\link[pbapply]{pblapply}} for a
live progress bar, when the \pkg{pbapply} package is installed). Whether the
\code{psoptim} (particle-swarm) optimizer stage is used on each refit is
model-dependent: it is enabled for the \code{"randeff"} family (matching this
function's original behavior, tested clean), but disabled for \code{"TFE"}
-- testing found \code{psoptim}'s wide exploration made \code{TFE}
repeatedly collapse onto a degenerate variance-boundary solution (reproduced
even refitting the SAME real data across different \code{psoptim} seeds), a
problem that disappeared entirely once that stage was skipped (matching
\code{psfm()}'s own default of \code{PSopt = FALSE}). As a defense-in-depth
safety net across every family, any replication whose refit still lands on a
degenerate variance-boundary mode (a scale parameter pinned near its
numerical floor -- the same well-understood boundary-mode phenomenon
documented for \code{ttsfm()}'s \code{TTNE}/\code{TTHN} models) is treated
as a per-replication FAILURE rather than silently corrupting the bootstrap
standard error with a meaningless draw.

Re-estimation failures (optimizer errors, or the degenerate-mode check just
described) are caught per replication, recorded with their message, and set
to \code{NA} in the output matrices rather than stopping the bootstrap; a
single \code{\link{warning}} at the end lists every failed replication index
together with its message. Bootstrap standard errors are computed with
\code{na.rm = TRUE}, so they remain well-defined even when some replications
failed.

Requires the \pkg{parallel} and \pkg{Formula} packages (and \pkg{truncnorm}
when bootstrapping an \code{"FD"} fit), and requires \pkg{sfa} (and any other
packages named in \code{pkgs}) to be installed and loadable on each cluster
worker.
}
\value{
A list with components:
  \item{boot_par}{A \code{BOOT} x (\code{nrow(psfm_object$out)} + 2) matrix of
    bootstrapped parameter draws, with the final two columns \code{"loglik"}
    and \code{"hours"}. Rows corresponding to failed replications are
    \code{NA}.}
  \item{boot_eff}{A \code{BOOT} x \code{nrow(psfm_object$data)} matrix of
    bootstrapped per-observation efficiency-score draws (\code{$U} for the
    \code{"randeff"} family, \code{$exp_u_hat} for \code{"TFE"}/\code{"FD"}
    -- the field is named differently across these models in \code{psfm.R},
    but always collected into this same \code{boot_eff} matrix here). Rows
    corresponding to failed replications are \code{NA}.}
  \item{boot_eff_h}{Only non-\code{NULL} for \code{model_name \%in\%
    c("GTRE", "GTRE_Z")} (the only two models that return \code{$H} at all);
    \code{NULL} otherwise. A \code{BOOT} x \code{length(psfm_object$H)}
    matrix of bootstrapped time-invariant, individual-specific
    inefficiency-score (\code{$H}) draws -- one column per unique individual
    (not one column per observation; \code{$H} does not vary over time
    within an individual). Columns are named by the individual ID values.
    Rows corresponding to failed replications are \code{NA}.}
  \item{se}{Named numeric vector of bootstrap standard errors for every
    parameter in \code{psfm_object$out}.}
  \item{tval}{Named numeric vector of bootstrap t-values
    (\code{estimate / boot SE}) for every parameter in \code{psfm_object$out}.}
  \item{failures}{Integer vector of replication indices \code{b} for which
    re-estimation failed, or \code{NULL} if none failed.}
  \item{model}{Only present when \code{write_back = TRUE}. A copy of
    \code{psfm_object} with \code{$out[,2]} and \code{$out[,3]} replaced by
    the bootstrap standard errors and t-values.}
}
\author{David H. Bernstein and Christopher F. Parmeter}
\seealso{
\code{\link[sfa]{psfm}}, \code{\link[parallel]{parLapply}},
\code{\link[pbapply]{pblapply}}, \code{\link[Formula]{Formula}}
}
\examples{
\dontrun{
## Not run: a meaningful bootstrap needs BOOT in the hundreds, each
## replication refitting the model on a PSOCK cluster, so this takes
## minutes rather than seconds.

data_trial <- data_gen_p(t = 5, N = 50, rand = 100,
                         sig_u = 1,   sig_v = 0.3,
                         sig_r = 0.2, sig_h = 0.4,
                         cons  = 0.5, beta1 = 0.5, beta2 = 0.5)

GTRE_ZZ <- psfm(y_gtre_zz ~ x1 + x2 | z_gtre | zp_gtre,
                model_name = "GTRE_Z",
                data       = data_trial,
                individual = "name")

boot_results <- psfm_bootstrap(
  psfm_object = GTRE_ZZ,
  numCores    = 2,
  BOOT        = 199,
  individual  = "name",
  h_type      = "auto",
  inefdec     = FALSE
)

## Bootstrap standard errors / t-values for every model parameter
boot_results$se
boot_results$tval

## Bootstrap draws of time-varying (U) and time-invariant (H)
## efficiency/inefficiency scores, one row per replication
boot_results$boot_eff
boot_results$boot_eff_h

## Fitted object with bootstrap SEs/t-values written into $out
boot_results$model$out
}
}

\keyword{models}
\keyword{htest}
