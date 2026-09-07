## extract() for "sfareg", so texreg renders these fits straight into LaTeX or
## HTML regression tables. See notes/code_history/texreg_methods.md.
##
## This is the last step of nearly every applied workflow: without it a user
## assembling a paper table formats `$out` by hand, which is where transcription
## errors come from.
##
## Registered with the `texreg::extract` form in NAMESPACE rather than a plain
## S3method(), because texreg is in Suggests and its generic does not exist when
## sfa loads. The sandwich methods were defined without that and were
## unreachable in an installed package for a full release; the same mistake is
## not worth making twice.

## Which rows of `out` are frontier coefficients rather than variance or shape
## parameters. A regression table wants the former in the body and the latter
## either omitted or clearly marked, not silently interleaved.
.sfa_scale_names <- c(
  "sigma_u", "sigma_v", "sigma", "lambda", "sigma_c", "gamma", "sigmaSq",
  "mu", "theta", "b", "nu", "alpha", "shape", "k", "p", "rho",
  "sigr", "sigh", "sigu", "sigv", "logit0", "gamma_uv", "gamma_hr",
  "sigmaSq_uv", "sigmaSq_hr"
)

.sfa_is_scale <- function(nm) {
  nm %in% .sfa_scale_names |
    grepl("^(sigma|sig|lambda|gamma|delta_|rho_|logit|theta|nu|shape)", nm)
}

extract.sfareg <- function(model,
                           include.scales = TRUE,
                           include.loglik = TRUE,
                           include.aic = TRUE,
                           include.bic = TRUE,
                           include.nobs = TRUE,
                           ...) {
  if (!requireNamespace("texreg", quietly = TRUE)) {
    stop("extract(): the texreg package is required.", call. = FALSE)
  }
  out <- model$out
  if (is.null(out) || !is.matrix(out)) {
    stop("extract(): this fit has no `out` matrix to render.", call. = FALSE)
  }

  nm <- rownames(out)
  keep <- if (isTRUE(include.scales)) rep(TRUE, length(nm)) else !.sfa_is_scale(nm)

  co <- as.numeric(out[keep, "par"])
  se <- as.numeric(out[keep, "st_err"])
  nms <- nm[keep]

  ## Two-sided normal p-values. The t-value column is already par/st_err, but
  ## it is recomputed here so a fit that stored NA standard errors yields NA
  ## p-values rather than whatever happens to sit in that column.
  pv <- 2 * stats::pnorm(-abs(co / se))

  gof <- numeric(0)
  gof.names <- character(0)
  gof.decimal <- logical(0)

  ## logLik() warns and returns NA for the moment-based estimators (C2SLS, the
  ## GTRE_SEQ pair, the classical panel models), which is correct behaviour --
  ## they maximise nothing. Suppress the warning here and simply omit the row,
  ## rather than printing an NA that invites a comparison that cannot be made.
  ll <- suppressWarnings(tryCatch(as.numeric(stats::logLik(model)),
    error = function(e) NA_real_
  ))
  if (isTRUE(include.loglik) && is.finite(ll)) {
    gof <- c(gof, ll)
    gof.names <- c(gof.names, "Log Likelihood")
    gof.decimal <- c(gof.decimal, TRUE)
  }
  if (isTRUE(include.aic) && is.finite(ll)) {
    a <- suppressWarnings(tryCatch(stats::AIC(model), error = function(e) NA_real_))
    if (is.finite(a)) {
      gof <- c(gof, a)
      gof.names <- c(gof.names, "AIC")
      gof.decimal <- c(gof.decimal, TRUE)
    }
  }
  if (isTRUE(include.bic) && is.finite(ll)) {
    b <- suppressWarnings(tryCatch(stats::BIC(model), error = function(e) NA_real_))
    if (is.finite(b)) {
      gof <- c(gof, b)
      gof.names <- c(gof.names, "BIC")
      gof.decimal <- c(gof.decimal, TRUE)
    }
  }
  if (isTRUE(include.nobs)) {
    n <- suppressWarnings(tryCatch(stats::nobs(model), error = function(e) NA_integer_))
    if (!is.na(n)) {
      gof <- c(gof, n)
      gof.names <- c(gof.names, "Num. obs.")
      gof.decimal <- c(gof.decimal, FALSE)
    }
  }

  texreg::createTexreg(
    coef.names = nms,
    coef = co,
    se = se,
    pvalues = pv,
    gof.names = gof.names,
    gof = gof,
    gof.decimal = gof.decimal,
    model.name = if (!is.null(model$model_name)) as.character(model$model_name) else ""
  )
}
