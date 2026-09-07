## Adapters between the tuning/weight functions and an "sfareg" NHN fit.
## Gap L17, supporting Bernstein, Parmeter and Wright (2026).
## See notes/code_history/robust_tuning.md.
##
## These are the only pieces that touch the fitted-object layout, so they are
## the only pieces to revisit if that layout changes. Everything else in
## robust_hscore.R, robust_tuning.R and robust_weights.R works from a residual
## vector and two scale parameters, and can be used without a fitted object.

## The package depends on R >= 4.0.0, where base has no `%||%`.
`%||%` <- function(a, b) if (is.null(a)) b else a

## sfm() reports NHN as (lambda, sigma), NOT as (sigma_v, sigma_u) -- see
## start_cs(). Inverting the reparameterisation here, with the same two lines
## sfm()'s own robust branch uses, keeps the accessor and the estimator in
## agreement by construction. Direct sigma_v/sigma_u names are still honoured
## first, for a hand-built list or a future layout that carries them.
.NHN_SV <- c("sigma_v", "sigmav", "s_v")
.NHN_SU <- c("sigma_u", "sigmau", "s_u")

.robust_scales <- function(cf) {
  hv <- intersect(.NHN_SV, names(cf))
  hu <- intersect(.NHN_SU, names(cf))
  if (length(hv) && length(hu)) {
    return(c(sigma_v = as.numeric(cf[[hv[1]]]), sigma_u = as.numeric(cf[[hu[1]]])))
  }
  if (all(c("lambda", "sigma") %in% names(cf))) {
    lam <- as.numeric(cf[["lambda"]])
    sig <- as.numeric(cf[["sigma"]])
    su <- (lam * sig) / sqrt(1 + lam^2)
    return(c(sigma_v = su / lam, sigma_u = su))
  }
  NULL
}

.robust_get <- function(object, what) {
  if (inherits(object, "sfareg")) {
    if (what == "c") {
      cc <- object$robust_c
      return(if (is.null(cc) || !is.finite(cc)) 0 else as.numeric(cc))
    }
    s <- .robust_scales(object$coefficients)
    if (is.null(s)) {
      stop("Could not recover ", what, " from the fit: its coefficients are ",
        "neither (sigma_v, sigma_u) nor (lambda, sigma). Names present: ",
        paste(names(object$coefficients), collapse = ", "), call. = FALSE)
    }
    return(unname(s[[what]]))
  }
  object[[what]]
}

## Composed residuals in the v - u orientation, matching what sfm()'s own
## likelihood forms: eps = inefdec_n * (y - x'beta), with inefdec_n = -1 for a
## cost frontier. sfm() stores no model frame, so the data is re-fetched from
## the call. It is evaluated in the formula's own environment rather than in
## the caller's: the formula was created where the data lives, whereas the
## caller here is whatever function happened to invoke the accessor.
.robust_residuals <- function(object) {
  if (!inherits(object, "sfareg")) {
    stop("Not an 'sfareg' fit.", call. = FALSE)
  }
  if (!is.null(object$residuals)) return(as.numeric(object$residuals))
  cl <- object$call
  if (is.null(cl)) {
    stop("The fit carries no call; cannot rebuild residuals.", call. = FALSE)
  }
  env <- environment(object$formula) %||% parent.frame()
  dat <- tryCatch(eval(cl$data, envir = env), error = function(e) NULL)
  if (is.null(dat)) dat <- tryCatch(eval(cl$data, envir = parent.frame()),
    error = function(e) NULL)
  if (is.null(dat)) {
    stop("Could not re-fetch the data named in the fit's call (`",
      deparse(cl$data), "`). Supply the residuals directly.", call. = FALSE)
  }
  mf <- stats::model.frame(object$formula, data = dat)
  X <- stats::model.matrix(object$formula, mf)
  y <- as.numeric(stats::model.response(mf))
  cf <- object$coefficients
  beta <- cf[intersect(colnames(X), names(cf))]
  if (length(beta) != ncol(X)) {
    stop("Could not match the frontier coefficients to the design matrix.",
      call. = FALSE)
  }
  sgn <- if (isFALSE(eval(cl$inefdec, envir = env))) -1 else 1
  as.numeric(sgn * (y - X[, names(beta), drop = FALSE] %*% beta))
}

## Reference start for the multistart guard. `start_v` is reassigned through
## every optimiser stage in sfm(), so on a returned fit it holds the parameters
## the fit ENDED at, not the ones it began from -- which is what is wanted here
## anyway: the MLE is the fixed point to fall back to when a warm start has
## wandered into the degenerate basin. The fit itself is returned, because
## sfm(start_from = ) takes an "sfareg" object and matches by parameter name.
.robust_ref_start <- function(object) object

## A closure that refits the same model at a given tuning value from a given
## starting fit, and returns the pieces hscore_select() needs. Built by
## re-issuing the original call with the robust arguments replaced, so it
## inherits every option the user set.
.robust_refit_fn <- function(object, method) {
  if (!inherits(object, "sfareg")) {
    stop("hscore_select() needs a fit from sfm().", call. = FALSE)
  }
  if (!identical(object$model_name, "NHN")) {
    stop("Tuning selection is implemented for the Normal--Half-Normal model; ",
      "this fit is '", object$model_name, "'.", call. = FALSE)
  }
  cl <- object$call
  env <- environment(object$formula) %||% parent.frame(2)
  tune_arg <- switch(method, mlqe = "c_mlqe", psi = "eta", mdpd = "alpha")

  function(cc, start) {
    cl2 <- cl
    cl2$robust <- if (cc <= 1e-10) "mle" else method
    cl2[[tune_arg]] <- cc
    cl2$start_from <- start
    cl2$verbose <- FALSE
    fit <- tryCatch(suppressMessages(suppressWarnings(eval(cl2, envir = env))),
      error = function(e) NULL)
    if (is.null(fit)) return(list(ok = FALSE))
    s <- tryCatch(.robust_scales(fit$coefficients), error = function(e) NULL)
    sv <- if (is.null(s)) NA_real_ else unname(s[["sigma_v"]])
    su <- if (is.null(s)) NA_real_ else unname(s[["sigma_u"]])
    e <- tryCatch(.robust_residuals(fit), error = function(e) NULL)
    ok <- is.finite(sv) && is.finite(su) && sv > 1e-8 && su > 1e-8 && !is.null(e)
    list(ok = ok, fit = fit,
      sigma_v = sv, sigma_u = su, residuals = e,
      ## The optimiser's attained objective, negated so larger is better. Only
      ## ever compared between two starts at the SAME tuning value, which is
      ## the only comparison it supports: across c the objective changes.
      objective = tryCatch(-as.numeric(fit$opt[["value"]]),
        error = function(e) NA_real_))
  }
}
