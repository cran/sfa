## ------------------------------------------------------------
## S3 methods for class "sfareg"
##
## psfm(), sfm(), zsfm(), and ttsfm() all return objects of class
## "sfareg", but previously only print.sfareg() and summary.sfareg()
## existed (see data.processing.R). Users had to reach into
## object$coefficients / object$opt by hand for anything else. These
## add the standard R modeling generics so sfareg objects behave like
## other model objects (coef(fit), logLik(fit), AIC(fit), BIC(fit),
## vcov(fit)).
##
## Sign convention: every estimator in this package minimizes its
## objective function (bobyqa/psoptim/optim all minimize by default;
## see opts.R, which does not flip the sign of the user-supplied fn),
## and each fn is written to return the NEGATIVE summed
## log-likelihood -- so logLik = -object$opt$value. This matches
## print.sfareg()/summary.sfareg(), which already display
## -object$opt$value as "log likelihood".
##
## psfm()'s GTRE_SEQ1/GTRE_SEQ2 branches are moment-based (not MLE)
## and carry no $opt component at all; logLik/AIC/BIC are undefined
## for those fits and will return NA with a warning rather than
## erroring.
## ------------------------------------------------------------

coef.sfareg <- function(object, ...) {
  object$coefficients
}

vcov.sfareg <- function(object, ...) {
  p <- length(object$coefficients)
  nm <- names(object$coefficients)

  if (!is.null(object$opt) && !is.null(object$opt$hessian)) {
    V <- tryCatch(solve(object$opt$hessian), error = function(e) NULL)
    if (!is.null(V)) {
      dimnames(V) <- list(nm, nm)
      return(V)
    }
    warning("Hessian stored on this fit could not be inverted; falling back to a diagonal approximation built from the reported standard errors.", call. = FALSE)
  }

  if (!is.null(object$std.errors) && !all(is.na(object$std.errors))) {
    V <- diag(object$std.errors^2, nrow = p)
    dimnames(V) <- list(nm, nm)
    return(V)
  }

  warning("No Hessian or standard errors are available on this fit (was optHessian = FALSE?); returning a matrix of NAs.", call. = FALSE)
  matrix(NA_real_, p, p, dimnames = list(nm, nm))
}

logLik.sfareg <- function(object, ...) {
  if (is.null(object$opt) || is.null(object$opt$value)) {
    warning("This fit has no stored optimizer output (e.g. psfm()'s GTRE_SEQ1/GTRE_SEQ2 are moment-based, not maximum likelihood), so logLik() is not defined for it.", call. = FALSE)
    ## Return a properly classed logLik carrying NA rather than a bare
    ## NA_real_. stats::AIC()/BIC() read the "df" and "nobs" attributes off
    ## whatever logLik() returns; given an unclassed NA they find no "df" and
    ## silently produce numeric(0) instead of the documented NA -- a missing
    ## value that disappears rather than propagating.
    val <- NA_real_
    attr(val, "df") <- length(object$coefficients)
    attr(val, "nobs") <- nobs.sfareg(object)
    class(val) <- "logLik"
    return(val)
  }
  val <- -object$opt$value
  attr(val, "df") <- length(object$coefficients)
  attr(val, "nobs") <- nobs.sfareg(object)
  class(val) <- "logLik"
  val
}

nobs.sfareg <- function(object, ...) {
  if (!is.null(object$data)) {
    return(nrow(object$data))
  }
  ## sfm()/zsfm()/ttsfm() do not store the data on the fitted object, so fall
  ## back to re-evaluating the `data` argument of the recorded call.
  ##
  ## It must be evaluated in the environment the model was FIT in, which is
  ## captured on the formula, not in parent.frame(). parent.frame() here is
  ## whoever happened to call nobs() -- fine from the console, where the data
  ## is usually a global, but NA from inside any function whose caller cannot
  ## see it. That silently propagated: BIC() needs the "nobs" attribute, so
  ## BIC() on an sfm() fit returned NA whenever it was called from inside a
  ## function, while AIC() (which needs only "df") kept working.
  dcall <- object$call$data
  if (is.null(dcall)) {
    return(NA_integer_)
  }
  envs <- list(
    if (inherits(object$formula, "formula")) environment(object$formula) else NULL,
    parent.frame(), globalenv()
  )
  for (e in envs) {
    if (is.null(e)) next
    n <- tryCatch(nrow(eval(dcall, envir = e)), error = function(err) NULL)
    if (!is.null(n)) {
      return(as.integer(n))
    }
  }
  NA_integer_
}


## ---------------------------------------------------------------------------
## predict() / fitted() / residuals() for "sfareg".
##
## These are the standard modelling generics `frontier` and `sfaR` provide and
## sfa previously lacked. Two details make a generic implementation possible
## across every model in the package:
##
##   * The frontier coefficients are always NAMED after the regressors, while
##     the auxiliary parameters are not (lambda/sigma for NHN, sigv/sigu for
##     NE, sigmaSq/gamma for PL80, and so on). Selecting coefficients by
##     matching against the design-matrix column names therefore isolates
##     beta without needing a per-model index table, which would otherwise
##     have to be maintained in lock step with every new model.
##   * Only the FIRST part of a pipe formula (y ~ x | z | zp) describes the
##     frontier; the later parts parameterize variances. Prediction uses
##     rhs = 1 only.
## ---------------------------------------------------------------------------

## Recover the data a fit was built from: explicit newdata wins, then the
## copy some models store on the object, then the `data` argument of the
## original call. Fails loudly rather than silently predicting from the
## wrong frame.
.sfa_data <- function(object, newdata = NULL) {
  if (!is.null(newdata)) {
    return(as.data.frame(newdata))
  }
  if (!is.null(object$data)) {
    return(as.data.frame(object$data))
  }
  dat <- tryCatch(eval(object$call$data, parent.frame(3)), error = function(e) NULL)
  if (is.null(dat)) {
    stop("Cannot recover the data used to fit this model. Pass `newdata` explicitly.",
      call. = FALSE
    )
  }
  as.data.frame(dat)
}

## Frontier design matrix and the matching beta, for any sfareg model.
.sfa_xb <- function(object, newdata = NULL) {
  dat <- .sfa_data(object, newdata)
  f1 <- stats::formula(Formula::Formula(object$formula), lhs = 1, rhs = 1)
  mm <- stats::model.matrix(stats::delete.response(stats::terms(f1, data = dat)), data = dat)
  cf <- object$coefficients
  nm <- intersect(colnames(mm), names(cf))
  if (!length(nm)) {
    stop("None of this fit's coefficients match the frontier design matrix; ",
      "cannot form a prediction.",
      call. = FALSE
    )
  }
  list(xb = as.numeric(mm[, nm, drop = FALSE] %*% cf[nm]), data = dat, formula = f1)
}

predict.sfareg <- function(object, newdata = NULL,
                           type = c("frontier", "response", "efficiency"), ...) {
  type <- match.arg(type)
  if (identical(type, "efficiency")) {
    if (!is.null(newdata)) {
      stop("type = \"efficiency\" is only available for the estimation sample: ",
        "predicting efficiency requires the composed residual, which needs ",
        "the response.",
        call. = FALSE
      )
    }
    te <- object$exp_u_hat
    if (is.null(te)) te <- object$U
    if (is.null(te)) {
      stop("This model (", object$model_name, ") does not return an efficiency ",
        "prediction; see ?sfm for which models do.",
        call. = FALSE
      )
    }
    return(as.numeric(te))
  }
  z <- .sfa_xb(object, newdata)
  if (identical(type, "frontier")) {
    return(z$xb)
  }

  ## "response": the frontier shifted by predicted inefficiency, i.e. an
  ## estimate of E[y | x] rather than of the frontier itself. Sign follows
  ## the production/cost convention the model was fitted under.
  u <- object$u_hat
  if (is.null(u) && !is.null(object$exp_u_hat)) u <- -log(pmax(object$exp_u_hat, .Machine$double.xmin))
  if (is.null(u)) {
    stop("type = \"response\" needs an inefficiency prediction, which model ",
      object$model_name, " does not return.",
      call. = FALSE
    )
  }
  if (!is.null(newdata)) {
    stop("type = \"response\" is only available for the estimation sample.", call. = FALSE)
  }
  s <- if (isTRUE(object$call$inefdec == FALSE)) -1 else 1
  z$xb - s * as.numeric(u)
}

fitted.sfareg <- function(object, ...) .sfa_xb(object)$xb

residuals.sfareg <- function(object, ...) {
  z <- .sfa_xb(object)
  yv <- all.vars(z$formula)[1]
  if (!(yv %in% names(z$data))) {
    stop("Response '", yv, "' not found in the fit's data; cannot form residuals.",
      call. = FALSE
    )
  }
  as.numeric(z$data[[yv]]) - z$xb
}
