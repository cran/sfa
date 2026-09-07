## S3 methods for class "sfareg"

coef.sfareg <- function(object, ...) {
  object$coefficients
}

vcov.sfareg <- function(object, type = c("hessian", "bhhh"), ...) {
  type <- match.arg(type)
  p <- length(object$coefficients)
  nm <- names(object$coefficients)

  ## BHHH: the inverse of the outer product of the per-observation scores.
  ## Worth having because it needs no Hessian at all -- it is defined whenever
  ## the scores are, which is exactly the case the default path fails on. When
  ## the Hessian is singular vcov() currently falls back to a DIAGONAL
  ## approximation, discarding every covariance; BHHH keeps them.
  ##
  ## It is only as good as the information-matrix equality, so it disagrees
  ## with the Hessian under misspecification -- which is a feature when the two
  ## are compared deliberately and a trap when they are not. Not the default.
  if (identical(type, "bhhh")) {
    G <- tryCatch(estfun.sfareg(object), error = function(e) NULL)
    if (is.null(G)) {
      stop("vcov(type = \"bhhh\"): the score matrix is unavailable for this ",
        "fit. Refit with keep_objective = TRUE.",
        call. = FALSE
      )
    }
    V <- tryCatch(solve(crossprod(G)), error = function(e) NULL)
    if (is.null(V)) {
      stop("vcov(type = \"bhhh\"): the outer product of gradients is singular.",
        call. = FALSE
      )
    }
    dimnames(V) <- list(nm, nm)
    return(V)
  }

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
    ## NA_real_.
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
  ## Rows USED, not rows supplied. Re-evaluating the call (below) counts the
  ## latter, so a fit on data with any missing value reported too many
  ## observations and BIC() was computed against the wrong n.
  if (!is.null(object$nobs) && is.finite(object$nobs)) {
    return(as.integer(object$nobs))
  }
  if (!is.null(object$data)) {
    return(nrow(object$data))
  }
  ## Failing an explicit count, any vector the fit stores one element of per
  ## observation is exact where re-evaluating the call is not.
  for (nm in c("exp_u_hat", "u_hat", "residuals", "med_u_hat")) {
    v <- object[[nm]]
    if (is.numeric(v) && length(v) > 0L) {
      return(length(v))
    }
  }
  if (is.numeric(object$u_posterior$mu_star)) {
    return(length(object$u_posterior$mu_star))
  }
  ## sfm()/zsfm()/ttsfm() do not store the data on the fitted object, so fall
  ## back to re-evaluating the `data` argument of the recorded call.
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


## predict() / fitted() / residuals() for "sfareg".

## Recover the data a fit was built from: explicit newdata wins, then the copy
## some models store on the object.
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

  ## "response": the frontier shifted by predicted inefficiency, i.e.
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
