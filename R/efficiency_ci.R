## Horrace and Schmidt (1996) intervals for individual inefficiency


## Conditional on the fitted parameters, u_i | e_i is normal with mean mu_star
## and standard deviation sigma_star, truncated below at zero.
efficiency_ci <- function(object, level = 0.95, type = c("both", "u", "te")) {
  if (!inherits(object, "sfareg")) {
    stop("`object` must be an \"sfareg\" fit, as returned by sfm().", call. = FALSE)
  }
  type <- match.arg(type)

  post <- object$u_posterior
  if (is.null(post)) {
    stop(
      "efficiency_ci() needs the posterior of u given the composed residual, ",
      "which is a truncated normal only for model_name \"NHN\", \"NHN_Z\", ",
      "\"NE\" and \"NTN\". This fit is model_name \"",
      if (is.null(object$model_name)) "unknown" else object$model_name,
      "\", whose posterior has a different shape, so the Horrace-Schmidt ",
      "bounds do not apply to it.",
      call. = FALSE
    )
  }

  ci <- .horrace_schmidt_ci(post$mu_star, post$sigma_star, level = level)

  ## The point predictors are taken from the fit rather than recomputed, so
  ## what this function reports and what print()/summary() report cannot drift.
  u_hat <- object$u_hat
  if (is.null(u_hat)) {
    u_hat <- .jlms_u(post$mu_star, post$sigma_star)
  }
  te_hat <- object$exp_u_hat

  out <- data.frame(row.names = seq_along(ci$lower))
  if (type %in% c("both", "u")) {
    out$u_lower <- ci$lower
    out$u_hat <- as.numeric(u_hat)
    out$u_upper <- ci$upper
  }
  if (type %in% c("both", "te")) {
    ## Monotone decreasing transform, so the bounds swap.
    out$te_lower <- exp(-ci$upper)
    out$te_hat <- if (is.null(te_hat)) exp(-as.numeric(u_hat)) else as.numeric(te_hat)
    out$te_upper <- exp(-ci$lower)
  }

  attr(out, "level") <- level
  out
}
