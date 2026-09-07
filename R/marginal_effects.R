## Marginal effects of the variance determinants on inefficiency


## Constant columns of the z design are dropped: the derivative with respect
## to an intercept is not a marginal effect.
marginal_effects <- function(object, average = FALSE, component = c("u", "h")) {
  if (!inherits(object, "sfareg")) {
    stop("`object` must be an \"sfareg\" fit, as returned by sfm() or psfm().",
      call. = FALSE
    )
  }
  ## "h" selects the PERSISTENT block of psfm("GTRE_Z"); every other model has
  ## only the one.
  component <- match.arg(component)
  if (identical(component, "h")) {
    if (is.null(object$z_spec_h)) {
      stop(
        "component = \"h\" is only available for a psfm(model_name = \"GTRE_Z\") ",
        "fit with a third formula segment, y ~ x | z | zp, which is what ",
        "parameterizes the persistent sigma_h. This fit is model_name \"",
        if (is.null(object$model_name)) "unknown" else object$model_name,
        "\".",
        call. = FALSE
      )
    }
    zs <- object$z_spec_h
  } else {
    zs <- object$z_spec
  }
  if (is.null(zs)) {
    stop(
      "marginal_effects() needs a model whose inefficiency scale depends on ",
      "covariates. This fit is model_name \"",
      if (is.null(object$model_name)) "unknown" else object$model_name,
      "\", which has a single homoskedastic sigma_u and therefore no ",
      "marginal effect to report. Use one of the _Z models with a formula ",
      "such as y ~ x1 + x2 | z.",
      call. = FALSE
    )
  }

  Z <- zs$Z
  delta <- zs$delta
  eta <- as.numeric(Z %*% delta)

  ## sigma_u on the scale each family actually uses.
  sigma_u <- if (identical(zs$link, "sd")) exp(eta) else sqrt(exp(eta))

  ## Moments of u given sigma_u.
  if (identical(zs$family, "halfnormal")) {
    e_u <- sigma_u * sqrt(2 / pi)
    var_u <- sigma_u^2 * (1 - 2 / pi)
  } else if (identical(zs$family, "exponential")) {
    e_u <- sigma_u
    var_u <- sigma_u^2
  } else {
    stop("unsupported inefficiency family: ", zs$family, call. = FALSE)
  }

  ## d sigma_u / d z_k = scale_k * sigma_u, with scale_k = delta_k under the
  ## SD link and delta_k/2 under the variance link.
  half <- if (identical(zs$link, "sd")) 1 else 0.5

  keep <- vapply(seq_len(ncol(Z)), function(j) length(unique(Z[, j])) > 1L, logical(1))
  if (!any(keep)) {
    stop("every column of the variance-determinant design is constant, so ",
      "there is no covariate to differentiate with respect to.",
      call. = FALSE
    )
  }
  nm <- colnames(Z)[keep]
  d <- delta[keep]

  me_e <- outer(e_u, half * d)
  me_v <- outer(var_u, 2 * half * d)
  ## Name the columns after the component actually differentiated, so a u table
  ## and an h table cannot be confused once separated from their call.
  colnames(me_e) <- paste0("dE_", component, ".d", nm)
  colnames(me_v) <- paste0("dVar_", component, ".d", nm)

  out <- data.frame(sigma_u, e_u, var_u, me_e, me_v, check.names = FALSE)
  names(out)[1:3] <- paste0(c("sigma_", "E_", "Var_"), component)

  avg <- colMeans(cbind(me_e, me_v))
  attr(out, "average") <- avg
  attr(out, "link") <- zs$link
  attr(out, "family") <- zs$family
  attr(out, "component") <- component

  if (isTRUE(average)) {
    return(avg)
  }
  out
}
