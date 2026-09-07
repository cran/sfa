## efficiency(): the three point predictors on either scale of the dependent
## variable. See notes/code_history/efficiency.md.
##
## Two gaps closed at once, because they are the same decision seen twice: what
## exactly is being reported as "efficiency".
##
##  - H7. The package assumed throughout that the dependent variable is logged,
##    without saying so anywhere in a signature. `logDepVar` makes the
##    assumption visible and gives the level-scale alternative.
##  - G7. Only two of the three usual predictors were available. `type = "mode"`
##    adds the third.

## Mode of the posterior of u given the composed residual. For the models whose
## posterior is a truncated normal N+(mu_star, sigma_star^2) this is simply
## mu_star censored at zero: the untruncated normal's mode is its mean, and
## truncation to the positive half-line moves it to the boundary whenever the
## mean is negative.
##
## The mode is the ONLY one of the three predictors that can be exactly zero,
## and it is zero for every observation whose posterior mean is negative --
## which is most of the efficient tail. That is a property of the estimator,
## not a defect: the most likely value of u really is the boundary there.
.u_mode <- function(mu_star) pmax(mu_star, 0)

efficiency <- function(object,
                       type = c("bc", "jlms", "mode"),
                       logDepVar = TRUE,
                       newdata = NULL) {
  if (!inherits(object, "sfareg")) {
    stop("`object` must be an \"sfareg\" fit.", call. = FALSE)
  }
  type <- match.arg(type)
  if (length(logDepVar) != 1L || !is.logical(logDepVar) || is.na(logDepVar)) {
    stop("`logDepVar` must be TRUE or FALSE.", call. = FALSE)
  }

  ## The point estimate of u implied by each predictor. `bc` is a certainty
  ## equivalent rather than a moment of u, which is exactly why it differs from
  ## `jlms` by Jensen's inequality and not by noise.
  u <- switch(type,
    jlms = {
      ## Prefer the stored value; fall back to rebuilding it from the posterior
      ## with the SAME helper efficiency_ci() uses, so the two cannot drift.
      ## Most cross-sectional fits store exp_u_hat and the posterior but not
      ## E[u | e] itself.
      v <- object$u_hat
      if (is.null(v)) v <- object$jlms
      if (is.null(v)) {
        post <- object$u_posterior
        if (!is.null(post) && !is.null(post$mu_star)) {
          v <- .jlms_u(post$mu_star, post$sigma_star)
        }
      }
      if (is.null(v)) {
        stop("efficiency(type = \"jlms\"): this fit stores neither E[u | e] ",
          "nor the posterior needed to rebuild it.",
          call. = FALSE
        )
      }
      as.numeric(v)
    },
    bc = {
      v <- object$exp_u_hat
      if (is.null(v)) v <- object$efficiency
      if (is.null(v)) {
        stop("efficiency(type = \"bc\"): this fit stores no E[exp(-u) | e]. ",
          "Several models report only the JLMS predictor; use type = \"jlms\".",
          call. = FALSE
        )
      }
      -log(pmax(as.numeric(v), .SFA_CONSTANTS$MIN_POSITIVE))
    },
    mode = {
      post <- object$u_posterior
      if (is.null(post) || is.null(post$mu_star)) {
        stop("efficiency(type = \"mode\"): the mode needs the posterior of u ",
          "given the composed residual, which is a truncated normal only for ",
          "model_name \"NHN\", \"NHN_Z\", \"NE\" and \"NTN\". This fit is ",
          "model_name \"",
          if (is.null(object$model_name)) "unknown" else object$model_name,
          "\", whose posterior has a different shape.",
          call. = FALSE
        )
      }
      .u_mode(as.numeric(post$mu_star))
    }
  )

  if (isTRUE(logDepVar)) {
    ## y in logs: u is a proportional shortfall and TE = exp(-u). For `bc` this
    ## returns E[exp(-u) | e] unchanged, since the log and the exp cancel --
    ## deliberately, so the default path reports the stored number rather than
    ## a round-tripped copy of it.
    return(exp(-u))
  }

  ## y in levels: the frontier is f(x) and TE = (f - u)/f = 1 - u/f. This is a
  ## RATIO of two estimated quantities, so it is less well behaved than the log
  ## version -- it is undefined where the fitted frontier is zero and negative
  ## where u exceeds it, both of which happen near zero output.
  ## .sfa_xb() rebuilds the frontier by re-evaluating the call's `data`
  ## argument, which fails whenever the fit was made inside a function whose
  ## environment has since gone -- a limitation shared by fitted(), residuals()
  ## and predict(). `newdata` is the escape hatch, and the error says so rather
  ## than reporting that the frontier does not exist.
  f <- tryCatch(as.numeric(.sfa_xb(object, newdata)$xb), error = function(e) NULL)
  if (is.null(f) || length(f) != length(u)) {
    stop("efficiency(logDepVar = FALSE): could not rebuild the fitted frontier ",
      "for this fit. This happens when the data used to fit it can no longer be ",
      "recovered from the call -- pass it explicitly as `newdata`.",
      call. = FALSE
    )
  }
  te <- 1 - u / f
  if (any(!is.finite(te)) || any(te < 0, na.rm = TRUE)) {
    warning("efficiency(logDepVar = FALSE): ",
      sum(!is.finite(te) | te < 0, na.rm = TRUE),
      " of ", length(te), " scores are negative or non-finite, which happens ",
      "where the fitted frontier is at or below zero. Check that the dependent ",
      "variable really is on a level scale.",
      call. = FALSE
    )
  }
  te
}
