## Goodness of fit for the assumed inefficiency distribution.
## Wang, Amsler and Schmidt (2011), J Prod Anal 35:95-118.
## See notes/code_history/gof_test.md.
##
## The distribution of u is the assumption in this model that is least often
## defended and least often tested. It is testable: hold the normality of v
## fixed, and the assumed u implies a distribution for the composed error, so a
## rejection of that is a rejection of the assumed u.
##
## Wang, Amsler and Schmidt make the case for testing on eps rather than on
## u-hat: u-hat = E[u | eps] is a MONOTONIC function of eps, so the KS test is
## identical either way and the chi-square test is identical when the cells are
## defined conformably -- but eps is far easier, and the distribution of u-hat
## is not the distribution of u, which is the trap the paper opens by warning
## about. Comparing the observed spread of u-hat with the assumed density of u
## is a mistake, not a diagnostic.
##
## Both statistics are computed from the composed-error CDF, which is why this
## needed pcomposed_model() first.

## Two-sided KS distance between the PIT values and uniform.
.gof_ks <- function(pit) {
  n <- length(pit)
  p <- sort(pit)
  max(max((1:n) / n - p), max(p - (0:(n - 1)) / n))
}

## Pearson statistic on EQUIPROBABLE cells: bin the PIT on (0,1) into k equal
## bins, so every expected count is n/k by construction and no cell can be
## sparse enough to invalidate the approximation.
.gof_chisq <- function(pit, cells) {
  n <- length(pit)
  O <- tabulate(pmin(cells, pmax(1L, ceiling(pit * cells))), nbins = cells)
  E <- n / cells
  sum((O - E)^2 / E)
}

#' Goodness-of-fit test for the assumed inefficiency distribution
#'
#' @param object An `"sfareg"` fit from [sfm()].
#' @param data The data the model was fitted to (`sfm()` does not retain it).
#' @param test `"ks"`, `"chisq"`, or both.
#' @param null How to get the null distribution: `"bootstrap"` (the parametric
#'   bootstrap of Wang, Amsler and Schmidt) or `"asymptotic"`.
#' @param B Bootstrap replications.
#' @param cells Number of equiprobable cells for the chi-square test.
#' @param seed Optional seed for the bootstrap.
#' @return A data frame with one row per test.
#' @export
gof_test <- function(object, data = NULL, test = c("ks", "chisq"),
                     null = c("bootstrap", "asymptotic"), B = 199,
                     cells = 10, seed = NULL) {
  if (!inherits(object, "sfareg")) {
    stop("`object` must be an \"sfareg\" fit.", call. = FALSE)
  }
  test <- match.arg(test, several.ok = TRUE)
  null <- match.arg(null)
  mn <- object$model_name
  par <- object$out[, "par"]
  inefdec <- .sfa_inefdec(object)
  if (is.null(data)) {
    data <- tryCatch(eval(object$call$data, envir = parent.frame()),
      error = function(e) NULL
    )
  }
  if (is.null(data)) {
    stop("gof_test() needs the data the model was fitted to, and sfm() does ",
      "not retain it. Pass `data =`.",
      call. = FALSE
    )
  }
  data <- as.data.frame(data)
  eps <- .sfa_eps_hat(object, data)
  n <- length(eps)
  if (length(cells) != 1L || !is.finite(cells) || cells < 3 || cells > n / 5) {
    stop("gof_test(): `cells` must be a single number in [3, n/5].", call. = FALSE)
  }

  pit <- pcomposed_model(eps, mn, par, inefdec = inefdec)
  stat <- c(ks = .gof_ks(pit), chisq = .gof_chisq(pit, cells))[test]

  ## Number of estimated parameters that shape the composed error: everything
  ## except the frontier slopes. The intercept counts, because it locates eps.
  m <- .n_shape_par(mn) + 1L

  if (identical(null, "asymptotic")) {
    p <- rep(NA_real_, length(stat))
    names(p) <- names(stat)
    if ("chisq" %in% test) {
      df <- cells - 1L - m
      if (df < 1L) {
        stop("gof_test(): ", cells, " cells leaves ", df, " degrees of freedom ",
          "after ", m, " estimated parameters. Use more cells.",
          call. = FALSE
        )
      }
      p[["chisq"]] <- stats::pchisq(stat[["chisq"]], df = df, lower.tail = FALSE)
    }
    if ("ks" %in% test) {
      ## Deliberately NOT the Kolmogorov distribution: with theta estimated the
      ## statistic is stochastically smaller, so those critical values give a
      ## conservative test of unknown size. Bai (2003) fixes it properly; the
      ## bootstrap fixes it here.
      warning("gof_test(): there is no valid asymptotic null for the KS ",
        "statistic once the parameters are estimated (Bai 2003 is not ",
        "implemented). Returning NA; use null = \"bootstrap\".",
        call. = FALSE
      )
    }
    out <- data.frame(
      test = names(stat), statistic = unname(stat),
      null = ifelse(names(stat) == "chisq", paste0("chi2(", cells - 1L - m, ")"), "none"),
      p.value = unname(p), B = NA_integer_, stringsAsFactors = FALSE
    )
    attr(out, "note") <- paste0(
      "The chi-square p-value uses the MLE with (k-1-m) degrees of freedom, ",
      "which Wang, Amsler and Schmidt note is CONSERVATIVE: that reference ",
      "distribution belongs to the minimum-chi-square estimator, not the MLE. ",
      "Prefer null = \"bootstrap\"."
    )
    rownames(out) <- NULL
    return(out)
  }

  ## Parametric bootstrap. The point is to copy the estimation step exactly:
  ## each replication re-estimates and forms ITS OWN residuals, so the effect
  ## of using eps-hat rather than eps is reproduced rather than assumed away.
  if (!is.null(seed)) {
    st <- .rng_snapshot()
    on.exit(.rng_restore(st), add = TRUE)
    set.seed(seed)
  }
  fx <- stats::formula(Formula::Formula(object$formula), lhs = 1, rhs = 1)
  yname <- all.vars(fx)[1]
  X <- stats::model.matrix(fx, data = data)
  bn <- intersect(colnames(X), names(par))
  fitted_part <- as.numeric(X[, bn, drop = FALSE] %*% par[bn])
  sgn <- if (inefdec) 1 else -1

  boot <- matrix(NA_real_, B, length(stat), dimnames = list(NULL, names(stat)))
  for (b in seq_len(B)) {
    db <- data
    db[[yname]] <- fitted_part + sgn * .composed_rgen(n, mn, par)
    fb <- tryCatch(suppressWarnings(sfm(object$formula,
      model_name = mn, data = db, inefdec = inefdec
    )), error = function(e) NULL)
    if (is.null(fb)) next
    eb <- tryCatch(.sfa_eps_hat(fb, db), error = function(e) NULL)
    if (is.null(eb)) next
    pb <- tryCatch(pcomposed_model(eb, mn, fb$out[, "par"], inefdec = inefdec),
      error = function(e) NULL
    )
    if (is.null(pb) || !all(is.finite(pb))) next
    boot[b, ] <- c(ks = .gof_ks(pb), chisq = .gof_chisq(pb, cells))[names(stat)]
  }
  ok <- colSums(!is.na(boot))
  ## The +1s make the p-value the standard bootstrap one, which can never be
  ## exactly zero -- a p-value of 0 from B draws is a statement B cannot support.
  p <- vapply(seq_along(stat), function(j) {
    bj <- boot[, j][!is.na(boot[, j])]
    if (!length(bj)) NA_real_ else (1 + sum(bj >= stat[j])) / (1 + length(bj))
  }, numeric(1))

  out <- data.frame(
    test = names(stat), statistic = unname(stat),
    null = "parametric bootstrap", p.value = p, B = unname(ok),
    stringsAsFactors = FALSE
  )
  attr(out, "boot") <- boot
  rownames(out) <- NULL
  out
}

## How many parameters shape the composed error, apart from the frontier.
.n_shape_par <- function(model_name) {
  switch(model_name,
    "NHN" = , "NE" = , "NR" = , "NU" = , "NGE" = 2L,
    "NTN" = , "NLN" = , "NW" = , "NG" = , "NNAK" = , "TSL" = , "THT" = , "tHN" = 3L,
    2L
  )
}
