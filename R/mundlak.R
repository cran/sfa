## Mundlak / Mundlak-Maddala adjustment terms for the true effects panel
## models. Gap L10, from Karagiannis and Kellermann (2019), J Prod Anal
## 51:175-187. See notes/code_history/mundlak.md.
##
## The true random effects model of Greene (2005) assumes the firm effect is
## INDEPENDENT of the regressors. When it is not -- and in production data it
## usually is not, because more capable managers choose different input levels
## -- the effect is partly absorbed into the slopes and partly into
## inefficiency, so the technology parameters are biased and inefficiency is
## understated.
##
## Mundlak's (1978) device is to admit the correlation and model it rather than
## assume it away:
##
##   alpha_i = pi' xbar_i + delta_i
##
## with xbar_i the firm's own mean of the regressors. Substituting into the
## frontier leaves a model that is ALGEBRAICALLY the same likelihood with a
## wider design matrix -- which is why this is a covariate-block change and not
## a new model_name.
##
## Karagiannis and Kellermann widen it further. Their Mundlak-Maddala form
##
##   alpha_i = pi' xbar*_i + gamma' z_i + delta_i,   x*_it = (x_it, z_it)
##
## adds ENVIRONMENTAL variables: things outside the producer's control that
## shape the conditions it works in -- agro-ecological conditions, topography,
## population density, the geography an airline serves. Those enter as
## ordinary time-invariant regressors, so `mundlak` here covers the piece a
## user cannot easily build by hand: the within-firm means, computed against
## the same panel index the fit uses.
##
## Their third variant (Hausman-Taylor) lets only a SUBSET of the regressors
## correlate with the effect, on the argument that the effect is plausibly
## uncorrelated with, say, a neutral time trend. That is why `mundlak` takes a
## formula naming variables rather than a single TRUE/FALSE: `~ .` gives the
## Mundlak-Maddala case, and naming a subset gives the Hausman-Taylor one.

## Frontier regressors, i.e. the variables in the first pipe segment. The later
## segments parameterize variances and are not part of f(x), so averaging them
## would be meaningless.
.frontier_vars <- function(formula) {
  s <- paste(deparse(formula), collapse = " ")
  rhs <- sub("^[^~]*~", "", s)
  rhs1 <- strsplit(rhs, "|", fixed = TRUE)[[1L]][1L]
  all.vars(stats::as.formula(paste("~", rhs1)))
}

## Insert new terms into the FIRST pipe segment, leaving any variance segments
## untouched. Done on the deparsed string because `Formula` has no insertion
## operator and rebuilding a multi-part formula term by term loses the
## environment; the string form keeps the pipes exactly where they were.
.add_to_frontier <- function(formula, new_terms) {
  if (!length(new_terms)) return(formula)
  s <- paste(deparse(formula), collapse = " ")
  lhs <- sub("~.*$", "", s)
  rhs <- sub("^[^~]*~", "", s)
  parts <- strsplit(rhs, "|", fixed = TRUE)[[1L]]
  parts[1L] <- paste(parts[1L], "+", paste(new_terms, collapse = " + "))
  out <- stats::as.formula(
    paste0(lhs, "~", paste(parts, collapse = "|")),
    env = environment(formula)
  )
  out
}

## Group means of `vars` within `id`, appended to `data` under a prefix.
##
## A variable that is ALREADY constant within firm gets no mean column: its
## own mean is itself, so including both would be exactly collinear and the
## design matrix would drop rank. Karagiannis and Kellermann's environmental
## z_i are exactly this case -- they enter the auxiliary equation directly, not
## through a mean -- so this is the specification, not a convenience.
.mundlak_augment <- function(data, id, vars, prefix = "mbar_") {
  keep <- character(0)
  dropped <- character(0)
  for (v in vars) {
    x <- data[[v]]
    if (is.null(x) || !is.numeric(x)) {
      dropped <- c(dropped, v)
      next
    }
    m <- stats::ave(x, id, FUN = function(z) mean(z, na.rm = TRUE))
    if (isTRUE(all.equal(as.numeric(m), as.numeric(x)))) {
      dropped <- c(dropped, v)
      next
    }
    nm <- paste0(prefix, v)
    data[[nm]] <- m
    keep <- c(keep, nm)
  }
  list(data = data, added = keep, dropped = dropped)
}

## The whole transform: returns the widened formula and data, or the originals
## when `mundlak` is NULL.
.apply_mundlak <- function(formula, data, id, mundlak, verbose = FALSE) {
  if (is.null(mundlak)) {
    return(list(formula = formula, data = data, added = character(0)))
  }
  if (!inherits(mundlak, "formula")) {
    stop("`mundlak` must be a formula naming the variables whose within-firm ",
      "means enter the frontier, such as ~ x1 + x2, or ~ . for all of them.",
      call. = FALSE
    )
  }
  fv <- .frontier_vars(formula)
  rhs <- paste(deparse(mundlak[[length(mundlak)]]), collapse = " ")
  vars <- if (identical(trimws(rhs), ".")) fv else all.vars(mundlak)

  unknown <- setdiff(vars, names(data))
  if (length(unknown)) {
    stop("`mundlak` names variable(s) not in `data`: ",
      paste(unknown, collapse = ", "), ".",
      call. = FALSE
    )
  }
  not_frontier <- setdiff(vars, fv)
  if (length(not_frontier)) {
    warning("`mundlak` names variable(s) that are not frontier regressors: ",
      paste(not_frontier, collapse = ", "),
      ". Their within-firm means are still added, but Mundlak's device is ",
      "about the regressors the firm effect is correlated WITH.",
      call. = FALSE
    )
  }

  aug <- .mundlak_augment(data, id, vars)
  if (!length(aug$added)) {
    stop("`mundlak`: none of the named variables vary within firm, so their ",
      "group means are the variables themselves and would be exactly ",
      "collinear with them. Time-invariant environmental factors belong in ",
      "the formula directly.",
      call. = FALSE
    )
  }
  if (length(aug$dropped) && isTRUE(verbose)) {
    message("mundlak: no mean added for ", paste(aug$dropped, collapse = ", "),
      " (constant within firm, or not numeric)."
    )
  }
  list(
    formula = .add_to_frontier(formula, aug$added),
    data = aug$data,
    added = aug$added,
    dropped = aug$dropped
  )
}
