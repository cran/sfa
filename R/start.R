.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "sfa version ", utils::packageVersion("sfa"), "\n",
    "Type citation('sfa') for citing this package in publications."
  )
}

## sandwich is in Suggests, so its generics cannot be imported. Register the
## methods at load time when it is present -- the standard pattern for making
## an optional package's generics work without depending on it.
.onLoad <- function(libname, pkgname) {
  if (requireNamespace("sandwich", quietly = TRUE)) {
    registerS3method("bread", "sfareg", bread.sfareg,
                     envir = asNamespace("sandwich"))
    registerS3method("estfun", "sfareg", estfun.sfareg,
                     envir = asNamespace("sandwich"))
  }
  invisible()
}
