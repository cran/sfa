.onAttach <- function(libname, pkgname) {
  packageStartupMessage(
    "sfa version ", utils::packageVersion("sfa"), "\n",
    "Type citation('sfa') for citing this package in publications."
  )
}
