#' @importFrom Rcpp evalCpp
#' @useDynLib o2geosocial, .registration = TRUE
#'
.onAttach <- function(libname, pkgname) {
  ## pkg_version <- packageDescription("o2geosocial", fields = "Version")
  ## startup_txt <- paste("\n   === o2geosocial", pkg_version, "is loaded ===\n")
  ## packageStartupMessage(startup_txt)
}
