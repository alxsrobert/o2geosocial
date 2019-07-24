#' @importFrom Rcpp evalCpp
#' @useDynLib measlesoutbreaker, .registration = TRUE
#'
.onAttach <- function(libname, pkgname) {
  ## pkg_version <- packageDescription("measlesoutbreaker", fields = "Version")
  ## startup_txt <- paste("\n   === measlesoutbreaker", pkg_version, "is loaded ===\n")
  ## packageStartupMessage(startup_txt)
}
