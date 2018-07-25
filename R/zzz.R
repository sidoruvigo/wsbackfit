.onAttach <- function(libname, pkgname) {
  packageStartupMessage("Weighted Smooth Backfitting for Structured Models")
}


# Opciones
#===================
.onLoad <- function(libname, pkgname) {

  requireNamespace("wsbackfit", quietly = TRUE)
}
