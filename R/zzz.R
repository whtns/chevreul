# global reference to scvelo (will be initialized in .onLoad)
scvelo <- NULL
matplotlib <- NULL
pyplot <- NULL
.onLoad <- function(libname, pkgname) {
  # reticulate::configure_environment(pkgname, force = TRUE)
  # use superassignment to update global reference to scvelo
  # scvelo <<- reticulate::import("scvelo", delay_load = TRUE)
  # matplotlib <<- reticulate::import("matplotlib", convert = TRUE)
  # matplotlib$use("Agg", force = TRUE)
  # pyplot <<- reticulate::import("matplotlib.pyplot")
}

install_scvelo <- function(method = "auto", conda = "auto") {
  reticulate::py_install("scvelo", method = method, conda = conda, pip = TRUE)
  reticulate::py_install("matplotlib", method = method, conda = conda)
}
