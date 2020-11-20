# global reference to scvelo (will be initialized in .onLoad)
scvelo <- NULL
matplotlib <- NULL
pyplot <- NULL
.onLoad <- function(libname, pkgname) {
  # use superassignment to update global reference to scvelo
  scvelo <<- reticulate::import("scvelo", delay_load = TRUE)
  matplotlib <<- reticulate::import("matplotlib", convert = TRUE)
  matplotlib$use("Agg", force = TRUE)
  pyplot <<- reticulate::import("matplotlib.pyplot")
}
