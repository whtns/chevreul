# global reference to scipy (will be initialized in .onLoad)
scvelo <- NULL

# .onLoad <- function(libname, pkgname) {
#   # use superassignment to update global reference to scipy
#   scvelo <<- reticulate::import("scvelo", delay_load = TRUE)
# }
