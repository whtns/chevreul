# global reference to scvelo (will be initialized in .onLoad)
#

# .onAttach <- function(libname, pkgname){
#   suppressPackageStartupMessages()
# }

# scvelo <- NULL
# matplotlib <- NULL
# pyplot <- NULL
.onLoad <- function(libname, pkgname) {
    loadRData <- function(fileName) {
        # loads an RData file, and returns it
        load(fileName)
        get(ls()[ls() != "fileName"])
    }
}
