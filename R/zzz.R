# global reference to scvelo (will be initialized in .onLoad)
#

# .onAttach <- function(libname, pkgname){
#   suppressPackageStartupMessages()
# }

scvelo <- NULL
matplotlib <- NULL
pyplot <- NULL
.onLoad <- function(libname, pkgname) {
    # reticulate::configure_environment(pkgname, force = TRUE)
    # use superassignment to update global reference to scvelo
    scvelo <<- reticulate::import("scvelo", delay_load = TRUE)
    matplotlib <<- reticulate::import("matplotlib", convert = TRUE)
    matplotlib$use("Agg", force = TRUE)
    pyplot <<- reticulate::import("matplotlib.pyplot", delay_load = TRUE)

    loadRData <- function(fileName){
      #loads an RData file, and returns it
      load(fileName)
      get(ls()[ls() != "fileName"])
    }

    human_gene_transcript_seu <<- loadRData(url("http://cobrinik-1.saban-chla.usc.edu/human_gene_transcript_seu.rda"))
    human_gene_transcript_sce <<- loadRData(url("http://cobrinik-1.saban-chla.usc.edu/human_gene_transcript_sce.rda"))
    human_fetal_retina_sce <<- readRDS(url("http://cobrinik-1.saban-chla.usc.edu/chevreuldata/human_fetal_retina_sce.rds"))
    human_fetal_retina_seu <<- readRDS(url("http://cobrinik-1.saban-chla.usc.edu/chevreuldata/human_fetal_retina_seu.rds"))
    baron2016singlecell <<- loadRData(url("http://cobrinik-1.saban-chla.usc.edu/chevreuldata/baron2016singlecell.rda"))
    panc8 <<- loadRData(url("http://cobrinik-1.saban-chla.usc.edu/chevreuldata/panc8.rda"))
    human_count <<- loadRData(url("http://cobrinik-1.saban-chla.usc.edu/human_count.rda"))
    human_meta <<- loadRData(url("http://cobrinik-1.saban-chla.usc.edu/human_meta.rda"))
}

#' Install scvelo
#'
#' @param method deterministic, stochastic, or dynamical
#' @param conda conda env
#'
#' @return
#' @export
#'
#' @examples
#' \donttest{install_scvelo()}
install_scvelo <- function(method = "auto", conda = "auto") {
    reticulate::py_install("scvelo", method = method, conda = conda, pip = TRUE)
    reticulate::py_install("matplotlib", method = method, conda = conda)
}
