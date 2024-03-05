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

    # chevreuldata::human_gene_transcript_sce <<- loadRData(url("http://cobrinik-1.saban-chla.usc.edu/chevreuldata::human_gene_transcript_sce.rda"))
    # chevreul_sce <<- chevreuldata::chevreuldata::human_gene_transcript_sce()
    # human_gene_transcript_loom_path <<- "http://cobrinik-1.saban-chla.usc.edu/human_gene_transcript.loom"
    # human_fetal_retina_sce <<- readRDS(url("http://cobrinik-1.saban-chla.usc.edu/chevreuldata/human_fetal_retina_sce.rds"))
    # baron2016singlecell <<- loadRData(url("http://cobrinik-1.saban-chla.usc.edu/chevreuldata/baron2016singlecell.rda"))
}
