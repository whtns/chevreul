#' Regress SingleCellExperiment Object by Given Set of Genes
#'
#' @param object A object
#'
#' @return a SingleCellExperiment object with features regressed
#' @export
#' @examples
#' chevreul_sce <- chevreuldata::human_gene_transcript_sce()
#' data(cc.genes.cyclone)
#' regressed_object <- regress_cell_cycle(chevreul_sce)
#'
regress_cell_cycle <- function(object) {
    message("regressing objects by cell cycle")

    if (dim(object)[2] < 100) {
        ctrl <- dim(object)[2] / 10
    }
    if (!"Phase" %in% colnames(colData(object))) {
        object <- annotate_cell_cycle(object)
    }
    if (!any(str_detect(c(mainExpName(object), altExpNames(object)), pattern = ".*regress.*"))) {
        dec.nocycle <- modelGeneVar(object, block = colData(object)[["Phase"]])
        reg.nocycle <- regressBatches(object, batch = colData(object)[["Phase"]])

        reg.nocycle <- runPCA(reg.nocycle,
            exprs_values = "corrected",
            subset_row = getTopHVGs(dec.nocycle, prop = 0.1)
        )

        original_experiment <- mainExpName(object)

        altExp(object, glue("{original_experiment}_regressed")) <- reg.nocycle
        # browser()
        object <- swapAltExp(object, as.character(glue("{original_experiment}_regressed")))

        reductions <- reducedDimNames(object)
        resolutions <- str_extract(colnames(colData(object))[grepl(glue("{original_experiment}_snn_res."), colnames(colData(object)))], "[0-9].*$")
        object <- runTSNE(x = object, dimred = "PCA", n_dimred = seq(30))
        object <- runUMAP(x = object, dimred = "PCA", n_dimred = seq(30))
        object <- object_cluster(object, resolution = resolutions)
    }
    return(object)
}
