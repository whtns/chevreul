#' Regress SingleCellExperiment Object by Given Set of Genes
#'
#' @param object A object
#' @param feature_set a set of features
#' @param set_name as a string
#' @param regress whether to regress
#'
#' @return a single cell object with features regressed
#' @export
#'
#' @examples
#'
#' regressed_object <- regress_cell_cycle(human_gene_transcript_sce)
#'
regress_cell_cycle <- function (object, regress = TRUE)
          {
            message("regressing objects by cell cycle")

            if (dim(object)[2] < 100) {
              ctrl <- dim(object)[2]/10
            }
  if(!"Phase" %in% colnames(colData(object))){
    object <- annotate_cell_cycle(object)
  }

            if (regress) {
              dec.nocycle <- modelGeneVar(object, block=colData(object)[["Phase"]])
              reg.nocycle <- regressBatches(object, batch=colData(object)[["Phase"]])

              reg.nocycle <- runPCA(reg.nocycle, exprs_values="corrected",
                                    subset_row=getTopHVGs(dec.nocycle, prop=0.1))
              mainExpName(object) <- "original"
              altExp(object, "gene") <- reg.nocycle
              original_experiment = mainExpName(object)
              swapAltExp(object, "gene")

              reductions <- reducedDimNames(object)
              resolutions <- str_extract(colnames(colData(object))[grepl(glue("{original_experiment}_snn_res."), colnames(colData(object)))], "[0-9].*$")
              object <- runTSNE(x = object, dimred = "PCA", n_dimred = 1:30)
              object <- runUMAP(x = object, dimred = "PCA", n_dimred = 1:30)
              object <- object_cluster(object, resolution = resolutions)
            }
            return(object)
          }
