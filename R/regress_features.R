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
#' regressed_object <- regress_by_features(human_gene_transcript_sce, feature_set = cc.genes$s.genes, set_name = "s_genes")
#'
#'@importFrom scran modelGeneVar
#'@importFrom batchelor regressBatches
regress_by_features <- function (object, feature_set, set_name, regress = TRUE)
          {
            message(paste0("regressing objects by ", set_name))
            if (!is.list(feature_set))
              feature_set <- list(feature_set)
            # ctrl <- 100
            if (dim(object)[2] < 100) {
              ctrl <- dim(object)[2]/10
            }

            # object <- AddModuleScore(object, feature_set, name = set_name, ctrl = ctrl)
            set_name <- paste0(set_name, length(feature_set))
            message(paste0("Module score stored as ", set_name))
            # if ("integrated" %in% names(object@experiments)) {
            #   default_experiment <- "integrated"
            # }
            # else {
            #   default_experiment <- "gene"
            # }
            # SingleCellExperiment::DefaultAssay(object) <- default_experiment
            if (regress) {
              dec.nocycle <- modelGeneVar(object, block=colData(object)[["Phase"]])
              reg.nocycle <- regressBatches(object, batch=colData(object)[["Phase"]])
              # object <- ScaleData(object, vars.to.regress = set_name)
              reductions <- reducedDimNames(object)
              resolutions <- stringr::str_extract(names(get_cell_metadata(object))[grepl("snn", names(get_cell_metadata(object)))], "[0-9].*$")
              resolutions <- sort(unique(as.numeric(resolutions)))
              object <- object_reduce_dimensions(object)
              object <- object_cluster(object, resolution = resolutions)
            }
            return(object)
          }
