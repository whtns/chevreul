#' Regress Seurat Object by Given Set of Genes
#'
#' @param object A object
#' @param feature_set
#' @param set_name as a string
#' @param regress
#'
#' @return
#' @export
#'
#' @examples
#'
#' regressed_object <- regress_by_features(panc8, feature_set = cc.genes$s.genes, set_name = "s_genes")
#'
regress_by_features <- function(object, feature_set, set_name, regress = TRUE, ...) {
    message(paste0("regressing objects by ", set_name))

    if (!is.list(feature_set)) feature_set <- list(feature_set)

    ctrl <- 100
    if (dim(object)[2] < 100) {
        ctrl <- dim(object)[2] / 10
    }

    object <- AddModuleScore(object, feature_set, name = set_name, ctrl = ctrl)

    set_name <- paste0(set_name, length(feature_set))
    message(paste0("Module score stored as ", set_name))

    if ("integrated" %in% names(object@assays)) {
        default_assay <- "integrated"
    } else {
        default_assay <- "gene"
    }

    # revert default assay to "integrated"
    Seurat::DefaultAssay(object) <- default_assay

    # regress out feature set
    if (regress) {
        object <- ScaleData(object, vars.to.regress = set_name)
        # rerun dimensional reduction
        reductions <- names(object@reductions)
        resolutions <- stringr::str_extract(names(pull_metadata(object))[grepl("snn", names(pull_metadata(object)))], "[0-9].*$")
        resolutions <- sort(as.numeric(resolutions))

        object <- object_reduce_dimensions(object)
        object <- object_cluster(object, resolution = resolutions)
    }

    return(object)
}
