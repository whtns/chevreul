#' Regress Seurat Object by Given Set of Genes
#'
#' @param seu_list
#' @param gene_set as a length 1 list containing a character vector of gene symbols
#' @param set_name as a string
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
regress_by_features <- function(seu, feature_set, set_name, regress = TRUE, ...) {
  message(paste0("regressing seurat objects by ", set_name))
  message(paste0("Module score stored as ", set_name, "1"))

  if (!is.list(feature_set)) feature_set <- list(feature_set)

  #set default assay to "RNA"
  DefaultAssay(seu) <- "RNA"

  ctrl <- 100
  if (dim(seu)[2] < 100){
    ctrl <- dim(seu)[2]/10
  }

  seu <- AddModuleScore(seu, feature_set, name = set_name, ctrl = ctrl)

  if (any(grepl("integrated", names(seu[[]])))){
    default_assay = "integrated"
  } else {
    default_assay = "RNA"
  }

  #revert default assay to "integrated"
  DefaultAssay(seu) <- default_assay

  #regress out feature set
  if (regress){
    seu <- ScaleData(seu, vars.to.regress = paste0(set_name, "1"))
    # rerun dimensional reduction
    reductions <- names(seu@reductions)
    resolutions <- stringr::str_extract(names(seu[[]])[grepl("snn", names(seu[[]]))], "[0-9].*$")
    resolutions <- sort(as.numeric(resolutions))

    seu <- seurat_reduce_dimensions(seu)
    seu <- seurat_cluster(seu, resolution = resolutions)
  }

  return(seu)
}
