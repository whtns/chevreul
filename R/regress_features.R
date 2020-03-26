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
regress_by_features <- function(seu, feature_set, set_name, ...) {
  message(paste0("regressing seurat objects by ", set_name))
  message(paste0("Module score stored as ", set_name, "1"))

  feature_set <- list(feature_set)

  #set default assay to "RNA"
  DefaultAssay(seu) <- "RNA"

  seu <- AddModuleScore(seu, feature_set, name = set_name)

  if (any(grepl("integrated", names(seu[[]])))){
    default_assay = "integrated"
  } else {
    default_assay = "RNA"
  }

  #revert default assay to "integrated"
  DefaultAssay(seu) <- default_assay

  #regress out feature set
  seu <- ScaleData(seu, vars.to.regress = paste0(set_name, "1"))

  # rerun dimensional reduction
  seu <- seurat_reduce_dimensions(seu)

  seu <- seurat_cluster(seu, resolution = seq(0.2, 2.0, 0.2))
}
