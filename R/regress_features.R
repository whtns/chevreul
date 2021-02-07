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
#'
#' regressed_seu <- regress_by_features(panc8, feature_set  = cc.genes$s.genes, set_name = "s_genes")
#'
regress_by_features <- function(seu, feature_set, set_name, regress = TRUE, ...) {
  message(paste0("regressing seurat objects by ", set_name))

  if (!is.list(feature_set)) feature_set <- list(feature_set)

  ctrl <- 100
  if (dim(seu)[2] < 100) {
    ctrl <- dim(seu)[2] / 10
  }

  seu <- AddModuleScore(seu, feature_set, name = set_name, ctrl = ctrl)

  set_name = paste0(set_name, length(feature_set))
  message(paste0("Module score stored as ", set_name))

  if ("integrated" %in% names(seu@assays)) {
    default_assay <- "integrated"
  } else {
    default_assay <- "gene"
  }

  # revert default assay to "integrated"
  Seurat::DefaultAssay(seu) <- default_assay

  # regress out feature set
  if (regress) {
    seu <- ScaleData(seu, vars.to.regress = set_name)
    # rerun dimensional reduction
    reductions <- names(seu@reductions)
    resolutions <- stringr::str_extract(names(seu[[]])[grepl("snn", names(seu[[]]))], "[0-9].*$")
    resolutions <- sort(as.numeric(resolutions))

    seu <- seurat_reduce_dimensions(seu)
    seu <- seurat_cluster(seu, resolution = resolutions)
  }

  return(seu)
}
