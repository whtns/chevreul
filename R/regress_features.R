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
regress_by_gene_set <- function(seu_list, gene_set, set_name, ...) {
  message(paste0("regressing seurat objects by ", set_name))
  message(paste0("Module score stored as ", set_name, "1"))

  transcript_set <- dplyr::left_join(tibble::enframe(unlist(gene_set)), annotables::grch38, by = c("value" = "symbol")) %>%
    dplyr::left_join(annotables::grch38_tx2gene, by = "ensgene") %>%
    dplyr::pull(enstxp)

  transcript_set <- list(transcript_set)

  #set default assay to "RNA"
  seu_list <- purrr::map(seu_list, SetDefaultAssay, "RNA")

  seu_list <- purrr::map2(seu_list, list(gene_set, transcript_set),  ~AddModuleScore(.x, .y, name = set_name))

  #revert default assay to "integrated"
  seu_list <- map(seu_list, SetDefaultAssay, "integrated")

  #regress out stress score
  seu_list <- purrr::map(seu_list, ScaleData, vars.to.regress = paste0(set_name, "1"))

  # rerun dimensional reduction
  seu_list <- purrr::map(seu_list, seurat_reduce_dimensions)

  seu_list <- purrr::map(seu_list, seurat_cluster, resolution = seq(0.2, 2.0, 0.2))
}
