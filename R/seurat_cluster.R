

#' Preprocess Seurat Object
#'
#' @param seu
#' @param scale
#'
#' @return
#' @export
#'
#' @examplesse
seurat_preprocess <- function(seu, scale=TRUE){
  # Normalize data
  seu <- NormalizeData(object = seu, verbose = FALSE)

  # Filter out only variable genes
  seu <- FindVariableFeatures(object = seu, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

  # Regress out unwanted sources of variation
  if (scale){
    seu <- ScaleData(object = seu, features = rownames(x = seu))

  }

  return(seu)
}

#' Find All Markers at a range of resolutions
#'
#' @param seu
#'
#' @return
#' @export
#'
#' @examples
find_all_markers <- function(seu){
  clusters <- paste0("clusters_", seq(0.2, 2.0, by = 0.2)) %>%
    set_names(.)

  marker_features <- purrr::map(clusters, stash_marker_features, seu)

  seu@misc$markers <- marker_features
  return(seu)

}

#' Stash Marker Genes in a Seurat Object
#'
#' Marker Genes will be stored in slot `@misc$markers`
#'
#' @param resolution
#' @param seu
#'
#' @return
#' @export
#'
#' @examples
stash_marker_features <- function(resolution, seu){

  Idents(seu) <- resolution

  markers <- Seurat::FindAllMarkers(object = seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25) %>%
    dplyr::group_by(cluster) %>%
    dplyr::top_n(n = 2, wt = avg_logFC) %>%
    dplyr::pull(gene)

  return(markers)

}


