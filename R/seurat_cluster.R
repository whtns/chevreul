

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
  seu <- Seurat::NormalizeData(object = seu, verbose = FALSE)

  # Filter out only variable genes
  seu <- Seurat::FindVariableFeatures(object = seu, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

  # Regress out unwanted sources of variation
  if (scale){
    seu <- Seurat::ScaleData(object = seu, features = rownames(x = seu))

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
find_all_markers <- function(seu, resolution = 0.6, ...){
  # browser()
  clusters <- paste0("clusters_", resolution) %>%
    set_names(.)

  clusters <- seu[[clusters]]

  cluster_levels <- purrr::map_int(clusters, ~length(levels(.x)))
  cluster_levels <- cluster_levels[cluster_levels > 1]

  clusters <- dplyr::select(clusters, dplyr::one_of(names(cluster_levels)))

  # marker_features <- purrr::map(clusters, mod_stash_marker_features, seu)
  marker_features <- purrr::map(names(clusters), stash_marker_features, seu)
  names(marker_features) <- names(clusters)

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

  markers <- presto::wilcoxauc(seu, resolution) %>%
    dplyr::group_by(group) %>%
    dplyr::top_n(n = 5, wt = logFC) %>%
    dplyr::pull(feature)

  return(markers)

}

