

#' Preprocess Seurat Object
#'
#' @param seu
#' @param scale
#'
#' @return
#' @export
#'
#' @examples

seurat_preprocess <- function(seu, scale=TRUE, normalize = TRUE, ...){
  # Normalize data

  if (normalize){
    seu <- Seurat::NormalizeData(object = seu, verbose = FALSE)
  }

  # Filter out only variable genes
  seu <- Seurat::FindVariableFeatures(object = seu, selection.method = "vst", nfeatures = 2000, verbose = FALSE)

  # Regress out unwanted sources of variation
  if (scale){
    seu <- Seurat::ScaleData(object = seu, features = rownames(x = seu), ...)

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
find_all_markers <- function(seu, ...){
  # browser()

  if (any(grepl("integrated", colnames(seu[[]])))){
    default_assay = "integrated"
  } else {
    default_assay = "RNA"
  }


  resolutions <- colnames(seu[[]])[grepl(paste0(default_assay, "_snn_res."), colnames(seu[[]]))]

  cluster_index <- grepl(paste0(default_assay, "_snn_res."), colnames(seu[[]]))

  if(!any(cluster_index)) {
    stop("no clusters found in metadata. Please run seurat_cluster")
  }

  clusters <- seu[[]][,cluster_index]

  cluster_levels <- purrr::map_int(clusters, ~length(unique(.x)))
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
    dplyr::top_n(n = 200, wt = logFC) %>%
    # dplyr::pull(feature) %>%
    identity()

  return(markers)

}

#' Convenience Function to transition existing seus
#'
#' @param seu
#'
#' @return
#' @export
#'
#' @examples
rename_markers <- function(seu){
  default_assay <- DefaultAssay(seu)
  test <- names(seu@misc$markers)
  names(seu@misc$markers) <- gsub("clusters_", paste0(default_assay, "_snn_res."), test)
  return(seu)
}

#' Convenience function to save all seurat objects of a given project by feature and suffix
#'
#' @param seu_list list of names seurat objects by suffix then feature; can generate with load_seurat_from_proj and purrr
#' @param proj_dir
#'
#' @return
#' @export
#'
#' @examples
save_proj_feature_seus <- function(seu_list, proj_dir){
  # browser()
  suffixes <- names(seu_list)
  features <- names(seu_list[[1]])

  args <- expand.grid(features, suffixes)
  seus <- unlist(seu_list)

  args <- dplyr::mutate(args, seu = seus) %>%
    dplyr::rename(feature = "Var1", suffix = "Var2")

  purrr::pmap(args, save_seurat, proj_dir = proj_dir)

  return(args)

}

#' Convenience function to replace markers in seurat objects
#'
#' @param suffixes
#' @param features
#' @param proj_dir
#'
#' @return
#' @export
#'
#' @examples
replace_markers_in_proj <- function(suffixes, features, proj_dir){
  names(suffixes) <- suffixes

  names(suffixes)[1] <- "unfiltered"

  seus <- purrr::map(suffixes, ~load_seurat_from_proj(proj_dir, features = features, suffix = .x))

  seus <- purrr::map(seus, ~purrr::map(.x, find_all_markers))

  save_proj_feature_seus(seus, proj_dir)
}

