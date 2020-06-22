

#' Preprocess Seurat Object
#'
#' @param seu
#' @param scale
#'
#' @return
#' @export
#'
#' @examples

seurat_preprocess <- function(seu, scale=TRUE, normalize = TRUE, features = NULL, ...){
  # Normalize data

  if (normalize){
    seu <- Seurat::NormalizeData(object = seu, verbose = FALSE, ...)
  }

  # Filter out only variable genes
  seu <- Seurat::FindVariableFeatures(object = seu, selection.method = "vst", verbose = FALSE, ...)

  # Regress out unwanted sources of variation
  if (scale){
    seu <- Seurat::ScaleData(object = seu, features = rownames(x = seu), ...)

  }

  return(seu)
}

#' Find All Markers at a range of resolutions
#'
#' @param seu
#' @param metavar
#'
#' @return
#' @export
#'
#' @examples
find_all_markers <- function(seu, metavar = NULL, ...){
  if (is.null(metavar)){
    if ("integrated" %in% names(seu@assays)) {
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
    metavar = names(clusters)
  }

  marker_features <- purrr::map(metavar, stash_marker_features, seu)
  names(marker_features) <- metavar

  new_markers <- marker_features[setdiff(names(marker_features), names(seu@misc$markers))]

  seu@misc$markers <- c(seu@misc$markers, new_markers)

  return(seu)

}


#' Stash Marker Genes in a Seurat Object
#'
#' Marker Genes will be stored in slot `@misc$markers`
#'
#' @param metavar
#' @param seu
#' @param top_n
#'
#' @return
#' @export
#'
#' @examples
stash_marker_features <- function(metavar, seu, top_n = 200){
  markers <- list()
    markers$presto <- presto::wilcoxauc(seu, metavar) %>%
      dplyr::group_by(group) %>%
      dplyr::filter(padj < 0.05) %>%
      dplyr::top_n(n = top_n, wt = logFC) %>%
      dplyr::arrange(group, desc(logFC)) %>%
      dplyr::select(feature, group) %>%
      dplyr::mutate(rn = row_number()) %>%
      tidyr::pivot_wider(names_from = group, values_from = feature) %>%
      dplyr::select(-rn)

    markers$genesorteR <- genesorteR::sortGenes(
      seu@assays$RNA@data,
      seu[[]][[metavar]]
    )

    markers$genesorteR <-
      apply(markers$genesorteR$specScore, 2, function(x) names(head(sort(x, decreasing = TRUE), n = top_n))) %>%
      as.data.frame()

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

