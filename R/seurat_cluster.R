
#' Preprocess Seurat Object
#'
#' @param seu
#' @param assay
#' @param scale
#' @param normalize
#' @param features
#' @param legacy_settings
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#'
#' panc8[["gene"]] <- seurat_preprocess(panc8[["gene"]])
seurat_preprocess <- function(assay, scale = TRUE, normalize = TRUE, features = NULL, legacy_settings = FALSE, ...) {

  # Normalize data

  if (legacy_settings) {
    message("using legacy settings")

    logtransform_exp <- as.matrix(log1p(Seurat::GetAssayData(assay)))

    seu <- Seurat::SetAssayData(assay, slot = "data", logtransform_exp) %>%
      Seurat::ScaleData(features = rownames(.))

    return(assay)
  }

  if (normalize) {
    assay <- Seurat::NormalizeData(assay, verbose = FALSE, ...)
  }

  # Filter out only variable genes
  assay <- Seurat::FindVariableFeatures(assay, selection.method = "vst", verbose = FALSE, ...)

  # Regress out unwanted sources of variation
  if (scale) {
    assay <- Seurat::ScaleData(assay, features = rownames(assay), ...)
  }

  return(assay)
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
#' markers_stashed_seu <- find_all_markers(panc8)
#' marker_genes <- Misc(markers_stashed_seu, "markers")
#' str(marker_genes)
find_all_markers <- function(seu, metavar = NULL, seurat_assay = "gene", ...) {

  if (is.null(metavar)) {
    resolutions <- colnames(seu[[]])[grepl(paste0(seurat_assay, "_snn_res."), colnames(seu[[]]))]

    cluster_index <- grepl(paste0(seurat_assay, "_snn_res."), colnames(seu[[]]))

    if (!any(cluster_index)) {
      warning("no clusters found in metadata. runnings seurat_cluster")
      seu <- seurat_cluster(seu, resolution = seq(0.2, 2.0, by = 0.2))
    }

    clusters <- seu[[]][, cluster_index]

    cluster_levels <- purrr::map_int(clusters, ~ length(unique(.x)))
    cluster_levels <- cluster_levels[cluster_levels > 1]

    clusters <- dplyr::select(clusters, dplyr::one_of(names(cluster_levels)))
    metavar <- names(clusters)
  }

  new_markers <- purrr::map(metavar, stash_marker_features, seu, seurat_assay = seurat_assay)
  names(new_markers) <- metavar

  old_markers <- seu@misc$markers[!names(seu@misc$markers) %in% names(new_markers)]

  seu@misc$markers <- c(old_markers, new_markers)

  return(seu)
}


#' Stash Marker Genes in a Seurat Object
#'
#' Marker Genes will be stored in slot `@misc$markers`
#'
#' @param metavar
#' @param seu
#' @param seurat_assay
#' @param top_n
#'
#' @return
#'
#' @examples
#'
#' seu <- stash_marker_features(metavar = "batch", seu, seurat_assay = "gene")
#'
stash_marker_features <- function(metavar, seu, seurat_assay, top_n = 200) {

  markers <- list()
  markers$presto <- presto::wilcoxauc(seu, metavar, seurat_assay = seurat_assay) %>%
    dplyr::group_by(group) %>%
    dplyr::filter(padj < 0.05) %>%
    dplyr::top_n(n = top_n, wt = logFC) %>%
    dplyr::arrange(group, desc(logFC)) %>%
    dplyr::select(feature, group) %>%
    dplyr::mutate(rn = row_number()) %>%
    tidyr::pivot_wider(names_from = group, values_from = feature) %>%
    dplyr::select(-rn)

  markers$genesorteR <- tryCatch(
    {
      genesorter_results <- genesorteR::sortGenes(
        Seurat::GetAssayData(seu, assay = seurat_assay, slot = "data"),
        seu[[]][[metavar]]
      )

      genesorter_results <- apply(genesorter_results$specScore, 2, function(x) names(head(sort(x, decreasing = TRUE), n = top_n))) %>%
        as.data.frame()
    },
    error = function(e) {
      message(sprintf("Error in %s: %s", deparse(e[["call"]]), e[["message"]]))
      NULL
    },
    finally = {
    }
  )


  return(markers)
}
