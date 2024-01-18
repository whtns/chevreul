#' Transfer Labels Between Seurat Objects
#'
#' @param ref_seu
#' @param query_seu
#' @param ref_name
#' @param query_name
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
label_transfer <- function(ref_seu, query_seu, ref_name = NULL, query_name = NULL, ...) {
    label_transerred <- purrr::map2(organoid_seu, seu, reference_integrate, query_name = query_name, ref_name = ref_name)
}

#' Transfer labels between gene or transcript objects
#'
#' @param ref_seu
#' @param query_seu
#' @param query_name
#' @param ref_name
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
reference_integrate <- function(ref_seu, query_seu, query_name = "fetal", ref_name = "organoid", ...) {
    seu.anchors <- FindTransferAnchors(reference = ref_seu, query = query_seu, dims = 1:30, project.query = TRUE)

    reference_clusters <- colnames(ref_seu[[]])[grepl(paste0("gene", "_snn_res."), colnames(ref_seu[[]]))]

    refdata <- t(ref_seu[[reference_clusters]])

    cellids <- colnames(refdata)
    refdata <- setNames(split(refdata, seq(nrow(refdata))), rownames(refdata)) %>%
        purrr::map(purrr::set_names, cellids)

    ref_names <- stringr::str_replace(reference_clusters, ".*snn_res", ref_name)

    predictions <- purrr::map(refdata, ~ TransferData(
        anchorset = seu.anchors, refdata = .x,
        dims = 1:30
    ))
    predictions0 <-
        predictions %>%
        purrr::map(as_tibble, rownames = "sample_id") %>%
        purrr::map(~ dplyr::select(.x, sample_id, predicted.id)) %>%
        dplyr::bind_rows(.id = "resolution") %>%
        tidyr::spread(resolution, predicted.id) %>%
        purrr::set_names(c("sample_id", ref_names)) %>%
        tibble::column_to_rownames("sample_id") %>%
        identity()

    query_seu <- AddMetaData(query_seu, predictions0)
}


## ------------------------------------------------------------------------
## find markers for every cluster compared to all remaining cells, report only the positive ones
#' Find Cell Type Markers in a Seurat Object
#'
#' @param seu A seurat object
#'
#' @return
#' @export
#'
#' @examples
seurat_find_markers <- function(seu, num_features) {
    markers <- Seurat::FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    markers %>%
        group_by(cluster) %>%
        top_n(n = 2, wt = avg_logFC)

    cluster_markers <- reference.markers$gene %>%
        group_by(cluster) %>%
        top_n(n = num_features, wt = avg_logFC) %>%
        dplyr::pull(gene)

    DoHeatmap(reference_integrated[[1]], features = cluster_markers)
    return(cluster_markers)
}
