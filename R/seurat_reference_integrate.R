#' Transfer Labels Between Seurat Objects
#'
#' @param ref_object
#' @param query_object
#' @param ref_name
#' @param query_name
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
label_transfer <- function(ref_object, query_object, ref_name = NULL, query_name = NULL, ...) {
    label_transerred <- purrr::map2(organoid_object, object, reference_integrate, query_name = query_name, ref_name = ref_name)
}

#' Transfer labels between gene or transcript objects
#'
#' @param ref_object
#' @param query_object
#' @param query_name
#' @param ref_name
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
reference_integrate <- function(ref_object, query_object, query_name = "fetal", ref_name = "organoid", ...) {
    object.anchors <- FindTransferAnchors(reference = ref_object, query = query_object, dims = 1:30, project.query = TRUE)

    reference_clusters <- colnames(ref_object[[]])[grepl(paste0("gene", "_snn_res."), colnames(ref_object[[]]))]

    refdata <- t(ref_object[[reference_clusters]])

    cellids <- colnames(refdata)
    refdata <- setNames(split(refdata, seq(nrow(refdata))), rownames(refdata)) %>%
        map(purrr::set_names, cellids)

    ref_names <- stringr::str_replace(reference_clusters, ".*snn_res", ref_name)

    predictions <- map(refdata, ~ TransferData(
        anchorset = object.anchors, refdata = .x,
        dims = 1:30
    ))
    predictions0 <-
        predictions %>%
        map(as_tibble, rownames = "sample_id") %>%
        map(~ dplyr::select(.x, sample_id, predicted.id)) %>%
        dplyr::bind_rows(.id = "resolution") %>%
        tidyr::spread(resolution, predicted.id) %>%
        purrr::set_names(c("sample_id", ref_names)) %>%
        tibble::column_to_rownames("sample_id") %>%
        identity()

    query_object <- AddMetaData(query_object, predictions0)
}


## ------------------------------------------------------------------------
## find markers for every cluster compared to all remaining cells, report only the positive ones
#' Find Cell Type Markers in a Seurat Object
#'
#' @param object A object
#'
#' @return
#' @export
#'
#' @examples
object_find_markers <- function(object, num_features) {
    markers <- Seurat::FindAllMarkers(object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
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
