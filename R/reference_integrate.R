#' Transfer Labels Between SingleCellExperiment Objects
#'
#' @param ref_object reference object
#' @param query_object query object
#' @param ref_name reference name
#' @param query_name query name
#'
#' @return a single cell object
#' @export
label_transfer <- function(ref_object, query_object, ref_name = NULL, query_name = NULL) {
    label_transerred <- purrr::map2(ref_object, query_object, reference_integrate, query_name = query_name, ref_name = ref_name)
}

#' Transfer labels between gene or transcript objects
#'
#' @param ref_object reference object
#' @param query_object query object
#' @param query_name query name
#' @param ref_name reference name
#'
#' @return a single cell object
#' @export
reference_integrate <- function(ref_object, query_object, query_name = "fetal", ref_name = "organoid") {
    object.anchors <- FindTransferAnchors(reference = ref_object, query = query_object, dims = 1:30, project.query = TRUE)

    reference_clusters <- colnames(ref_object[[]])[grepl(paste0("gene", "_snn_res."), colnames(ref_object[[]]))]

    refdata <- t(ref_object[[reference_clusters]])

    cellids <- colnames(refdata)
    refdata <- setNames(split(refdata, seq(nrow(refdata))), rownames(refdata)) %>%
        map(purrr::set_names, cellids)

    ref_names <- str_replace(reference_clusters, ".*snn_res", ref_name)

    predictions <- map(refdata, ~ TransferData(
        anchorset = object.anchors, refdata = .x,
        dims = 1:30
    ))
    predictions0 <-
        predictions %>%
        map(as_tibble, rownames = "sample_id") %>%
        map(~ select(.x, sample_id, predicted.id)) %>%
        bind_rows(.id = "resolution") %>%
        tidyr::spread(resolution, predicted.id) %>%
        purrr::set_names(c("sample_id", ref_names)) %>%
        column_to_rownames("sample_id") %>%
        identity()

    query_object <- AddMetaData(query_object, predictions0)
}



## find markers for every cluster compared to all remaining cells, report only the positive ones
#' Find Cell Type Markers in a SingleCellExperiment Object
#'
#' @param object A object
#' @param num_features number of features to retrieve for marker genes
#'
#' @return a single cell object
#' @export
find_markers <- function(object, num_features) {
    markers <- SingleCellExperiment::FindAllMarkers(object, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
    markers %>%
        group_by(cluster) %>%
        top_n(n = 2, wt = avg_logFC)

    cluster_markers <- reference.markers$gene %>%
        group_by(cluster) %>%
        top_n(n = num_features, wt = avg_logFC) %>%
        pull(gene)

    DoHeatmap(reference_integrated[[1]], features = cluster_markers)
    return(cluster_markers)
}
