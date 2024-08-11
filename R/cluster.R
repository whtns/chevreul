#' Preprocess Single Cell Object
#'
#' Performs standard pre-processing workflow for scRNA-seq data
#'
#' @param object Assay to use
#' @param scale Perform linear transformation 'Scaling'
#' @param normalize Perform normalization
#' @param features Identify highly variable features
#' @param legacy_settings Use legacy settings
#' @param ... extra args passed to scaling functions
#'
#' @return a preprocessed SingleCellExperiment object
object_preprocess <- function(object, scale = TRUE, normalize = TRUE, features = NULL, legacy_settings = FALSE, ...) {
    clusters <- quickCluster(object)
    object <- computeSumFactors(object, clusters = clusters)
    # summary(sizeFactors(object))

    object <- logNormCounts(object)

    dec <- modelGeneVar(object)

    # Get the top 10% of genes.
    top.hvgs <- getTopHVGs(dec, prop = 0.1)
    object <- runPCA(object, subset_row = top.hvgs)

    output <- getClusteredPCs(reducedDim(object))

    g1 <- buildSNNGraph(object, use.dimred = "PCA")


    return(object)
}


#' Find All Markers
#'
#' Find all markers at a range of resolutions
#'
#' @param object An object.
#' @param group_by A metadata variable to group by.
#' @param experiment Assay to use, Default "gene".
#' @param ... extra args passed to stash_marker_features
#'
#' @return a SingleCellExperiment object containing marker genes
#' @export
#'
#' @examples
#' chevreul_sce <- chevreuldata::human_gene_transcript_sce()
#' markers_stashed_object <- find_all_markers(chevreul_sce, group_by = "Age")
find_all_markers <- function(object, group_by = NULL, experiment = "gene", ...) {
    if (is.null(group_by)) {
        resolutions <- colnames(get_cell_metadata(object))[grepl(paste0(experiment, "_snn_res."), colnames(get_cell_metadata(object)))]
        cluster_index <- grepl(paste0(experiment, "_snn_res."), colnames(get_cell_metadata(object)))
        if (!any(cluster_index)) {
            warning("no clusters found in metadata. runnings object_cluster")
            object <- object_cluster(object, resolution = seq(0.2, 1, by = 0.2))
        }
        clusters <- get_cell_metadata(object)[, cluster_index]
        cluster_levels <- map_int(clusters, ~ length(unique(.x)))
        cluster_levels <- cluster_levels[cluster_levels > 1]
        clusters <- select(clusters, one_of(names(cluster_levels)))
        group_by <- names(clusters)
    }
    new_markers <- map(group_by, ~ stash_marker_features(object, .x, experiment = experiment, ...))
    names(new_markers) <- group_by
    old_markers <- metadata(object)$markers[!names(metadata(object)[["markers"]]) %in% names(new_markers)]
    metadata(object)[["markers"]] <- c(old_markers, new_markers)
    return(object)
}


#' Enframe Markers
#'
#' @param marker_table a table of marker genes
#'
#' @return a table of marker genes
#' @export
#' @examples
#' chevreul_sce <- chevreuldata::human_gene_transcript_sce()
#' marker_table <- metadata(chevreul_sce)$markers[["batch"]]
#' enframe_markers(marker_table)
#'
enframe_markers <- function(marker_table) {
    marker_table %>%
        select(Gene.Name, Cluster) %>%
        mutate(rn = row_number()) %>%
        pivot_wider(names_from = Cluster, values_from = Gene.Name) %>%
        select(-rn)
}

#' Stash Marker Genes in a SingleCellExperiment Object
#'
#' Marker Genes will be stored in object metadata as `markers`
#'
#' @param group_by A metadata variable to group by
#' @param object A object
#' @param experiment An experiment to use
#' @param top_n Use top n genes, Default 200
#' @param p_val_cutoff p value cut-off, Default value is "0.5"
#'
#' @return a SingleCellExperiment object with marker genes
#'
stash_marker_features <- function(object, group_by, experiment = "gene", top_n = 200, p_val_cutoff = 0.5) {
    message("stashing markers for ", group_by)
    markers <- list()
    markers <-
        findMarkers(object, test.type = "t", groups = colData(object)[[group_by]]) %>%
        map(as.data.frame) %>%
        map(rownames_to_column, "feature") %>%
        bind_rows(.id = "group") %>%
        group_by(group) %>%
        filter(FDR < p_val_cutoff) %>%
        top_n(n = top_n, wt = summary.logFC) %>%
        arrange(group, desc(summary.logFC)) %>%
        select(Gene.Name = feature, Average.Log.Fold.Change = summary.logFC, Adjusted.pvalue = FDR, Cluster = group) %>%
        identity()
    return(markers)
}
