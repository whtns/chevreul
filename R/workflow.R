#' Integration Workflow
#'
#' Integrate multiple objects and save to file
#'
#' @param batches objects for all batches provided as a list. If named, the
#' resulting integrated object will be identified with corresponding values in
#' 'batch' metadata
#' @param excluded_cells named list of cells to exclude
#' @param resolution value(s) to control the clustering resolution
#' via `scran::findMarkers`
#' @param experiment_name arbitrary name to identify experiment
#' @param organism either "human" or "mouse"
#' @param ... extra args passed to object_integration_pipeline
#'
#' @return an integrated SingleCellExperiment object
#' @export
#' @examples
#' chevreul_sce <- chevreuldata::human_gene_transcript_sce()
#' batches <- splitByCol(chevreul_sce, "batch")
#' integration_workflow(batches)
integration_workflow <- function(batches, excluded_cells = NULL,
                                 resolution = seq(0.2, 2.0, by = 0.2),
                                 experiment_name = "default_experiment",
                                 organism = "human", ...) {
    # organisms <- map(batches, Misc, c("experiment", "organism"))

    organisms <- map(batches, list("meta.data", "organism", 1))

    if (any(map_lgl(organisms, is.null))) {
        organisms <- case_when(
            grepl("Hs", names(batches)) ~ "human",
            grepl("Mm", names(batches)) ~ "mouse",
            TRUE ~ "human"
        )
        names(organisms) <- names(batches)
    }

    experiment_names <- names(batches)

    batches <- pmap(list(batches, experiment_names, organisms),
                    record_experiment_data)

    merged_batches <- object_integration_pipeline(batches,
                                                  resolution = resolution,
                                                  organism = "human", ...)
    metadata(merged_batches)$batches <- names(batches)

    merged_batches <- record_experiment_data(merged_batches, experiment_name,
                                             organism)

    return(merged_batches)
}


#' Clustering Workflow
#'
#' Cluster and Reduce Dimensions of a object
#'
#' @param object a SingleCellExperiment object
#' @param excluded_cells named list of cells to exclude
#' @param resolution resolution(s) to use for clustering cells
#' @param organism Organism
#' @param experiment_name name of the experiment
#' @param ... extra args passed to object_pipeline
#'
#' @return a clustered SingleCellExperiment object
#' @export
#' @examples
#' chevreul_sce <- chevreuldata::human_gene_transcript_sce()
#' clustering_workflow(chevreul_sce)
clustering_workflow <- function(object, excluded_cells,
                                resolution = seq(0.2, 2.0, by = 0.2),
                                organism = "human",
                                experiment_name = "default_experiment", ...) {
    object <- object_pipeline(object, resolution = resolution, ...)

    object <- record_experiment_data(object, experiment_name, organism)
}
