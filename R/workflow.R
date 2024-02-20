#' Integration Workflow
#'
#' Integrate multiple objects and save to file
#'
#' @param batches objects for all batches provided as a list. If named, the resulting integrated object will be identified with corresponding values in 'batch' metadata
#' @param excluded_cells named list of cells to exclude
#' @param resolution value(s) to control the clustering resolution via `Seurat::FindMarkers`
#' @param experiment_name arbitrary name to identify experiment
#' @param organism either "human" or "mouse"
#' @param ... extra args passed to object_integration_pipeline
#'
#' @return an integrated single cell object
#' @export
#'
setGeneric("integration_workflow", function(batches, excluded_cells = NULL, resolution = seq(0.2, 2.0, by = 0.2), experiment_name = "default_experiment", organism = "human", ...) {
    standardGeneric("integration_workflow")
})

setMethod(
    "integration_workflow", "Seurat",
    function(batches, excluded_cells = NULL, resolution = seq(0.2, 2.0, by = 0.2), experiment_name = "default_experiment", organism = "human", ...) {

        # organisms <- map(batches, Misc, c("experiment", "organism"))

        organisms <- map(batches, list("meta.data", "organism", 1))

        if (any(purrr::map_lgl(organisms, is.null))) {
            organisms <- case_when(
                grepl("Hs", names(batches)) ~ "human",
                grepl("Mm", names(batches)) ~ "mouse",
                TRUE ~ "human"
            )
            names(organisms) <- names(batches)
        }

        experiment_names <- names(batches)

        batches <- purrr::pmap(list(batches, experiment_names, organisms), record_experiment_data)

        batch_organisms <- map_chr(batches, list("misc", "experiment", "organism"))

        if (all(batch_organisms == "mouse")) {
            merged_batches <- purrr::imap(batches, object_integration_pipeline, resolution = resolution, organism = "mouse", ...)
            Misc(merged_batches)$batches <- names(batches)
        } else if (all(batch_organisms == "human")) {
            merged_batches <- object_integration_pipeline(batches, resolution = resolution, organism = "human", ...)
            Misc(merged_batches)$batches <- names(batches)
        } else {
            mouse_object_list <- batches[names(organisms[organisms == "mouse"])]
            human_object_list <- batches[names(organisms[organisms == "human"])]
            merged_batches <- cross_species_integrate(mouse_object_list = mouse_object_list, human_object_list = human_object_list)
            Misc(merged_batches)$batches <- names(batches)
        }

        merged_batches <- record_experiment_data(merged_batches, experiment_name, organism)

        return(merged_batches)
    }
)

setMethod(
    "integration_workflow", "SingleCellExperiment",
    function(batches, excluded_cells = NULL, resolution = seq(0.2, 2.0, by = 0.2), experiment_name = "default_experiment", organism = "human", ...) {

        # organisms <- map(batches, Misc, c("experiment", "organism"))

        organisms <- map(batches, list("meta.data", "organism", 1))

        if (any(purrr::map_lgl(organisms, is.null))) {
            organisms <- case_when(
                grepl("Hs", names(batches)) ~ "human",
                grepl("Mm", names(batches)) ~ "mouse",
                TRUE ~ "human"
            )
            names(organisms) <- names(batches)
        }

        experiment_names <- names(batches)

        batches <- purrr::pmap(list(batches, experiment_names, organisms), record_experiment_data)

        batch_organisms <- map_chr(batches, list("misc", "experiment", "organism"))

        if (all(batch_organisms == "mouse")) {
            merged_batches <- purrr::imap(batches, object_integration_pipeline, resolution = resolution, organism = "mouse", ...)
            Misc(merged_batches)$batches <- names(batches)
        } else if (all(batch_organisms == "human")) {
            merged_batches <- object_integration_pipeline(batches, resolution = resolution, organism = "human", ...)
            Misc(merged_batches)$batches <- names(batches)
        } else {
            mouse_object_list <- batches[names(organisms[organisms == "mouse"])]
            human_object_list <- batches[names(organisms[organisms == "human"])]
            merged_batches <- cross_species_integrate(mouse_object_list = mouse_object_list, human_object_list = human_object_list)
            Misc(merged_batches)$batches <- names(batches)
        }

        merged_batches <- record_experiment_data(merged_batches, experiment_name, organism)

        return(merged_batches)
    }
)

#' Clustering Workflow
#'
#' Cluster and Reduce Dimensions of a object
#'
#' @param object a single cell object
#' @param excluded_cells named list of cells to exclude
#' @param resolution resolution(s) to use for clustering cells
#' @param organism Organism
#' @param experiment_name name of the experiment
#' @param ... extra args passed to object_pipeline
#'
#' @return a clustered single cell object
#' @export
setGeneric("clustering_workflow", function(object, excluded_cells, resolution = seq(0.2, 2.0, by = 0.2), organism = "human", experiment_name = "default_experiment", ...) {
    standardGeneric("clustering_workflow")
})

setMethod(
    "clustering_workflow", "Seurat",
    function(object, excluded_cells, resolution = seq(0.2, 2.0, by = 0.2), organism = "human", experiment_name = "default_experiment", ...) {
        object <- object_pipeline(object, resolution = resolution, organism = organism, ...)

        object <- record_experiment_data(object, experiment_name, organism)

        # save_object(feature_objects, proj_dir = proj_dir, ...)
    }
)

setMethod(
    "clustering_workflow", "SingleCellExperiment",
    function(object, excluded_cells, resolution = seq(0.2, 2.0, by = 0.2), organism = "human", experiment_name = "default_experiment", ...) {
        object <- object_pipeline(object, resolution = resolution, organism = organism, ...)

        object <- record_experiment_data(object, experiment_name, organism)

        # save_object(feature_objects, proj_dir = proj_dir, ...)
    }
)
