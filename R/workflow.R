
# integration workflow ------------------------------

setGeneric("integration_workflow", function(batches, excluded_cells = NULL, resolution = seq(0.2, 2.0, by = 0.2), experiment_name = "default_experiment", organism = "human", ...)
  standardGeneric("integration_workflow"))

#' Integration Workflow
#'
#' Integrate multiple objects and save to file
#'
#' @param batches objects for all batches provided as a list. If named, the resulting integrated object will be identified with corresponding values in 'batch' metadata
#' @param excluded_cells named list of cells to exclude
#' @param resolution value(s) to control the clustering resolution via `Seurat::FindMarkers`
#' @param experiment_name arbitrary name to identify experiment
#' @param organism either "human" or "mouse"
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
setMethod("integration_workflow", "Seurat",
          function(batches, excluded_cells = NULL, resolution = seq(0.2, 2.0, by = 0.2), experiment_name = "default_experiment", organism = "human", ...) {
            checkmate::check_list(batches)

            checkmate::check_character(excluded_cells)

            # organisms <- purrr::map(batches, Misc, c("experiment", "organism"))

            organisms <- purrr::map(batches, list("meta.data", "organism", 1))

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
              merged_batches@misc$batches <- names(batches)

            } else if (all(batch_organisms == "human")) {
              merged_batches <- object_integration_pipeline(batches, resolution = resolution, organism = "human", ...)
              merged_batches@misc$batches <- names(batches)

            } else {
              mouse_object_list <- batches[names(organisms[organisms == "mouse"])]
              human_object_list <- batches[names(organisms[organisms == "human"])]
              merged_batches <- cross_species_integrate(mouse_object_list = mouse_object_list, human_object_list = human_object_list)
              merged_batches@misc$batches <- names(batches)
            }

            merged_batches <- record_experiment_data(merged_batches, experiment_name, organism)

            return(merged_batches)
          }
          )

#' Integration Workflow
#'
#' Integrate multiple objects and save to file
#'
#' @param batches objects for all batches provided as a list. If named, the resulting integrated object will be identified with corresponding values in 'batch' metadata
#' @param excluded_cells named list of cells to exclude
#' @param resolution value(s) to control the clustering resolution via `Seurat::FindMarkers`
#' @param experiment_name arbitrary name to identify experiment
#' @param organism either "human" or "mouse"
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
setMethod("integration_workflow", "SingleCellExperiment",
          function(batches, excluded_cells = NULL, resolution = seq(0.2, 2.0, by = 0.2), experiment_name = "default_experiment", organism = "human", ...) {
            checkmate::check_list(batches)

            checkmate::check_character(excluded_cells)

            # organisms <- purrr::map(batches, Misc, c("experiment", "organism"))

            organisms <- purrr::map(batches, list("meta.data", "organism", 1))

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
              merged_batches@misc$batches <- names(batches)

            } else if (all(batch_organisms == "human")) {
              merged_batches <- object_integration_pipeline(batches, resolution = resolution, organism = "human", ...)
              merged_batches@misc$batches <- names(batches)

            } else {
              mouse_object_list <- batches[names(organisms[organisms == "mouse"])]
              human_object_list <- batches[names(organisms[organisms == "human"])]
              merged_batches <- cross_species_integrate(mouse_object_list = mouse_object_list, human_object_list = human_object_list)
              merged_batches@misc$batches <- names(batches)
            }

            merged_batches <- record_experiment_data(merged_batches, experiment_name, organism)

            return(merged_batches)
          }
)


# clustering workflow ------------------------------


setGeneric("clustering_workflow", function(object, excluded_cells, resolution = seq(0.2, 2.0, by = 0.2), organism = "human", experiment_name = "default_experiment", ...)
  standardGeneric("clustering_workflow"))

#' Clustering Workflow
#'
#' Cluster and Reduce Dimensions of a object
#'
#' @param feature_objects list of objects named according to feature of interest ("gene" or "transcript")
#' @param excluded_cells named list of cells to exclude
#' @param resolution resolution(s) to use for clustering cells
#' @param organism Organism
#' @param experiment_name
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
setMethod("clustering_workflow", "Seurat",
          function(object, excluded_cells, resolution = seq(0.2, 2.0, by = 0.2), organism = "human", experiment_name = "default_experiment", ...) {
            object <- object_pipeline(object, resolution = resolution, organism = organism, ...)

            object <- record_experiment_data(object, experiment_name, organism)

            # save_object(feature_objects, proj_dir = proj_dir, ...)
          }
          )

#' Clustering Workflow
#'
#' Cluster and Reduce Dimensions of a object
#'
#' @param feature_objects list of objects named according to feature of interest ("gene" or "transcript")
#' @param excluded_cells named list of cells to exclude
#' @param resolution resolution(s) to use for clustering cells
#' @param organism Organism
#' @param experiment_name
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
setMethod("clustering_workflow", "SingleCellExperiment",
          function(object, excluded_cells, resolution = seq(0.2, 2.0, by = 0.2), organism = "human", experiment_name = "default_experiment", ...) {
            object <- object_pipeline(object, resolution = resolution, organism = organism, ...)

            object <- record_experiment_data(object, experiment_name, organism)

            # save_object(feature_objects, proj_dir = proj_dir, ...)
          }
          )
