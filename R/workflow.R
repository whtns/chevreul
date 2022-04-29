
#' Integration Workflow
#'
#' Integrate multiple seurat objects and save to file
#'
#' @param batches seurat objects for all batches provided as a list. If named, the resulting integrated object will be identified with corresponding values in 'batch' metadata
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
#' batches <- panc8 %>%
#'   Seurat::SplitObject(split.by = "tech")
#'
#' integrated_seu <- integration_workflow(batches)
integration_workflow <- function(batches, excluded_cells = NULL, resolution = seq(0.2, 2.0, by = 0.2), experiment_name = "default_experiment", organism = "human", ...) {
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
    merged_batches <- purrr::imap(batches, seuratTools::seurat_integration_pipeline, resolution = resolution, organism = "mouse", ...)
    merged_batches@misc$batches <- names(batches)

  } else if (all(batch_organisms == "human")) {
    merged_batches <- seurat_integration_pipeline(batches, resolution = resolution, organism = "human", ...)
    merged_batches@misc$batches <- names(batches)

  } else {
    mouse_seu_list <- batches[names(organisms[organisms == "mouse"])]
    human_seu_list <- batches[names(organisms[organisms == "human"])]
    merged_batches <- cross_species_integrate(mouse_seu_list = mouse_seu_list, human_seu_list = human_seu_list)
    merged_batches@misc$batches <- names(batches)
  }

  merged_batches <- record_experiment_data(merged_batches, experiment_name, organism)

  return(merged_batches)
}


#' Clustering Workflow
#'
#' Cluster and Reduce Dimensions of a seurat object
#'
#' @param feature_seus list of seurat objects named according to feature of interest ("gene" or "transcript")
#' @param excluded_cells named list of cells to exclude
#' @param resolution resolution(s) to use for clustering cells
#' @param organism
#' @param experiment_name
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' clustered_human_seu <- clustering_workflow(panc8)
#' clustered_mouse_seu <- clustering_workflow(baron2016singlecell)
clustering_workflow <- function(seu, excluded_cells, resolution = seq(0.2, 2.0, by = 0.2), organism = "human", experiment_name = "default_experiment", ...) {
  seu <- seurat_pipeline(seu, resolution = resolution, organism = organism, ...)

  seu <- record_experiment_data(seu, experiment_name, organism)

  # save_seurat(feature_seus, proj_dir = proj_dir, ...)
}
