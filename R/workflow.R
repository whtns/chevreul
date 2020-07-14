
#' Integration Workflow
#'
#' Integrate multiple seurat objects and save to file
#'
#' @param batches seurat objects for each all batches provided as a list. If named, the resulting integrated object will be identified with corresponding values in 'batch' metadata
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
#' batches <- seurat_pancreas_reduced %>%
#'   purrr::map(Seurat::SplitObject, split.by = "dataset") %>%
#'   purrr::transpose()
#'
#' inegrated_seu <- integration_workflow(batches)

integration_workflow <- function(batches, excluded_cells, resolution = seq(0.2, 2.0, by = 0.2), experiment_name = "default_experiment", organism = "human", ...) {

  # names(child_proj_dirs) <- gsub("_proj", "", fs::path_file(child_proj_dirs))

  # load seurat objects from 'child' projects
  # seus <- purrr::map(child_proj_dirs, load_seurat_from_proj)

  organisms <- purrr::map(batches, list(1, "meta.data", "organism", 1))

  if (any(map_lgl(organisms, is.null))){
    organisms <- case_when(
        grepl("Hs", names(batches)) ~ "human",
        grepl("Mm", names(batches)) ~ "mouse"
      )
    names(organisms) <- names(batches)
  }

  experiment_names <- names(batches)

  batches <- purrr::transpose(batches)
  for (i in names(batches)){
    batches[[i]] <- purrr::pmap(list(batches[[i]], experiment_names, organisms), record_experiment_data)
  }
  batches <- purrr::transpose(batches)

  if (all(purrr::map(batches, list(1, "misc", "experiment", "organism")) == "human")){
    batches <- purrr::transpose(batches)
    merged_batches <- purrr::imap(batches, seuratTools::seurat_integration_pipeline, resolution = resolution, organism = "human", ...)
    for (i in names(merged_batches)){
      merged_batches[[i]]@misc$batches <- names(batches)
    }

  } else if (all(purrr::map(batches, list(1, "misc", "experiment", "organism")) == "mouse")){
    batches <- purrr::transpose(batches)
    merged_batches <- purrr::imap(batches, seuratTools::seurat_integration_pipeline, resolution = resolution, organism = "mouse", ...)
    for (i in names(merged_batches)){
      merged_batches[[i]]@misc$batches <- names(batches)
    }

  }  else {

    # mouse_seu_list <- batches[grepl("Mm", names(batches))]
    mouse_seu_list <- batches[names(organisms[organisms == "mouse"])]
    # human_seu_list <- batches[grepl("Hs", names(batches))]
    human_seu_list <- batches[names(organisms[organisms == "human"])]
    merged_batches <- cross_species_integrate(mouse_seu_list = mouse_seu_list, human_seu_list = human_seu_list)
    for (i in names(merged_batches)){
      merged_batches[[i]]@misc$batches <- names(batches)
    }
  }

  merged_batches <- purrr::map(merged_batches, record_experiment_data, experiment_name, organism)

  return(merged_batches)

}


#' Clustering Workflow
#'
#' Cluster and Reduce Dimensions of a seurat object
#'
#' @param feature_seus list of seurat objedevelcts named according to feature of interest ("gene" or "transcript")
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
#' clustered_seu <- clustering_workflow(seurat_pancrease_reduced)
clustering_workflow <- function(feature_seus = NULL, excluded_cells, resolution = seq(0.2, 2.0, by = 0.2), organism = "human", experiment_name = "default_experiment", ...){

  feature_seus <- purrr::imap(feature_seus, seuratTools::seurat_pipeline, resolution = resolution, organism = organism, ...)

  feature_seus <- purrr::map(feature_seus, record_experiment_data, experiment_name, organism)

  # save_seurat(feature_seus, proj_dir = proj_dir, ...)

}


