
#' Integration Workflow
#'
#' Integrate multiple seurat objects and save to file
#'
#' @param child_proj_dirs child projects to be integrated
#' @param excluded_cells named list of cells to exclude
#' @param ...
#' @param resolution
#'
#' @return
#' @export
#'
#' @examples
integration_workflow <- function(child_proj_dirs, excluded_cells, resolution = seq(0.2, 2.0, by = 0.2), ...) {

  # return(list(gene = mtcars))

  names(child_proj_dirs) <- gsub("_proj", "", fs::path_file(child_proj_dirs))

  # load seurat objects from 'child' projects
  seus <- purrr::map(child_proj_dirs, load_seurat_from_proj)

  # check species of child projects
  # project_names <- purrr::map(seus, ~.x[[1]]@project.name)

  if (all(grepl("Hs", names(seus)))){
    seus <- purrr::transpose(seus)
    merged_seus <- purrr::imap(seus, seuratTools::seurat_integration_pipeline, resolution = resolution, organism = "human", ...)
    for (i in names(merged_seus)){
      merged_seus[[i]]@misc$child_projs <- names(child_proj_dirs)
    }

  } else if (all(grepl("Mm", names(seus)))){
    seus <- purrr::transpose(seus)
    merged_seus <- purrr::imap(seus, seuratTools::seurat_integration_pipeline, resolution = resolution, organism = "mouse", ...)
    for (i in names(merged_seus)){
      merged_seus[[i]]@misc$child_projs <- names(child_proj_dirs)
    }

  }  else {

    mouse_seu_list <- seus[grepl("Mm", names(seus))]
    human_seu_list <- seus[grepl("Hs", names(seus))]
    merged_seus <- cross_species_integrate(mouse_seu_list = mouse_seu_list, human_seu_list = human_seu_list)
    for (i in names(merged_seus)){
      merged_seus[[i]]@misc$child_projs <- names(child_proj_dirs)
    }
  }

  return(merged_seus)


}


#' Clustering Workflow
#'
#' Integrate multiple seurat objects and save to file
#'
#' @param feature_seus list of seurat objects named according to feature of interest ("gene" or "transcript")
#' @param excluded_cells named list of cells to exclude
#' @param cell_cycle whether to score and regress cell cycle related features
#' @param resolution resolution(s) to use for clustering cells
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
clustering_workflow <- function(feature_seus = NULL, excluded_cells, cell_cycle = FALSE, resolution = seq(0.2, 2.0, by = 0.2), ...) {

  feature_seus <- purrr::map(feature_seus, seuratTools::seurat_pipeline, resolution = resolution, ...)

  if(cell_cycle){
    # add cell cycle scoring to seurat objects
    feature_seus <- purrr::imap(feature_seus, seuratTools::annotate_cell_cycle, ...)
  }

  # annotate low read count category in seurat metadata
  feature_seus <- purrr::map(feature_seus, seuratTools::add_read_count_col)

  # annotate mitochondrial percentage in seurat metadata
  feature_seus <- purrr::imap(feature_seus, seuratTools::add_percent_mito)

  # save_seurat(feature_seus, proj_dir = proj_dir, ...)

}


