
#' Integration Workflow
#'
#' Integrate multiple seurat objects and save to file
#'
#' @param child_proj_dirs child projects to be integrated
#' @param excluded_cells named list of cells to exclude
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
integration_workflow <- function(child_proj_dirs, excluded_cells, ...) {

  # return(list(gene = mtcars))

  names(child_proj_dirs) <- gsub("_proj", "", fs::path_file(child_proj_dirs))

  # load seurat objects from 'child' projects
  seus <- purrr::map(child_proj_dirs, load_seurat_from_proj)

  # check species of child projects
  # project_names <- purrr::map(seus, ~.x[[1]]@project.name)

  if (all(grepl("Hs", names(seus)))){
    seus <- purrr::transpose(seus)
    merged_seus <- purrr::map(seus, seuratTools::seurat_integration_pipeline, ...)
    for (i in names(merged_seus)){
      merged_seus[[i]]@misc$child_projs <- names(child_proj_dirs)
    }

  } else if (any(grepl("Mm", names(seus)))){

    mouse_seu_list <- seus[grepl("Mm", names(seus))]
    human_seu_list <- seus[grepl("Hs", names(seus))]
    merged_seus <- cross_species_integrate(mouse_seu_list = mouse_seu_list, human_seu_list = human_seu_list)
    for (i in names(merged_seus)){
      merged_seus[[i]]@misc$child_projs <- names(child_proj_dirs)
    }
  }

  return(merged_seus)


}



seurat_integration_pipeline <- function(seus, resolution = seq(0.2, 2.0, by = 0.2), suffix = '', ...) {

  corrected_seu <- seurat_integrate(seus, ...)

  # cluster merged seurat objects
  corrected_seu <- seurat_cluster(corrected_seu, resolution = resolution, ...)

  corrected_seu <- find_all_markers(corrected_seu)

  # corrected_seu <- save_seurat(corrected_seu, feature = feature, suffix = suffix, ...)

}

#' Clustering Workflow
#'
#' Integrate multiple seurat objects and save to file
#'
#' @param proj_dir home directory of current project
#' @param feature_seus
#' @param excluded_cells named list of cells to exclude
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
clustering_workflow <- function(proj_dir, feature_seus = NULL, excluded_cells, cell_cycle = T, ...) {

  feature_seus <- map(feature_seus, seuratTools::seurat_pipeline, resolution = seq(0.2, 2.0, by = 0.2), ...)

  if(cell_cycle){
    # add cell cycle scoring to seurat objects
    feature_seus <- purrr::imap(feature_seus, seuratTools::annotate_cell_cycle, ...)
  }

  # add marker genes to seurat objects
  feature_seus <- map(feature_seus, seuratTools::find_all_markers)

  # annotate low read count category in seurat metadata
  feature_seus <- map(feature_seus, seuratTools::add_read_count_col)

  save_seurat(feature_seus, proj_dir = proj_dir)



}


