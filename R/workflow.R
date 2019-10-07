<<<<<<< HEAD

#' Integration Workflow
#'
#' Integrate multiple seurat objects and save to file
#'
#' @param proj_dir home directory of current project
#' @param child_proj_dirs child projects to be integrated
#' @param excluded_cells named list of cells to exclude
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
integration_workflow <- function(proj_dir, child_proj_dirs, excluded_cells, ...) {

  names(child_proj_dirs) <- gsub("_proj", "", fs::path_file(child_proj_dirs))

  # load seurat objects from 'child' projects
  seus <- purrr::map(child_proj_dirs, load_seurat_from_proj)

  seus <- purrr::transpose(seus)

  # run seurat batch correction on 'child' projects
  corrected_seus <- purrr::map(seus, seuratTools::seurat_batch_correct)

  # cluster merged seurat objects
  corrected_seus <- map(corrected_seus, seuratTools::seurat_cluster, resolution = seq(0.2, 2.0, by = 0.2))

  # add read count column
  corrected_seus <- map(corrected_seus, add_read_count_col)

  purrr::imap(corrected_seus, ~save_seurat(proj_dir, .x, .y))

  # annotate cell cycle scoring to seurat objects

  corrected_seus <- annotate_cell_cycle(corrected_seus)

  #annotate excluded cells

  corrected_seus <- map(corrected_seus, annotate_excluded, excluded_cells)

  # add marker genes to seurat objects

  corrected_seus <- map(corrected_seus, find_all_markers)

  # load cell cycle annotated, batch corrected seurat objects
  purrr::imap(corrected_seus, save_seurat)

  # annotated photoreceptors: filter seurat objects by metadata criteria

  excluded_PRs_seus <- filter_merged_seus(corrected_seus, filter_var = "excluded_because", filter_val = NA)

  excluded_PRs_seus <- reintegrate_seus(excluded_PRs_seus, prefix = "remove_annotated")

  #  low read count: filter seurat objects by metadata criteria

  high_rc_seus <- filter_merged_seus(corrected_seus, filter_var = "read_count", filter_val = NA)

  high_rc_seus <- reintegrate_seus(high_rc_seus, prefix = "remove_lowrc")

  ## filter seurat objects by both read count and exclusion annotation

  excluded_PRs_seus <- filter_merged_seus(corrected_seus, filter_var = "excluded_because", filter_val = NA)

  high_rc_seus <- filter_merged_seus(excluded_PRs_seus, filter_var = "read_count", filter_val = NA)

  rc_and_excl_seus <- reintegrate_seus(high_rc_seus, prefix = "remove_annotated_and_lowrc")

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


||||||| merged common ancestors
=======

#' Integration Workflow
#'
#' Integrate multiple seurat objects and save to file
#'
#' @param proj_dir home directory of current project
#' @param child_proj_dirs child projects to be integrated
#' @param excluded_cells named list of cells to exclude
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
integration_workflow <- function(proj_dir, child_proj_dirs, excluded_cells, ...) {

  # # load packages
  #
  # library(tidyverse)
  # library(Seurat)
  # library(fs)
  # library(annotables)
  # library(rprojroot)
  # library(seuratTools)
  #
  # # future_max_size = 850*1024^2
  # # options(future.globals.maxSize= future_max_size)
  # # library(future)
  # # plan(multiprocess)
  #
  # features <- c("gene", "transcript")
  # proj_dir = rprojroot::find_root(criterion = has_file_pattern("*.Rproj"))

  names(child_proj_dirs) <- gsub("_proj", "", fs::path_file(child_proj_dirs))

  # load seurat objects from 'child' projects
  seus <- purrr::map(child_proj_dirs, load_seurat_from_proj)

  seus <- purrr::transpose(seus)

  # run seurat batch correction on 'child' projects
  corrected_seus <- purrr::map(seus, seuratTools::seurat_batch_correct)

  # cluster merged seurat objects
  corrected_seus <- map(corrected_seus, seuratTools::seurat_cluster, resolution = seq(0.2, 2.0, by = 0.2))

  # add read count column
  corrected_seus <- map(corrected_seus, add_read_count_col)

  purrr::imap(corrected_seus, ~save_seurat(proj_dir, .x, .y))

  # annotate cell cycle scoring to seurat objects

  corrected_seus <- annotate_cell_cycle(corrected_seus)

  #annotate excluded cells

  corrected_seus <- map(corrected_seus, annotate_excluded, excluded_cells)

  # add marker genes to seurat objects

  corrected_seus <- map(corrected_seus, find_all_markers)

  # load cell cycle annotated, batch corrected seurat objects
  purrr::imap(corrected_seus, save_seurat)

  # annotated photoreceptors: filter seurat objects by metadata criteria

  excluded_PRs_seus <- filter_merged_seus(corrected_seus, filter_var = "excluded_because", filter_val = NA)

  excluded_PRs_seus <- reintegrate_seus(excluded_PRs_seus, prefix = "remove_annotated")

  #  low read count: filter seurat objects by metadata criteria

  high_rc_seus <- filter_merged_seus(corrected_seus, filter_var = "read_count", filter_val = NA)

  high_rc_seus <- reintegrate_seus(high_rc_seus, prefix = "remove_lowrc")

  ## filter seurat objects by both read count and exclusion annotation

  excluded_PRs_seus <- filter_merged_seus(corrected_seus, filter_var = "excluded_because", filter_val = NA)

  high_rc_seus <- filter_merged_seus(excluded_PRs_seus, filter_var = "read_count", filter_val = NA)

  rc_and_excl_seus <- reintegrate_seus(high_rc_seus, prefix = "remove_annotated_and_lowrc")

}



>>>>>>> 2d43e2725a02764750212f5c9acf866ac6e6a483
