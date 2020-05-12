#' Create a Table of single Cell Projects
#'
#' @param projects_dir
#' @param sub_dirs
#'
#' @return
#' @export
#'
#' @examples
create_proj_matrix <- function(proj_list){

  proj_list <- unlist(proj_list)

  proj_tbl <- tibble::tibble(project_path = proj_list, project_name = fs::path_file(proj_list))

  patterns <- c("{date}-{user}-{note}-{species}_proj")

  proj_matrix <- unglue::unglue_data(proj_list, patterns) %>%
    dplyr::mutate(date = fs::path_file(date)) %>%
    dplyr::bind_cols(proj_tbl) %>%
    identity()

  primary_projects <-
    proj_matrix %>%
    dplyr::filter(!grepl("integrated_projects", project_path)) %>%
    dplyr::filter(stringr::str_count(project_name, "_") == 1) %>%
  identity()

  integrated_projects <-
    proj_matrix %>%
    dplyr::anti_join(primary_projects) %>%
    identity()

  proj_matrices <- list(primary_projects = primary_projects, integrated_projects = integrated_projects)

  return(proj_matrices)
}


#' Create Project database
#'
#' @return
#' @export
#'
#' @examples
create_proj_db <- function(projects_dir = "/dataVolume/storage/single_cell_projects/"){

  system_command <- "updatedb -l 0 -U /dataVolume/storage/single_cell_projects/ -o /dataVolume/storage/single_cell_projects/single_cell_projects.db"
  system(system_command, wait = TRUE)
  print(system_command)

}


#' subset by new metadata
#'
#' @param meta_path
#' @param seu
#'
#' @return
#' @export
#'
#' @examples
subset_by_meta <- function(meta_path, seu){
  upload_meta <- read.csv(meta_path, row.names = 1, header = TRUE)

  upload_cells <- rownames(upload_meta)

  seu <- seu[, colnames(seu) %in% upload_cells]

  seu <- AddMetaData(seu, upload_meta)

  # seu@meta.data <- upload_meta
  return(seu)

}

#' Title
#'
#' @param projectPaths
#' @param newProjectPath
#'
#' @return
#' @export
#'
#' @examples
combine_looms <- function(projectPaths, newProjectPath){
  #loom combine
  loompy <- reticulate::import("loompy")

  loom_filenames <- stringr::str_replace(fs::path_file(projectPaths), "_proj", ".loom")

  selected_looms <- fs::path(projectPaths, "output", "velocyto", loom_filenames)
  loompy$combine(selected_looms, newProjectPath)
}

