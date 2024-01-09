#' Create a Table of single Cell Projects
#'
#' Uses a list of projects to create a matrix of single cell projects
#'
#' @param proj_list List of projects
#'
#' @return
#' @export
#'
#' @examples
create_proj_matrix <- function(proj_list) {
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


#' Subset by new metadata
#'
#' Subset the object using new metadata
#'
#' @param meta_path Path to new metadata
#' @param object A object
#'
#' @return
#' @export
#'
#' @examples
subset_by_meta <- function(meta_path, object) {

  upload_meta <- readr::read_csv(meta_path, col_names = "sample_id") %>%
    dplyr::filter(!is.na(sample_id) & !sample_id == "sample_id") %>%
    dplyr::mutate(name = sample_id) %>%
    tibble::column_to_rownames("sample_id") %>%
    identity()

  upload_cells <- rownames(upload_meta)

  object <- object[, colnames(object) %in% upload_cells]

  object <- Seurat::AddMetaData(object, upload_meta)

  return(object)
}

#' Combine Loom Files
#'
#' @param projectPaths
#' @param newProjectPath
#'
#' @return
#' @export
#'
#' @examples
combine_looms <- function(projectPaths, newProjectPath) {
  # loom combine
  loompy <- reticulate::import("loompy")

  loom_filenames <- stringr::str_replace(fs::path_file(projectPaths), "_proj", ".loom")

  selected_looms <- fs::path(projectPaths, "output", "velocyto", loom_filenames)

  if (all(fs::is_file(selected_looms))) loompy$combine(selected_looms, newProjectPath)
}
