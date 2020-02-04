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
    dplyr::filter(str_count(project_name, "_") == 1) %>%
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

  projects_db <- paste0(projects_dir, "single_cell_projects.db")

  system(paste0("updatedb -l 0 -U ", projects_dir, " -o ", projects_db))

}




