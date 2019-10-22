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

  patterns <- c("{date}-{user}-{note}-{species}_proj")
  proj_matrix <- unglue::unglue_data(proj_list, patterns) %>%
    tidyr::unite(project_path, date, user, note, species, sep = "-", remove = F) %>%
    dplyr::mutate(project_path = paste0(project_path, "_proj")) %>%
    dplyr::mutate(project_name = fs::path_file(project_path)) %>%
    dplyr::mutate(date = fs::path_file(date)) %>%
    identity()

  primary_projects <-
    proj_matrix %>%
    dplyr::filter(species %in% c("Hs", "Mm"))

  integrated_projects <-
    proj_matrix %>%
    dplyr::filter(!species %in% c("Hs", "Mm"))

  proj_matrices <- list(primary_projects = primary_projects, integrated_projects = integrated_projects)

  return(proj_matrices)
}

#' Create List of Projects
#'
#' @param proj_matrices
#'
#' @return
#' @export
#'
#' @examples
create_proj_list <- function(projects_dir = "/dataVolume/storage/single_cell_projects", sub_dirs = c("sc_cone_devel", "sc_RB_devel", "integrated_projects", "resources")){

  project_list <- fs::dir_ls(fs::path(projects_dir, sub_dirs), recurse = T, glob = "*_proj") %>%
    # fs::path_filter("*_proj") %>%
    tibble::enframe("name", "path") %>%
    dplyr::mutate(sub_dir = dplyr::case_when(grepl("sc_cone_devel", path) ~ "sc_cone_devel",
                                      grepl("sc_RB_devel", path) ~ "sc_RB_devel",
                                      grepl("integrated_projects", path) ~ "integrated_projects",
                                      grepl("resources", path) ~ "resources")) %>%
    dplyr::select(-name) %>%
    split(.$sub_dir) %>%
    purrr::map(~dplyr::pull(.x, path)) %>%
    identity()

  project_list <- purrr::map(project_list, ~purrr::set_names(.x, fs::path_file(.x)))



  # names(proj_list) <- fs::path_file(fs::path_dir(proj_list))



# proj_list <- fs::dir_ls(fs::path(projects_dir, sub_dirs), recurse = T) %>%
#     fs::path_filter("*_proj") %>%
#     identity()
#
#   primary_project_list <- proj_matrices$primary_projects %>%
#     dplyr::pull(project_path) %>%
#     purrr::set_names(fs::path_file(.))
#
#   integrated_project_list <- proj_matrices$integrated_projects %>%
#     dplyr::pull(project_path) %>%
#     purrr::set_names(fs::path_file(.))
#
#   project_list <- list(primary = primary_project_list,
#                     integrated = integrated_project_list)

  return(project_list)
}

