#' Create a Table of single Cell Projects
#'
#' Uses a list of projects to create a matrix of single cell projects
#'
#' @param proj_list List of projects
#'
#' @return a tibble of single cell projects
#'
create_proj_matrix <- function(proj_list) {
    proj_list <- unlist(proj_list)

    proj_tbl <- tibble(project_path = proj_list, project_name = path_file(proj_list))

    patterns <- c("{date}-{user}-{note}-{species}_proj")

    proj_matrix <- unglue::unglue_data(proj_list, patterns) %>%
        mutate(date = path_file(date)) %>%
        bind_cols(proj_tbl) %>%
        identity()

    primary_projects <-
        proj_matrix %>%
        filter(!grepl("integrated_projects", project_path)) %>%
        filter(str_count(project_name, "_") == 1) %>%
        identity()

    integrated_projects <-
        proj_matrix %>%
        anti_join(primary_projects) %>%
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
#' @return a single cell object
#'
subset_by_meta <- function(meta_path, object) {
    upload_meta <- read_csv(meta_path, col_names = "sample_id") %>%
        filter(!is.na(sample_id) & !sample_id == "sample_id") %>%
        mutate(name = sample_id) %>%
        column_to_rownames("sample_id") %>%
        identity()

    upload_cells <- rownames(upload_meta)

    object <- object[, colnames(object) %in% upload_cells]

    object <- SingleCellExperiment::AddMetaData(object, upload_meta)

    return(object)
}

#' Combine Loom Files
#'
#' @param projectPaths project paths
#' @param newProjectPath new project path
#'
#' @return a combined loom file
#' @export
combine_looms <- function(projectPaths, newProjectPath) {
    # loom combine
    loompy <- reticulate::import("loompy")

    loom_filenames <- str_replace(path_file(projectPaths), "_proj", ".loom")

    selected_looms <- path(projectPaths, "output", "velocyto", loom_filenames)

    file_not_dir <- function(f){file.exists(f) && !dir.exists(f)}

    if (all(file_not_dir(selected_looms))) loompy$combine(selected_looms, newProjectPath)
}

#' Read in Gene and Transcript SingleCellExperiment Objects
#'
#' @param proj_dir path to project directory
#' @param prefix default "unfiltered"
#'
#' @return a single cell object
#' @export
#' @examples
load_object_path <- function(proj_dir = getwd(), prefix = "unfiltered") {
  object_regex <- paste0(paste0(".*/", prefix, "_object.rds"))

  object_path <- path(proj_dir, "output", "seurat") %>%
    dir_ls(regexp = object_regex)

  if (!rlang::is_empty(object_path)) {
    return(object_path)
  }

  stop("'", object_path, "' does not exist",
       paste0(" in current working directory ('", getwd(), "')"),
       ".",
       call. = FALSE
  )
}


#' Load SingleCellExperiment Files from a single project path
#'
#' @param proj_dir project directory
#' @param ... extra args passed to load_object_path
#'
#' @return a single cell object
#' @examples
load_object_from_proj <- function(proj_dir, ...) {
  object_file <- load_object_path(proj_dir, ...)
  object_file <- readRDS(object_file)
}
