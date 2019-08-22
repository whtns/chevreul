#' Collate list of variables to be plotted
#'
#' @param seu
#'
#' @return plot_types
#'
#' @examples
list_plot_types <- function(seu){
  meta_types <- unlist(lapply(seu[[]], class))

  meta_types <- meta_types[!grepl(paste0("^", DefaultAssay(seu), "_snn_res"), names(meta_types))]

  meta_types <- split(meta_types, factor(meta_types))
  meta_types <- lapply(meta_types, names)


  continuous_vars <- c("custom", unlist(meta_types[c("integer", "numeric")]))

  continuous_names <- stringr::str_to_title(gsub("[[:punct:]]", " ", continuous_vars))

  names(continuous_vars) <- continuous_names


  category_vars <- c("seurat", unlist(meta_types[c("character", "factor")]))
  category_vars <- category_vars[!grepl("snn_res", category_vars)]

  category_names <- stringr::str_to_title(gsub("[[:punct:]]", " ", category_vars))

  names(category_vars) <- category_names

  plot_types <- list(category_vars = category_vars, continuous_vars = continuous_vars)

  return(plot_types)


}

#' Reorganize seurat objects
#'
#' @param proj_dir
#'
#' @return
#' @export
#'
#' @examples
reorg_seurat_files <- function(proj_dir = NULL){
  seurat_dir <- fs::path(proj_dir, "output", "seurat")

  features <- c("gene", "transcript")

  rds_files <- fs::dir_ls(seurat_dir) %>%
    fs::path_filter(paste0("*", features, "_seu", "*", ".rds", collapse = "|")) %>%
    purrr::set_names(gsub(".*seu|.rds", "", fs::path_file(.)))

  names(rds_files) <- gsub("^_", "", names(rds_files))
  names(rds_files)[names(rds_files) == ""] <- "unfiltered"

  names(rds_files) <- paste0(names(rds_files), "_seu.rds")

  rds_files <- rds_files[order(names(rds_files))]

  rds_files <- split(rds_files, names(rds_files))

  rds_files <- purrr::map(rds_files, ~purrr::set_names(.x, c("gene", "transcript")))

  message(paste0("reading in files: ", unlist(rds_files)))
  rds_files <- purrr::map(rds_files, ~purrr::map(.x, readRDS))

  new_rds_paths <- fs::path(seurat_dir, names(rds_files))

  purrr::map2(rds_files, new_rds_paths, saveRDS)

  old_files <- fs::dir_ls(seurat_dir) %>%
    fs::file_info() %>%
    dplyr::filter(modification_time < lubridate::today()) %>%
    dplyr::pull(path)

  message(paste0("deleting old files: ", old_files))
  fs::file_delete(old_files)

  return(rds_files)
}
