#' Collate list of variables to be plotted
#'
#' @param seu
#'
#' @return plot_types
#'
#' @examples
list_plot_types <- function(seu){

  meta_types <- tibble(
    vars = colnames(seu[[]]),
    var_type = unlist(purrr::map(seu[[]], class)),
    num_levels = unlist(purrr::map(seu[[]], ~length(unique(.x))))
  )

  meta_types <-
    meta_types %>%
    dplyr::filter(!grepl("_snn_res", vars)) %>%
    dplyr::mutate(meta_type = case_when(var_type %in% c("integer", "numeric") ~ "continuous",
                                        var_type %in% c("character", "factor") ~ "category")) %>%
    dplyr::mutate(meta_type = ifelse(meta_type == "continuous" & num_levels < 30, "category", meta_type)) %>%
    dplyr::filter(num_levels > 1) %>%
    identity()

  continuous_vars <- meta_types %>%
    dplyr::filter(meta_type == "continuous") %>%
    dplyr::pull(vars)

  continuous_vars <- c("custom", continuous_vars)

  continuous_names <- stringr::str_to_title(gsub("[[:punct:]]", " ", continuous_vars))

  names(continuous_vars) <- continuous_names


  category_vars <- meta_types %>%
    dplyr::filter(meta_type == "category") %>%
    dplyr::pull(vars)

  category_vars <- c("seurat", category_vars)

  category_names <- stringr::str_to_title(gsub("[[:punct:]]", " ", category_vars))

  names(category_vars) <- category_names

  plot_types <- list(category_vars = category_vars, continuous_vars = continuous_vars)

  return(plot_types)


}
<<<<<<< HEAD

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

#' Get Transcripts in Seurat Object
#'
#' @param seu
#' @param gene
#' @param organism
#'
#' @return
#' @export
#'
#' @examples
get_transcripts_from_seu <- function(seu, gene, organism = "human") {
  transcripts <- genes_to_transcripts(gene, organism)

  transcripts <- transcripts[transcripts %in%
                               rownames(GetAssay(seu$transcript, "RNA"))]
}
||||||| merged common ancestors
=======

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
>>>>>>> 2d43e2725a02764750212f5c9acf866ac6e6a483
