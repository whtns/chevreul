
#' Rename cell ids from annoying old notation
#'
#' @param cell_ids
#' @param batch_id
#'
#' @return
#' @export
#'
#' @examples
rename_from_x_notation <- function(cell_ids, batch_id){
  cell_ids <- str_replace(cell_ids, "X", "")
  cell_ids <- paste0(batch_id, str_pad(cell_ids, width = max(nchar(cell_ids)), pad = "0"))
}

#' delete seurat path
#'
#' @param seu_path
#'
#' @return
#' @export
#'
#' @examples
delete_seurat <- function(seu_path){
  print(seu_path)
  fs::file_delete(seu_path)
}

#' Reformat Seurat Object Metadata
#'
#'
#' @param seu
#' @param cols
#' @param new_col
#'
#' @return
#' @export
#'
#' @examples
combine_cols <- function(seu, cols, new_col) {

  new_col <- janitor::make_clean_names(new_col)

  # ensure that new colname will not be dropped
  drop_cols <- cols[!cols == new_col]

  # make sure that none of the columns to be coalesced are entirely NA
  na_cols <- purrr::map_lgl(cols, ~all(is.na(seu$gene[[.x]])))
  cols <- cols[!na_cols]

  # check class of cols to be coalesced


  meta <- tibble::rownames_to_column(seu$gene[[]]) %>%
    dplyr::mutate_at(vars(one_of(cols)), as.character) %>%
    dplyr::mutate(!!new_col := dplyr::coalesce(!!!syms(cols))) %>%
    dplyr::select(-drop_cols) %>%
    tibble::column_to_rownames(var = "rowname") %>%
    identity()
}

#' Filter Rows to Top
#'
#' @param df
#' @param column
#' @param values
#'
#' @return
#' @export
#'
#' @examples
filter_rows_to_top <- function(df, column, values){
  matched_df <- df[df[[column]] %in% values,]

  matched_df <- matched_df[match(values, matched_df[[column]]),]

  unmatched_df <- df[!(df[[column]] %in% values),]

  total_df <- list(matched_df = matched_df, unmatched_df = unmatched_df)
  total_df <- dplyr::bind_rows(total_df)

  return(total_df)

}

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

#' pad sample numbers to prep for armor
#'
#' @param proj_dir
#'
#' @return
#' @export
#'
#' @examples
pad_sample_files <- function(proj_dir) {
  fastq_files <-
    fs::dir_ls(fs::path(proj_dir, "data", "FASTQ"), glob = "*.fastq.gz") %>%
    purrr::set_names(path_file(.)) %>%
    tibble::enframe("file", "path") %>%
    tidyr::separate(file, into = c("proj_id", "sample_number", "pair"), sep = "-|_") %>%
    dplyr::mutate(sample_number = str_pad(sample_number, max(nchar(sample_number)), pad = "0")) %>%
    dplyr::mutate(file = fs::path(proj_dir, "data", "FASTQ", paste0(proj_id, "-", sample_number, "_", pair))) %>%
    # dplyr::mutate() %>%
    identity()

  purrr::map2(fastq_files$path, fastq_files$file, file_move)

}

#' prep armor metadata file
#'
#' @param proj_dir
#'
#' @return
#' @export
#'
#' @examples
prep_armor_meta <- function(proj_dir){
  seu <- load_seurat_from_proj(proj_dir)
  meta <- seu[[1]]@meta.data %>%
    dplyr::mutate(names = sample_id) %>%
    dplyr::mutate(type = "PE")

  readr::write_tsv(meta, fs::path(proj_dir, "data", "metadata.txt"))
}

