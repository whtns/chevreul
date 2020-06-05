
#' Title
#'
#' @param datapath
#' @param header
#' @param row.names
#'
#' @return
#' @export
#'
#' @examples
format_new_metadata <- function(datapath, header, row.names = 1){
  new_meta <- read.csv(datapath, header = header, row.names = row.names) %>%
    dplyr::mutate(across(contains("snn"), as.factor) )
}

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

  meta_types <- tibble::tibble(
    vars = colnames(seu[[]]),
    var_type = unlist(purrr::map(seu[[]], class)),
    num_levels = unlist(purrr::map(seu[[]], ~length(unique(.x))))
  )

  meta_types <-
    meta_types %>%
    dplyr::filter(!grepl("_snn_res", vars)) %>%
    dplyr::mutate(meta_type = dplyr::case_when(var_type %in% c("integer", "numeric") ~ "continuous",
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
                               rownames(GetAssay(seu, "RNA"))]
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

#' Title
#'
#' @param cds
#' @param mygenes
#' @param resolution
#'
#' @return
#' @export
#'
#' @examples
prep_plot_genes_in_pseudotime <- function(cds, mygenes, resolution, partition = FALSE) {
  if (partition){
    partition_cells <- monocle3::partitions(cds)
    # partition_cells <-  split(names(partition_cells), partition_cells)[[input$partitions]]
    partition_cells <-  split(names(partition_cells), partition_cells)[[1]]

    cds = cds[, colnames(cds) %in% partition_cells]

  }

  cds <- cds[rownames(cds) %in% mygenes,]

  if (any(grepl("integrated", colnames(colData(cds))))){
    default_assay = "integrated"
  } else {
    default_assay = "RNA"
  }

  color_cells_by = paste0(default_assay, "_snn_res.", resolution)

  gene_ptime_plot <- monocle3::plot_genes_in_pseudotime(cds,
                                                        color_cells_by=color_cells_by,
                                                        min_expr=0.5)

  return(gene_ptime_plot)
}

#' Record Experiment Metadata
#'
#' @param object
#' @param experiment_name
#' @param organism
#' @param column_sample
#' @param column_cluster
#' @param column_nUMI
#' @param column_nGene
#' @param column_cell_cycle_seurat
#' @param column_cell_cycle_cyclone
#' @param add_all_meta_data
#'
#' @return
#' @export
#'
#' @examples
record_experiment_data <- function(object, experiment_name, organism){
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' needed for this function to work. Please install it.",
         call. = FALSE)
  }


  message(paste0("[", format(Sys.time(), "%H:%M:%S"), "] Initializing Cerebro object..."))
  export <- list(experiment = list(experiment_name = experiment_name,
                                   organism = organism))
  export$experiment$date_of_analysis <- object@misc$experiment$date_of_analysis
  export$experiment$date_of_export <- Sys.Date()

  object@misc$experiment <- list(
    experiment_name = experiment_name,
    organism = organism,
    date_of_analysis = Sys.Date()
  )
  object@misc$parameters <- list(
    gene_nomenclature = 'gene_symbol',
    discard_genes_expressed_in_fewer_cells_than = 10,
    keep_mitochondrial_genes = TRUE,
    variables_to_regress_out = 'ncount_RNA',
    number_PCs = 30,
    tSNE_perplexity = 30,
    cluster_resolution = seq(0.2, 2.0, by = 0.2)
  )
  object@misc$parameters$filtering <- list(
    UMI_min = 50,
    UMI_max = Inf,
    genes_min = 10,
    genes_max = Inf
  )
  object@misc$technical_info <- list(
    'R' = capture.output(session_info())
  )

  if (!is.null(object@version)) {
    export$technical_info$seurat_version <- object@version
  }

  export$technical_info$seuratTools_version <- utils::packageVersion("seuratTools")



  object@misc <- c(object@misc, export)

  return(object)
}

#' Title
#'
#' @param seu
#' @param transcripts
#'
#' @return
#' @export
#'
#' @examples
plot_all_transcripts <- function(seu_transcript, seu_gene, features, embedding){

  # browser()
  transcript_cols <- as.data.frame(t(as.matrix(seu_transcript[["RNA"]][features,])))

  cells <- rownames(transcript_cols)
  transcript_cols <- as.list(transcript_cols) %>%
    purrr::map(~purrr::set_names(.x, cells))

  seu_gene[[features]] <- transcript_cols

  pList <- purrr::map(features, ~plot_feature(seu_gene,
                                                   embedding = embedding, features = .x))
  names(pList) <- features

  return(pList)

}
