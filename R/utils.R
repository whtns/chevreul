
#' Title
#'
#' @param seu
#' @param datapath
#'
#' @return
#' @export
#'
#' @examples
format_new_metadata <- function(seu, datapath){
  new_meta <- read_csv(datapath) %>%
    dplyr::mutate(across(contains("snn"), as.factor) )

  rowname_col <- colnames(new_meta)[1]

  new_meta <- tibble::column_to_rownames(new_meta, rowname_col)

  seu <- Seurat::AddMetaData(seu, new_meta)
  DefaultAssay(seu) <- "RNA"
  ncalc <- Seurat:::CalcN(seu)
  seu$nFeature_RNA <- ncalc$nFeature
  seu$nCount_RNA <- ncalc$nCount

  return(seu)
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
  cell_ids <- stringr::str_replace(cell_ids, "X", "")
  cell_ids <- paste0(batch_id, stringr::str_pad(cell_ids, width = max(nchar(cell_ids)), pad = "0"))
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
#' @export
#' @examples
list_plot_types <- function(seu){


  meta_types <- tibble::tibble(vars = colnames(seu[[]]),
                               var_type = purrr::map_chr(purrr::map(seu[[]], pillar::new_pillar_type), "type"),
                               num_levels = unlist(purrr::map(seu[[]], ~length(unique(.x)))))

  meta_types <- meta_types %>%
    dplyr::filter(!grepl("_snn_res", vars)) %>%
    dplyr::mutate(meta_type = dplyr::case_when(var_type %in% c("int", "dbl") ~ "continuous",
                                               var_type %in% c("chr", "fct", "ord") ~ "category")) %>%
    dplyr::mutate(meta_type = ifelse(meta_type == "continuous" & num_levels < 30, "category", meta_type)) %>%
    dplyr::filter(num_levels > 1) %>% identity()

  continuous_vars <- meta_types %>%
    dplyr::filter(meta_type =="continuous") %>%
    dplyr::pull(vars)

  continuous_vars <- c("custom", continuous_vars) %>%
    purrr::set_names(stringr::str_to_title(stringr::str_replace(., "[[:punct:]]", " ")))


  category_vars <- meta_types %>%
    dplyr::filter(meta_type =="category") %>%
    dplyr::pull(vars)

  category_vars <- c("seurat", category_vars) %>%
    purrr::set_names(stringr::str_to_title(stringr::str_replace(., "[[:punct:]]", " ")))

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
    dplyr::mutate(sample_number = stringr::str_pad(sample_number, max(nchar(sample_number)), pad = "0")) %>%
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
#' @importFrom purrr %||%
record_experiment_data <- function(object, experiment_name = "default_experiment", organism = "human"){
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' needed for this function to work. Please install it.",
         call. = FALSE)
  }

  organism <- object@misc$experiment$organism %||% organism

  experiment_name <- object@misc$experiment$experiment_name %||% experiment_name

  message(paste0("[", format(Sys.time(), "%H:%M:%S"), "] Logging Technical Details..."))
  export <- list(experiment = list(experiment_name = experiment_name,
                                   organism = organism))
  export$experiment$date_of_analysis <- object@misc$experiment$date_of_analysis
  export$experiment$date_of_export <- Sys.Date()
  export$experiment$date_of_analysis <- Sys.Date()

  export$experiment$parameters <- list(
    gene_nomenclature = 'gene_symbol',
    discard_genes_expressed_in_fewer_cells_than = 10,
    keep_mitochondrial_genes = TRUE,
    variables_to_regress_out = 'nCount_RNA',
    number_PCs = 30,
    tSNE_perplexity = 30,
    cluster_resolution = seq(0.2, 2.0, by = 0.2)
  )
  export$experiment$filtering <- list(
    UMI_min = 50,
    UMI_max = Inf,
    genes_min = 10,
    genes_max = Inf
  )
  export$experiment$technical_info <- list(
    # capture.output(sessioninfo::session_info())
  )

  if (!is.null(object@version)) {
    export$experiment$technical_info$seurat_version <- object@version
  }

  export$experiment$technical_info$seuratTools_version <- utils::packageVersion("seuratTools")

  object@misc <- c(object@misc, export)

  return(object)
}


#' Title
#'
#' @param seu
#' @param feature
#' @param resolution
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
update_seuratTools_object <- function(seu, feature, resolution = seq(0.2, 2.0, by = 0.2), ...){

  # set appropriate assay
  if ("integrated" %in% names(seu@assays)) {
    default_assay = "integrated"
  } else {
    default_assay = "RNA"
  }

  DefaultAssay(seu) <- default_assay

  seuratTools_version <- seu@misc$experiment$technical_info$seuratTools_version

  seuratTools_version <- ifelse(is.null(seuratTools_version), 0.1, seuratTools_version)

  if(seuratTools_version < getNamespaceVersion("seuratTools")){

    if (!any(grepl("_snn_res", colnames(seu@meta.data)))){
      seu <- seurat_cluster(seu = seu, resolution = resolution, reduction = "pca", ...)

    }
      seu <- find_all_markers(seu, resolution = resolution)
      seu <- record_experiment_data(seu, ...)
      seu <- seu_calcn(seu)
  }

  return(seu)

}


#' recalculate counts/features per cell for a seurat object
#'
#' @param seu
#' @param assay
#' @param slot
#'
#' @return
#' @export
#'
#' @examples
seu_calcn <- function(seu, assay = "RNA", slot = "counts"){

  n.calc <- Seurat:::CalcN(object = GetAssay(seu, assay))
  if (!is.null(x = n.calc)) {
    names(x = n.calc) <- paste(names(x = n.calc), assay, sep = '_')
    seu[[names(x = n.calc)]] <- n.calc
  }

  return(seu)
}

#' Propagate Metadata Changes
#'
#' @param meta
#' @param changes
#'
#' @return
#' @export
#'
#' @examples
propagate_spreadsheet_changes <- function(changes){
  meta <- rhandsontable::hot_to_r(changes)

  sample_ids <- rownames(meta)

  meta <- meta %>%
    dplyr::mutate(meta, across(contains("snn"), as.factor)) %>%
    mutate(across(where(is.ordered), ~as.factor(as.character(.x))))

  rownames(meta) <- sample_ids

  return(meta)
}

#' Create a database of seuratTools projects
#'
#' @param destdir
#' @param destfile
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
create_project_db <- function(destdir = "~/.cache/seuratTools",
          destfile = "single-cell-projects.db", verbose = TRUE){

  if (!dir.exists(destdir)) {
    dir.create(destdir)
  }

    con <- DBI::dbConnect(RSQLite::SQLite(), fs::path(destdir, destfile))

    projects_tbl <- tibble::tibble(
        project_name = character(),
        project_path = character(),
        project_slug = character(),
        project_type = character(),
      )

    message(paste0("building table of seuratTools projects at ", fs::path(destdir, destfile)))

    DBI::dbWriteTable(con, "projects_tbl", projects_tbl)

    DBI::dbDisconnect(con)

}

#' Update a database of seuratTools projects
#'
#' @param projects_dir
#' @param destdir
#' @param destfile
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
update_project_db <- function(projects_dir = NULL,
                              destdir = "~/.cache/seuratTools",
                              destfile = "single-cell-projects.db",
                              verbose = TRUE){


  if (!dir.exists(destdir)) {
    dir.create(destdir)
  }

  con <- DBI::dbConnect(RSQLite::SQLite(), fs::path(destdir, destfile))

  projects_tbl <-
    fs::dir_ls(projects_dir, glob = "*.here", recurse = TRUE, fail = FALSE, all = TRUE) %>%
    fs::path_dir(.) %>%
    purrr::set_names(fs::path_file(.)) %>%
    tibble::enframe("project_name", "project_path") %>%
    dplyr::mutate(project_slug = stringr::str_remove(project_name, "_proj$")) %>%
    dplyr::mutate(project_type = fs::path_file(fs::path_dir(project_path))) %>%
    identity()

  current_projects_tbl <-
    DBI::dbReadTable(con, "projects_tbl") %>%
    dplyr::filter(fs::file_exists(project_path)) %>%
    dplyr::filter(!project_path %in% projects_tbl$project_path) %>%
    dplyr::bind_rows(projects_tbl) %>%
    dplyr::distinct(project_path, .keep_all = TRUE)

  DBI::dbWriteTable(con, "projects_tbl", projects_tbl, overwrite = TRUE)

  DBI::dbDisconnect(con)
}

#' Make Bigwig Database
#'
#' @param destdir
#' @param destfile
#'
#' @return
#' @export
#'
#' @examples
make_bigwig_db <- function(destdir = "~/.cach/seuratTools/", destfile = "bw-files.db") {
  bigwigfiles <- dir_ls(destdir, glob = "*.bw", recurse = TRUE) %>%
    set_names(path_file(.)) %>%
    enframe("name", "bigWig") %>%
    dplyr::mutate(sample_id = str_remove(name, "_Aligned.sortedByCoord.out.bw")) %>%
    dplyr::filter(!str_detect(path, "integrated")) %>%
    dplyr::distinct(sample_id, .keep_all = TRUE) %>%
    identity()

  con <- dbConnect(RSQLite::SQLite(), dbname = fs::path(destdir, destfile))

  DBI::dbWriteTable(con, "bigwigfiles", bigwigfiles)
}

#' Retrieve Metadata from Batch
#'
#' @param batch
#' @param projects_dir
#' @param db_path
#'
#' @return
#'
#' @examples
metadata_from_batch <- function(batch, projects_dir = "/dataVolume/storage/single_cell_projects",
                                db_path = "single-cell-projects.db"){

  mydb <- DBI::dbConnect(RSQLite::SQLite(), fs::path(projects_dir, db_path))

  projects_tbl <- DBI::dbReadTable(mydb, "projects_tbl") %>%
    dplyr::filter(!project_type %in% c("integrated_projects", "resources"))

  DBI::dbDisconnect(mydb)

  metadata <-
    projects_tbl %>%
    dplyr::filter(project_slug == batch) %>%
    dplyr::pull(project_path) %>%
    fs::path("data") %>%
    fs::dir_ls(glob = "*.csv") %>%
    identity()

}


