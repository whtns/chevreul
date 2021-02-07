
#' Title
#'
#' @param seu
#' @param datapath
#'
#' @return
#' @export
#'
#' @examples
format_new_metadata <- function(seu, datapath) {
  new_meta <- read_csv(datapath) %>%
    dplyr::mutate(across(contains("snn"), as.factor))

  rowname_col <- colnames(new_meta)[1]

  new_meta <- tibble::column_to_rownames(new_meta, rowname_col)

  seu <- Seurat::AddMetaData(seu, new_meta)
  DefaultAssay(seu) <- "gene"
  ncalc <- Seurat:::CalcN(seu)
  seu$nFeature_RNA <- ncalc$nFeature
  seu$nCount_RNA <- ncalc$nCount

  return(seu)
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
  na_cols <- purrr::map_lgl(cols, ~ all(is.na(seu[[.x]])))
  cols <- cols[!na_cols]

  # check class of cols to be coalesced


  meta <- tibble::rownames_to_column(seu[[]]) %>%
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
filter_rows_to_top <- function(df, column, values) {
  matched_df <- df[df[[column]] %in% values, ]

  matched_df <- matched_df[match(values, matched_df[[column]]), ]

  unmatched_df <- df[!(df[[column]] %in% values), ]

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
list_plot_types <- function(seu) {
  meta_types <- tibble::tibble(
    vars = colnames(seu[[]]),
    var_type = purrr::map_chr(purrr::map(seu[[]], pillar::new_pillar_type), "type"),
    num_levels = unlist(purrr::map(seu[[]], ~ length(unique(.x))))
  )

  meta_types <- meta_types %>%
    dplyr::filter(!grepl("_snn_res", vars)) %>%
    dplyr::mutate(meta_type = dplyr::case_when(
      var_type %in% c("int", "dbl") ~ "continuous",
      var_type %in% c("chr", "fct", "ord") ~ "category"
    )) %>%
    dplyr::mutate(meta_type = ifelse(meta_type == "continuous" & num_levels < 30, "category", meta_type)) %>%
    dplyr::filter(num_levels > 1) %>%
    identity()

  continuous_vars <- meta_types %>%
    dplyr::filter(meta_type == "continuous") %>%
    dplyr::pull(vars)

  continuous_vars <- c("custom", continuous_vars) %>%
    purrr::set_names(stringr::str_to_title(stringr::str_replace(., "[[:punct:]]", " ")))


  category_vars <- meta_types %>%
    dplyr::filter(meta_type == "category") %>%
    dplyr::pull(vars)

  category_vars <- c("seurat", category_vars) %>%
    purrr::set_names(stringr::str_to_title(stringr::str_replace(., "[[:punct:]]", " ")))

  plot_types <- list(category_vars = category_vars, continuous_vars = continuous_vars)



  return(plot_types)
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
#' RXRG_transcripts <- get_transcripts_from_seu(human_gene_transcript_seu, "RXRG")
get_transcripts_from_seu <- function(seu, gene, organism = "human") {
  transcripts <- genes_to_transcripts(gene, organism)

  transcripts <- transcripts[transcripts %in% rownames(GetAssay(seu, "transcript"))]
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
  if (partition) {
    partition_cells <- monocle3::partitions(cds)
    # partition_cells <-  split(names(partition_cells), partition_cells)[[input$partitions]]
    partition_cells <- split(names(partition_cells), partition_cells)[[1]]

    cds <- cds[, colnames(cds) %in% partition_cells]
  }

  cds <- cds[rownames(cds) %in% mygenes, ]

  if (any(grepl("integrated", colnames(colData(cds))))) {
    default_assay <- "integrated"
  } else {
    default_assay <- "gene"
  }

  color_cells_by <- paste0(default_assay, "_snn_res.", resolution)

  gene_ptime_plot <- monocle3::plot_genes_in_pseudotime(cds,
    color_cells_by = color_cells_by,
    min_expr = 0.5
  )

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
#' logged_seu <- record_experiment_data(human_gene_transcript_seu, experiment_name = "human_gene_transcript", organism = "mouse")
#' Misc(logged_seu, "experiment")
#' @importFrom purrr %||%
record_experiment_data <- function(object, experiment_name = "default_experiment", organism = "human") {
  if (!requireNamespace("Seurat", quietly = TRUE)) {
    stop("Package 'Seurat' needed for this function to work. Please install it.",
      call. = FALSE
    )
  }

  organism <- Seurat::Misc(object, "experiment")[["organism"]] %||% organism

  experiment_name <- Seurat::Misc(object, "experiment")[["experiment_name"]] %||% experiment_name

  message(paste0("[", format(Sys.time(), "%H:%M:%S"), "] Logging Technical Details..."))
  experiment <- list(
    experiment_name = experiment_name,
    organism = organism
  )
  experiment$date_of_export <- Sys.Date()
  experiment$date_of_analysis <- Sys.Date()

  experiment$parameters <- list(
    gene_nomenclature = "gene_symbol",
    discard_genes_expressed_in_fewer_cells_than = 10,
    keep_mitochondrial_genes = TRUE,
    variables_to_regress_out = "nCount_RNA",
    number_PCs = 30,
    tSNE_perplexity = 30,
    cluster_resolution = seq(0.2, 2.0, by = 0.2)
  )
  experiment$filtering <- list(
    UMI_min = 50,
    genes_min = 10
  )
  experiment$session_info <- list(
    capture.output(sessioninfo::session_info())
  )

  if (!is.null(object@version)) {
    experiment$seurat_version <- object@version
  }

  experiment$seuratTools_version <- utils::packageVersion("seuratTools")

  Seurat::Misc(object, "experiment") <- NULL
  Seurat::Misc(object, "experiment") <- experiment

  return(object)
}


#' Update a SeuratTools Object
#'
#' @param seu_path
#' @param feature
#' @param resolution
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
update_seuratTools_object <- function(seu_path, feature, resolution = seq(0.2, 2.0, by = 0.2), return_seu = TRUE, ...) {
  message(seu_path)
  seu <- readRDS(seu_path)

  if (is.list(seu)) {
    seu <- convert_seu_list_to_multimodal(seu)
    # seu <- Seurat::UpdateSeuratObject(seu)
  } else if (all(names(seu@assays) == "RNA")) {
    seu <- RenameAssays(seu, RNA = "gene")
  } else if (identical(names(seu@assays), c("RNA", "integrated"))) {
    seu <- RenameAssays(seu, RNA = "gene")
  }

  # set appropriate assay
  if ("integrated" %in% names(seu@assays)) {
    default_assay <- "integrated"
  } else {
    default_assay <- "gene"
  }

  DefaultAssay(seu) <- default_assay

  seuratTools_version <- seu@misc$experiment$seuratTools_version

  seuratTools_version <- ifelse(is.null(seuratTools_version), 0.1, seuratTools_version)

  if (seuratTools_version < getNamespaceVersion("seuratTools")) {
    if (!any(grepl("_snn_res", colnames(seu@meta.data)))) {
      seu <- seurat_cluster(seu = seu, resolution = resolution, reduction = "pca", ...)
    }

    seu <- find_all_markers(seu, resolution = resolution)
    seu <- record_experiment_data(seu, ...)
    seu <- seu_calcn(seu)
  }

  if (return_seu) {
    return(seu)
  } else {
    message(paste0("saving ", seu_path))
    # saveRDS(seu, gsub(".rds", "_multimodal.rds", seu_path))
    saveRDS(seu, seu_path)
  }
}


#' Calculate Read Count Metrics for a Seurat object
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
seu_calcn <- function(seu, assay = "gene", slot = "counts") {
  n.calc <- Seurat:::CalcN(object = GetAssay(seu, assay))
  if (!is.null(x = n.calc)) {
    names(x = n.calc) <- paste(names(x = n.calc), assay, sep = "_")
    seu[[names(x = n.calc)]] <- n.calc
  }

  return(seu)
}

#' Propagate Metadata Changes
#'
#' @param updated_table
#' @param seu
#'
#' @return
#' @export
#'
#' @examples
propagate_spreadsheet_changes <- function(updated_table, seu) {
  meta <- updated_table

  sample_ids <- rownames(meta)

  meta <- meta %>%
    dplyr::mutate(meta, across(contains("snn"), as.factor)) %>%
    mutate(across(where(is.ordered), ~ as.factor(as.character(.x))))

  rownames(meta) <- sample_ids

  seu@meta.data <- meta

  return(seu)
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
                              destfile = "single-cell-projects.db", verbose = TRUE) {
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
                              verbose = TRUE) {
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
make_bigwig_db <- function(destdir = "~/.cache/seuratTools/", destfile = "bw-files.db") {
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
                                db_path = "single-cell-projects.db") {
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

#' Swap counts from Feature
#'
#' @param cds
#' @param featureType
#'
#' @return
#' @export
#'
#' @examples
swap_counts_from_feature <- function(cds, featureType) {
  print(featureType)
  #
  #   if (featureType == "transcript"){
  #     rowData(cds[[featureType]])$gene_short_name <- rownames(cds[[featureType]])
  #   }

  assay(cds$traj, withDimnames = FALSE) <- assay(cds[[featureType]])
  rowData(cds$traj) <- rowData(cds[[featureType]])
  rownames(cds$traj) <- rownames(cds[[featureType]])
  cds$traj@preprocess_aux$gene_loadings <- cds[[featureType]]@preprocess_aux$gene_loadings
  # counts(cds$traj) <- counts(cds[[featureType]])
  cds$traj
}

#' convert seurat list to multimodal object
#'
#' @param seu_list
#'
#' @return
#' @export
#'
#' @examples
convert_seu_list_to_multimodal <- function(seu_list) {
  colnames(seu_list[["gene"]]@meta.data) <- gsub("RNA_", "gene_", colnames(seu_list[["gene"]]@meta.data))

  multimodal_seu <- seu_list$gene
  multimodal_seu <- RenameAssays(multimodal_seu, RNA = "gene")

  if ("transcript" %in% names(seu_list)) {
    if (identical(length(Cells(seu_list$gene)), length(Cells(seu_list$transcript)))) {
      colnames(seu_list[["transcript"]]@meta.data) <- gsub("RNA_", "transcript_", colnames(seu_list[["transcript"]]@meta.data))
      multimodal_seu[["transcript"]] <- seu_list$transcript$RNA
      transcript_markers <- grepl("transcript_", names(seu_list$transcript@meta.data))
      transcript_cluster_cols <- seu_list[["transcript"]]@meta.data[transcript_markers]
      if (length(transcript_cluster_cols) > 0) {
        multimodal_seu <- AddMetaData(multimodal_seu, transcript_cluster_cols)
      }
    }
  }

  marker_names <- names(Misc(multimodal_seu)[["markers"]])

  if (!is.null(multimodal_seu@misc$markers)) {
    names(multimodal_seu@misc$markers) <- gsub("RNA", "gene", marker_names)
  }


  return(multimodal_seu)
}
