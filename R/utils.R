
#' Add object Metadata
#'
#' Adds data to the object to produce a object with metadata added.
#'
#' @param object A object
#' @param datapath Path to file containing metadata
#'
#' @return
#' @export
#'
#' @examples
#'
setGeneric("format_new_metadata", function(object, datapath)
  standardGeneric("format_new_metadata"))

setMethod("format_new_metadata", "Seurat",
          function(object, datapath) {
            new_meta <- read_csv(datapath) %>%
              dplyr::mutate(across(contains("snn"), as.factor))

            rowname_col <- colnames(new_meta)[1]

            new_meta <- tibble::column_to_rownames(new_meta, rowname_col)

            object <- object::AddMetaData(object, new_meta)
            DefaultAssay(object) <- "gene"
            ncalc <- object:::CalcN(object)
            object$nFeature_RNA <- ncalc$nFeature
            object$nCount_RNA <- ncalc$nCount

            return(object)
          }
          )

setMethod("format_new_metadata", "SingleCellExperiment",
          function(object, datapath) {
            new_meta <- read_csv(datapath) %>%
              dplyr::mutate(across(contains("snn"), as.factor))

            rowname_col <- colnames(new_meta)[1]

            new_meta <- tibble::column_to_rownames(new_meta, rowname_col)

            object <- object::AddMetaData(object, new_meta)
            DefaultAssay(object) <- "gene"
            ncalc <- object:::CalcN(object)
            object$nFeature_RNA <- ncalc$nFeature
            object$nCount_RNA <- ncalc$nCount

            return(object)
          }
          )

#' Reformat object Metadata
#'
#' Reformat object Metadata by Coalesced Columns
#'
#' @param object A object
#' @param cols Columns
#' @param new_col New columns
#'
#' @return
#' @export
#'
#' @examples
setGeneric("combine_cols", function(object, cols, new_col) {
  standardGeneric("combine_cols")
})

setMethod(
  "combine_cols", "Seurat",
  function(object, cols, new_col) {
    new_col <- janitor::make_clean_names(new_col)
    drop_cols <- cols[!cols == new_col]
    na_cols <- purrr::map_lgl(cols, ~ all(is.na(object[[.x]])))
    cols <- cols[!na_cols]
    meta <- tibble::rownames_to_column(pull_metadata(object)) %>%
      dplyr::mutate_at(vars(one_of(cols)), as.character) %>%
      dplyr::mutate(`:=`(!!new_col, dplyr::coalesce(!!!syms(cols)))) %>%
      dplyr::select(-drop_cols) %>%
      tibble::column_to_rownames(var = "rowname") %>%
      identity()
  }
)

setMethod(
  "combine_cols", "SingleCellExperiment",
  function(object, cols, new_col) {
    new_col <- janitor::make_clean_names(new_col)
    drop_cols <- cols[!cols == new_col]
    na_cols <- purrr::map_lgl(cols, ~ all(is.na(object[[.x]])))
    cols <- cols[!na_cols]
    meta <- tibble::rownames_to_column(pull_metadata(object)) %>%
      dplyr::mutate_at(vars(one_of(cols)), as.character) %>%
      dplyr::mutate(`:=`(!!new_col, dplyr::coalesce(!!!syms(cols)))) %>%
      dplyr::select(-drop_cols) %>%
      tibble::column_to_rownames(var = "rowname") %>%
      identity()
  }
)



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


#' Get Transcripts in object
#'
#' Get transcript ids in objects for one or more gene of interest
#'
#' @param object A object
#' @param gene Gene of intrest
#' @param organism Organism
#'
#' @return
#' @export
#'
#' @examples
#' RXRG_transcripts <- get_transcripts_from_object(human_gene_transcript_object, "RXRG")
#'
get_transcripts_from_object <- function(object, gene, organism = "human") {
  transcripts <- genes_to_transcripts(gene, organism)

  transcripts <- transcripts[transcripts %in% rownames(GetAssay(object, "transcript"))]
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
#' @param object A object objet
#' @param experiment_name
#' @param organism
#'
#' @return
#' @export
#'
#' @examples
#' logged_object <- record_experiment_data(human_gene_transcript_object, experiment_name = "human_gene_transcript", organism = "mouse")
#' Misc(logged_object, "experiment")
#' @importFrom purrr %||%
setGeneric("record_experiment_data", function(object, experiment_name = "default_experiment", organism = "human")
  standardGeneric("record_experiment_data"))

setMethod("record_experiment_data", "Seurat",
          function(object, experiment_name = "default_experiment", organism = "human") {
            if (!requireNamespace("Seurat", quietly = TRUE)) {
              stop("Package 'object' needed for this function to work. Please install it.",
                   call. = FALSE
              )
            }

            organism <- object::Misc(object, "experiment")[["organism"]] %||% organism

            experiment_name <- object::Misc(object, "experiment")[["experiment_name"]] %||% experiment_name

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

            experiment$chevreul_version <- utils::packageVersion("chevreul")

            Misc(object)[["experiment"]] <- NULL
            Misc(object)[["experiment"]] <- experiment

            return(object)
          }
          )

setMethod("record_experiment_data", "SingleCellExperiment",
          function(object, experiment_name = "default_experiment", organism = "human") {
            if (!requireNamespace("Seurat", quietly = TRUE)) {
              stop("Package 'object' needed for this function to work. Please install it.",
                   call. = FALSE
              )
            }

            organism <- object::Misc(object, "experiment")[["organism"]] %||% organism

            experiment_name <- object::Misc(object, "experiment")[["experiment_name"]] %||% experiment_name

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

            experiment$chevreul_version <- utils::packageVersion("chevreul")

            Misc(object)[["experiment"]] <- NULL
            Misc(object)[["experiment"]] <- experiment

            return(object)
          }
          )

#' Update a chevreul Object
#'
#' @param object_path Path to a object
#' @param feature
#' @param resolution
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
update_chevreul_object <- function(object_path, feature, resolution = seq(0.2, 2.0, by = 0.2), return_object = TRUE, ...) {
  message(object_path)
  object <- readRDS(object_path)

  if (is.list(object)) {
    object <- convert_object_list_to_multimodal(object)
    # object <- object::UpdateobjectObject(object)
  } else if (all(names(object@assays) == "RNA")) {
    object <- RenameAssays(object, RNA = "gene")
  } else if (identical(names(object@assays), c("RNA", "integrated"))) {
    object <- RenameAssays(object, RNA = "gene")
  }

  seurat_version <- Misc(object)$experiment$seurat_version

  if(packageVersion("Seurat") == '5.0.0' & (seurat_version < 5 || is.null(seurat_version))){
    object <- convert_v3_to_v5(object)
  }

  object <- propagate_spreadsheet_changes(pull_metadata(object), object)

  # set appropriate assay
  if ("integrated" %in% names(object@assays)) {
    default_assay <- "integrated"
  } else {
    default_assay <- "gene"
  }

  DefaultAssay(object) <- default_assay

  cluster_tag <- glue::glue("{DefaultAssay(object)}_snn_res\\.")

  cluster_names <- str_subset(names(pull_metadata(object)), cluster_tag)
  new_cluster_names <- str_replace(cluster_names, cluster_tag, "cluster_resolution_")

  new_cluster_cols <- pull_metadata(object)[cluster_names]
  names(new_cluster_cols) <- new_cluster_names

  new_meta <- cbind(pull_metadata(object), new_cluster_cols)

  object <- set_metadata(object, new_meta)

  chevreul_version <- Misc(object)$experiment$chevreul_version

  chevreul_version <- ifelse(is.null(chevreul_version), 0.1, chevreul_version)

  # update human gene symbols to grch38
  old_symbol <- "CTC-378H22.2"
  if (old_symbol %in% rownames(object[["gene"]])){
    for (i in names(object@assays)[names(object@assays) %in% c("gene", "integrated")]){
      # object <- update_human_gene_symbols(object, assay = i)
    }
  }

  if (chevreul_version < getNamespaceVersion("chevreul")) {
    message(paste0(object_path, " is out of date! updating..."))
    if (!any(grepl("_snn_res", colnames(pull_metadata(object))))) {
      object <- object_cluster(object = object, resolution = resolution, reduction = "pca", ...)
    }

    for (i in names(object@assays)[names(object@assays) %in% c("gene", "integrated")]){
      object <- find_all_markers(object, object_assay = i)
    }

    object <- record_experiment_data(object, ...)
    object <- object_calcn(object)
  }


  if (return_object) {
    return(object)
  } else {
    message(paste0("saving ", object_path))
    # saveRDS(object, gsub(".rds", "_multimodal.rds", object_path))
    saveRDS(object, object_path)
  }
}


#' Calculate Read Count Metrics for a object
#'
#' Recalculate counts/features per cell for a object
#'
#' @param object A object
#' @param assay Assay to use, Default = "gene"
#' @param slot
#'
#' @return
#' @export
#'
#' @examples
object_calcn <- function(object, assay = "gene", slot = "counts") {
  n.calc <- object:::CalcN(object = GetAssay(object, assay))
  if (!is.null(x = n.calc)) {
    names(x = n.calc) <- paste(names(x = n.calc), assay, sep = "_")
    object[[names(x = n.calc)]] <- n.calc
  }

  return(object)
}

#' Propagate Metadata Changes
#'
#' @param updated_table
#' @param object
#'
#' @return
#' @export
#'
#' @examples
propagate_spreadsheet_changes <- function(updated_table, object) {
  meta <- updated_table

  sample_ids <- rownames(meta)

  meta <- meta %>%
    dplyr::mutate(meta, across(contains("snn"), as.factor)) %>%
    mutate(across(where(is.ordered), ~ as.factor(as.character(.x))))

  rownames(meta) <- sample_ids

  object <- set_metadata(object, meta)

  return(object)
}

#' Create a database of chevreul projects
#'
#' Create a database containing chevreul projects
#'
#' @param cache_location Path to cache "~/.cache/chevreul"
#' @param sqlite_db Database to be created
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
create_project_db <- function(cache_location = "~/.cache/chevreul",
                              sqlite_db = "single-cell-projects.db", verbose = TRUE) {
  if (!dir.exists(cache_location)) {
    dir.create(cache_location)
  }

  con <- DBI::dbConnect(RSQLite::SQLite(), fs::path(cache_location, sqlite_db))

  projects_tbl <- tibble::tibble(
    project_name = character(),
    project_path = character(),
    project_slug = character(),
    project_type = character(),
  )

  message(paste0("building table of chevreul projects at ", fs::path(cache_location, sqlite_db)))
  # DBI::dbWriteTable(con, "projects_tbl", projects_tbl)

  tryCatch({
    DBI::dbWriteTable(con, "projects_tbl", projects_tbl)
  }, warning = function(w) {
      message(sprintf("Warning in %s: %s", deparse(w[["call"]]), w[["message"]]))

  }, error = function(e) {
      message("projects db already exists!")

  }, finally = {
  })

  DBI::dbDisconnect(con)
}

#' Update a database of chevreul projects
#'
#' Add new/update existing projects to the database by recursing fully
#'
#' @param projects_dir The project directory to be updated
#' @param cache_location Path to cache "~/.cache/chevreul"
#' @param sqlite_db sqlite db
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
update_project_db <- function(projects_dir = NULL,
                              cache_location = "~/.cache/chevreul",
                              sqlite_db = "single-cell-projects.db",
                              verbose = TRUE) {
  if (!dir.exists(cache_location)) {
    dir.create(cache_location)
  }

  con <- DBI::dbConnect(RSQLite::SQLite(), fs::path(cache_location, sqlite_db))

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

#' Update a database of chevreul projects
#'
#' Append projects to datatbase
#'
#' @param new_project_path
#' @param projects_dir
#' @param cache_location Path to cache "~/.cache/chevreul"
#' @param sqlite_db sqlite db
#' @param verbose
#'
#'
#' @return
#' @export
#'
#' @examples
append_to_project_db <- function(new_project_path, projects_dir = NULL,
                              cache_location = "~/.cache/chevreul",
                              sqlite_db = "single-cell-projects.db",
                              verbose = TRUE) {
  if (!dir.exists(cache_location)) {
    dir.create(cache_location)
  }

  con <- DBI::dbConnect(RSQLite::SQLite(), fs::path(cache_location, sqlite_db))

  projects_tbl <-
    new_project_path %>%
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

  DBI::dbWriteTable(con, "projects_tbl", current_projects_tbl, overwrite = TRUE)

  DBI::dbDisconnect(con)
}

#' Read a database of chevreul projects
#'
#' Reads database of chevreul projects to a data frame
#'
#' @param projects_dir
#' @param cache_location Path to cache "~/.cache/chevreul"
#' @param sqlite_db sqlite db
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
read_project_db <- function(projects_dir = NULL,
                              cache_location = "~/.cache/chevreul",
                              sqlite_db = "single-cell-projects.db",
                              verbose = TRUE) {
  if (!dir.exists(cache_location)) {
    dir.create(cache_location)
  }

  con <- DBI::dbConnect(RSQLite::SQLite(), fs::path(cache_location, sqlite_db))

  current_projects_tbl <-
    DBI::dbReadTable(con, "projects_tbl")

  DBI::dbDisconnect(con)

  return(current_projects_tbl)
}

#' Make Bigwig Database
#'
#'
#' @param new_project Project directory
#' @param cache_location Path to cache "~/.cache/chevreul"
#' @param sqlite_db sqlite db containing bw files
#'
#' @return
#' @export
#'
#' @examples
make_bigwig_db <- function(new_project = NULL, cache_location = "~/.cache/chevreul/", sqlite_db = "bw-files.db") {
  new_bigwigfiles <- fs::dir_ls(new_project, glob = "*.bw", recurse = TRUE) %>%
    purrr::set_names(fs::path_file(.)) %>%
    tibble::enframe("name", "bigWig") %>%
    dplyr::mutate(sample_id = stringr::str_remove(name, "_Aligned.sortedByCoord.out.*bw$")) %>%
    dplyr::filter(!stringr::str_detect(name, "integrated")) %>%
    dplyr::distinct(sample_id, .keep_all = TRUE) %>%
    identity()

  con <- DBI::dbConnect(RSQLite::SQLite(), dbname = fs::path(cache_location, sqlite_db))

  all_bigwigfiles <-
    dbReadTable(con, "bigwigfiles") %>%
    dplyr::bind_rows(new_bigwigfiles)

  DBI::dbWriteTable(con, "bigwigfiles", all_bigwigfiles, overwrite = TRUE)

  return(all_bigwigfiles)
}

#' Retrieve Metadata from Batch
#'
#' @param batch
#' @param projects_dir path to project dir
#' @param db_path path to .db file
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

#' convert object list to multimodal object
#'
#' @param object_list
#'
#' @return
#' @export
#'
#' @examples
convert_object_list_to_multimodal <- function(object_list) {
  colnames(object_list[["gene"]]@meta.data) <- gsub("RNA_", "gene_", colnames(object_list[["gene"]]@meta.data))

  multimodal_object <- object_list$gene
  multimodal_object <- RenameAssays(multimodal_object, RNA = "gene")

  if ("transcript" %in% names(object_list)) {
    if (identical(length(Cells(object_list$gene)), length(Cells(object_list$transcript)))) {
      colnames(object_list[["transcript"]]@meta.data) <- gsub("RNA_", "transcript_", colnames(object_list[["transcript"]]@meta.data))
      multimodal_object[["transcript"]] <- object_list$transcript$RNA
      transcript_markers <- grepl("transcript_", names(object_list$transcript@meta.data))
      transcript_cluster_cols <- object_list[["transcript"]]@meta.data[transcript_markers]
      if (length(transcript_cluster_cols) > 0) {
        multimodal_object <- AddMetaData(multimodal_object, transcript_cluster_cols)
      }
    }
  }

  marker_names <- names(Misc(multimodal_object)[["markers"]])

  if (!is.null(Misc(multimodal_object)$markers)) {
    names(Misc(multimodal_object)$markers) <- gsub("RNA", "gene", marker_names)
  }


  return(multimodal_object)
}

#' Clean Vector of Chevreul Names
#'
#' Cleans names of objects provided in a vector form
#'
#' @param myvec A vector of object names
#'
#' @return
#' @export
#'
#' @examples
make_chevreul_clean_names <- function(myvec){
  myvec %>%
    purrr::set_names(stringr::str_to_title(stringr::str_replace_all(., "[^[:alnum:][:space:]\\.]", " ")))
}


#' Convert seurat to seurat V5 format
#'
#' Convert seurat from v3 to v5 format
#'
#' @param seurat_v3 a version 3 seurat
#'
#' @return
#' @export
#'
#' @examples
#' convert_v3_to_v5(human_gene_transcript_seurat)
convert_v3_to_v5 <- function(seurat_v3){

  seurat_version <- Misc(seurat_v3)$experiment$seurat_version

  if(seurat_version < 5 || is.null(seurat_version)){
  meta <- seurat_v3@meta.data

  seurat_v5 <- CreateSeuratObject(counts = seurat_v3$gene@counts, data = seurat_v3$gene@data, assay = "gene", meta.data = meta)

  transcript_assay.v5 <- CreateAssay5Object(counts = seurat_v3$transcript@counts, data = seurat_v3$transcript@data)
  seurat_v5$transcript <- transcript_assay.v5

  seurat_v5$gene <- seurat_preprocess(seurat_v5$gene, normalize = FALSE)
  seurat_v5$transcript <- seurat_preprocess(seurat_v5$transcript, normalize = FALSE)

  # seurat_v5 <- clustering_workflow(seurat_v5)
  seurat_v5@reductions <- seurat_v3@reductions
  seurat_v5@graphs <- seurat_v3@graphs
  seurat_v5@neighbors <- seurat_v3@neighbors

  Misc(seurat_v5) <- Misc(seurat_v3)

  Idents(seurat_v5) <- Idents(seurat_v3)

  Misc(seurat_v5)$experiment$seurat_version <- packageVersion("seurat")
  } else{
   seurat_v5 <- seurat_v3
  }

  return(seurat_v5)

}
