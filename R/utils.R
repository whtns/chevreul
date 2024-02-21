#' Add object Metadata
#'
#' Adds data to the object to produce a object with metadata added.
#'
#' @param object A object
#' @param datapath Path to file containing metadata
#'
#' @return a single cell object
#' @export
#' @examples
format_new_metadata <- function(object, datapath) {
        new_meta <- read_csv(datapath) %>%
            mutate(across(contains("snn"), as.factor))

        rowname_col <- colnames(new_meta)[1]

        new_meta <- column_to_rownames(new_meta, rowname_col)

        colData(object) <- new_meta

        return(object)
    }

#' Reformat object Metadata
#'
#' Reformat object Metadata by Coalesced Columns
#'
#' @param object A object
#' @param cols Columns
#' @param new_col New columns
#'
#' @return updated cell level metadata
#' @export
#' @examples
combine_cols <- function(object, cols, new_col) {
        new_col <- make_clean_names(new_col)
        drop_cols <- cols[!cols == new_col]
        na_cols <- map_lgl(cols, ~ all(is.na(object[[.x]])))
        cols <- cols[!na_cols]
        meta <- rownames_to_column(get_cell_metadata(object)) %>%
            mutate_at(vars(one_of(cols)), as.character) %>%
            mutate(`:=`(!!new_col, coalesce(!!!syms(cols)))) %>%
            select(-drop_cols) %>%
            column_to_rownames(var = "rowname") %>%
            identity()
    }

#' Get Transcripts in object
#'
#' Get transcript ids in objects for one or more gene of interest
#'
#' @param object A object
#' @param gene Gene of intrest
#' @param organism Organism
#'
#' @return transcripts constituting a gene of interest in a single cell object
#' @export
#'
#' @examples
#'
#' NRL_transcripts <- get_transcripts_from_object(human_gene_transcript_sce, "NRL")
#'
get_transcripts_from_object <- function(object, gene, organism = "human") {
        transcripts <- genes_to_transcripts(gene, organism)

        transcripts <- transcripts[transcripts %in% rownames(retrieve_experiment(object, "transcript"))]
    }


#' Record Experiment Metadata
#'
#' Records miscellaneous data
#' @param object A object
#' @param experiment_name name of the experiment
#' @param organism human or mouse
#'
#' @return a single cell object
#' @export
#' @examples
#' record_experiment_data(human_gene_transcript_sce)
#'
record_experiment_data <- function(object, experiment_name = "default_experiment", organism = "human") {
        if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
            stop("Package 'object' needed for this function to work. Please install it.",
                call. = FALSE
            )
        }

        organism <- metadata(human_gene_transcript_sce)[["experiment"]][["organism"]] %||% organism

        experiment_name <- metadata(object)[["experiment"]][["experiment_name"]] %||% experiment_name

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
        experiment$sessionInfo <- list(
            capture.output(sessionInfo())
        )

        if (!is.null(objectVersion(object))) {
            experiment$SingleCellExperiment_version <- objectVersion(object)
        }

        experiment$chevreul_version <- packageVersion("chevreul")

        metadata(object)[["experiment"]] <- experiment

        return(object)
    }

#' Update a chevreul Object
#'
#' @param object_path Path to a object
#' @param feature gene or transcript of interest
#' @param resolution resolution for louvain clustering
#' #' @param return_object whether to return the single cell object
#' @param ... extra args passed to object_cluster
#'
#' @return a single cell object
#' @export
#' @importFrom glue glue
#' @examples
update_chevreul_object <- function(object_path, feature, resolution = seq(0.2, 2.0, by = 0.2), return_object = TRUE, ...) {
        message(object_path)
        object <- readRDS(object_path)

        if (is.list(object)) {
            object <- convert_object_list_to_multimodal(object)
            # object <- SingleCellExperiment::UpdateSingleCellExperimentObject(object)
        } else if (all(names(object@experiments) == "RNA")) {
            object <- RenameAssays(object, RNA = "gene")
        } else if (identical(names(object@experiments), c("RNA", "integrated"))) {
            object <- RenameAssays(object, RNA = "gene")
        }

        seurat_version <- metadata(object)$experiment$seurat_version

        if (packageVersion("SingleCellExperiment") == "5.0.0" & (seurat_version < 5 || is.null(seurat_version))) {
            object <- convert_v3_to_v5(object)
        }

        object <- propagate_spreadsheet_changes(get_cell_metadata(object), object)

        # set appropriate experiment
        if ("integrated" %in% names(object@experiments)) {
            default_experiment <- "integrated"
        } else {
            default_experiment <- "gene"
        }

        DefaultAssay(object) <- default_experiment

        cluster_tag <- glue("{DefaultAssay(object)}_snn_res\\.")

        cluster_names <- str_subset(names(get_cell_metadata(object)), cluster_tag)
        new_cluster_names <- str_replace(cluster_names, cluster_tag, "cluster_resolution_")

        new_cluster_cols <- get_cell_metadata(object)[cluster_names]
        names(new_cluster_cols) <- new_cluster_names

        new_meta <- cbind(get_cell_metadata(object), new_cluster_cols)

        object <- set_metadata(object, new_meta)

        chevreul_version <- metadata(object)$experiment$chevreul_version

        chevreul_version <- ifelse(is.null(chevreul_version), 0.1, chevreul_version)

        # update human gene symbols to grch38
        old_symbol <- "CTC-378H22.2"
        if (old_symbol %in% rownames(object[["gene"]])) {
            for (i in names(object@experiments)[names(object@experiments) %in% c("gene", "integrated")]) {
                # object <- update_human_gene_symbols(object, experiment = i)
            }
        }

        if (chevreul_version < getNamespaceVersion("chevreul")) {
            message(paste0(object_path, " is out of date! updating..."))
            if (!any(grepl("_snn_res", colnames(get_cell_metadata(object))))) {
                object <- object_cluster(object = object, resolution = resolution, reduction = "PCA", ...)
            }

            for (i in names(object@experiments)[names(object@experiments) %in% c("gene", "integrated")]) {
                object <- find_all_markers(object, experiment = i)
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
#' @param object A single cell object
#' @param experiment Assay to use, Default = "gene"
#' @param slot slot for data
#'
#' @return a single cell object with nfeatures and ngenes stored in metadata
#' @export
#' @examples
object_calcn <- function(object, experiment = "gene", slot = "counts") {

  object <- addPerCellQC(object)
  object$nFeature_gene <- object$detected
  object$nCount_gene <- object$sum

        return(object)
    }


#' Propagate Metadata Changes
#'
#' @param updated_table updated metadata
#' @param object a single cell object
#'
#' @return a single cell object
#' @export
#' @examples
propagate_spreadsheet_changes <- function(updated_table, object) {
    meta <- updated_table

    sample_ids <- rownames(meta)

    meta <- meta %>%
        mutate(meta, across(contains("snn"), as.factor)) %>%
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
#' @param verbose print messages
#'
#' @return a sqlite database with single cell objects
#' @export
#' @examples
create_project_db <- function(cache_location = "~/.cache/chevreul", sqlite_db = "single-cell-projects.db", verbose = TRUE) {
        if (!dir.exists(cache_location)) {
            dir.create(cache_location)
        }
        con <- dbConnect(SQLite(), path(cache_location, sqlite_db))
        projects_tbl <- tibble(project_name = character(), project_path = character(), project_slug = character(), project_type = character(), )
        message(paste0("building table of chevreul projects at ", path(cache_location, sqlite_db)))
        tryCatch({
            dbWriteTable(con, "projects_tbl", projects_tbl)
        }, warning = function(w) {
            message(sprintf("Warning in %s: %s", deparse(w[["call"]]), w[["message"]]))
        }, error = function(e) {
            message("projects db already exists!")
        }, finally = {
        })
        dbDisconnect(con)
}

#' Update a database of chevreul projects
#'
#' Add new/update existing projects to the database by recursing fully
#'
#' @param projects_dir The project directory to be updated
#' @param cache_location Path to cache "~/.cache/chevreul"
#' @param sqlite_db sqlite db
#' @param verbose print messages
#'
#' @return a sqlite database with single cell objects
#' @export
#' @examples
update_project_db <- function(projects_dir = NULL,
    cache_location = "~/.cache/chevreul",
    sqlite_db = "single-cell-projects.db",
    verbose = TRUE) {
    if (!dir.exists(cache_location)) {
        dir.create(cache_location)
    }

    con <- dbConnect(SQLite(), path(cache_location, sqlite_db))

    projects_tbl <-
        dir_ls(projects_dir, glob = "*.here", recurse = TRUE, fail = FALSE, all = TRUE) %>%
        path_dir(.) %>%
        set_names(path_file(.)) %>%
        enframe("project_name", "project_path") %>%
        mutate(project_slug = str_remove(project_name, "_proj$")) %>%
        mutate(project_type = path_file(path_dir(project_path))) %>%
        identity()

    current_projects_tbl <-
        dbReadTable(con, "projects_tbl") %>%
        filter(file.exists(project_path)) %>%
        filter(!project_path %in% projects_tbl$project_path) %>%
        bind_rows(projects_tbl) %>%
        distinct(project_path, .keep_all = TRUE)

    dbWriteTable(con, "projects_tbl", projects_tbl, overwrite = TRUE)

    dbDisconnect(con)
}

#' Update a database of chevreul projects
#'
#' Append projects to datatbase
#'
#' @param new_project_path new project path
#' @param projects_dir project directory
#' @param cache_location Path to cache "~/.cache/chevreul"
#' @param sqlite_db sqlite db
#' @param verbose print messages
#'
#'
#' @return a sqlite database with single cell objects
#' @export
#' @examples
append_to_project_db <- function(new_project_path, projects_dir = NULL,
    cache_location = "~/.cache/chevreul",
    sqlite_db = "single-cell-projects.db",
    verbose = TRUE) {
    if (!dir.exists(cache_location)) {
        dir.create(cache_location)
    }

    con <- dbConnect(SQLite(), path(cache_location, sqlite_db))

    projects_tbl <-
        new_project_path %>%
        set_names(path_file(.)) %>%
        enframe("project_name", "project_path") %>%
        mutate(project_slug = str_remove(project_name, "_proj$")) %>%
        mutate(project_type = path_file(path_dir(project_path))) %>%
        identity()

    current_projects_tbl <-
        dbReadTable(con, "projects_tbl") %>%
        filter(file.exists(project_path)) %>%
        filter(!project_path %in% projects_tbl$project_path) %>%
        bind_rows(projects_tbl) %>%
        distinct(project_path, .keep_all = TRUE)

    dbWriteTable(con, "projects_tbl", current_projects_tbl, overwrite = TRUE)

    dbDisconnect(con)
}

#' Read a database of chevreul projects
#'
#' Reads database of chevreul projects to a data frame
#'
#' @param projects_dir project directory
#' @param cache_location Path to cache "~/.cache/chevreul"
#' @param sqlite_db sqlite db
#' @param verbose print messages
#'
#' @return a tibble with single cell objects
#' @export
#' @examples
read_project_db <- function(projects_dir = NULL,
    cache_location = "~/.cache/chevreul",
    sqlite_db = "single-cell-projects.db",
    verbose = TRUE) {
    if (!dir.exists(cache_location)) {
        dir.create(cache_location)
    }

    con <- dbConnect(SQLite(), path(cache_location, sqlite_db))

    current_projects_tbl <-
        dbReadTable(con, "projects_tbl")

    dbDisconnect(con)

    return(current_projects_tbl)
}

#' Make Bigwig Database
#'
#'
#' @param new_project Project directory
#' @param cache_location Path to cache "~/.cache/chevreul"
#' @param sqlite_db sqlite db containing bw files
#'
#' @return a sqlite database of bigwig files for cells in a single cell object
#' @export
#' @examples
make_bigwig_db <- function(new_project = NULL, cache_location = "~/.cache/chevreul/", sqlite_db = "bw-files.db") {
    new_bigwigfiles <- dir_ls(new_project, glob = "*.bw", recurse = TRUE) %>%
        set_names(path_file(.)) %>%
        enframe("name", "bigWig") %>%
        mutate(sample_id = str_remove(name, "_Aligned.sortedByCoord.out.*bw$")) %>%
        filter(!str_detect(name, "integrated")) %>%
        distinct(sample_id, .keep_all = TRUE) %>%
        identity()

    con <- dbConnect(SQLite(), dbname = path(cache_location, sqlite_db))

    all_bigwigfiles <-
        dbReadTable(con, "bigwigfiles") %>%
        bind_rows(new_bigwigfiles)

    dbWriteTable(con, "bigwigfiles", all_bigwigfiles, overwrite = TRUE)

    return(all_bigwigfiles)
}

#' Retrieve Metadata from Batch
#'
#' @param batch batch
#' @param projects_dir path to project dir
#' @param db_path path to .db file
#'
#' @return a tibble with cell level metadata from a single cell object
#' @examples
metadata_from_batch <- function(batch, projects_dir = "/dataVolume/storage/single_cell_projects",
    db_path = "single-cell-projects.db") {
    mydb <- dbConnect(SQLite(), path(projects_dir, db_path))

    projects_tbl <- dbReadTable(mydb, "projects_tbl") %>%
        filter(!project_type %in% c("integrated_projects", "resources"))

    dbDisconnect(mydb)

    metadata <-
        projects_tbl %>%
        filter(project_slug == batch) %>%
        pull(project_path) %>%
        path("data") %>%
        dir_ls(glob = "*.csv") %>%
        identity()
}

#' convert object list to multimodal object
#'
#' @param object_list a list of objects
#'
#' @return a list of single cell objects
#' @export
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

    marker_names <- names(metadata(multimodal_object)[["markers"]])

    if (!is.null(metadata(multimodal_object)$markers)) {
        names(metadata(multimodal_object)$markers) <- gsub("RNA", "gene", marker_names)
    }


    return(multimodal_object)
}

#' Clean Vector of Chevreul Names
#'
#' Cleans names of objects provided in a vector form
#'
#' @param myvec A vector of object names
#'
#' @return a clean vector of object names
#' @export
#' @examples
make_chevreul_clean_names <- function(myvec) {
    myvec %>%
        set_names(str_to_title(str_replace_all(., "[^[:alnum:][:space:]\\.]", " ")))
}


#' Convert seurat to seurat V5 format
#'
#' Convert seurat from v3 to v5 format
#'
#' @param seurat_v3 a version 3 seurat
#'
#' @return a verstion 5 seurat
#' @export
#'
#' @examples
#' convert_v3_to_v5(human_gene_transcript_seurat)
convert_v3_to_v5 <- function(seurat_v3) {
    seurat_version <- metadata(seurat_v3)$experiment$seurat_version

    if (seurat_version < 5 || is.null(seurat_version)) {
        meta <- seurat_v3@meta.data

        seurat_v5 <- CreateSingleCellExperimentObject(counts = seurat_v3$gene@counts, data = seurat_v3$gene@data, experiment = "gene", meta.data = meta)

        transcript_experiment.v5 <- CreateAssay5Object(counts = seurat_v3$transcript@counts, data = seurat_v3$transcript@data)
        seurat_v5$transcript <- transcript_experiment.v5

        seurat_v5$gene <- seurat_preprocess(seurat_v5$gene, normalize = FALSE)
        seurat_v5$transcript <- seurat_preprocess(seurat_v5$transcript, normalize = FALSE)

        # seurat_v5 <- clustering_workflow(seurat_v5)
        seurat_v5@reductions <- seurat_v3@reductions
        seurat_v5@graphs <- seurat_v3@graphs
        seurat_v5@neighbors <- seurat_v3@neighbors

        metadata(seurat_v5) <- metadata(seurat_v3)

        Idents(seurat_v5) <- Idents(seurat_v3)

        metadata(seurat_v5)$experiment$seurat_version <- packageVersion("seurat")
    } else {
        seurat_v5 <- seurat_v3
    }

    return(seurat_v5)
}

#' Get genes from Object
#'
#' Get genes from the object
#'
#' @param object a single cell object
#'
#' @return a vector of genes in a single cell object
#' @export
#' @examples
genes_from_object <- function(object) {
    rownames(object)
}

#' Get metadata from object
#'
#' Get metadata from the given object
#'
#' @param object a single cell object
#'
#' @return a tibble with metadata from a single cell object
#' @export
#' @examples
metadata_from_object <- function(object) {
    colnames(colData(object))
}

convert_seurat_to_sce <- function(seu) {
    sce <- as.SingleCellExperiment(seu, experiment = DefaultAssay(seu))

    alt_exp_names <- SingleCellExperiment::Assays(seu)[!SingleCellExperiment::Assays(seu) == DefaultAssay(seu)]

    for (i in alt_exp_names) {
        altExp(sce, i) <- as.SingleCellExperiment(seu, experiment = i)
    }

    sce@metadata <- seu@misc

    return(sce)
}

#' Save object to <project>/output/sce/<feature>_object.rds
#'
#' @param object a single cell object
#' @param prefix a prefix for saving
#' @param proj_dir path to a project directory
#'
#' @return a path to an rds file containing a single cell object
#' @export
#'
#' @examples
#' \dontrun{
#' save_object(gene = feature_objects$gene, transcript = feature_objects$transcript, proj_dir = proj_dir)
#'
#' save_object(gene = feature_objects$gene, transcript = feature_objects$transcript, prefix = "remove_nonPRs", proj_dir = proj_dir)
#' }
save_object <- function(object, prefix = "unfiltered", proj_dir = getwd()) {
  object_dir <- path(proj_dir, "output", "seurat")

  dir.create(object_dir)

  object_path <- path(object_dir, paste0(prefix, "_object.rds"))

  message(paste0("saving to ", object_path))
  saveRDS(object, object_path)

  return(object)
}
