#' Replace object Metadata
#'
#' Replace object metadata.
#'
#' @param object A object
#' @param datapath Path to csv file containing metadata
#'
#' @return a single cell object
#' @export
#' @examples
#' new_meta <- system.file("extdata", "new_meta.csv", package="chevreul")
#' replace_object_metadata(human_gene_transcript_sce, new_meta)
replace_object_metadata <- function(object, datapath) {
        new_meta <- read_csv(datapath) %>%
            mutate(across(contains("snn"), as.factor))

        new_meta <- DataFrame(column_to_rownames(new_meta, colnames(new_meta)[1]))

        colData(object) <- new_meta

        return(object)
    }

#' Reformat object Metadata
#'
#' Reformat object Metadata by Coalesced Columns
#'
#' @param object A SingleCellExperiment object
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

        return(meta)
    }

#' Get Transcripts in object
#'
#' Get transcript ids in objects for one or more gene of interest
#'
#' @param object A SingleCellExperiment object
#' @param gene Gene of interest
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

#' Calculate Read Count Metrics for a object
#'
#' Recalculate counts/features per cell for a object
#'
#' @param object A single cell object
#'
#' @return a single cell object with nfeatures and ngenes stored in metadata
#' @export
#' @examples
#' object_calcn(human_gene_transcript_sce)
object_calcn <- function(object) {

  object <- addPerCellQC(object)
  object[[glue("nFeature_{mainExpName(object)}")]] <- object$detected
  object[[glue("nCount_{mainExpName(object)}")]] <- object$sum

  for (alt_exp_name in altExpNames(object)){
    altExp(object, alt_exp_name) <- addPerCellQC(altExp(object, alt_exp_name))
    altExp(object, alt_exp_name)[[glue("nFeature_{alt_exp_name}")]] <- altExp(object, alt_exp_name)[["detected"]]
    altExp(object, alt_exp_name)[[glue("nCount_{alt_exp_name}")]] <- altExp(object, alt_exp_name)[["sum"]]
  }

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

    colData(object) <- meta

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
#' \donttest{create_project_db()}
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
#' \donttest{update_project_db()}
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
#' Append projects to database
#'
#' @param new_project_path new project path
#' @param cache_location Path to cache "~/.cache/chevreul"
#' @param sqlite_db sqlite db
#' @param verbose print messages
#'
#'
#' @return a sqlite database with single cell objects
#' @export
#' @examples
#' \donttest{append_to_project_db("example_project_path")}
append_to_project_db <- function(new_project_path,
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
#' @param cache_location Path to cache "~/.cache/chevreul"
#' @param sqlite_db sqlite db
#' @param verbose print messages
#'
#' @return a tibble with single cell objects
#' @export
#' @examples
#' \donttest{read_project_db()}
read_project_db <- function(
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
#' \donttest{make_bigwig_db("example_project")}
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

#' Clean Vector of Chevreul Names
#'
#' Cleans names of objects provided in a vector form
#'
#' @param myvec A vector of object names
#'
#' @return a clean vector of object names
#' @export
#' @examples
#' make_chevreul_clean_names(colnames(get_cell_metadata(human_gene_transcript_sce)))
make_chevreul_clean_names <- function(myvec) {
    myvec %>%
        set_names(str_to_title(str_replace_all(., "[^[:alnum:][:space:]\\.]", " ")))
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
#' metadata_from_object(human_gene_transcript_sce)
#'
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
