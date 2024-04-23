#' Get Transcripts in object
#'
#' Get transcript ids in objects for one or more gene of interest
#'
#' @param object A SingleCellExperiment object
#' @param gene Gene of interest
#' @param organism Organism
#'
#' @return transcripts constituting a gene of interest in a SingleCellExperiment object
#' @export
#'
#' @examples
#' chevreul_sce <- chevreuldata::human_gene_transcript_sce()
#'
#' NRL_transcripts <- get_transcripts_from_object(chevreul_sce, "NRL")
#'
get_transcripts_from_object <- function(object, gene, organism = "human") {
    transcripts <- genes_to_transcripts(gene, organism)

    transcripts <- transcripts[transcripts %in% get_features(object, "transcript")]
}


#' Record Experiment Metadata
#'
#' Records miscellaneous data
#' @param object A object
#' @param experiment_name name of the experiment
#' @param organism human or mouse
#'
#' @return a SingleCellExperiment object
#' @export
#' @examples
#' chevreul_sce <- chevreuldata::human_gene_transcript_sce()
#' record_experiment_data(chevreul_sce)
#'
record_experiment_data <- function(object, experiment_name = "default_experiment", organism = "human") {
    if (!requireNamespace("SingleCellExperiment", quietly = TRUE)) {
        stop("Package 'object' needed for this function to work. Please install it.",
            call. = FALSE
        )
    }

    organism <- metadata(object)[["experiment"]][["organism"]] %||% organism

    experiment_name <- metadata(object)[["experiment"]][["experiment_name"]] %||% experiment_name

    message(glue("[{format(Sys.time(), '%H:%M:%S')} Logging Technical Details..."))
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
#' @param object A SingleCellExperiment object
#'
#' @return a SingleCellExperiment object with nfeatures and ngenes stored in metadata
#' @export
#' @examples
#' chevreul_sce <- chevreuldata::human_gene_transcript_sce()
#' object_calcn(chevreul_sce)
object_calcn <- function(object) {
    object <- addPerCellQC(object)
    object[[glue("nFeature_{mainExpName(object)}")]] <- object$detected
    object[[glue("nCount_{mainExpName(object)}")]] <- object$sum

    for (alt_exp_name in altExpNames(object)[!altExpNames(object) == "velocity"]) {
        altExp(object, alt_exp_name) <- addPerCellQC(altExp(object, alt_exp_name))
        altExp(object, alt_exp_name)[[glue("nFeature_{alt_exp_name}")]] <- altExp(object, alt_exp_name)[["detected"]]
        altExp(object, alt_exp_name)[[glue("nCount_{alt_exp_name}")]] <- altExp(object, alt_exp_name)[["sum"]]
    }

    return(object)
}


#' Propagate Metadata Changes
#'
#' @param meta updated metadata
#' @param object a SingleCellExperiment object
#'
#' @return a SingleCellExperiment object
#' @export
#' @examples
#'
#' chevreul_sce <- chevreuldata::human_gene_transcript_sce()
#' new_meta <- data.frame(row.names = colnames(chevreul_sce))
#' new_meta$example <- "example"
#'
#' propagate_spreadsheet_changes(new_meta, chevreul_sce)
propagate_spreadsheet_changes <- function(meta, object) {
    meta <- meta %>%
        tibble::rownames_to_column("cell") %>%
        mutate(across(contains("snn"), as.factor)) %>%
        mutate(across(where(is.ordered), ~ as.factor(as.character(.x)))) %>%
        tibble::column_to_rownames("cell") %>%
        identity()

    colData(object) <- DataFrame(meta)

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
#' @return a sqlite database with SingleCellExperiment objects
#' @export
#' @examples
#' \donttest{
#' create_project_db()
#' }
create_project_db <- function(cache_location = "~/.cache/chevreul", sqlite_db = "single-cell-projects.db", verbose = TRUE) {
    if (!dir.exists(cache_location)) {
        dir.create(cache_location)
    }
    con <- dbConnect(SQLite(), path(cache_location, sqlite_db))
    projects_tbl <- tibble(project_name = character(), project_path = character(), project_slug = character(), project_type = character(), )
    message(glue("building table of chevreul projects at {path(cache_location, sqlite_db)}"))
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
#' @return a sqlite database with SingleCellExperiment objects
#' @export
#' @examples
#' \donttest{
#' update_project_db()
#' }
update_project_db <- function(
        projects_dir = NULL,
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
#' @return a sqlite database with SingleCellExperiment objects

append_to_project_db <- function(
        new_project_path,
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
#' @return a tibble with SingleCellExperiment objects
#

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
#' @return a sqlite database of bigwig files for cells in a SingleCellExperiment object
#' @export
#' @examples
#' \donttest{
#' make_bigwig_db("example_project")
#' }
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
#' @return a tibble with cell level metadata from a SingleCellExperiment object
#' @examples
#' \donttest{
#' chevreul_sce <- chevreuldata::human_gene_transcript_sce()
#' metadata_from_batch(chevreul_sce)
#' }
#'
metadata_from_batch <- function(
        batch, projects_dir = "/dataVolume/storage/single_cell_projects",
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
#' chevreul_sce <- chevreuldata::human_gene_transcript_sce()
#' make_chevreul_clean_names(colnames(get_cell_metadata(chevreul_sce)))
make_chevreul_clean_names <- function(myvec) {
    myvec %>%
        set_names(str_to_title(str_replace_all(., "[^[:alnum:][:space:]\\.]", " ")))
}

#' Get metadata from object
#'
#' Get metadata from the given object
#'
#' @param object a SingleCellExperiment object
#'
#' @return a tibble with metadata from a SingleCellExperiment object
#' @export
#' @examples
#' chevreul_sce <- chevreuldata::human_gene_transcript_sce()
#' metadata_from_object(chevreul_sce)
#'
metadata_from_object <- function(object) {
    colnames(colData(object))
}

#' Save object to <project>/output/sce/<feature>_object.rds
#'
#' @param object a SingleCellExperiment object
#' @param prefix a prefix for saving
#' @param proj_dir path to a project directory
#'
#' @return a path to an rds file containing a SingleCellExperiment object
#'
#'
#'
save_object <- function(object, prefix = "unfiltered", proj_dir = getwd()) {
    object_dir <- path(proj_dir, "output", "singlecellexperiment")

    dir.create(object_dir)

    object_path <- path(object_dir, paste0(prefix, "_sce.rds"))

    message(glue("saving to {object_path}"))
    saveRDS(object, object_path)

    return(object)
}
