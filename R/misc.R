#' Rename cell ids from annoying old notation
#'
#' @param cell_ids cell ids
#' @param batch_id bath id
#'
#' @return a vector of cell ids
rename_from_x_notation <- function(cell_ids, batch_id) {
    cell_ids <- str_replace(cell_ids, "X", "")
    cell_ids <- paste0(batch_id, stringr::str_pad(cell_ids, width = max(nchar(cell_ids)), pad = "0"))
}

#' Reorganize objects to a multimodal format
#'
#' @param proj_dir project directory
#'
#' @return a list of single cell objects
#' @export
reorg_object_files <- function(projects_db = "~/.cache/chevreul/single-cell-projects.db") {
    conn <- dbConnect(RSQLite::SQLite(), projects_db)

    project_paths <- tbl(conn, "projects_tbl") %>%
        select(project_name, project_path) %>%
        collect() %>%
        tibble::deframe() %>%
        identity()

    dbDisconnect(conn)

    object_dirs <- map(project_paths, path, "output", "seurat")

    get_objects <- function(object_dir) {
        rds_files <- dir_ls(object_dir) %>%
            path_filter("*.rds")
    }

    objects <- map(object_dirs, possibly(get_objects, NA))

    objects <- objects[!is.na(objects)]

    objects <- objects[lapply(objects, length) > 0]

    objects <- objects[!names(objects) %in% c("20170407-SHL-FACS-Hs_proj", "20171031-SHL-FACS-Hs_proj")]

    # objects after conversion
    new_objects <- enframe(objects) %>%
        tidyr::unnest() %>%
        filter(!stringr::str_detect(value, "_multimodal.rds")) %>%
        mutate(multi_copy = str_replace(value, ".rds", "_multimodal.rds")) %>%
        mutate(exists = fs::file_exists(multi_copy)) %>%
        filter(exists) %>%
        mutate(readable = fs::file_access(multi_copy, "read")) %>%
        mutate(writable = fs::file_access(value, "write")) %>%
        identity()

    safe_update <- purrr::safely(update_chevreul_object)

    map(objects, safe_update, return_object = FALSE)
}

#' pad sample numbers to prep for armor
#'
#' @param proj_dir project directory
#'
#' @return a list of fastq files
pad_sample_files <- function(proj_dir) {
    fastq_files <-
        dir_ls(path(proj_dir, "data", "FASTQ"), glob = "*.fastq.gz") %>%
        purrr::set_names(path_file(.)) %>%
        enframe("file", "path") %>%
        tidyr::separate(file, into = c("proj_id", "sample_number", "pair"), sep = "-|_") %>%
        mutate(sample_number = stringr::str_pad(sample_number, max(nchar(sample_number)), pad = "0")) %>%
        mutate(file = path(proj_dir, "data", "FASTQ", paste0(proj_id, "-", sample_number, "_", pair))) %>%
        # mutate() %>%
        identity()

    purrr::map2(fastq_files$path, fastq_files$file, file_move)
}

#' prep armor metadata file
#'
#' @param proj_dir project directory
#'
#' @return a path to a tsv
prep_armor_meta <- function(proj_dir) {
    object <- load_object_from_proj(proj_dir)
    meta <- pull_metadata(object) %>%
        mutate(names = sample_id) %>%
        mutate(type = "PE")

    readr::write_tsv(meta, path(proj_dir, "data", "metadata.txt"))
}

