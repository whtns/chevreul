#' Rename cell ids from annoying old notation
#'
#' @param cell_ids
#' @param batch_id
#'
#' @return
#'
#' @examples
rename_from_x_notation <- function(cell_ids, batch_id) {
    cell_ids <- stringr::str_replace(cell_ids, "X", "")
    cell_ids <- paste0(batch_id, stringr::str_pad(cell_ids, width = max(nchar(cell_ids)), pad = "0"))
}

#' Reorganize objects to a multimodal format
#'
#' @param proj_dir
#'
#' @return
#' @export
#'
#' @examples
reorg_object_files <- function(projects_db = "~/.cache/chevreul/single-cell-projects.db") {
    conn <- DBI::dbConnect(RSQLite::SQLite(), projects_db)

    project_paths <- tbl(conn, "projects_tbl") %>%
        select(project_name, project_path) %>%
        collect() %>%
        tibble::deframe() %>%
        identity()

    DBI::dbDisconnect(conn)

    object_dirs <- map(project_paths, fs::path, "output", "seurat")

    get_object_objects <- function(object_dir) {
        rds_files <- fs::dir_ls(object_dir) %>%
            fs::path_filter("*.rds")
    }

    object_objects <- map(object_dirs, possibly(get_object_objects, NA))

    object_objects <- object_objects[!is.na(object_objects)]

    object_objects <- object_objects[lapply(object_objects, length) > 0]

    object_objects <- object_objects[!names(object_objects) %in% c("20170407-SHL-FACS-Hs_proj", "20171031-SHL-FACS-Hs_proj")]

    # # objects before conversion
    # object_objects <- tibble::enframe(object_objects) %>%
    #   tidyr::unnest() %>%
    #   dplyr::filter(!stringr::str_detect(value, "_multimodal.rds")) %>%
    #   dplyr::mutate(multi_copy = stringr::str_replace(value, ".rds", "_multimodal.rds")) %>%
    #   dplyr::mutate(exists = fs::file_exists(multi_copy)) %>%
    #   dplyr::filter(!exists) %>%
    #   dplyr::filter(!str_detect(value, "sce")) %>%
    #   dplyr::mutate(readable = fs::file_access(value, "read")) %>%
    #   dplyr::mutate(writable = fs::file_access(value, "write")) %>%
    #   dplyr::select(name, value) %>%
    #   tibble::deframe() %>%
    #   identity()

    # objects after conversion
    new_object_objects <- tibble::enframe(object_objects) %>%
        tidyr::unnest() %>%
        dplyr::filter(!stringr::str_detect(value, "_multimodal.rds")) %>%
        dplyr::mutate(multi_copy = stringr::str_replace(value, ".rds", "_multimodal.rds")) %>%
        dplyr::mutate(exists = fs::file_exists(multi_copy)) %>%
        dplyr::filter(exists) %>%
        dplyr::mutate(readable = fs::file_access(multi_copy, "read")) %>%
        dplyr::mutate(writable = fs::file_access(value, "write")) %>%
        # dplyr::filter(!str_detect(value, "sce")) %>%
        # dplyr::mutate(readable = fs::file_access(value, "read")) %>%
        # dplyr::mutate(writable = fs::file_access(value, "write")) %>%
        # dplyr::select(name, value) %>%
        # tibble::deframe() %>%
        identity()

    safe_update <- purrr::safely(update_chevreul_object)

    map(object_objects, safe_update, return_object = FALSE)

    # message(paste0("deleting old files: ", old_files))
    # fs::file_delete(old_files)

    # return(rds_files)
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
prep_armor_meta <- function(proj_dir) {
    object <- load_object_from_proj(proj_dir)
    meta <- pull_metadata(object) %>%
        dplyr::mutate(names = sample_id) %>%
        dplyr::mutate(type = "PE")

    readr::write_tsv(meta, fs::path(proj_dir, "data", "metadata.txt"))
}

#' Profile the memory of an s4object
#'
#' @param s4object
#'
#' @return
#' @export
#'
#' @examples
mem_str <- function(s4object) {
    slot_names <- slotNames(s4object) %>%
        set_names(.)

    map(slot_names, ~ slot(s4object, .x)) %>%
        map(pryr::object_size) %>%
        tibble::enframe("slot", "hr_size") %>%
        dplyr::mutate(hr_size = fs::as_fs_bytes(hr_size)) %>%
        dplyr::arrange(desc(hr_size)) %>%
        identity()
}
