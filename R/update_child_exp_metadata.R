#' Get Metadata for a project
#'
#' @param proj_path a project path
#'
#' @return a filepath to a project metadata file
#' @examples
get_meta <- function(proj_path) {
    meta_path <- path(proj_path, "data", gsub("_proj", "_metadata.csv", path_file(proj_path)))
}


#' Update experiment Metadata
#'
#' @param original_meta original metadata
#' @param corrected_meta corrected metadata
#'
#' @return a cell-level metadata tibble
#' @export
#' @examples
update_exp_meta <- function(original_meta, corrected_meta) {
    coltypes <- sapply(corrected_meta, is.numeric)
    numcols <- colnames(corrected_meta)[coltypes]

    common_cols <- intersect(colnames(original_meta), colnames(corrected_meta))
    original_meta <- mutate_at(original_meta, .vars = vars(one_of(numcols)), .funs = funs(as.numeric))
    updated_meta <- left_join(original_meta, corrected_meta, by = "sample_id")

    left_side_common <- paste0(common_cols, ".x")
    right_side_common <- paste0(common_cols, ".y")

    right_side_columns <- select(updated_meta, one_of(right_side_common))
    colnames(right_side_columns) <- gsub("\\.y$", "", colnames(right_side_columns))

    updated_meta <- cbind(updated_meta, right_side_columns) %>%
        select(-one_of(left_side_common)) %>%
        select(-one_of(right_side_common)) %>%
        select(sample_id, everything()) %>%
        # select(common_cols) %>%
        identity()
}
