#' Update Seurat Metadata from Project
#'
#' Given a project metadata file located in <proj_dir>/data/<meta_file>, update an existing object with the project data
#'
#' @param object A object
#' @param proj_dir Project directory
#' @param numcols number of columns
#'
#' @return updated object metadata
#' @export
update_object_meta <- function(object, proj_dir, numcols) {
    object_meta <- as.data.frame(get_cell_metadata(object))

    project_meta <- readr::read_csv(get_meta(proj_dir))

    common_cols <- intersect(colnames(object_meta), colnames(project_meta))
    object_meta <- mutate_at(object_meta, .vars = vars(one_of(numcols)), .funs = funs(as.numeric))
    updated_object_meta <- left_join(object_meta, project_meta, by = "sample_id")

    left_side_common <- paste0(common_cols, ".y")
    right_side_common <- paste0(common_cols, ".x")

    left_side_columns <- select(updated_object_meta, one_of(left_side_common))
    colnames(left_side_columns) <- gsub("\\.y$", "", colnames(left_side_columns))

    updated_object_meta <- cbind(updated_object_meta, left_side_columns) %>%
        select(-one_of(right_side_common)) %>%
        select(-one_of(left_side_common)) %>%
        select(Sample_ID, everything())
}


#' Reset Seurat Metadata
#'
#' @param object A object
#' @param new_meta new object
#'
#' @return a single cell object with new metadata
#' @export
reset_object_meta <- function(object, new_meta) {
    object@meta.data <- as.data.frame(new_meta, row.names = new_meta$sample_id)
    return(object)
}
