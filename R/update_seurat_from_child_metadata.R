
#' Update Seurat Metadata from Project
#'
#' Given a project metadata file located in <proj_dir>/data/<meta_file>, update an existing seurat object with the project data
#'
#' @param seu
#' @param proj_dir
#' @param numcols
#'
#' @return
#' @export
#'
#' @examples
update_seu_meta <- function(seu, proj_dir, numcols) {
	# browser()

  seu_meta <- as.data.frame(seu[[]])

  project_meta <- readr::read_csv(get_meta(proj_dir))

	common_cols <- intersect(colnames(seu_meta), colnames(project_meta))
  seu_meta <- mutate_at(seu_meta, .vars = vars(one_of(numcols)), .funs = funs(as.numeric))
  updated_seu_meta <- dplyr::left_join(seu_meta, project_meta, by = "sample_id")

  left_side_common <- paste0(common_cols, ".y")
  right_side_common <- paste0(common_cols, ".x")

  left_side_columns <- dplyr::select(updated_seu_meta, one_of(left_side_common))
  colnames(left_side_columns) <- gsub("\\.y$", "", colnames(left_side_columns))

  updated_seu_meta <- cbind(updated_seu_meta, left_side_columns) %>%
  	dplyr::select(-one_of(right_side_common)) %>%
  	dplyr::select(-one_of(left_side_common)) %>%
  	dplyr::select(Sample_ID, everything())
}



#' Reset Seurat Metadata
#'
#' @param seu
#' @param new_meta
#'
#' @return
#' @export
#'
#' @examples
reset_seu_meta <- function(seu, new_meta){
	# browser()
	seu@meta.data <- as.data.frame(new_meta, row.names = new_meta$sample_id)
	return(seu)
}

