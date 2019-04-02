
#' Update Seurat Metadata
#'
#' @param seu_meta
#' @param merged_new_meta
#' @param numcols
#'
#' @return
#' @export
#'
#' @examples
update_seu_meta <- function(seu_meta, merged_new_meta, numcols) {
	# browser()
	common_cols <- intersect(colnames(seu_meta), colnames(merged_new_meta))
  seu_meta <- mutate_at(seu_meta, .vars = vars(one_of(numcols)), .funs = funs(as.numeric))
  updated_seu_meta <- dplyr::left_join(seu_meta, merged_new_meta, by = "Sample_ID")

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
	seu@meta.data <- as.data.frame(new_meta, row.names = new_meta$Sample_ID)
	return(seu)
}

