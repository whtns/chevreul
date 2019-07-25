#' Collate list of variables to be plotted
#'
#' @param seu
#'
#' @return plot_types
#'
#' @examples
list_plot_types <- function(seu){
  meta_types <- unlist(lapply(seu[[]], class))

  meta_types <- meta_types[!grepl(paste0("^", DefaultAssay(seu), "_snn_res"), names(meta_types))]

  meta_types <- split(meta_types, factor(meta_types))
  meta_types <- lapply(meta_types, names)


  continuous_vars <- c("custom", unlist(meta_types[c("integer", "numeric")]))

  continuous_names <- stringr::str_to_title(gsub("[[:punct:]]", " ", continuous_vars))

  names(continuous_vars) <- continuous_names


  category_vars <- c("seurat", unlist(meta_types[c("character", "factor")]))
  category_vars <- category_vars[!grepl("snn_res", category_vars)]

  category_names <- stringr::str_to_title(gsub("[[:punct:]]", " ", category_vars))

  names(category_vars) <- category_names

  plot_types <- list(category_vars = category_vars, continuous_vars = continuous_vars)

  return(plot_types)


}
