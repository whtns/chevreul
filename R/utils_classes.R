setGeneric("list_plot_types", function(object)
  standardGeneric("list_plot_types"))

#' Collate list of variables to be plotted
#'
#' @param object
#'
#' @return plot_types a list of category_vars or continuous_vars
#' @export
#' @examples
setMethod("list_plot_types", "Seurat",
          function(object) {
  meta_types <- tibble::tibble(
    vars = colnames(pull_metadata(object)),
    var_type = purrr::map_chr(purrr::map(pull_metadata(object), pillar::new_pillar_type), 1),
    num_levels = unlist(purrr::map(pull_metadata(object), ~ length(unique(.x))))
  )

  meta_types <- meta_types %>%
    dplyr::filter(!grepl("_snn_res", vars)) %>%
    dplyr::mutate(meta_type = dplyr::case_when(
      var_type %in% c("int", "dbl") ~ "continuous",
      var_type %in% c("chr", "fct", "ord", "lgl") ~ "category"
    )) %>%
    dplyr::mutate(meta_type = ifelse(meta_type == "continuous" & num_levels < 30, "category", meta_type)) %>%
    dplyr::filter(num_levels > 1) %>%
    identity()

  continuous_vars <- meta_types %>%
    dplyr::filter(meta_type == "continuous") %>%
    dplyr::pull(vars)

  continuous_vars <- c("feature", continuous_vars) %>%
    purrr::set_names(stringr::str_to_title(stringr::str_replace_all(., "[[:punct:]]", " ")))


  category_vars <- meta_types %>%
    dplyr::filter(meta_type == "category") %>%
    dplyr::pull(vars) %>%
    purrr::set_names(stringr::str_to_title(stringr::str_replace_all(., "[^[:alnum:][:space:]\\.]", " ")))

  plot_types <- list(category_vars = category_vars, continuous_vars = continuous_vars)



  return(plot_types)
})

#' Collate list of variables to be plotted
#'
#' @param object
#'
#' @return plot_types a list of category_vars or continuous_vars
#' @export
#' @examples
setMethod("list_plot_types", "SingleCellExperiment",
          function(object) {
            meta_types <- tibble::tibble(
              vars = colnames(colData(object)),
              var_type = purrr::map_chr(purrr::map(colData(object), pillar::new_pillar_type), 1),
              num_levels = unlist(purrr::map(colData(object), ~ length(unique(.x))))
            )

            meta_types <- meta_types %>%
              dplyr::filter(!grepl("_snn_res", vars)) %>%
              dplyr::mutate(meta_type = dplyr::case_when(
                var_type %in% c("int", "dbl") ~ "continuous",
                var_type %in% c("chr", "fct", "ord", "lgl") ~ "category"
              )) %>%
              dplyr::mutate(meta_type = ifelse(meta_type == "continuous" & num_levels < 30, "category", meta_type)) %>%
              dplyr::filter(num_levels > 1) %>%
              identity()

            continuous_vars <- meta_types %>%
              dplyr::filter(meta_type == "continuous") %>%
              dplyr::pull(vars)

            continuous_vars <- c("feature", continuous_vars) %>%
              purrr::set_names(stringr::str_to_title(stringr::str_replace_all(., "[[:punct:]]", " ")))


            category_vars <- meta_types %>%
              dplyr::filter(meta_type == "category") %>%
              dplyr::pull(vars) %>%
              purrr::set_names(stringr::str_to_title(stringr::str_replace_all(., "[^[:alnum:][:space:]\\.]", " ")))

            plot_types <- list(category_vars = category_vars, continuous_vars = continuous_vars)



            return(plot_types)
          })

# pull_metadata ------------------------------

setGeneric("pull_metadata", function(object)
  standardGeneric("pull_metadata"))

#' Pull object metadata
#'
#' @param object
#'
#' @return object metadata
#' @export
#' @examples
setMethod("pull_metadata", "Seurat",
          function(object) {
            object[[]]
          })

#' Pull object metadata
#'
#' @param object
#'
#' @return object metadata
#' @export
#' @examples
setMethod("pull_metadata", "SingleCellExperiment",
          function(object) {
            colData(object) %>%
              as.data.frame()
          })

# set_metadata ------------------------------

setGeneric("set_metadata", function(object, meta.data)
  standardGeneric("set_metadata"))

#' Pull object metadata
#'
#' @param object
#'
#' @return object metadata
#' @export
#' @examples
setMethod("set_metadata", "Seurat",
          function(object, meta.data) {
            object@meta.data <- meta.data
            return(object)
          })

#' Pull object metadata
#'
#' @param object
#'
#' @return object metadata
#' @export
#' @examples
setMethod("set_metadata", "SingleCellExperiment",
          function(object, meta.data) {
            colData(object) <- meta.data
            return(object)
          })
