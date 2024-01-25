# get_feature_types ------------------------------


# list_plot_types ------------------------------

#' Collate list of variables to be plotted
#'
#' @param object
#'
#' @return plot_types a list of category_vars or continuous_vars
#' @export
#' @examples
setGeneric("list_plot_types", function(object)
  standardGeneric("list_plot_types"))

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

#' Pull object metadata
#'
#' @param object
#'
#' @return object metadata
#' @export
#' @examples
setGeneric("pull_metadata", function(object)
  standardGeneric("pull_metadata"))

setMethod("pull_metadata", "Seurat",
          function(object) {
            object[[]]
          })

setMethod("pull_metadata", "SingleCellExperiment",
          function(object) {
            colData(object) %>%
              as.data.frame()
          })

# get variable features ------------------------------

#' Pull object metadata
#'
#' @param object
#'
#' @return variable features from a single cell object
#' @export
#' @examples
setGeneric("get_variable_features", function(object, ...)
  standardGeneric("get_variable_features"))

setMethod("get_variable_features", "Seurat",
          function(object, ...) {
            VariableFeatures(object, ...)
          })

setMethod("get_variable_features", "SingleCellExperiment",
          function(object, ...) {
            scran::getTopHVGs(object, ...)
          })


# set_metadata ------------------------------

#' Pull object metadata
#'
#' @param object
#'
#' @return object metadata
#' @export
#' @examples
setGeneric("set_metadata", function(object, meta.data)
  standardGeneric("set_metadata"))

setMethod("set_metadata", "Seurat",
          function(object, meta.data) {
            object@meta.data <- meta.data
            return(object)
          })

setMethod("set_metadata", "SingleCellExperiment",
          function(object, meta.data) {
            colData(object) <- meta.data
            return(object)
          })

setGeneric("get_feature_types", function(object) {
  standardGeneric("get_feature_types")
})

setMethod("get_feature_types", "Seurat", function(object) {
  names(object()@assays)
})

setMethod("get_feature_types", "SingleCellExperiment", function(object) {
  c(mainExpName(human_gene_transcript_sce), altExpNames(human_gene_transcript_sce))
})
