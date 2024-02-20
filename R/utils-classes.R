#' Collate list of variables to be plotted
#'
#' @param object a single cell object
#'
#' @return plot_types a list of category_vars or continuous_vars
#' @export
#' @importFrom pillar new_pillar_type
list_plot_types <- function(object) {
        meta_types <- tibble::tibble(
            vars = colnames(colData(object)),
            var_type = purrr::map_chr(map(colData(object), new_pillar_type), 1),
            num_levels = unlist(map(colData(object), ~ length(unique(.x))))
        )

        meta_types <- meta_types %>%
            filter(!grepl("_snn_res", vars)) %>%
            mutate(meta_type = case_when(
                var_type %in% c("int", "dbl") ~ "continuous",
                var_type %in% c("chr", "fct", "ord", "lgl") ~ "category"
            )) %>%
            mutate(meta_type = ifelse(meta_type == "continuous" & num_levels < 30, "category", meta_type)) %>%
            filter(num_levels > 1) %>%
            identity()

        continuous_vars <- meta_types %>%
            filter(meta_type == "continuous") %>%
            pull(vars)

        continuous_vars <- c("feature", continuous_vars) %>%
            purrr::set_names(stringr::str_to_title(stringr::str_replace_all(., "[[:punct:]]", " ")))


        category_vars <- meta_types %>%
            filter(meta_type == "category") %>%
            pull(vars) %>%
            purrr::set_names(stringr::str_to_title(stringr::str_replace_all(., "[^[:alnum:][:space:]\\.]", " ")))

        plot_types <- list(category_vars = category_vars, continuous_vars = continuous_vars)



        return(plot_types)
    }

# Pull metadata

#' Pull the metadata from a given object
#'
#' @param object a single cell object
#'
#' @return dataframe containing object metadata
get_cell_metadata <- function(object) {
        colData(object) %>%
            as.data.frame()
    }

#' Get object metadata
#'
#' @param object a single cell object
#' @param ... extra args passed to seurat class method
#'
#' @return variable features from a single cell object
#' @export
get_object_metadata <- function(object, ...) {
    metadata(object)
  }

#' Get variable features
#'
#' @param object a single cell object
#' @param ... extra args passed to seurat class method
#'
#' @return variable features from a single cell object
#' @export
get_variable_features <- function(object, ...) {
        scran::getTopHVGs(object)
    }

# get feature names ------------------------------

#' Get feature names
#'
#' @param object a single cell object
#' @param ... extra args passed to seurat class method
#'
#' @return variable features from a single cell object
#' @export
get_features <- function(object, ...) {
        rownames(object)
    }


#' Add new metadata to the object metadata
#'
#' @param object a single cell object
#' @param meta.data new metadata
#'
#' @return object metadata
#' @export
set_metadata <- function(object, meta.data) {
        colData(object) <- meta.data
        return(object)
    }

# get_feature_types ------------------------------

get_feature_types <- function(object) {
    sort(c(mainExpName(object), altExpNames(object)))
}
set_feature_types <- function(object, feature_type) {
  if(feature_type %in% altExpNames(object)){
    object <- swapAltExp(object, feature_type, saved = mainExpName(object), withColData = TRUE)
  }
  return(object)
}

# check_integrated ------------------------------

check_integrated <- function(object) {
    assay <- str_subset(assayNames(scMerge_unsupervised), "scMerge")

    return(assay)
}

# retrieve_assays ------------------------------
retrieve_assay <- function(object, assay) {
    if (assay %in% mainExpName(object)) {
        return(object)
    } else if (assay %in% altExpNames(object)) {
        return(altExp(object, assay))
    }
}

# query_assay ------------------------------
query_assay <- function(object, assay) {
    return(assay %in% c(mainExpName(object), altExpNames(object)))
}
