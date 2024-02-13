# get_feature_types ------------------------------


# list_plot_types ------------------------------

#' Collate list of variables to be plotted
#'
#' @param object
#'
#' @return plot_types a list of category_vars or continuous_vars
#' @export
#' @examples
setGeneric("list_plot_types", function(object) {
    standardGeneric("list_plot_types")
})

setMethod(
    "list_plot_types", "Seurat",
    function(object) {
        meta_types <- tibble::tibble(
            vars = colnames(pull_metadata(object)),
            var_type = purrr::map_chr(map(pull_metadata(object), pillar::new_pillar_type), 1),
            num_levels = unlist(map(pull_metadata(object), ~ length(unique(.x))))
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
)

setMethod(
    "list_plot_types", "SingleCellExperiment",
    function(object) {
        meta_types <- tibble::tibble(
            vars = colnames(colData(object)),
            var_type = purrr::map_chr(map(colData(object), pillar::new_pillar_type), 1),
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
)

# pull_metadata ------------------------------

#' Pull object metadata
#'
#' @param object
#'
#' @return object metadata
#' @export
#' @examples
setGeneric("pull_metadata", function(object) {
    standardGeneric("pull_metadata")
})

setMethod(
    "pull_metadata", "Seurat",
    function(object) {
        object[[]]
    }
)

setMethod(
    "pull_metadata", "SingleCellExperiment",
    function(object) {
        colData(object) %>%
            as.data.frame()
    }
)

# get variable features ------------------------------

#' Get variable features
#'
#' @param object
#'
#' @return variable features from a single cell object
#' @export
#' @examples
setGeneric("get_variable_features", function(object, ...) {
    standardGeneric("get_variable_features")
})

setMethod(
    "get_variable_features", "Seurat",
    function(object, ...) {
        VariableFeatures(object, ...)
    }
)

setMethod(
    "get_variable_features", "SingleCellExperiment",
    function(object, ...) {
        scran::getTopHVGs(object)
    }
)

# get feature names ------------------------------

#' Get feature names
#'
#' @param object
#'
#' @return variable features from a single cell object
#' @export
#' @examples
setGeneric("get_features", function(object, ...) {
    standardGeneric("get_features")
})

setMethod(
    "get_features", "Seurat",
    function(object, ...) {
        Features(object, ...)
    }
)

setMethod(
    "get_features", "SingleCellExperiment",
    function(object, ...) {
        rownames(object)
    }
)


# set_metadata ------------------------------

#' Pull object metadata
#'
#' @param object
#'
#' @return object metadata
#' @export
#' @examples
setGeneric("set_metadata", function(object, meta.data) {
    standardGeneric("set_metadata")
})

setMethod(
    "set_metadata", "Seurat",
    function(object, meta.data) {
        object@meta.data <- meta.data
        return(object)
    }
)

setMethod(
    "set_metadata", "SingleCellExperiment",
    function(object, meta.data) {
        colData(object) <- meta.data
        return(object)
    }
)

# get_feature_types ------------------------------

setGeneric("get_feature_types", function(object) {
    standardGeneric("get_feature_types")
})

setMethod("get_feature_types", "Seurat", function(object) {
    sort(names(object@assays))
})

setMethod("get_feature_types", "SingleCellExperiment", function(object) {
    sort(c(mainExpName(object), altExpNames(object)))
})

setGeneric("set_feature_types", function(object, feature_type) {
  standardGeneric("set_feature_types")
})

setMethod("set_feature_types", "Seurat", function(object, feature_type) {
  DefaultAssay(object) <- feature_type
})

setMethod("set_feature_types", "SingleCellExperiment", function(object, feature_type) {
  if(feature_type %in% altExpNames(object)){
    object <- swapAltExp(object, feature_type, saved = mainExpName(object), withColData = TRUE)
  }
  return(object)
})

# check_integrated ------------------------------

setGeneric("check_integrated", function(object) {
    standardGeneric("check_integrated")
})

setMethod("check_integrated", "Seurat", function(object) {
    if ("integrated" %in% names(object@assays)) {
        assay <- "integrated"
    } else {
        assay <- "gene"
    }

    return(assay)
})

setMethod("check_integrated", "SingleCellExperiment", function(object) {
    assay <- str_subset(assayNames(scMerge_unsupervised), "scMerge")

    return(assay)
})

# retrieve_assays ------------------------------
setGeneric("retrieve_assay", function(object, assay) {
    standardGeneric("retrieve_assay")
})

setMethod("retrieve_assay", "Seurat", function(object, assay) {
    assay <- GetAssay(object, assay = assay)

    return(assay)
})

setMethod("retrieve_assay", "SingleCellExperiment", function(object, assay) {
    if (assay %in% mainExpName(object)) {
        return(object)
    } else if (assay %in% altExpNames(object)) {
        return(altExp(object, assay))
    }
})

# query_assay ------------------------------
setGeneric("query_assay", function(object, assay) {
    standardGeneric("query_assay")
})

setMethod("query_assay", "Seurat", function(object, assay) {
    return(assay %in% Seurat::Assays(object))
})

setMethod("query_assay", "SingleCellExperiment", function(object, assay) {
    return(assay %in% c(mainExpName(object), altExpNames(object)))
})
