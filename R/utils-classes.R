#' Collate list of variables to be plotted
#'
#' @param object a SingleCellExperiment object
#'
#' @return plot_types a list of category_vars or continuous_vars
#' @export
#' @examples
#' 
#' 
#' chevreul_sce <- mockSCE(ncells=200, ngenes=1000)
#' list_plot_types(chevreul_sce)
list_plot_types <- function(object) {
    meta_types <- tibble(
        vars = colnames(colData(object)),
        var_type = map_chr(map(colData(object), new_pillar_type), 1),
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
        set_names(str_to_title(str_replace_all(., "[[:punct:]]", " ")))


    category_vars <- meta_types %>%
        filter(meta_type == "category") %>%
        pull(vars) %>%
        set_names(str_to_title(str_replace_all(., "[^[:alnum:][:space:]\\.]", " ")))

    plot_types <- list(category_vars = category_vars, continuous_vars = continuous_vars)



    return(plot_types)
}

# Get cell metadata

#' Get cell metadata from a given object
#'
#' @param object a SingleCellExperiment object
#'
#' @return dataframe containing object metadata
#' @export
#' @examples
#' 
#' 
#' chevreul_sce <- mockSCE(ncells=200, ngenes=1000)
#' get_cell_metadata(chevreul_sce)
get_cell_metadata <- function(object) {
    colData(object) %>%
        as.data.frame()
}

# set cell metadata

#' Set cell metadata from a given object
#'
#' @param object a SingleCellExperiment object
#' @param meta a dataframe containing object metadata
#'
#' @return a SingleCellExperiment object with new colData
#' @export
#' @examples
#' 
#' 
#' chevreul_sce <- mockSCE(ncells=200, ngenes=1000)
#' new_meta <- data.frame(row.names = colnames(chevreul_sce))
#' new_meta$example <- "example"
#' set_cell_metadata(chevreul_sce, new_meta)
set_cell_metadata <- function(object, meta) {
    colData(object) <- DataFrame(meta)
    return(object)
}

#' Get object metadata
#'
#' @param object a SingleCellExperiment object
#'
#' @return variable features from a SingleCellExperiment object
#' @export
#' @importFrom S4Vectors metadata
#' @examples
#' 
#' 
#' chevreul_sce <- mockSCE(ncells=200, ngenes=1000)
#' get_object_metadata(chevreul_sce)
get_object_metadata <- function(object) {
    metadata(object)
}

#' Get variable features
#'
#' @param object a SingleCellExperiment object
#' @param experiment "gene" or "transcript"
#'
#' @return variable features from a SingleCellExperiment object
#' @export
#' @examples
#' 
#' chevreul_sce <- chevreuldata::human_gene_transcript_sce()
#' get_variable_features(chevreul_sce)
get_variable_features <- function(object, experiment = "gene") {
    if (experiment == mainExpName(object)) {
        getTopHVGs(object)
    } else if (experiment %in% altExpNames(object)) {
        getTopHVGs(altExp(object, experiment))
    }
}

# get feature names ------------------------------

#' Get feature names
#'
#' @param object a SingleCellExperiment object
#' @param experiment "gene" or "transcript"
#'
#' @return variable features from a SingleCellExperiment object
#' @export
#' @examples
#' 
#' 
#' chevreul_sce <- chevreuldata::human_gene_transcript_sce()
#' get_features(chevreul_sce)
get_features <- function(object, experiment = "gene") {
    if (experiment == mainExpName(object)) {
        rownames(object)
    } else if (experiment %in% altExpNames(object)) {
        rownames(altExp(object, experiment))
    }
}


#' Get Feature Types
#'
#' @param object a SingleCellExperiment object
#'
#' @return vector of feature types in an object
#' @export
#'
#' @examples
#' 
#' 
#' chevreul_sce <- mockSCE(ncells=200, ngenes=1000)
#' get_feature_types(chevreul_sce)
get_feature_types <- function(object) {
    sort(c(mainExpName(object), altExpNames(object)))
}

#' Set Feature Types
#'
#' @param object a SingleCellExperiment object
#' @param feature_type feature type
#' @return an object with assigned feature type
#' @export
#'
#' @examples
#' 
#' 
#' chevreul_sce <- mockSCE(ncells=200, ngenes=1000)
#' set_feature_type(chevreul_sce, "transcript")
set_feature_type <- function(object, feature_type) {
    if (feature_type %in% altExpNames(object)) {
        object <- swapAltExp(object, feature_type, saved = mainExpName(object), withColData = TRUE)
    }
    return(object)
}

#' Retrieve Assay
#'
#' @param object a SingleCellExperiment object
#' @param experiment an experiment name
#'
#' @return Main or alt experiment in an object
#' @export
#'
#' @examples
#' 
#' 
#' chevreul_sce <- mockSCE(ncells=200, ngenes=1000)
#' mainExpName(chevreul_sce) <- "gene"
#' retrieve_experiment(chevreul_sce, experiment = "gene")
retrieve_experiment <- function(object, experiment) {
    if (experiment %in% mainExpName(object)) {
        return(object)
    } else if (experiment %in% altExpNames(object)) {
        return(altExp(object, experiment))
    }
}

#' Query Assay
#'
#' @param object a SingleCellExperiment object
#' @param experiment an experiment name
#'
#' @return logical scalar indicating if experiment is present in object
#' @export
#'
#' @examples
#'
#'chevreul_sce <- mockSCE(ncells=200, ngenes=1000)
#' query_experiment(chevreul_sce, "gene")
query_experiment <- function(object, experiment) {
    return(experiment %in% c(mainExpName(object), altExpNames(object)))
}
