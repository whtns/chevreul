#' Convert SingleCellExperiment Objects from Mouse to Human
#'
#' @param object Mouse object
#'
#' @return a single cell object
#' @export
#'
#' @examples
#'
#' convert_mouse_object_to_human(baron2016singlecell)
convert_mouse_object_to_human <- function(object) {
    # transfer default species expression data to a species-specific experiment
    object[["mouse"]] <- object[["gene"]]

    new_rownames <- convert_symbols_by_species(src_genes = rownames(object), src_species = "mouse")

    object_slots <- c("counts", "data", "scale.data", "meta.features")

    for (i in object_slots) {
        current_slot <- slot(object@experiments[["gene"]], i)
        if (!(dim(current_slot) == c(0, 0))) {
            rownames(slot(object@experiments[["gene"]], i)) <- new_rownames
        }
    }

    return(object)
}

#' Convert SingleCellExperiment Objects from Human to Mouse
#' @param object Human SingleCellExperiment object
#' @param ... to be passed to \code{convert_symbols_by_species}
#'
#' @return a single cell object
#' @export
convert_human_object_to_mouse <- function(object, ...) {
    new_rownames <- convert_symbols_by_species(src_genes = rownames(object), src_species = "human")

    object_slots <- c("counts", "data", "scale.data", "meta.features")



    for (i in object_slots) {
        rownames(slot(object@experiments[["gene"]], i)) <- new_rownames
    }

    return(object)
}

#' Convert gene symbols between mouse and human
#'
#' @param src_genes Source gene symbol to be converted
#' @param src_species Source species
#'
#' @return a single cell object
#' @export
#'
#' @examples
#'
#' convert_symbols_by_species("NRL", "human")
#'
#' convert_symbols_by_species("NRL", "mouse")
#'
convert_symbols_by_species <- function(src_genes, src_species) {
    if (src_species == "human") {
        dest_species <- "mouse"

        dest_symbols <- src_genes %>%
            enframe("gene_index", "HGNC.symbol") %>%
            left_join(human_to_mouse_homologs, by = "HGNC.symbol") %>%
            distinct(HGNC.symbol, .keep_all = TRUE) %>%
            mutate(MGI.symbol = case_when(
                is.na(MGI.symbol) ~ str_to_sentence(HGNC.symbol),
                TRUE ~ MGI.symbol
            )) %>%
            select(-gene_index) %>%
            identity()
    } else if (src_species == "mouse") {
        dest_species <- "human"

        dest_symbols <- src_genes %>%
            enframe("gene_index", "MGI.symbol") %>%
            left_join(human_to_mouse_homologs, by = "MGI.symbol") %>%
            distinct(MGI.symbol, .keep_all = TRUE) %>%
            mutate(HGNC.symbol = case_when(
                is.na(HGNC.symbol) ~ str_to_upper(MGI.symbol),
                TRUE ~ HGNC.symbol
            )) %>%
            select(-gene_index) %>%
            identity()
    }

    return(make.unique(dest_symbols[[2]]))
}

#' Integrate SingleCellExperiment Objects from Mouse to Human
#'
#' @param mouse_object_list List of mouse single cell objects
#' @param human_object_list List of human single cell objects
#'
#' @return a single cell object
#' @export
#'
#' @examples
#'
#' cross_species_integrate(list(baron2016singlecell = baron2016singlecell), list(human_gene_transcript_sce = human_gene_transcript_sce))
#'
cross_species_integrate <- function(mouse_object_list, human_object_list){
    mouse_object_list <- map(mouse_object_list, convert_mouse_object_to_human)

    object_list <- c(mouse_object_list, human_object_list)


    integrated_object <- object_integrate(object_list)

    integrated_object <- object_reduce_dimensions(integrated_object)

    # cluster merged objects
    integrated_object <- object_cluster(integrated_object, resolution = seq(0.2, 2.0, by = 0.2))

    # annotate cell cycle scoring to objects

    integrated_object <- annotate_cell_cycle(integrated_object)

    # add marker genes to objects

    integrated_object <- find_all_markers(integrated_object)

    return(integrated_object)
}

#' Update human gene symbols in object
#'
#' @param object A SingleCellExperiment object
#' @param experiment Assay to use, Default = "gene"
#'
#' @return a single cell object
#' @export
#' @importFrom EnsDb.Hsapiens.v86 EnsDb.Hsapiens.v86
#'
#' @examples
#'
#'
update_human_gene_symbols <- function(object, experiment = "gene") {

    ensdb <- EnsDb.Hsapiens.v86

    symbols <- rownames(object[[experiment]])

    new_rownames <-
        mapIds(ensdb, symbols, keytype = "SYMBOL", columns = c("SYMBOL", "GENEID")) %>%
        enframe("old_symbol", "ensgene")

    rownames(new_rownames) <- new_rownames$old_symbol

    object <- SingleCellExperiment::AddMetaData(object, new_rownames)

    new_rownames <-
        new_rownames %>%
        left_join(annotables::grch38, by = "ensgene") %>%
        distinct(old_symbol, .keep_all = TRUE) %>%
        mutate(new_symbol = symbol) %>%
        mutate(symbol = coalesce(new_symbol, old_symbol)) %>%
        identity()

    object_slots <- c("counts", "data", "scale.data", "meta.features")

    for (i in object_slots) {
        if (length(slot(object@experiments[[experiment]], i)) > 0) {
            rownames(slot(object@experiments[[experiment]], i)) <- make.unique(new_rownames$symbol)
        }
    }

    variable_features <- get_variable_features(object[[experiment]])
    if (length(variable_features) > 1) {
        new_variable_features <-
            filter(new_rownames, old_symbol %in% variable_features) %>%
            pull(symbol)

        get_variable_features(object[[experiment]]) <- new_variable_features
    }

    return(object)
}
