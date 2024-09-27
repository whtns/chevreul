#' Convert SingleCellExperiment Objects from Human to Mouse
#' @param object Human SingleCellExperiment object
#' @param ... to be passed to \code{convert_symbols_by_species}
#'
#' @return a SingleCellExperiment object
#'
convert_human_object_to_mouse <- function(object, ...) {
    new_rownames <- convert_symbols_by_species(
        src_genes = rownames(object), src_species = "human")

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
#' @return a SingleCellExperiment object
#'
#'
convert_symbols_by_species <- function(src_genes, src_species) {
    if (src_species == "human") {
        dest_species <- "mouse"

        dest_symbols <- src_genes |>
            enframe("gene_index", "HGNC.symbol") |>
            left_join(human_to_mouse_homologs, by = "HGNC.symbol") |>
            distinct(HGNC.symbol, .keep_all = TRUE) |>
            mutate(MGI.symbol = case_when(
                is.na(MGI.symbol) ~ str_to_sentence(HGNC.symbol),
                TRUE ~ MGI.symbol
            )) |>
            select(-gene_index) |>
            identity()
    } else if (src_species == "mouse") {
        dest_species <- "human"

        dest_symbols <- src_genes |>
            enframe("gene_index", "MGI.symbol") |>
            left_join(human_to_mouse_homologs, by = "MGI.symbol") |>
            distinct(MGI.symbol, .keep_all = TRUE) |>
            mutate(HGNC.symbol = case_when(
                is.na(HGNC.symbol) ~ str_to_upper(MGI.symbol),
                TRUE ~ HGNC.symbol
            )) |>
            select(-gene_index) |>
            identity()
    }

    return(make.unique(dest_symbols[[2]]))
}
