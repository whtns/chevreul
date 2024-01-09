
# you must set the rownames of the raw.data, data, and scale.data slots directly,
# as well as the the rownames of any gene loadings for PCAs you may have calculated.

#' Convert Seurat Objects from Mouse to Human
#'
#' @param object Mouse object
#'
#' @return
#' @export
#'
#' @examples
#'
#' convert_mouse_object_to_human(baron2016singlecell)
convert_mouse_object_to_human <- function(object) {

  # transfer default species expression data to a species-specific assay
  object[["mouse"]] <- object[["gene"]]

  new_rownames <- convert_symbols_by_species(src_genes = rownames(object), src_species = "mouse")

  object_slots <- c("counts", "data", "scale.data", "meta.features")

  for (i in object_slots) {
    current_slot <- slot(object@assays[["gene"]], i)
    if (!(dim(current_slot) == c(0, 0))) {
      rownames(slot(object@assays[["gene"]], i)) <- new_rownames
    }
  }

  return(object)
}

#' Convert Seurat Objects from Human to Mouse
#' @param object Human Seurat object
#' @param ... to be passed to \code{convert_symbols_by_species}
#'
#' @return
#' @export
#'
#' @examples
convert_human_object_to_mouse <- function(object, ...) {
  new_rownames <- convert_symbols_by_species(src_genes = rownames(object), src_species = "human")

  object_slots <- c("counts", "data", "scale.data", "meta.features")



  for (i in object_slots) {
    rownames(slot(object@assays[["gene"]], i)) <- new_rownames
  }

  return(object)
}

#' Convert gene symbols between mouse and human
#'
#' @param src_genes Source gene symbol to be converted
#' @param src_species Source species
#'
#' @return
#' @export
#'
#' @examples
#'
#' convert_symbols_by_species("RXRG", "human")
#'
#' convert_symbols_by_species("Rxrg", "mouse")
#'
convert_symbols_by_species <- function(src_genes, src_species) {

  if (src_species == "human") {
    dest_species <- "mouse"

    dest_symbols <- src_genes %>%
      tibble::enframe("gene_index", "HGNC.symbol") %>%
      dplyr::left_join(human_to_mouse_homologs, by = "HGNC.symbol") %>%
      dplyr::distinct(HGNC.symbol, .keep_all = TRUE) %>%
      dplyr::mutate(MGI.symbol = dplyr::case_when(
        is.na(MGI.symbol) ~ stringr::str_to_sentence(HGNC.symbol),
        TRUE ~ MGI.symbol
      )) %>%
      dplyr::select(-gene_index) %>%
      identity()
  } else if (src_species == "mouse") {
    dest_species <- "human"

    dest_symbols <- src_genes %>%
      tibble::enframe("gene_index", "MGI.symbol") %>%
      dplyr::left_join(human_to_mouse_homologs, by = "MGI.symbol") %>%
      dplyr::distinct(MGI.symbol, .keep_all = TRUE) %>%
      dplyr::mutate(HGNC.symbol = dplyr::case_when(
        is.na(HGNC.symbol) ~ stringr::str_to_upper(MGI.symbol),
        TRUE ~ HGNC.symbol
      )) %>%
      dplyr::select(-gene_index) %>%
      # dplyr::mutate(HGNC.symbol = make.unique(HGNC.symbol)) %>%
      identity()
  }

  return(make.unique(dest_symbols[[2]]))
}

#' Integrate Seurat Objects from Mouse to Human
#'
#' @param mouse_object_list Mouse Seurat object
#' @param human_object_list Human Seurat object
#'
#' @return
#' @export
#'
#' @examples
#'
#' cross_species_integrate(list(baron2016singlecell = baron2016singlecell), list(panc8 = panc8))
#'
cross_species_integrate <- function(mouse_object_list, human_object_list, excluded_cells = NULL, ...) {
  mouse_object_list <- purrr::map(mouse_object_list, convert_mouse_object_to_human)

  object_list <- c(mouse_object_list, human_object_list)


  integrated_object <- object_integrate(object_list)

  # cluster merged objects
  integrated_object <- object_cluster(integrated_object, resolution = seq(0.2, 2.0, by = 0.2))

  # add read count column
  integrated_object <- add_read_count_col(integrated_object)

  # annotate cell cycle scoring to objects

  integrated_object <- annotate_cell_cycle(integrated_object, feature = "gene")

  # annotate excluded cells

  if (!is.null(excluded_cells)) {
    integrated_object <- annotate_excluded(integrated_object, excluded_cells)
  }

  # add marker genes to objects

  integrated_object <- find_all_markers(integrated_object)

  return(integrated_object)
}

#' Update human gene symbols in object
#'
#' @param object A Seurat object
#' @param assay Assay to use, Default = "gene"
#'
#' @return
#' @export
#'
#' @examples
update_human_gene_symbols <- function(object, assay = "gene") {
  # browser()

  ensdb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86

  symbols <- rownames(object[[assay]])

  new_rownames <-
    AnnotationDbi::mapIds(ensdb, symbols, keytype = "SYMBOL", columns = c("SYMBOL", "GENEID")) %>%
    tibble::enframe("old_symbol", "ensgene")

  rownames(new_rownames) <- new_rownames$old_symbol

  object[[assay]] <- Seurat::AddMetaData(object[[assay]], new_rownames)

  new_rownames <-
    new_rownames %>%
    dplyr::left_join(annotables::grch38, by = "ensgene") %>%
    dplyr::distinct(old_symbol, .keep_all = TRUE) %>%
    dplyr::mutate(new_symbol = symbol) %>%
    dplyr::mutate(symbol = dplyr::coalesce(new_symbol, old_symbol)) %>%
    # tidyr::drop_na(symbol) %>%
    # dplyr::pull(symbol) %>%
    identity()

  object_slots <- c("counts", "data", "scale.data", "meta.features")

  for (i in object_slots) {
    if (length(slot(object@assays[[assay]], i)) > 0) {
      rownames(slot(object@assays[[assay]], i)) <- make.unique(new_rownames$symbol)
    }
  }

  variable_features <- VariableFeatures(object[[assay]])
  if(length(variable_features) > 1){
    new_variable_features <-
      dplyr::filter(new_rownames, old_symbol %in% variable_features) %>%
      dplyr::pull(symbol)

    VariableFeatures(object[[assay]]) <- new_variable_features

  }

  return(object)
}
