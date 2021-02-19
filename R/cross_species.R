
# you must set the rownames of the raw.data, data, and scale.data slots directly,
# as well as the the rownames of any gene loadings for PCAs you may have calculated.

#' Convert Seurat Objects from Mouse to Human
#'
#' @param seu
#'
#' @return
#' @export
#'
#' @examples
#'
#' convert_mouse_seu_to_human(baron2016singlecell)
convert_mouse_seu_to_human <- function(seu) {

  # transfer default species expression data to a species-specific assay
  seu[["mouse"]] <- seu[["gene"]]

  new_rownames <- convert_symbols_by_species(src_genes = rownames(seu), src_species = "mouse")

  seu_slots <- c("counts", "data", "scale.data", "meta.features")

  for (i in seu_slots) {
    current_slot <- slot(seu@assays[["gene"]], i)
    if (!(dim(current_slot) == c(0, 0))) {
      rownames(slot(seu@assays[["gene"]], i)) <- new_rownames
    }
  }

  return(seu)
}

#' Convert Seurat Objects from Human to Mouse
#' @param seu
#' @param ... to be passed to \code{convert_symbols_by_species}
#'
#' @return
#' @export
#'
#' @examples
convert_human_seu_to_mouse <- function(seu, ...) {
  new_rownames <- convert_symbols_by_species(src_genes = rownames(seu), src_species = "human")

  seu_slots <- c("counts", "data", "scale.data", "meta.features")



  for (i in seu_slots) {
    rownames(slot(seu@assays[["gene"]], i)) <- new_rownames
  }

  return(seu)
}

#' Convert gene symbols between mouse and human
#'
#' @param src_genes
#' @param src_species
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
#' @param mouse_seu_list
#' @param human_seu_list
#'
#' @return
#' @export
#'
#' @examples
#'
#' cross_species_integrate(list(baron2016singlecell = baron2016singlecell), list(panc8 = panc8))
#'
cross_species_integrate <- function(mouse_seu_list, human_seu_list, excluded_cells = NULL, ...) {
  mouse_seu_list <- purrr::map(mouse_seu_list, convert_mouse_seu_to_human)

  seu_list <- c(mouse_seu_list, human_seu_list)


  integrated_seu <- seurat_integrate(seu_list)

  # cluster merged seurat objects
  integrated_seu <- seurat_cluster(integrated_seu, resolution = seq(0.2, 2.0, by = 0.2))

  # add read count column
  integrated_seu <- add_read_count_col(integrated_seu)

  # annotate cell cycle scoring to seurat objects

  integrated_seu <- annotate_cell_cycle(integrated_seu, feature = "gene")

  # annotate excluded cells

  if (!is.null(excluded_cells)) {
    integrated_seu <- annotate_excluded(integrated_seu, excluded_cells)
  }

  # add marker genes to seurat objects

  integrated_seu <- find_all_markers(integrated_seu)

  return(integrated_seu)
}

#' Update human gene symbols in seurat object
#'
#' @param seu
#'
#' @return
#' @export
#'
#' @examples
update_human_gene_symbols <- function(seu, assay = "gene") {
  ensdb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86

  symbols <- rownames(seu[[assay]])

  new_rownames <-
    AnnotationDbi::mapIds(ensdb, symbols, keytype = "SYMBOL", columns = c("SYMBOL", "GENEID")) %>%
    tibble::enframe("old_symbol", "ensgene") %>%
    dplyr::left_join(annotables::grch38, by = "ensgene") %>%
    dplyr::distinct(old_symbol, .keep_all = TRUE) %>%
    dplyr::mutate(symbol = dplyr::coalesce(symbol, old_symbol)) %>%
    # tidyr::drop_na(symbol) %>%
    dplyr::pull(symbol) %>%
    identity()

  seu_slots <- c("counts", "data", "scale.data", "meta.features")

  for (i in seu_slots) {
    if (length(slot(seu@assays[[assay]], i)) > 0) {
      rownames(slot(seu@assays[[assay]], i)) <- new_rownames
    }
  }

  return(seu)
}
