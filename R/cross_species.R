
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
convert_mouse_seu_to_human <- function(seu){
  new_rownames <- convert_symbols_by_species(src_genes = rownames(seu$gene), src_species = "mouse")

  new_rownames <- new_rownames %>%
    dplyr::group_by(human) %>%
    dplyr::mutate(duplicate_symbol = dplyr::row_number()) %>%
    dplyr::filter()

  keep_rows <- new_rownames$duplicate_symbol == 1 & !is.na(new_rownames$human)

  seu$gene <- seu$gene[keep_rows,]

  new_rownames <- new_rownames[keep_rows,] %>%
    dplyr::pull(human)

  seu_slots <- c("counts", "data", "scale.data", "meta.features")

  for (i in seu_slots){
    current_slot <- slot(seu$gene@assays$RNA, i)
    if(!(dim(current_slot) == c(0,0))){
      rownames(slot(seu$gene@assays$RNA, i)) <- new_rownames
    }
  }

  return(seu)
}

#' Convert Seurat Objects from Human to Mouse
#' @param seu
#'
#' @return
#' @export
#'
#' @examples
convert_human_seu_to_mouse <- function(seu){
  new_rownames <- convert_symbols_by_species(src_genes = rownames(seu$gene), src_species = "human")

  new_rownames <- new_rownames %>%
    dplyr::group_by(mouse) %>%
    dplyr::mutate(duplicate_symbol = dplyr::row_number()) %>%
    dplyr::filter()

  keep_rows <- new_rownames$duplicate_symbol == 1 & !is.na(new_rownames$mouse)

  seu$gene <- seu$gene[keep_rows,]

  new_rownames <- new_rownames[keep_rows,] %>%
    dplyr::pull(mouse)

  seu_slots <- c("counts", "data", "scale.data", "meta.features")



  for (i in seu_slots){
    rownames(slot(seu$gene@assays$RNA, i)) <- new_rownames
  }

  return(seu)
}


#' convert gene symbols between species
#'
#' @param src_genes
#' @param src_species
#'
#' @return
#' @export
#'
#' @examples
convert_symbols_by_species <- function(src_genes, src_species){
  # browser()
  if(src_species == "human"){
    dest_species = "mouse"
    src_species_mart = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    dest_species_mart = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    src_attribute = "hgnc_symbol"
    dest_attribute = "mgi_symbol"
  } else if (src_species == "mouse"){
    dest_species = "human"
    src_species_mart = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl")
    dest_species_mart = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    src_attribute = "mgi_symbol"
    dest_attribute = "hgnc_symbol"
  }

  genesV2 = biomaRt::getLDS(attributes = src_attribute, filters = src_attribute, values = src_genes , mart = src_species_mart, attributesL = dest_attribute, martL = dest_species_mart, uniqueRows=T)

  dest_symbols <- genesV2[match(src_genes, genesV2[,1]),]
  names(dest_symbols) <- c(src_species, dest_species)

  return(dest_symbols)
}

#' Integrate Seurat Objects from Moust to Human
#'
#' @param mouse_seu_list
#' @param human_seu_list
#'
#' @return
#' @export
#'
#' @examples
cross_species_integrate <- function(mouse_seu_list, human_seu_list, excluded_cells = NULL, ...){
  mouse_seu_list <- purrr::map(mouse_seu_list, convert_mouse_seu_to_human)

  seu_list <- c(mouse_seu_list, human_seu_list)

  seu_list <- purrr::transpose(seu_list)
  seu_list <- seu_list["gene"]

  integrated_seu <- purrr::map(seu_list, seurat_integrate)

  # cluster merged seurat objects
  integrated_seu <- purrr::map(integrated_seu, seurat_cluster, resolution = seq(0.2, 2.0, by = 0.2))

  # add read count column
  integrated_seu <- purrr::map(integrated_seu, add_read_count_col)

  # annotate cell cycle scoring to seurat objects

  integrated_seu <- purrr::map(integrated_seu, annotate_cell_cycle, feature = "gene")

  #annotate excluded cells

  if (!is.null(excluded_cells)){
    integrated_seu <- purrr::map(integrated_seu, annotate_excluded, excluded_cells)
  }

  # add marker genes to seurat objects

  integrated_seu <- purrr::map(integrated_seu, find_all_markers)

  return(integrated_seu)

}


