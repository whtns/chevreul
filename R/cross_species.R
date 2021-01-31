
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
#' convert_mouse_seu_to_human(seurat_mouse_data)
#'

convert_mouse_seu_to_human <- function(seu){

  # transfer default species expression data to a species-specific assay
  seu[["mouse"]] <- seu[["gene"]]

  new_rownames <- convert_symbols_by_species(src_genes = rownames(seu), src_species = "mouse")

  # new_rownames <- new_rownames %>%
  #   dplyr::group_by(human) %>%
  #   dplyr::mutate(duplicate_symbol = dplyr::row_number()) %>%
  #   dplyr::filter()
  #
  # keep_rows <- new_rownames$duplicate_symbol == 1 & !is.na(new_rownames$human)
  #
  # subset_seu <- seu[keep_rows,]
  #
  # seu[["gene"]] <- subset_seu[["gene"]]

  seu_slots <- c("counts", "data", "scale.data", "meta.features")

  for (i in seu_slots){
    current_slot <- slot(seu@assays[["gene"]], i)
    if(!(dim(current_slot) == c(0,0))){
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
convert_human_seu_to_mouse <- function(seu, ...){

  new_rownames <- convert_symbols_by_species(src_genes = rownames(seu), src_species = "human")

  seu_slots <- c("counts", "data", "scale.data", "meta.features")



  for (i in seu_slots){
    rownames(slot(seu@assays[["gene"]], i)) <- new_rownames
  }

  return(seu)
}


#' convert gene symbols between species
#'
#' @param src_genes a vector of hgnc/mgi symbols
#' @param src_species "mouse" or "human"
#' @param host the biomart host to use try also useast.ensembl.org
#'
#' @return
#' @export
#'
#' @examples
convert_symbols_by_biomart <- function(src_genes, src_species, host = "uswest.ensembl.org"){

  if(src_species == "human"){
    dest_species = "mouse"
    src_species_mart = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = host)
    dest_species_mart = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = host)
    src_attribute = "hgnc_symbol"
    dest_attribute = "mgi_symbol"
  } else if (src_species == "mouse"){
    dest_species = "human"
    src_species_mart = biomaRt::useMart("ensembl", dataset = "mmusculus_gene_ensembl", host = host)
    dest_species_mart = biomaRt::useMart("ensembl", dataset = "hsapiens_gene_ensembl", host = host)
    src_attribute = "mgi_symbol"
    dest_attribute = "hgnc_symbol"
  }

  chunk2 <- function(x,n) split(x, cut(seq_along(x), n, labels = FALSE))

  src_genes <- chunk2(src_genes, 4)

  genesv2 <- purrr::map(src_genes, ~biomaRt::getLDS(attributes = src_attribute,
                                                    filters = src_attribute,
                                                    values = .x ,
                                                    mart = src_species_mart,
                                                    attributesL = dest_attribute,
                                                    martL = dest_species_mart,
                                                    uniqueRows=T)
  )

  genesv2 <- dplyr::bind_rows(genesv2)

  dest_symbols <- genesV2[match(src_genes, genesV2[,1]),]
  names(dest_symbols) <- c(src_species, dest_species)

  return(dest_symbols)
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
convert_symbols_by_species <- function(src_genes, src_species){

  # library(biomaRt)
  #
  # seu <- readRDS("~/single_cell_projects/integrated_projects/7-seq_050120/output/seurat/Final_dataset_duplicate_070320.rds")
  #
  # hs_genes <- rownames(seu)
  #
  # ## Not run:
  # ## The code to prepare the .Rda file file from the marker file is:
  # listMarts(host="www.ensembl.org")
  # human <- useMart(host="www.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")
  # mouse <- useMart(host="www.ensembl.org", "ENSEMBL_MART_ENSEMBL", dataset="mmusculus_gene_ensembl")
  #
  # hs_entrez <- "378938"
  #
  # human_to_mouse_homologs_manual = getLDS(attributes = c("hgnc_symbol","entrezgene_id","ensembl_gene_id"),
  #                                         filters = "hgnc_symbol", values = hs_genes,
  #                                         # filters = "entrezgene_id", values = hs_entrez,
  #                                         mart = human,
  #                                         attributesL = c("mgi_symbol","ensembl_gene_id","entrezgene_id"),
  #                                         martL = mouse)


  if(src_species == "human"){
    dest_species = "mouse"

    dest_symbols <- src_genes %>%
      tibble::enframe("gene_index", "HGNC.symbol") %>%
      dplyr::left_join(human_to_mouse_homologs, by = "HGNC.symbol") %>%
      dplyr::distinct(HGNC.symbol, .keep_all = TRUE) %>%
      dplyr::mutate(MGI.symbol = dplyr::case_when(is.na(MGI.symbol) ~ stringr::str_to_sentence(HGNC.symbol),
                                                   TRUE ~ MGI.symbol)) %>%
      dplyr::select(-gene_index) %>%
      identity()

  } else if (src_species == "mouse"){
    dest_species = "human"

    dest_symbols <- src_genes %>%
      tibble::enframe("gene_index", "MGI.symbol") %>%
      dplyr::left_join(mouse_to_human_homologs, by = "MGI.symbol") %>%
      dplyr::distinct(MGI.symbol, .keep_all = TRUE) %>%
      dplyr::mutate(HGNC.symbol = dplyr::case_when(is.na(HGNC.symbol) ~ stringr::str_to_upper(MGI.symbol),
                                                   TRUE ~ HGNC.symbol)) %>%
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
#'
#'
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

#' Update human gene symbols in seurat object
#'
#' @param seu
#'
#' @return
#' @export
#'
#' @examples
update_human_gene_symbols <- function(seu, assay = "gene"){
  browser()
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

  for (i in seu_slots){
    if (length(slot(seu@assays[[assay]], i)) > 0){
      rownames(slot(seu@assays[[assay]], i)) <- new_rownames
    }
  }

  return(seu)

}


