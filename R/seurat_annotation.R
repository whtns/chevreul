#' Annotate Cell Cycle
#'
#' Annotate Cell Cycle for Gene and Transcript Seurat Objects
#'
#' @param seu_list
#'
#' @return
#' @export
#'
#' @examples
annotate_cell_cycle <- function(seu_list){



  # S phase genes
  s.genes <- cc.genes[1:43]

  # G2/M phase genes
  g2m.genes <- cc.genes[44:97]

  # S phase transcripts
  s.transcripts <- genes_to_transcripts(s.genes)

  # G2/M phase transcripts
  g2m.transcripts <- genes_to_transcripts(g2m.genes)

  cc.features <- list(
    "gene" = list("s" = s.genes, "g2m" = g2m.genes),
    "transcript" = list("s" = s.transcripts, "g2m" = g2m.transcripts)
  )

  # setdefaultassay to "RNA"
  seu_list <- furrr::future_map(seu_list, SetDefaultAssay, "RNA")

  seu_list <- furrr::future_map2(seu_list, cc.features, ~CellCycleScoring(.x, s.features = .y$s, g2m.features = .y$g2m, set.ident = FALSE))

}

#' Ensembl Genes to Transcripts
#'
#' convert genes to transcripts
#'
#' @param genes
#'
#' @return
#' @export
#'
#' @examples
genes_to_transcripts <- function(genes) {
  tibble::tibble(symbol = genes) %>%
    dplyr::left_join(annotables::grch38, by = "symbol") %>%
    dplyr::left_join(annotables::grch38_tx2gene, by = "ensgene") %>%
    dplyr::pull("enstxp") %>%
    identity()
}

#' Annotate Low Read Count Category
#'
#'  Add a Read Count Categorical Variable to Seurat Object (based on nCount_RNA)
#'
#' @param seu
#' @param thresh
#'
#' @return
#' @export
#'
#' @examples
add_read_count_col <- function(seu, thresh = 1e5){
  rc <- as_tibble(seu[["nCount_RNA"]], rownames = "Sample_ID") %>%
    mutate(read_count = ifelse(nCount_RNA > thresh, "keep", "low_read_count")) %>%
    select(-nCount_RNA) %>%
    tibble::deframe()

  seu <- AddMetaData(
    object = seu,
    metadata = rc,
    col.name = "read_count"
  )
}


#' Annotate Exclusion Criteria
#'
#' @param seu
#' @param excluded_cells a named list of cells to be excluded of the form list(set1 = c(cell1, celll2), set2 = c(cell3, cell4))
#' all other cells will be marked "keep" in a column titled "excluded_because"
#'
#' @return
#' @export
#'
#' @examples
annotate_excluded <- function(seu, excluded_cells){

  excluded_cells <- furrr::future_map2(excluded_cells, names(excluded_cells), ~rep(.y, length(.x))) %>%
    unlist() %>%
    set_names(unlist(excluded_cells))

  excluded_because <- as_tibble(seu[["nCount_RNA"]], rownames = "Sample_ID") %>%
    mutate(excluded_because = ifelse(Sample_ID %in% names(excluded_cells), excluded_cells, "keep")) %>%
    select(-nCount_RNA) %>%
    tibble::deframe()

  seu <- AddMetaData(
    object = seu,
    metadata = excluded_because,
    col.name = "excluded_because"
  )

  return(seu)

}


