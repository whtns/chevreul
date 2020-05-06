#' Annotate Cell Cycle
#'
#' Annotate Cell Cycle for Gene and Transcript Seurat Objects
#'
#' @param seu_list A list of seurat objects
#'
#' @return
#' @export
#'
#' @examples

annotate_cell_cycle <- function(seu, feature, organism = "human", ...){


  # setdefaultassay to "RNA"
  DefaultAssay(seu) <- "RNA"

  s_genes <- cc.genes$s.genes
  g2m_genes <- cc.genes$g2m.genes

  s_transcripts <- cc.transcripts$s.genes
  g2m_transcripts <- cc.transcripts$g2m.genes

  if(organism == "mouse"){
    s_genes <- stringr::str_to_title(s_genes)
    g2m_genes <- stringr::str_to_title(g2m_genes)
    s_transcripts <- genes_to_transcripts(s_genes, organism = organism)
    g2m_transcripts <- genes_to_transcripts(g2m_genes, organism = organism)
  }

  if (feature == "gene"){
    seu <- CellCycleScoring(seu, s.features = s_genes, g2m.features = g2m_genes, set.ident = FALSE)

  } else if (feature == "transcript"){
    seu <- CellCycleScoring(seu, s.features = s_transcripts, g2m.features = g2m_transcripts, set.ident = FALSE)
  }

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
genes_to_transcripts <- function(genes, organism = "human") {

  if(organism == "human"){
    gene_table = annotables::grch38
    transcript_table = annotables::grch38_tx2gene
  } else if (organism == "mouse"){
    gene_table = annotables::grcm38
    transcript_table = annotables::grcm38_tx2gene
  }

  tibble::tibble(symbol = genes) %>%
    dplyr::left_join(gene_table, by = "symbol") %>%
    dplyr::left_join(transcript_table, by = "ensgene") %>%
    dplyr::pull("enstxp") %>%
    identity()
}

#' Annotate Low Read Count Category
#'
#'  Add a Read Count Categorical Variable to Seurat Object (based on nCount_RNA)
#'
#' @param seu A seurat object
#' @param thresh
#'
#' @return
#' @export
#'
#' @examples
add_read_count_col <- function(seu, thresh = 1e5){
  rc <- tibble::as_tibble(seu[["nCount_RNA"]], rownames = "sample_id") %>%
    dplyr::mutate(read_count = ifelse(nCount_RNA > thresh, NA, "low_read_count")) %>%
    dplyr::select(-nCount_RNA) %>%
    tibble::deframe()

  seu <- Seurat::AddMetaData(
    object = seu,
    metadata = rc,
    col.name = "read_count"
  )
}

#' Annotate Low Read Count Category
#'
#'  Add a Read Count Categorical Variable to Seurat Object (based on nCount_RNA)
#'
#' @param seu A seurat object
#'
#' @return
#' @export
#'
#' @examples
add_percent_mito <- function(seu, feature, organism = "human"){

  # mito_features = list()
  # if (organism == "human"){
  #   mito_regex = c(mt = "^MT-", mt_pseudo = "^MT.*P.*")
  #   mito_features$gene <- annotables::grch38 %>%
  #     dplyr::filter(stringr::str_detect(symbol, paste(mito_regex, collapse = "|"))) %>%
  #     dplyr::pull(symbol)
  #
  #   mito_features$transcript <-
  #     annotables::grch38 %>%
  #     dplyr::filter(symbol %in% mito_features$gene) %>%
  #     dplyr::left_join(annotables::grch38_tx2gene, by = "ensgene") %>%
  #     dplyr::pull(enstxp)
  # } else if (organism == "mouse"){
  #   mito_regex = c(mt = "^Mt-", mt_pseudo = "^Mt.*p.*")
  #   mito_features$gene <- annotables::grcm38 %>%
  #     dplyr::filter(stringr::str_detect(symbol, paste(mito_regex, collapse = "|"))) %>%
  #     dplyr::pull(symbol)
  #
  #   mito_features$transcript <-
  #     annotables::grcm38 %>%
  #     dplyr::filter(symbol %in% mito_features$gene) %>%
  #     dplyr::left_join(annotables::grcm38_tx2gene, by = "ensgene") %>%
  #     dplyr::pull(enstxp)
  # }

  mito_features <- mito_features[[organism]][[feature]]

  mito_features <- mito_features[mito_features %in% rownames(seu)]

  seu[["percent.mt"]] <- PercentageFeatureSet(seu, features = mito_features)

  return(seu)

}


#' Annotate Exclusion Criteria
#'
#' @param seu A seurat object
#' @param excluded_cells a named list of cells to be excluded of the form list(set1 = c(cell1, celll2), set2 = c(cell3, cell4))
#' all other cells will be marked NA in a column titled "excluded_because"
#'
#' @return
#' @export
#'
#' @examples
annotate_excluded <- function(seu, ...){
  # consider replacing
  # mutate(coalesce_var = coalesce(!!! syms(vars_select(names(.), dplyr::starts_with("my")))))

  excluded_cells <- list(...)

  excluded_cells <- purrr::map2(excluded_cells, names(excluded_cells), ~rep(.y, length(.x))) %>%
    unlist() %>%
    purrr::set_names(unlist(excluded_cells)) %>%
    tibble::enframe("sample_id", "excluded_because")

  excluded_because <- tibble::as_tibble(seu[["nCount_RNA"]], rownames = "sample_id") %>%
    dplyr::full_join(excluded_cells, by = "sample_id")

  if ("excluded_because.x" %in% colnames(excluded_because)){
    excluded_because <- dplyr::coalesce(excluded_because, excluded_because.x, excluded_because.y) %>%
      dplyr::select(-nCount_RNA) %>%
      tibble::deframe() %>%
      identity()
  } else {
    excluded_because <- select(excluded_because, -nCount_RNA) %>%
      tibble::deframe() %>%
      identity()
  }

  seu <- Seurat::AddMetaData(
    object = seu,
    metadata = excluded_because,
    col.name = "excluded_because"
  )

  return(seu)

}


