#' Annotate Cell Cycle
#'
#' Annotate Cell Cycle for Gene and Transcript Seurat Objects
#'
#' @param object A object
#' @param organism organism. Default "human"
#'
#' @return
#' @export
#'
#' @examples
#' annotate_cell_cycle(panc8, organism = "human")
#' annotate_cell_cycle(baron2016singlecell, organism = "mouse")
setGeneric("annotate_cell_cycle", function (object, organism = "human", ...)  standardGeneric("annotate_cell_cycle"))

setMethod("annotate_cell_cycle", "Seurat",
          function (object, organism = "human", ...)
          {
            s_genes <- cc.genes$s.genes
            g2m_genes <- cc.genes$g2m.genes
            if (organism == "mouse") {
              s_genes <- dplyr::filter(human_to_mouse_homologs, HGNC.symbol %in% s_genes) %>% dplyr::pull(MGI.symbol)
              g2m_genes <- dplyr::filter(human_to_mouse_homologs, HGNC.symbol %in% g2m_genes) %>% dplyr::pull(MGI.symbol)
            }
            object <- CellCycleScoring(object, s.features = s_genes, g2m.features = g2m_genes, set.ident = FALSE)
            return(object)


          }
)

setMethod("annotate_cell_cycle", "SingleCellExperiment",
          function (object, organism = "human", ...)
          {
            hs_pairs0 <- cc.genes.cyclone
            if (organism == "mouse") {
              s_genes <- dplyr::filter(human_to_mouse_homologs, HGNC.symbol %in% s_genes) %>% dplyr::pull(MGI.symbol)
              g2m_genes <- dplyr::filter(human_to_mouse_homologs, HGNC.symbol %in% g2m_genes) %>% dplyr::pull(MGI.symbol)
            }

            assignments <- scran::cyclone(object , hs_pairs0, gene.names=rownames(object))
            colData(object)[colnames(assignments$scores)] <- assignments$scores
            colData(object)["phases"] <- assignments$phases
            return(object)

          }
)

#' Gene Symbols to Ensembl Transcript Ids
#'
#' convert hgnc gene symbols to ensembl transcript ids
#'
#' @param genes
#'
#' @return
#' @export
#'
#' @examples
#'
#' genes_to_transcripts("RXRG")
#'
genes_to_transcripts <- function(genes, organism = "human") {
  if (organism == "human") {
    feature_table <- ensembldb::transcripts(EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86,
                                            columns = c("gene_name", "gene_biotype", "gene_id"),
                                            return.type="DataFrame")
  } else if (organism == "mouse") {
    feature_table <- ensembldb::transcripts(EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79,
                                            columns = c("gene_name", "gene_biotype", "gene_id"),
                                            return.type="DataFrame")
  }

  feature_table %>%
    as_tibble() %>%
    dplyr::filter(gene_name %in% genes) %>%
    dplyr::pull(tx_id)

}

#' Ensembl Transcript Ids to Gene Symbols
#'
#' convert ensembl transcript ids to hgnc gene symbols
#'
#' @param transcripts
#'
#' @return
#' @export
#'
#' @examples
#'
#' RXRG_transcripts_hs <- c("ENST00000359842", "ENST00000470566", "ENST00000465764", "ENST00000619224")
#'
#' transcripts_to_genes(transcripts = RXRG_transcripts_hs)
#'
transcripts_to_genes <- function(transcripts, organism = "human") {
  if (organism == "human") {
    gene_table <- annotables::grch38
    transcript_table <- annotables::grch38_tx2gene
  } else if (organism == "mouse") {
    gene_table <- annotables::grcm38
    transcript_table <- annotables::grcm38_tx2gene
  }

  tibble::tibble(enstxp = transcripts) %>%
    dplyr::left_join(transcript_table, by = "enstxp") %>%
    dplyr::left_join(gene_table, by = "ensgene") %>%
    dplyr::pull("symbol") %>%
    identity()
}

#' Annotate Low Read Count Category
#'
#' Add a Read Count Categorical Variable to Seurat Object (based on nCount_RNA)
#'
#' @param object A object
#' @param thresh Set a threshold for low read count
#'
#' @return
#' @export
#'
#' @examples
setGeneric("add_read_count_col", function (object, thresh = 1e+05)  standardGeneric("add_read_count_col"))

setMethod("add_read_count_col", "Seurat",
          function (object, thresh = 1e+05)
          {
            rc <- object[["nCount_gene"]] < thresh
            object <- Seurat::AddMetaData(object = object, metadata = rc, col.name = "low_read_count")
          }
)

setMethod("add_read_count_col", "SingleCellExperiment",
          function (object, thresh = 1e+05)
          {
            rc <- object[["nCount_gene"]] < thresh


             #object <- Seurat::AddMetaData(object = object, metadata = rc, col.name = "low_read_count")
          }
)

#' Annotate percent mitochondrial reads per cell
#'
#'  Add a Read Count Categorical Variable to Seurat Object (based on nCount_RNA)
#'
#' @param object A object
#' @param organism mouse
#' @param object_assay gene
#'
#' @return
#' @export
#'
#' @examples
#' add_percent_mito(panc8)
#' add_percent_mito(baron2016singlecell, organism = "mouse")
setGeneric("add_percent_mito", function (object, organism = "human", object_assay = "gene")  standardGeneric("add_percent_mito"))

setMethod("add_percent_mito", "Seurat",
          function (object, object_assay = "gene")
          {
            object[["percent.mt"]] <- PercentageFeatureSet(object = object, pattern = "^MT-")
            return(object)
          }
)

setMethod("add_percent_mito", "SingleCellExperiment",
          function (object, organism = "human", object_assay = "gene")
          {
            # mito_features <- mito_features[[organism]][["gene"]]
            # mito_features <- mito_features[mito_features %in% rownames(object[[object_assay]])]
            # object[["percent.mt"]] <- PercentageFeatureSet(object, features = mito_features)
            # return(object)
            is.mito <-grepl("^MT-*", rownames(object))
            object <- scuttle::addPerCellQCMetrics(object, subsets=list(Mito=is.mito))
            return(object)
          }
)

#' Annotate Exclusion Criteria
#'
#' @param object A object
#' @param excluded_cells a named list of cells to be excluded of the form list(set1 = c(cell1, celll2), set2 = c(cell3, cell4))
#' all other cells will be marked NA in a column titled "excluded_because"
#'
#' @return
#' @export
#'
#' @examples
annotate_excluded <- function(object, ...) {
  # consider replacing
  # mutate(coalesce_var = coalesce(!!! syms(vars_select(names(.), dplyr::starts_with("my")))))

  excluded_cells <- list(...)

  excluded_cells <- purrr::map2(excluded_cells, names(excluded_cells), ~ rep(.y, length(.x))) %>%
    unlist() %>%
    purrr::set_names(unlist(excluded_cells)) %>%
    tibble::enframe("sample_id", "excluded_because")

  excluded_because <- tibble::as_tibble(object[["nCount_RNA"]], rownames = "sample_id") %>%
    dplyr::full_join(excluded_cells, by = "sample_id")

  if ("excluded_because.x" %in% colnames(excluded_because)) {
    excluded_because <- dplyr::coalesce(excluded_because, excluded_because.x, excluded_because.y) %>%
      dplyr::select(-nCount_RNA) %>%
      tibble::deframe() %>%
      identity()
  } else {
    excluded_because <- select(excluded_because, -nCount_RNA) %>%
      tibble::deframe() %>%
      identity()
  }

  object <- Seurat::AddMetaData(
    object = object,
    metadata = excluded_because,
    col.name = "excluded_because"
  )

  return(object)
}
