

#' Gene Symbols to Ensembl Transcript Ids
#'
#' convert hgnc gene symbols to ensembl transcript ids
#'
#' @param symbols character vector of gene symbols
#'
#' @return a vector of transcripts
#' @export
#'
#' @examples
#'
#' genes_to_transcripts("NRL")
genes_to_transcripts <- function(symbols, organism = "human") {
    if (organism == "human") {
        feature_table <- transcripts(EnsDb.Hsapiens.v86,
            columns = c("gene_name", "gene_biotype", "gene_id"),
            return.type = "DataFrame"
        )
    } else if (organism == "mouse") {
        feature_table <- transcripts(EnsDb.Mmusculus.v79,
            columns = c("gene_name", "gene_biotype", "gene_id"),
            return.type = "DataFrame"
        )
    }

  feature_table[(feature_table[["gene_name"]] %in% symbols),"tx_id"]

}

#' Ensembl Transcript Ids to Gene Symbols
#'
#' Convert ensembl transcript ids to hgnc gene symbols
#'
#' @param transcripts transcripts
#' @param organism human or mouse
#'
#' @return a vector of gene symbols
#' @export
#'
#' @examples
#'
#' NRL_transcripts_hs <- c("ENST00000359842", "ENST00000470566", "ENST00000465764", "ENST00000619224")
#'
#' transcripts_to_genes(transcripts = NRL_transcripts_hs)
#'
transcripts_to_genes <- function(transcripts, organism = "human") {
    if (organism == "human") {
        gene_table <- annotables::grch38
        transcript_table <- annotables::grch38_tx2gene
    } else if (organism == "mouse") {
        gene_table <- annotables::grcm38
        transcript_table <- annotables::grcm38_tx2gene
    }

    tibble(enstxp = transcripts) %>%
        left_join(transcript_table, by = "enstxp") %>%
        left_join(gene_table, by = "ensgene") %>%
        pull("symbol") %>%
        identity()
}

#' Annotate percent mitochondrial reads per cell
#'
#'  Add a Percentage of Mitochondrial Read Count Categorical Variable to the
#'  Object (based on nCount_RNA)
#'
#' @param object A object
#' @param organism mouse
#' @param experiment gene
#'
#' @return a single cell obejct with cell metadata column containing mitochondrial percentage
#' @export
#' @importFrom scuttle addPerCellQCMetrics
#' @examples
#' chevreul_sce <- chevreuldata::human_gene_transcript_sce()
#' add_percent_mito(chevreul_sce)
#'
add_percent_mito <- function(object, experiment = "gene") {
        is.mito <- grepl("^MT-*", rownames(object))
        object <- addPerCellQCMetrics(object, subsets = list(Mito = is.mito))
        return(object)
    }

