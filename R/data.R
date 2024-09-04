#' Small example SingleCellExperiment
#'
#' created with scuttle::mockSCE
#'
#' @format An SCE with 200 cells and 1000 genes
#' @source scuttle::mockSCE
"small_example_dataset"

#' Tiny example SingleCellExperiment
#'
#' subset to only NRL from chevreuldata::human_gene_transcript_sce()
#'
#' @format An SCE with only expression of NRL gene and NRL transripts
#' @source chevreuldata::human_gene_transcript_sce()
"tiny_sce"

#' Gene Homologs Between Human and Mouse
#'
#' Homologs drawn from Biomart
#'
#' @format A data frame with 23188 rows and 2 columns
#' \describe{
#'   \item{HGNC.symbol}{human gene symbols}
#'   \item{MGI.symbol}{mouse gene symbols}
#'   ...
#' }
#' @source bioMart
"human_to_mouse_homologs"

#' Cyclone cell cycle pairs by symbol
#'
#' cell cycle genes with paired expression represented by HGNC symbol
#'
#' @format a list of dataframes with G1, G2, and S gene expression
#' \describe{
#'   \item{G1}{G1 gene symbols}
#'   \item{G2}{G2 gene symbols}
#'   \item{S}{S gene symbols}
#'   ...
#' }
#' @source cyclone
"cc.genes.cyclone"

#' Human annotation data
#'
#' Human (*Homo sapiens*) annotations based on
#' genome assembly GRCH38 from Ensembl.
#'
#' @docType data
#' @keywords datasets
#'
#' @details
#' Variables:
#' 
#' - ensgene
#' - entrez
#' - symbol
#' - chr
#' - start
#' - end
#' - strand
#' - biotype
#' - description
#'
#' @source \url{http://ensembl.org/homo_sapiens}
#'
#' @examples
#' head(grch38)
"grch38"

#' Human transcripts to genes
#'
#' Lookup table for converting Human (*Homo sapiens*)
#' Ensembl transcript IDs to gene IDs based on genome assembly
#' GRCH38 from Ensembl.
#'
#' @docType data
#' @keywords datasets
#'
#' @details
#' Variables:
#' 
#' - enstxp
#' - ensgene
#'
#' @source \url{http://ensembl.org/homo_sapiens}
#'
#' @examples
#' head(grch38_tx2gene)
"grch38_tx2gene"