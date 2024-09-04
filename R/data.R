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
