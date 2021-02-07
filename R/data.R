#' Retinal Cell Type Marker Genes
#'
#' Retinal Cell Type Marker Genes
#'
#' @format A data frame with 99 rows and 3 variables:
#' \describe{
#'   \item{\code{gene}}{character COLUMN_DESCRIPTION}
#'   \item{\code{cell_type}}{character COLUMN_DESCRIPTION}
#'   \item{\code{species}}{character COLUMN_DESCRIPTION}
#' }
#' @source Cui et al. (2020). Transcriptomic Analysis of the Developmental Similarities and Differences Between the Native Retina and Retinal Organoids. Invest. Ophthalmol. Vis. Sci. 61, 6–6.
"celltype_markers"


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
