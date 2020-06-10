#' Run Seurat Integration
#'
#' run batch correction, followed by:
#' 1) stashing of batches in metadata 'batch'
#' 2) clustering with resolution 0.2 to 2.0 in increments of 0.2
#' 3) saving to <proj_dir>/output/sce/<feature>_seu_<suffix>.rds
#'
#' @param suffix a suffix to be appended to a file save in output dir
#' @param seu_list
#' @param resolution
#' @param feature
#' @param algorithm
#' @param organism
#' @param ...
#'
#' @return
#' @export
#'
#'
#' @examples
seurat_integration_pipeline <- function(seu_list, feature, resolution, suffix = '', algorithm = 1, organism, annotate_cell_cycle = FALSE, annotate_percent_mito = FALSE, ...) {

  integrated_seu <- seurat_integrate(seu_list, ...)

  # cluster merged seurat objects
  integrated_seu <- seurat_cluster(integrated_seu, resolution = resolution, algorithm = algorithm, ...)

  integrated_seu <- find_all_markers(integrated_seu)

  if (feature == "gene"){
    integrated_seu <- getEnrichedPathways(integrated_seu)
  }

  # add read count column
  integrated_seu <- add_read_count_col(integrated_seu)

  # annotate cell cycle scoring to seurat objects
  if (annotate_cell_cycle){
    integrated_seu <- annotate_cell_cycle(integrated_seu, feature, ...)
  }

  # annotate mitochondrial percentage in seurat metadata
  if (annotate_percent_mito){
    integrated_seu <- add_percent_mito(integrated_seu, feature, ...)
  }

  #annotate excluded cells
  # integrated_seu <- annotate_excluded(integrated_seu, excluded_cells)

  return(integrated_seu)

}

#' Run Seurat Pipeline
#'
#' Preprocess, Cluster and Reduce Dimensions for a single seurat object
#'
#' @param seu
#' @param resolution
#'
#' @return
#' @export
#'
#' @examples
seurat_pipeline <- function(seu, feature = "gene", resolution=0.6, reduction = "pca", annotate_cell_cycle = TRUE, annotate_percent_mito = TRUE, ...){

  seu <- seurat_preprocess(seu, scale = T, ...)

  # PCA
  seu <- seurat_reduce_dimensions(seu, check_duplicates = FALSE, reduction = reduction, ...)

  seu <- seurat_cluster(seu = seu, resolution = resolution, reduction = reduction, ...)

  seu <- find_all_markers(seu, resolution = resolution)

  if (feature == "gene"){
    seu <- getEnrichedPathways(seu)
  }

  # annotate low read count category in seurat metadata
  seu <- seuratTools::add_read_count_col(seu)

  # annotate cell cycle scoring to seurat objects
  if (annotate_cell_cycle){
    seu <- annotate_cell_cycle(seu, feature, ...)
  }

  # annotate mitochondrial percentage in seurat metadata
  if (annotate_percent_mito){
    seu <- add_percent_mito(seu, feature, ...)
  }

  return(seu)
}
