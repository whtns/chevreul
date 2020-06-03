#' Run Seurat Integration
#'
#' run batch correction, followed by:
#' 1) stashing of batches in metadata 'batch'
#' 2) clustering with resolution 0.2 to 2.0 in increments of 0.2
#' 3) saving to <proj_dir>/output/sce/<feature>_seu_<suffix>.rds
#'
#' @param suffix a suffix to be appended to a file save in output dir
#' @param seus
#' @param resolution
#' @param ...
#' @param feature
#'
#' @return
#' @export
#'
#'
#' @examples
seurat_integration_pipeline <- function(seus, feature, resolution, suffix = '', algorithm = 1, experiment_name, organism, ...) {

  corrected_seu <- seurat_integrate(seus, ...)

  # cluster merged seurat objects
  corrected_seu <- seurat_cluster(corrected_seu, resolution = resolution, algorithm = algorithm, ...)

  corrected_seu <- find_all_markers(corrected_seu)

  # add read count column
  corrected_seu <- add_read_count_col(corrected_seu)

  # annotate cell cycle scoring to seurat objects

  corrected_seu <- annotate_cell_cycle(corrected_seu, feature, ...)

  # annotate mitochondrial percentage in seurat metadata
  corrected_seu <- add_percent_mito(corrected_seu, feature, ...)

  #annotate excluded cells

  # corrected_seu <- annotate_excluded(corrected_seu, excluded_cells)

  # corrected_seu <- save_seurat(corrected_seu, feature = feature, suffix = suffix, ...)

  corrected_seu <- record_experiment_data(corrected_seu, experiment_name, organism)

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
seurat_pipeline <- function(seu = seu, resolution=0.6, reduction = "pca", ...){

  seu <- seurat_preprocess(seu, scale = T, ...)

  # PCA
  seu <- seurat_reduce_dimensions(seu, check_duplicates = FALSE, reduction = reduction, ...)

  seu <- seurat_cluster(seu = seu, resolution = resolution, reduction = reduction, ...)

  seu <- find_all_markers(seu, resolution = resolution, reduction = reduction)

  return(seu)
}
