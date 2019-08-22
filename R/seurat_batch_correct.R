
## ----load-funcs----------------------------------------------------------

#' Batch Correct Multiple Seurat Objects
#'
#' @param seu_list
#' @param method
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
seurat_batch_correct <- function(seu_list, method = "cca", ...) {
  #browser()
  # To construct a reference we will identify ‘anchors’ between the individual datasets. First, we split the combined object into a list, with each dataset as an element.

  # Prior to finding anchors, we perform standard preprocessing (log-normalization), and identify variable features individually for each. Note that Seurat v3 implements an improved method for variable feature selection based on a variance stabilizing transformation ("vst")

  for (i in 1:length(x = seu_list)) {
    seu_list[[i]] <- seurat_preprocess(seu_list[[i]], scale = TRUE, ...)
    seu_list[[i]]$batch <- names(seu_list)[[i]]
  }

  if (method == "rpca"){
    # scale and run pca for each separate batch in order to use reciprocal pca instead of cca
    features <- SelectIntegrationFeatures(object.list = seu_list)
    seu_list <- purrr::map(seu_list, Seurat::ScaleData, features = features)
    seu_list <- purrr::map(seu_list, Seurat::RunPCA, features = features)
    seu_list.anchors <- FindIntegrationAnchors(object.list = seu_list, reduction = "rpca", dims = 1:30)
  } else if (method == "cca"){
    # Next, we identify anchors using the FindIntegrationAnchors function, which takes a list of Seurat objects as input.
    seu_list.anchors <- Seurat::FindIntegrationAnchors(object.list = seu_list, dims = 1:30, k.filter = 50)
  }

  # proceed with integration
  seu_list.integrated  <- IntegrateData(anchorset = seu_list.anchors, dims = 1:30)

  # Next, we identify anchors using the FindIntegrationAnchors function, which takes a list of Seurat objects as input.

  # #stash batches
  Idents(seu_list.integrated) <- "batch"
  seu_list.integrated[["batch"]] <- Idents(seu_list.integrated)

  # switch to integrated assay. The variable features of this assay are
  # automatically set during IntegrateData
  DefaultAssay(object = seu_list.integrated) <- "integrated"

  # Run the standard workflow for visualization and clustering
  seu_list.integrated <- Seurat::ScaleData(object = seu_list.integrated, verbose = FALSE)
  seu_list.integrated <- seurat_reduce_dimensions(seu_list.integrated, ...)

  return(seu_list.integrated)
}

#' Run Louvain Clustering at Multiple Resolutions
#'
#' @param seu
#' @param resolution
#' @param custom_clust
#' @param reduction
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
seurat_cluster <- function(seu = seu, resolution = 0.6, custom_clust = NULL, reduction = "pca", ...){
  # browser()
  seu <- FindNeighbors(object = seu, dims = 1:10, reduction = reduction)

  if (length(resolution) > 1){
    for (i in resolution){
      # browser()
      seu <- Seurat::FindClusters(object = seu, resolution = i)
    }
  } else if (length(resolution) == 1){
    seu <- Seurat::FindClusters(object = seu, resolution = resolution, ...)
  }

  if (!is.null(custom_clust)){
    seu <- Seurat::StashIdent(object = seu, save.name = "old.ident")
    clusters <- tibble::tibble("Sample_ID" = rownames(seu[[]])) %>%
    tibbl::rownames_to_column("order") %>%
    dplyr::inner_join(custom_clust, by = "Sample_ID") %>%
    dplyr::pull(cluster) %>%
    identity()

    Idents(object = seu) <- clusters
    # browser()

    return(seu)
  }

  return(seu)
}

#' Read in Gene and Transcript Seurat Objects
#'
#' @param proj_dir
#' @param prefix
#'
#' @return
#' @export
#'
#'
#' @examples
load_seurat_path <- function(proj_dir = getwd(), prefix = "unfiltered"){
  # browser()

  seu_path <- paste0(paste0("*", prefix, "_seu.rds"))

  seu_path <- fs::path(proj_dir, "output", "seurat") %>%
    fs::dir_ls() %>%
    fs::path_filter(seu_path) %>%
    identity()
  return(seu_path)
}


#' Load Seurat Files from a signle project path
#'
#' @param proj_dir
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
load_seurat_from_proj <- function(proj_dir, ...){
  seu_file <- load_seurat_path(proj_dir, ...)

  seu_file <- readRDS(seu_file)
}


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
#'
#' @return
#' @export
#'
#'
#' @examples
seurat_integration_pipeline <- function(seus, resolution, suffix = '', ...) {

  corrected_seu <- seurat_batch_correct(seus, ...)

  # cluster merged seurat objects
  corrected_seu <- seurat_cluster(corrected_seu, resolution = resolution, ...)

  corrected_seu <- find_all_markers(corrected_seu)

  # corrected_seu <- save_seurat(corrected_seu, feature = feature, suffix = suffix, ...)

}

#' Dimensional Reduction
#'
#' Run PCA, TSNE and UMAP on a seurat object
#'
#' @param seu
#'
#' @return
#' @export
#'
#' @examples
seurat_reduce_dimensions <- function(seu, reduction = "pca", ...) {

  seu <- Seurat::RunPCA(object = seu, features = Seurat::VariableFeatures(object = seu), do.print = FALSE, ...)
  if (reduction == "harmony"){
    seu <- harmony::RunHarmony(seu, "batch")
  }
  seu <- Seurat::RunTSNE(object = seu, reduction = reduction, dims = 1:30, ...)
  seu <- Seurat::RunUMAP(object = seu, reduction = reduction, dims = 1:30)

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

  seu <- seurat_preprocess(seu, scale = T)

  # PCA
  seu <- seurat_reduce_dimensions(seu, check_duplicates = FALSE, reduction = reduction, ...)

  seu <- seurat_cluster(seu = seu, resolution = resolution, reduction = reduction, ...)

  seu <- find_all_markers(seu, resolution = resolution, reduction = reduction)

  return(seu)
}

#' Give a new project name to a seurat object
#'
#' @param seu
#' @param new_name
#'
#' @return
#' @export
#'
#' @examples
rename_seurat <- function(seu, new_name){
  seu@project.name <- new_name
  return(seu)
}

#' Reset the default assay of a seurat object
#'
#' @param seu
#' @param new_assay
#'
#' @return
#' @export
#'
#' @examples
SetDefaultAssay <- function(seu, new_assay){
  DefaultAssay(seu) <- new_assay
  return(seu)
}



#' Filter a List of Seurat Objects
#'
#' Filter Seurat Objects by custom variable and reset assay to uncorrected "RNA"
#'
#' @param seus
#' @param filter_var
#' @param filter_val
#' @param .drop
#'
#' @return
#' @export
#'
#' @examples
filter_merged_seus <- function(seus, filter_var, filter_val, .drop = F) {
  seus <- purrr::map(seus, ~filter_merged_seu(seu = .x, filter_var = filter_var, filter_val = filter_val, .drop = .drop))
}


#' Filter a Single Seurat Object
#'
#' @param seu
#' @param filter_var
#' @param filter_val
#' @param .drop
#'
#' @return
#' @export
#'
#' @examples
filter_merged_seu <- function(seu, filter_var, filter_val, .drop = .drop) {
  if(.drop){
    mycells <- seu[[]][[filter_var]] == filter_val
  } else {
    mycells <- seu[[]][[filter_var]] == filter_val | is.na(seu[[]][[filter_var]])
  }
  mycells <- colnames(seu)[mycells]
  seu <- seu[,mycells]
  return(seu)
}


#' Reintegrate (filtered) seurat objects
#'
#' 1) split by batch
#' 2) integrate
#' 3) run integration pipeline and save
#'
#' @param seu
#' @param suffix to be appended to file saved in output dir
#' @param reduction to use default is pca
#'
#' @return
#' @export
#'
#' @examples
reintegrate_seu <- function(seu, feature = "gene", suffix = "", reduction = "pca", ...){

  DefaultAssay(seu) <- "RNA"

  seus <- Seurat::SplitObject(seu, split.by = "batch")

  seu <- seurat_integration_pipeline(seus, suffix = suffix, ...)


}



