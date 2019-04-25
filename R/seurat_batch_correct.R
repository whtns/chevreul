
## ----load-funcs----------------------------------------------------------

#' Batch Correct Multiple Seurat Objects
#'
#' @param seu_list
#'
#' @return
#' @export
#'
#' @examples
seurat_batch_correct <- function(seu_list) {
  #browser()
  # To construct a reference, we will identify ‘anchors’ between the individual datasets. First, we split the combined object into a list, with each dataset as an element.

# Prior to finding anchors, we perform standard preprocessing (log-normalization), and identify variable features individually for each. Note that Seurat v3 implements an improved method for variable feature selection based on a variance stabilizing transformation ("vst")

  for (i in 1:length(x = seu_list)) {
  	seu_list[[i]] <- seurat_preprocess(seu_list[[i]], scale = F)
  }

  # Next, we identify anchors using the FindIntegrationAnchors function, which takes a list of Seurat objects as input.

  seu_list.anchors <- FindIntegrationAnchors(object.list = seu_list, dims = 1:30, k.filter = 50)


  # We then pass these anchors to the IntegrateData function, which returns a Seurat object.

  # The returned object will contain a new Assay, which holds an integrated (or ‘batch-corrected’) expression matrix for all cells, enabling them to be jointly analyzed.

  seu_list.integrated <- IntegrateData(anchorset = seu_list.anchors, dims = 1:30)

  # #stash batches
  # Idents(seu_list.integrated) <- "batch"
  # seu_list.integrated <- StashIdent(seu_list.integrated, save.name = "batch")

  # switch to integrated assay. The variable features of this assay are
  # automatically set during IntegrateData
  DefaultAssay(object = seu_list.integrated) <- "integrated"

  # Run the standard workflow for visualization and clustering
  seu_list.integrated <- ScaleData(object = seu_list.integrated, verbose = FALSE)
  seu_list.integrated <- RunPCA(object = seu_list.integrated, npcs = 30, verbose = FALSE)
  seu_list.integrated <- RunUMAP(object = seu_list.integrated, reduction = "pca",
  																dims = 1:30)
  seu_list.integrated <- RunTSNE(object = seu_list.integrated, reduction = "pca", dims = 1:30)

  return(seu_list.integrated)
}

#' Run Louvain Clustering at Multiple Resolutions
#'
#' @param seu
#' @param resolution
#' @param custom_clust
#'
#' @return
#' @export
#'
#' @examples
seurat_cluster <- function(seu, resolution = 0.6, custom_clust = NULL) {
  # browser()
  seu <- FindNeighbors(object = seu, dims = 1:10)
  seu <- FindClusters(object = seu, resolution = resolution)

  if (length(resolution) > 1){
    for (i in resolution){
      # browser()
      seu <- FindClusters(object = seu, resolution = i)
      seu <- StashIdent(object = seu, save.name = paste0("clusters_", i))
    }
  }

  if (!is.null(custom_clust)){
    seu <- StashIdent(object = seu, save.name = "old.ident")
    clusters <- tibble("Sample_ID" = rownames(seu[[]])) %>%
    rownames_to_column("order") %>%
    dplyr::inner_join(custom_clust, by = "Sample_ID") %>%
    pull(cluster) %>%
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
#'
#' @return
#' @export
#'
#' @examples
pull_gene_trx_seus <- function(proj_dir, suffix = ""){
  # browser()

  if(suffix != ""){
    suffix = paste0("_", suffix)
  }

  seu_path <- paste0("*", "_seu", suffix, ".rds")

  gene_trx_sces <- fs::path(proj_dir, "output", "sce") %>%
    dir_ls() %>%
    path_filter(seu_path) %>%
    set_names(c("gene", "transcript")) %>%
    identity()
  return(gene_trx_sces)
}


#' Load Seurat files from a vector of project paths
#'
#' @param proj_dirs
#'
#' @return
#' @export
#'
#' @examples
load_seurat_from_rds <- function(proj_dirs){
  seu_files <- map(proj_dirs, seuratTools::pull_gene_trx_seus)
  names(seu_files) <- gsub("_proj", "", basename(proj_dirs))

  seu_files <- purrr::transpose(seu_files)

  seu_files <- purrr::map(seu_files, ~map(.x, readRDS))
}

#' Run Seurat Integration
#'
#' run batch correction, followed by:
#' 1) stashing of batches in metadata 'batch'
#' 2) clustering with resolution 0.2 to 2.0 in increments of 0.2
#' 3) saving to <proj_dir>/output/sce/<feature>_seu_<suffix>.rds
#'
#' @param seus
#' @param suffix
#'
#' @return
#' @export
#'
#' @examples
seurat_integration_pipeline <- function(seus, res_low = 0.2, res_hi = 2.0, suffix = '') {
  corrected_seus <- purrr::map(seus, seuratTools::seurat_batch_correct)

  # cluster merged seurat objects
  corrected_seus <- map(corrected_seus, seuratTools::seurat_cluster, resolution = seq(res_low, res_hi, by = 0.2))

  corrected_seus <- map(corrected_seus, find_all_markers)

  corrected_seus <- purrr::imap(corrected_seus, save_seurat, suffix = suffix)
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
seurat_pipeline <- function(seu, resolution=0.6){

  seu <- seurat_preprocess(seu, scale = T)

  # PCA
  seu <- RunPCA(object = seu, features = VariableFeatures(object = seu), do.print = FALSE)
  seu <- seurat_cluster(seu, resolution = resolution)
  seu <- RunTSNE(object = seu, reduction = "pca", dims = 1:30)
  seu <- RunUMAP(object = seu, reduction = "pca", dims = 1:30)

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


#' Filter Seurat Objects
#'
#' Filter Seurat Objects by custom variable and reset assay to uncorrected "RNA"
#'
#' @param seus
#' @param filter_var
#' @param filter_val
#' @param tag
#'
#' @return
#' @export
#'
#' @examples
filter_merged_seus <- function(seus, filter_var, filter_val, tag = "") {
  # browser()

  if(tag != ""){
    tag = paste0(tag, "_")
  }

  seus <- map(seus, ~.x[, .x[[filter_var]] == filter_val])

  new_names <- paste0(tag, c("gene", "transcript"))

  seus <- map2(seus, new_names, rename_seurat) %>%
    set_names(new_names)

  seus <- map(seus, SetDefaultAssay, "RNA")

}

#' Reintegrate (filtered) seurat objects
#'
#' 1) split by batch
#' 2) integrate
#' 3) run integration pipeline and save
#'
#' @param seus
#'
#' @return
#' @export
#'
#' @examples
reintegrate_seus <- function(seus, suffix = "", ...){
  seus <- map(seus, SplitObject, split.by = "batch")

  seus <- seurat_integration_pipeline(seus, suffix = suffix)


}

#' Add Read Count Category
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
add_rc_cat_meta_to_seu <- function(seu, thresh = 1e5){
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

  excluded_cells <- map2(excluded_cells, names(excluded_cells), ~rep(.y, length(.x))) %>%
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


