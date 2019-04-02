
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
  	seu_list[[i]] <- NormalizeData(object = seu_list[[i]], verbose = FALSE)
  	seu_list[[i]] <- FindVariableFeatures(object = seu_list[[i]],
  																							selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  }

  # Next, we identify anchors using the FindIntegrationAnchors function, which takes a list of Seurat objects as input.

  seu_list.anchors <- FindIntegrationAnchors(object.list = seu_list, dims = 1:30, k.filter = 50)


  # We then pass these anchors to the IntegrateData function, which returns a Seurat object.

  # The returned object will contain a new Assay, which holds an integrated (or ‘batch-corrected’) expression matrix for all cells, enabling them to be jointly analyzed.


  seu_list.integrated <- IntegrateData(anchorset = seu_list.anchors, dims = 1:30)

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
pull_gene_trx_seus <- function(proj_dir){
  # browser()
  gene_trx_sces <- fs::path(proj_dir, "output", "sce") %>%
    dir_ls() %>%
    path_filter("*_seu.rds") %>%
    set_names(c("gene", "transcript")) %>%
    identity()
  return(gene_trx_sces)
}


#' Stash Batches in Integrated Seurat Object
#'
#' @param corrected_seus
#'
#' @return
#' @export
#'
#' @examples
stash_batches <- function(corrected_seus) {

  # stash batches in object metadata
  set_idents <- as_mapper(~ StashIdent(., save.name = "batch"))

  corrected_seus <- map(corrected_seus, set_idents)

}

