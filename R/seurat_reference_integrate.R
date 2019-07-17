


#' Predict One Seurat Object on another merged set
#'
#' @param ref_seu_list
#' @param query_seu
#'
#' @return
#' @export
#'
#' @examples
seurat_reference_integrate <- function(ref_seu_list, query_seu) {
  #browser()
  ## To construct a reference, we will identify ‘anchors’ between the individual datasets. First, we split the combined object into a list, with each dataset as an element.

  ## Prior to finding anchors, we perform standard preprocessing (log-normalization), and identify variable features individually for each. Note that Seurat v3 implements an improved method for variable feature selection based on a variance stabilizing transformation ("vst")

  for (i in 1:length(x = seu_list)) {
  	seu_list[[i]] <- NormalizeData(object = seu_list[[i]], verbose = FALSE)
  	seu_list[[i]] <- FindVariableFeatures(object = seu_list[[i]],
  																							selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  }

  ## Next, we identify anchors using the FindIntegrationAnchors function, which takes a list of Seurat objects as input.

  seu_list.anchors <- Seurat::FindIntegrationAnchors(object.list = seu_list, dims = 1:30, k.filter = 50)


  ## We then pass these anchors to the IntegrateData function, which returns a Seurat object.

  ## The returned object will contain a new Assay, which holds an integrated (or ‘batch-corrected’) expression matrix for all cells, enabling them to be jointly analyzed.


  seu_list.integrated <- IntegrateData(anchorset = seu_list.anchors, dims = 1:30)

  ## switch to integrated assay. The variable features of this assay are
  ## automatically set during IntegrateData
  DefaultAssay(object = seu_list.integrated) <- "integrated"

  ## Run the standard workflow for visualization and clustering
  seu_list.integrated <- ScaleData(object = seu_list.integrated, verbose = FALSE)
  seu_list.integrated <- RunPCA(object = seu_list.integrated, npcs = 30, verbose = FALSE)
  seu_list.integrated <- RunUMAP(object = seu_list.integrated, reduction = "pca",
  																dims = 1:30)
  seu_list.integrated <- RunTSNE(object = seu_list.integrated, reduction = "pca", dims = 1:30)

  return(seu_list.integrated)
}


## ------------------------------------------------------------------------
## find markers for every cluster compared to all remaining cells, report only the positive ones
#' Find Cell Type Markers in a Seurat Object
#'
#' @param seu A seurat object
#'
#' @return
#' @export
#'
#' @examples
seurat_find_markers <- function(seu, num_features){
  markers <- Seurat::FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
  markers %>%
    group_by(cluster) %>%
    top_n(n = 2, wt = avg_logFC)

  cluster_markers <- reference.markers$gene %>%
    group_by(cluster) %>%
    top_n(n = num_features, wt = avg_logFC) %>%
    dplyr::pull(gene)

  DoHeatmap(reference_integrated[[1]], features = cluster_markers)
  return(cluster_markers)
}

