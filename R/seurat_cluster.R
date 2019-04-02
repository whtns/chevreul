

#' Preprocess Seurat Object
#'
#' @param seu
#' @param scale
#'
#' @return
#' @export
#'
#' @examples
seurat_preprocess <- function(seu, scale=TRUE){
    # Normalize data
    seu <- NormalizeData(object = seu, normalization.method = "LogNormalize",
        scale.factor = 10000)

    # Filter out only variable genes
    seu <- FindVariableFeatures(object = seu, mean.function = ExpMean, dispersion.function = LogVMR,
        x.low.cutoff = 0.0125, x.high.cutoff = 3, y.cutoff = 0.5, do.plot = FALSE)

    FeatureScatter(object = seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")

    # Regress out unwanted sources of variation
    if (scale){
        seu <- ScaleData(object = seu, features = rownames(x = seu))

    }

    # seu <- CellCycleScoring(object = seu, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)
    #
    # # view cell cycle scores and phase assignments
    # head(x = seu@meta.data)
    #
    return(seu)
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
seurat_cluster2 <- function(seu, resolution = 0.6, custom_clust = NULL) {

  if (length(resolution) > 1){
    for (i in resolution){
      seu <- FindClusters(object = seu, resolution = i)
      seu <- StashIdent(object = seu, save.name = paste0("clusters_", i))
    }
  }

  if (!is.null(custom_clust)){
    seu <- StashIdent(object = seu, save.name = "old.ident")
    clusters <- tibble("Sample_ID" = rownames(seu[[]])) %>%
    rownames_to_column("order") %>%
    dplyr::left_join(custom_clust, by = "Sample_ID") %>%
    pull(cluster) %>%
    identity()

    Idents(object = seu) <- clusters
    #

    return(seu)
  }

  return(seu)
}

#' Run Seurat Pipeline
#'
#' @param seu
#' @param resolution
#'
#' @return
#' @export
#'
#' @examples
seurat_pipeline <- function(seu, resolution=0.6){

  seu <- seurat_preprocess(seu)

  # PCA
  seu <- RunPCA(object = seu, features = VariableFeatures(object = seu), do.print = FALSE)
  seu <- FindNeighbors(object = seu, dims = 1:10)
  seu <- seurat_cluster2(seu, resolution = resolution)
  seu <- RunTSNE(object = seu, reduction = "pca", dims = 1:30)
  seu <- RunUMAP(object = seu, reduction = "pca", dims = 1:30)
  seu <- ProjectDim(object = seu)

  return(seu)
}


seurat_batch_correct2 <- function(merge_sce) {
  # old function need to investigate why these the sequence of normalization alters output
  # To construct a reference, we will identify ‘anchors’ between the individual datasets. First, we split the combined object into a list, with each dataset as an element.

  merge_sce.list <- SplitObject(object = merge_sce, split.by = "batch")

  # Prior to finding anchors, we perform standard preprocessing (log-normalization), and identify variable features individually for each. Note that Seurat v3 implements an improved method for variable feature selection based on a variance stabilizing transformation ("vst")

  for (i in 1:length(x = merge_sce.list)) {
  	merge_sce.list[[i]] <- NormalizeData(object = merge_sce.list[[i]], verbose = FALSE)
  	merge_sce.list[[i]] <- FindVariableFeatures(object = merge_sce.list[[i]],
  																							selection.method = "vst", nfeatures = 2000, verbose = FALSE)
  }

  # Next, we identify anchors using the FindIntegrationAnchors function, which takes a list of Seurat objects as input.

  merge_sce.anchors <- FindIntegrationAnchors(object.list = merge_sce.list, dims = 1:30, k.filter = 150)


  # We then pass these anchors to the IntegrateData function, which returns a Seurat object.

  # The returned object will contain a new Assay, which holds an integrated (or ‘batch-corrected’) expression matrix for all cells, enabling them to be jointly analyzed.


  merge_sce.integrated <- IntegrateData(anchorset = merge_sce.anchors, dims = 1:30)

  # switch to integrated assay. The variable features of this assay are
  # automatically set during IntegrateData
  DefaultAssay(object = merge_sce.integrated) <- "integrated"

  # Run the standard workflow for visualization and clustering
  merge_sce.integrated <- ScaleData(object = merge_sce.integrated, verbose = FALSE)
  merge_sce.integrated <- RunPCA(object = merge_sce.integrated, npcs = 30, verbose = FALSE)
  merge_sce.integrated <- RunUMAP(object = merge_sce.integrated, reduction = "pca",
  																dims = 1:30)
  merge_sce.integrated <- RunTSNE(object = merge_sce.integrated, reduction = "pca", dims = 1:30)
}

