#' scvelo_assay
#'
#' run scvelo on a gene or transcript level seurat object
#'
#' @param seu a seurat object
#' @param loom_path path to matching loom file
#' @param fit.quantile
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
run_scvelo <- function(seu, loom_path, group.by = "batch", fit.quantile = 0.05, check_loom = FALSE, ...){
  # browser()

  h5seurat_path <- stringr::str_replace(loom_path, ".loom", ".h5Seurat")

  h5ad_path <- stringr::str_replace(loom_path, ".loom", ".h5ad")

  check_loom_dim <- function(seu){
    # browser()

    !is.null(seu@misc$vel)

    all(rownames(seu@misc$vel$cellKNN) %in% colnames(seu))
  }

  if (check_loom){
    if (check_loom_dim(seu)) return(seu)
  }

  ldat <- SeuratWrappers::ReadVelocity(file = loom_path)
  bm <- Seurat::as.Seurat(x = ldat)

  bm[["RNA"]] <- bm[["spliced"]]
  # bm <- seurat_preprocess(bm)
  # bm <- seurat_reduce_dimensions(bm)
  # bm <- FindNeighbors(bm, dims = 1:30)
  # bm <- FindClusters(bm)

  # subset bm by seurat object.size
  bm <- bm[,colnames(bm) %in% colnames(seu)]

  # subset seurat object by ldat
  sub_seu <- seu[,colnames(seu) %in% colnames(bm)]

  sub_seu@assays[names(bm@assays)] <- bm@assays
  DefaultAssay(sub_seu) <- "RNA"
  sub_seu@misc <- bm@misc

  SeuratDisk::SaveH5Seurat(sub_seu, filename = h5seurat_path, overwrite = TRUE)
  SeuratDisk::Convert(h5seurat_path, dest = "h5ad", overwrite = TRUE)

  # ## calculate velocity------------------------------------------------------------------------

  adata = scvelo$read(fs::path_expand(h5ad_path))

  scvelo$pp$moments(adata, n_pcs=30L, n_neighbors=30L)

  scvelo$tl$velocity(adata)

  scvelo$tl$velocity_graph(adata)

  num_cols <- length(unique(adata$obs[[group.by]]))

  mycols <- scales::hue_pal()(num_cols)
  scvelo$pl$velocity_embedding_stream(adata, basis="umap", palette = mycols, color=group.by)

  pyplot$show()

}
