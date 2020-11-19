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
run_scvelo <- function(seu, loom_path, fit.quantile = 0.05, check_loom = FALSE, ...){
  # browser()

  ldat <- SeuratWrappers::ReadVelocity(file = loom_path)
  bm <- Seurat::as.Seurat(x = ldat)

  bm[["RNA"]] <- bm[["spliced"]]

  # subset bm by seurat object.size
  bm <- bm[,colnames(bm) %in% colnames(seu)]

  # subset seurat object by ldat
  sub_seu <- seu[,colnames(seu) %in% colnames(bm)]

  sub_seu@assays[names(bm@assays)] <- bm@assays
  DefaultAssay(sub_seu) <- "RNA"
  sub_seu@misc$vel <- NULL
  sub_seu@misc[names(sub_seu@misc) == "experiment"] <- NULL

  convert_to_h5ad(sub_seu, file_path = loom_path)

  return(sub_seu)

}

convert_to_h5ad <- function(seu, file_path){
  h5seurat_path <- fs::path_ext_set(file_path, ".h5Seurat")
  SeuratDisk::SaveH5Seurat(seu, filename = h5seurat_path, overwrite = TRUE)
  SeuratDisk::Convert(h5seurat_path, dest = "h5ad", overwrite = TRUE)
}

#' scvelo_assay
#'
#' run scvelo on a gene or transcript level seurat object
#'
#' @param seu a seurat object
#' @param loom_path path to matching loom file
#' @param group.by metadata to color plot
#' @param plot_method plotting method to use from scvelo
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
prep_scvelo <- function(seu, loom_path, plot_method = c("stream", "arrow", "dynamics"), ...){
  # browser()

  h5ad_path <- fs::path_ext_set(loom_path, ".h5ad")

  adata_matches_seu <- function(seu, adata){
    all(adata$obs_names$values %in% colnames(seu))
  }

  if (fs::file_exists(h5ad_path)){
    adata = scvelo$read(fs::path_expand(h5ad_path))

    if(!adata_matches_seu(seu, adata)){
      run_scvelo(seu, loom_path)
    }

  } else {
    seu <- run_scvelo(seu, loom_path)
  }

  adata = scvelo$read(fs::path_expand(h5ad_path))

  scvelo$pp$moments(adata, n_pcs=30L, n_neighbors=30L)

  scvelo$tl$velocity(adata)

  scvelo$tl$velocity_graph(adata)

  if(plot_method == "dynamics"){
    scvelo$tl$recover_dynamics(adata)
    scvelo$tl$latent_time(adata)
  }

  return(adata)

}

#' Plot scvelo on embedding plot
#'
#' @param adata
#' @param group.by
#' @param plot_method
#'
#' @return
#' @export
#'
#' @examples
plot_scvelo <- function(adata, group.by = "batch", plot_method = c("stream", "arrow", "dynamics")){

  num_cols <- length(unique(adata$obs[[group.by]]))

  mycols <- scales::hue_pal()(num_cols)

  if(plot_method == "stream"){
    scvelo$pl$velocity_embedding_stream(adata, basis="umap", palette = mycols, color=group.by)
  } else if(plot_method == "arrow"){
    scvelo$pl$velocity_embedding(adata, basis="umap", palette = mycols, color=group.by, arrow_length=3, arrow_size=2, dpi=120)
  } else if(plot_method == "dynamics"){
    scvelo$pl$scatter(adata, color="latent_time", color_map="gnuplot")
  }

  pyplot$show()

}

scvelo_expression <- function(adata, features = c("RXRG")){
  scvelo$pl$velocity(adata, var_names = features)

  pyplot$show()
}
