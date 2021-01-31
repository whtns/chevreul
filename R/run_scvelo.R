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

  bm[["gene"]] <- bm[["spliced"]]

  # subset bm by seurat object.size
  bm <- bm[,colnames(bm) %in% colnames(seu)]

  # subset seurat object by ldat
  sub_seu <- seu[,colnames(seu) %in% colnames(bm)]

  sub_seu@assays[names(bm@assays)] <- bm@assays
  DefaultAssay(sub_seu) <- "gene"
  sub_seu@misc$vel <- NULL
  sub_seu@misc[names(sub_seu@misc) == "experiment"] <- NULL

  convert_to_h5ad(sub_seu, file_path = loom_path)

  return(sub_seu)

}

#' convert a seurat object to an on-disk anndata object
#'
#' @param seu
#' @param file_path
#'
#' @return
#' @export
#'
#' @examples
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
prep_scvelo <- function(seu, loom_path, velocity_mode = c("deterministic", "stochastic", "dynamical"), ...){
  browser()

  h5ad_path <- fs::path_ext_set(loom_path, ".h5ad")

    scvelo <- reticulate::import("scvelo")

    adata_matches_seu <- function(seu, adata){
      identical(sort(adata$obs_names$values), sort(colnames(seu)))
    }

    if (fs::file_exists(h5ad_path)){
      adata = scvelo$read(fs::path_expand(h5ad_path))

      if(!adata_matches_seu(seu, adata)){
        seu <- run_scvelo(seu, loom_path)
      }

    } else {
      seu <- run_scvelo(seu, loom_path)
    }

    adata = scvelo$read(fs::path_expand(h5ad_path))
    # reticulate::source_python("scripts/rename_raw.py")
    # adata$raw$var$rename(columns = list('_index' = 'symbol'), inplace = True)

    scvelo$pp$moments(adata, n_pcs=30L, n_neighbors=30L)

    if(velocity_mode == "dynamical"){
      if(!"recover_dynamics" %in% adata$uns_keys()){
        scvelo$tl$recover_dynamics(adata)
        reticulate::py_del_attr(adata, "raw")
        adata$write_h5ad(h5ad_path)
      }
      scvelo$tl$latent_time(adata)
    }

    scvelo$tl$velocity(adata, mode = velocity_mode)

    scvelo$tl$velocity_graph(adata)


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
    scvelo$pl$velocity_embedding_stream(adata, basis="umap", palette = mycols, color=group.by, dpi = 200, figsize=c(20,12))
  } else if(plot_method == "arrow"){
    scvelo$pl$velocity_embedding(adata, basis="umap", palette = mycols, color=group.by, arrow_length=3, arrow_size=2, dpi=200, figsize=c(20,12))
  } else if(plot_method == "dynamics"){
    scvelo$pl$scatter(adata, color="latent_time", color_map="gnuplot", figsize=c(20,12), dpi = 200)
  }

  # pyplot$show()

}

scvelo_expression <- function(adata, features = c("RXRG")){
  scvelo$pl$velocity(adata, var_names = features, figsize = c(10,10), dpi = 200)

  # pyplot$show()
}
