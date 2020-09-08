
#' run velocyto on a gene or transcript level seurat object
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
velocyto_assay <- function(seu, loom_path, fit.quantile = 0.05, check_loom = FALSE, ...){
  # browser()
  check_loom_dim <- function(seu){
    # browser()

    !is.null(seu@misc$vel)

    all(rownames(seu@misc$vel$cellKNN) %in% colnames(seu))
  }

  if (check_loom){
    if (check_loom_dim(seu)) return(seu)
  }

  ldat <- velocyto.R::read.loom.matrices(loom_path)

  # subset ldat by seurat object.size
  ldat <- purrr::map(ldat, ~.x[,colnames(.x) %in% colnames(seu)])

  # subset seurat object by ldat
  sub_seu <- seu[,colnames(seu) %in% colnames(ldat[[1]])]

  ## grab cell colors ------------------------------------------------------------------------

  p <- DimPlot(sub_seu, ...)
  col_vec <- unique(ggplot2::ggplot_build(p)$data[[1]]) %>%
    dplyr::arrange(group) %>%
    dplyr::pull(colour) %>%
    unique()

  ## format cell colors------------------------------------------------------------------------
  cell.colors <- Idents(sub_seu)

  # levels(cell.colors) <- as.character(as.numeric(levels(cell.colors)) + 1)
  levels(cell.colors) <- col_vec

  # filter expression matrices based on some minimum max-cluster averages
  ldat$spliced <- velocyto.R::filter.genes.by.cluster.expression(ldat$spliced,cell.colors,min.max.cluster.average = 0.5)
  ldat$unspliced <- velocyto.R::filter.genes.by.cluster.expression(ldat$unspliced,cell.colors,min.max.cluster.average = 1)
  ldat$spanning <- velocyto.R::filter.genes.by.cluster.expression(ldat$spanning,cell.colors,min.max.cluster.average = 0.5)
  # look at the resulting gene set
  length(intersect(rownames(ldat$spliced),rownames(ldat$unspliced)))


  ## ------------------------------------------------------------------------
  # and if we use spanning reads (ldat$spanning)
  length(intersect(intersect(rownames(ldat$spliced),rownames(ldat$unspliced)),rownames(ldat$spanning)))

  # plot_spliced_mag <- function(ldat) {
  #   hist(log10(rowSums(ldat$spliced)+1),col='wheat',xlab='log10[ number of reads + 1]',main='number of reads per gene')
  # }

  # ## calculate velocity------------------------------------------------------------------------
  rvel.qf <- velocyto.R::gene.relative.velocity.estimates(ldat$spliced, ldat$unspliced, deltaT=1, kCells = 5, fit.quantile = fit.quantile)

  cc_umap <- plot_velocity_arrows(seu, rvel.qf, reduction = "umap", cell.colors = cell.colors, plot_format = "arrow")

  # save graphical info to seurat object in slot @misc$cc
  seu@misc$cc <- cc_umap$cc

  # save velocity to seurat object in slot @misc$vel
  seu@misc$vel <- rvel.qf

  return(seu)
}

#' run velocyto on a seurat object
#'
#' @param seu
#' @param loom_path
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
velocyto_seu <- function(seu, loom_path, fit.quantile = 0.05, ...){
  seu <- purrr::map(seu, velocyto_assay, loom_path, fit.quantile = fit.quantile)
  return(seu)
}

#' Run RNA Velocity starting with only a loom File
#'
#' @param loom_path
#'
#' @return
#' @export
#'
#' @examples
velocyto_seurat_from_loom <- function(loom_path) {
  ldat <- ReadVelocity(file = loom_path)
  bm <- SeuratWrappers::as.Seurat(x = ldat)
  bm <- Seurat::SCTransform(object = bm, assay = "spliced")
  bm <- Seurat::RunPCA(object = bm, verbose = FALSE)
  bm <- Seurat::FindNeighbors(object = bm, dims = 1:20)
  bm <- Seurat::FindClusters(object = bm, algorithm = 3)
  bm <- Seurat::RunUMAP(object = bm, dims = 1:20)
  bm <-
    SeuratWrappers::RunVelocity(
      object = bm,
      deltaT = 1,
      kCells = 25,
      fit.quantile = 0.02
    )
  ident.colors <- (scales::hue_pal())(n = length(x = levels(x = bm)))
  names(x = ident.colors) <- levels(x = bm)
  cell.colors <- ident.colors[Idents(object = bm)]
  names(x = cell.colors) <- colnames(x = bm)
  velocyto.R::show.velocity.on.embedding.cor(
    emb = Embeddings(object = bm, reduction = "umap"),
    vel = Tool(object = bm,
               slot = "RunVelocity"),
    n = 200,
    scale = "sqrt",
    cell.colors = ac(x = cell.colors, alpha = 0.5),
    cex = 0.8,
    arrow.scale = 3,
    show.grid.flow = TRUE,
    min.grid.cell.mass = 0.5,
    grid.n = 40,
    arrow.lwd = 1,
    do.par = FALSE,
    cell.border.alpha = 0.1
  )
}
