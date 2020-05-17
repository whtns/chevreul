
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
velocyto_assay <- function(seu, loom_path, fit.quantile = 0.05, ...){

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

  ## organzie loom data ------------------------------------------------------------------------
  # exonic read (spliced) expression matrix
  emat <- ldat$spliced;
  # intronic read (unspliced) expression matrix
  nmat <- ldat$unspliced
  # spanning read (intron+exon) expression matrix
  smat <- ldat$spanning;
  # filter expression matrices based on some minimum max-cluster averages
  emat <- velocyto.R::filter.genes.by.cluster.expression(emat,cell.colors,min.max.cluster.average = 0.5)
  nmat <- velocyto.R::filter.genes.by.cluster.expression(nmat,cell.colors,min.max.cluster.average = 1)
  smat <- velocyto.R::filter.genes.by.cluster.expression(smat,cell.colors,min.max.cluster.average = 0.5)
  # look at the resulting gene set
  length(intersect(rownames(emat),rownames(nmat)))


  ## ------------------------------------------------------------------------
  # and if we use spanning reads (smat)
  length(intersect(intersect(rownames(emat),rownames(nmat)),rownames(smat)))

  # plot_spliced_mag <- function(ldat) {
  #   hist(log10(rowSums(ldat$spliced)+1),col='wheat',xlab='log10[ number of reads + 1]',main='number of reads per gene')
  # }
  #
  # pdf(fs::path(batch_dir, "reads_per_gene.pdf"))
  # plot_spliced_mag(ldat)
  # dev.off()

  # ## calculate velocity------------------------------------------------------------------------
  rvel.qf <- velocyto.R::gene.relative.velocity.estimates(emat, nmat, deltaT=1, kCells = 5, fit.quantile = fit.quantile)

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
