#' scvelo_assay
#'
#' run scvelo on a gene or transcript level object
#'
#' @param object a object
#' @param loom_path path to matching loom file
#' @param assay gene
#' @param fit.quantile how to fit velocity
#' @param check_loom FALSE
#'
#' @return a single cell object with RNA velocity calculated
#' @export
#' @importFrom LoomExperiment import
#' @examples \donttest{merge_loom(human_gene_transcript_sce, "my.loom")}
merge_loom <- function(object, loom_path, assay = "gene", fit.quantile = 0.05, check_loom = FALSE) {

  loom_object <- import(loom_path, type="SingleCellLoomExperiment", colnames_attr = "CellID", rownames_attr = "Gene")

  loom_object <- loom_object[,colnames(loom_object) %in% colnames(object)]

  object <- object[,colnames(object) %in% colnames(loom_object)]

  loom_object <- loom_object[,colnames(object)]

  altExp(object, "velocity") <- loom_object

  return(object)
}

#' plot scvelo
#'
#' run scvelo on a gene or transcript level object
#'
#' @param object a object
#' @param loom_path path to matching loom file
#' @param velocity_mode deterministic, stochastis, or dynamical
#' @param ... extra args passed to run_scvelo
#'
#' @return a single cell object with velocity calculated
#' @export
#' @examples plot_scvelo(human_gene_transcript_sce, embedding = "UMAP", color_by = "velocity_pseudotime")
plot_scvelo <- function(object, mode = c("steady_state", "deterministic", "stochastic", "dynamical"), embedding = c("UMAP", "PCA", "TSNE"), ...) {

  embedding <- match.arg(embedding)
  original_exp <- mainExpName(object)

  object <- swapAltExp(object, "velocity")

  object <- logNormCounts(object, assay.type=1)

  dec <- modelGeneVar(object)
  top.hvgs <- getTopHVGs(dec, n=2000)

  velo.out <- scvelo(object, subset.row=top.hvgs, assay.X="spliced", mode = mode)

  object$velocity_pseudotime <- velo.out$velocity_pseudotime

  object <- swapAltExp(object, original_exp)

  embedded <- embedVelocity(reducedDim(object, embedding), velo.out)
  grid.df <- gridVectors(reducedDim(object, embedding), embedded)

  colnames(grid.df) <- c("x", "y", "xend", "yend")

  umap_plot <- plotReducedDim(object, dimred = embedding, ...) +
    geom_segment(data=grid.df, mapping=aes(x=x, y=y,
                                           xend=xend, yend=yend), arrow=arrow(length=unit(0.05, "inches")))

  return(umap_plot)
}

#' scvelo expression
#'
#' run scvelo on a gene or transcript level object
#'
#' @param object a object
#' @param loom_path path to matching loom file
#' @param velocity_mode deterministic, stochastis, or dynamical
#' @param ... extra args passed to run_scvelo
#'
#' @return a single cell object with velocity calculated
#' @export
scvelo_expression <- function(object, mode = c("steady_state", "deterministic", "stochastic", "dynamical"), ...) {

  object <- swapAltExp(object, "velocity")

  object <- logNormCounts(object, assay.type=1)

  dec <- modelGeneVar(object)
  top.hvgs <- getTopHVGs(dec, n=2000)

  velo.out <- scvelo(object, subset.row=top.hvgs, assay.X="spliced", mode = mode)

  object$velocity_pseudotime <- velo.out$velocity_pseudotime

  object <- swapAltExp(object, "gene")

  embedded <- embedVelocity(reducedDim(object, "UMAP"), velo.out)
  grid.df <- gridVectors(reducedDim(object, "UMAP"), embedded)

  umap_plot <- plotUMAP(object, ...) +
    geom_segment(data=grid.df, mapping=aes(x=start.umap_1, y=start.umap_2,
                                           xend=end.umap_1, yend=end.umap_2), arrow=arrow(length=unit(0.05, "inches")))

  return(umap_plot)
}

# ## Tabulate important genes
# scvelo_genes <- function(adata){
#   scvelo$tl$rank_velocity_genes(adata, groupby='clusters_0.2', min_corr=.3)
#   df = scvelo$DataFrame(adata.uns['rank_velocity_genes']['names'])
#   df.head()
# }
