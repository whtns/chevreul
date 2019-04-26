#' Run RNA Velocity
#'
#' @param seu
#' @param loom_path
#'
#' @return
#' @export
#'
#' @examples
run_rna_velocity <- function(seu, loom_path) {

  ldat <- velocyto.R::read.loom.matrices(loom_path)

  ## load objects------------------------------------------------------------------------

  # subset ldat by seurat object.size
  ldat <- purrr::map(ldat, ~.x[,colnames(.x) %in% colnames(seu)])

  # subset seurat object by ldat
  seu <- seu[,colnames(seu) %in% colnames(ldat[[1]])]

  ## grab cell colors ------------------------------------------------------------------------

  clusters_meta <- seu[[]] %>%
    dplyr::select(starts_with("clusters_")) %>%
    dplyr::select(paste0("clusters_", seq(0.2, 2.0, by = 0.2))) %>%
    identity()

  colvecs <- purrr::map_int(clusters_meta, ~length(unique(.x))) %>%
    map(~scales::hue_pal()(.x))

  ## format cell colors------------------------------------------------------------------------

  df_col_to_nm_vec <- function(df_col, colvec, df_names){
    myvec <- unlist(df_col) %>%
      set_names(df_names) %>%
      forcats::lvls_revalue(colvec)
  }

  cell.colors <- purrr::map2(clusters_meta, colvecs, df_col_to_nm_vec, rownames(clusters_meta))


  ## organize loom data ------------------------------------------------------------------------
  # exonic read (spliced) expression matrix
  emat <- ldat$spliced;
  # intronic read (unspliced) expression matrix
  nmat <- ldat$unspliced
  # spanning read (intron+exon) expression matrix
  smat <- ldat$spanning;

  # pre-allocate a list for velocity estimates
  seu@misc$vel <- vector("list", length(cell.colors))

  for (i in names(cell.colors)){
    # filter expression matrices based on some minimum max-cluster averages
    emat <- filter.genes.by.cluster.expression(emat,cell.colors[[i]],min.max.cluster.average = 0.5)
    nmat <- filter.genes.by.cluster.expression(nmat,cell.colors[[i]],min.max.cluster.average = 1)
    smat <- filter.genes.by.cluster.expression(smat,cell.colors[[i]],min.max.cluster.average = 0.5)

    # ## calculate velocity------------------------------------------------------------------------
    fit.quantile <- 0.05;
    rvel.qf <- gene.relative.velocity.estimates(emat, nmat, deltaT=1, kCells = 5, fit.quantile = fit.quantile)

    # assign velocity to seurat object in slot @misc$vel$cluster_name

    seu@misc$vel[[i]] <- rvel.qf
  }

  return(seu)

}
