#' Compare two sets of cells with scde package
#'
#' @param data.use
#' @param cells.1
#' @param cells.2
#' @param verbose
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
scdeTest <- function(data.use, cells.1, cells.2, verbose = TRUE, ...)
{
  if (!Seurat:::PackageCheck("scde", error = FALSE)) {
    stop("Please install scde - learn more at https://bioconductor.org/packages/release/bioc/html/scde.html")
  }
  Seurat:::CheckDots(..., fxns = "DESeq2::results")
  group.info <- data.frame(row.names = c(cells.1, cells.2))
  group.info[cells.1, "group"] <- "Group1"
  group.info[cells.2, "group"] <- "Group2"
  group.info[, "group"] <- factor(x = group.info[, "group"])
  group.info$wellKey <- rownames(x = group.info)

  data.use <- as.matrix(data.use)
  data.use <-apply(data.use,2,function(x) {storage.mode(x) <- 'integer'; x})


  data.use <- scde::clean.counts(data.use, min.lib.size = 1000,
                           min.reads = 10, min.detected = 10)

  o.ifm <- scde::scde.error.models(counts = data.use, groups = group.info$group,
                                   n.cores = 1, threshold.segmentation = TRUE, save.crossfit.plots = FALSE,
                                   save.model.plots = FALSE, verbose = 1)
  valid.cells <- o.ifm$corr.a > 0
  table(valid.cells)
  o.ifm <- o.ifm[valid.cells, ]
  o.prior <- scde::scde.expression.prior(models = o.ifm, counts = data.use,
                                         length.out = 400, show.plot = FALSE)

  ediff <- scde::scde.expression.difference(o.ifm, data.use, o.prior,
                                            groups = group.info$group, n.randomizations = 100, n.cores = 1,
                                            verbose = 1)

  to.return <- data.frame(p_val = ediff$pvalue, row.names = rownames(ediff))
  return(to.return)

}


