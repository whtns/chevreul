#' Convert Seurat V3 Object to Monocle2 CellDataSet
#'
#' @param seu
#'
#' @return
#' @export
#'
#' @examples
convert_seuv3_to_monoclev2 <- function(seu){
  #Load Seurat object


  #Extract data, phenotype data, and feature data from the SeuratObject
  data <- as(as.matrix(seu@assays$RNA@data), 'sparseMatrix')

  pd <- new('AnnotatedDataFrame', data = seu@meta.data)

  fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
  fd <- new('AnnotatedDataFrame', data = fData)

  #Construct monocle cds
  monocle_cds <- newCellDataSet(data,
                                phenoData = pd,
                                featureData = fd,
                                lowerDetectionLimit = 0.5,
                                expressionFamily = negbinomial.size())

  monocle_cds <- estimateSizeFactors(monocle_cds)
  monocle_cds <- estimateDispersions(monocle_cds)

  # filter by gene expression (default monocle settings)
  monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)

  expressed_genes <- row.names(subset(featureData(monocle_cds)@data,
                                      num_cells_expressed >= 10))

  # look at distribution of mRNA totals across cells
  phenoData(monocle_cds)$Total_mRNAs <- Matrix::colSums(Biobase::exprs(monocle_cds))

  monocle_cds <- monocle_cds[,phenoData(monocle_cds)$Total_mRNAs < 1e6]

  upper_bound <- 10^(mean(log10(phenoData(monocle_cds)$Total_mRNAs)) +
                       2*sd(log10(phenoData(monocle_cds)$Total_mRNAs)))
  lower_bound <- 10^(mean(log10(phenoData(monocle_cds)$Total_mRNAs)) -
                       2*sd(log10(phenoData(monocle_cds)$Total_mRNAs)))

  # remove cells outside safe range of plot ---------------------------------

  monocle_cds <- monocle_cds[,phenoData(monocle_cds)$Total_mRNAs > lower_bound &
                               phenoData(monocle_cds)$Total_mRNAs < upper_bound]
  monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)

  # verify lognormal distribution of expression values ----------------------

  # Log-transform each value in the expression matrix.
  L <- log(Biobase::exprs(monocle_cds[expressed_genes,]))

  # Standardize each gene, so that they are all on the same scale,
  # Then melt the data with plyr so we can plot it easily
  # melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))

  # cluster cells without marker genes --------------------------------------

  # develop list of genes used for clustering
  disp_table <- dispersionTable(monocle_cds)
  unsup_clustering_genes <- subset(disp_table, mean_expression >= 0.1)
  monocle_cds <- setOrderingFilter(monocle_cds, unsup_clustering_genes$gene_id)

  # cluster cells by tsne
  monocle_cds_red <- reduceDimension(monocle_cds, max_components = 2, num_dim = 6,
                                     reduction_method = 'tSNE', verbose = T)
  monocle_cds_red <- clusterCells(monocle_cds_red, num_clusters = 2)

}

#' Run BEAM from Monocle2
#'
#' http://cole-trapnell-lab.github.io/monocle-release/docs/#differential-expression-analysis
#'
#' @param HSMM
#' @param branch_point
#' @param branches
#' @param pt_param
#' @param pt_paramval
#' @param colorval
#' @param top_genes
#'
#' @return
#' @export
#'
#' @examples
run_BEAM <- function(HSMM, branch_point, branches, pt_param, pt_paramval, colorval, top_genes){

  if(grepl(",", colorval)){
    colorval <- unlist(strsplit(colorval, ","))
  }

  # check if branches present in HSMM

  root_branch <- GM_state(HSMM, pt_param, pt_paramval, branches)
  branch_states = branches[which(branches != root_branch)]

  # test <- buildBranchCellDataSet(HSMM, progenitor_method = "sequential_split", branch_point = 2,
  #                                branch_labels = NULL, stretch = TRUE)
  #
  # test <- buildBranchCellDataSet(HSMM, progenitor_method = "duplicate", branch_point = 2,
  #                                branch_labels = NULL, stretch = TRUE)

  BEAM_res <- BEAM(HSMM, branch_point = branch_point, cores = 6,  progenitor_method = "sequential_split")

  BEAM_res <- BEAM_res[order(BEAM_res$qval),]
  BEAM_res <- BEAM_res[,c("pval", "qval")]

  return(BEAM_res)

}

#' Monocle2 differential expression
#'
#' @param seu
#'
#' @return
#' @export
#'
#' @examples
Monocle2_diffex <- function(seu){
  diff_test_res <- monocle::differentialGeneTest(cds_subset,
                                        fullModelFormulaStr = "~sm.ns(Pseudotime)")
}

#' Run heatmap from Monocle2
#'
#' http://cole-trapnell-lab.github.io/monocle-release/docs/#differential-expression-analysis
#'
#' @param HSMM
#' @param BEAM_res
#' @param branches
#' @param pt_param
#' @param pt_paramval
#' @param colorval
#' @param top_genes
#' @param out_pdf
#'
#' @return
#' @export
#'
#' @examples
run_hmap <- function(HSMM, BEAM_res, branches, pt_param, pt_paramval, colorval, top_genes, out_pdf){
  beam_pdf <- gsub(".pdf", "_beam.pdf", out_pdf)
  pdf(beam_pdf, height = 10)
  # browser()
  plotted_trx <- row.names(subset(BEAM_res, qval < 1e-4))
  gene_df <- data.frame(row.names = plotted_trx, gene_symbol = lookup_genes(plotted_trx))

  png(file = "test.png")
  mult_hmap <- plot_multiple_branches_heatmap(HSMM[plotted_trx,],
                                              branches = branches,
                                              num_clusters = 2,
                                              cores = 6,
                                              use_gene_short_name = F,
                                              show_rownames = T, return_heatmap = TRUE)
  dev.off()

  mult_hmap$gtable$grobs[[3]]$gp$fontsize <- 4

  return_trx <- mult_hmap$gtable$grobs[[3]]$label
  return_genes <- lookup_genes(return_trx)

  hmap_genes <- data.frame("ensembl_transcript_id" = return_trx, "gene_symbol" = return_genes, "pval" = BEAM_res[return_trx,]$pval, "qval" = BEAM_res[return_trx,]$qval)
  hmap_genes_csv <- gsub(".pdf", paste0("_beam_heatmap_genes.csv"), out_pdf)
  write.csv(hmap_genes, hmap_genes_csv)

  hmap_genes <- dplyr::mutate(hmap_genes, new_labels = paste(ensembl_transcript_id, gene_symbol, sep = "\t"))

  mult_hmap$gtable$grobs[[3]]$label <- hmap_genes$new_labels

  gridExtra::grid.arrange(mult_hmap[[4]])


  root_branch <- GM_state(HSMM, pt_param, pt_paramval)
  branch_states = branches[-root_branch]

  # plot individual genes ---------------------------------------------------
  # top_BEAM <- BEAM_res[1:top_genes,]
  top_BEAM <- BEAM_res[return_trx,]
  top_genes <- length(return_trx)



  batch_plot_multiple_branches_pseudotime(colorval, top_BEAM, top_genes, branches, HSMM)
  #
  # batch_plot_genes_branched_pseudotime(colorval, top_BEAM, top_genes, branch_states, HSMM)

  dev.off()
  return(mult_hmap)

}
