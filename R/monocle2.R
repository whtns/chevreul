#' Convert a Seurat V3 object to a Monocle v2 object
#'
#' @param seu
#'
#' @return
#' @export
#'
#' @examples
convert_seuv3_to_monoclev2 <- function(seu) {
  # browser()
  # Load Seurat object


  # Extract data, phenotype data, and feature data from the SeuratObject
  data <- as(as.matrix(seu@assays$RNA@data), "sparseMatrix")

  pd <- new("AnnotatedDataFrame", data = seu@meta.data)

  fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
  fd <- new("AnnotatedDataFrame", data = fData)

  # Construct monocle cds
  monocle_cds <- monocle::newCellDataSet(data,
                                         phenoData = pd,
                                         featureData = fd,
                                         lowerDetectionLimit = 0.5,
                                         expressionFamily = negbinomial.size()
  )

  # filter by gene expression
  # # keep genes that are expressed in at least 5 cells
  # min_expression <- 0.1
  # monocle_cds <- detectGenes(monocle_cds, min_expr=min_expression)
  # expressed_genes <- row.names(subset(featureData(monocle_cds), num_cells_expressed >= 5))

  # filter by gene expression (default monocle settings)
  monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)

  expressed_genes <- row.names(subset(
    Biobase::featureData(monocle_cds)@data,
    num_cells_expressed >= 10
  ))

  # look at distribution of mRNA totals across cells
  phenoData(monocle_cds)$Total_mRNAs <- Matrix::colSums(Biobase::exprs(monocle_cds))

  monocle_cds <- monocle_cds[, phenoData(monocle_cds)$Total_mRNAs < 1e6]

  upper_bound <- 10^(mean(log10(phenoData(monocle_cds)$Total_mRNAs)) +
                       2 * sd(log10(phenoData(monocle_cds)$Total_mRNAs)))
  lower_bound <- 10^(mean(log10(phenoData(monocle_cds)$Total_mRNAs)) -
                       2 * sd(log10(phenoData(monocle_cds)$Total_mRNAs)))

  # remove cells outside safe range of plot ---------------------------------

  monocle_cds <- monocle_cds[, phenoData(monocle_cds)$Total_mRNAs > lower_bound &
                               phenoData(monocle_cds)$Total_mRNAs < upper_bound]
  monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)

  monocle_cds <- estimateSizeFactors(monocle_cds)
  monocle_cds <- estimateDispersions(monocle_cds)

  # verify lognormal distribution of expression values ----------------------

  # Log-transform each value in the expression matrix.
  L <- log(Biobase::exprs(monocle_cds[expressed_genes, ]))

  # Standardize each gene, so that they are all on the same scale,
  # Then melt the data with plyr so we can plot it easily
  # melted_dens_df <- melt(Matrix::t(scale(Matrix::t(L))))


  # use dpFeature -----------------------------------------------------------
  # browser()
  monocle_cds <- detectGenes(monocle_cds, min_expr = 0.1)
  Biobase::featureData(monocle_cds)$use_for_ordering <- Biobase::featureData(monocle_cds)$num_cells_expressed > 0.05 * ncol(monocle_cds)

  # print(plot_pc_variance_explained(monocle_cds, return_all = F))

  monocle_cds_red <- reduceDimension(monocle_cds,
                                     max_components = 2,
                                     norm_method = "log",
                                     num_dim = 3,
                                     reduction_method = "tSNE",
                                     verbose = T
  )

  monocle_cds_red <- clusterCells(monocle_cds_red, verbose = F)

  # check clustering results
  print(plot_cell_clusters(monocle_cds_red, color_by = "as.factor(Cluster)"))

  # for (i in colorval){
  #   print(plot_cell_clusters(monocle_cds_red, color_by = paste0('as.factor(', i, ')')))
  # }


  # provide decision plot
  print(plot_rho_delta(monocle_cds_red, rho_threshold = 2, delta_threshold = 4))

  # rerun based on user-defined threshold
  monocle_cds_red <- clusterCells(monocle_cds_red,
                                  rho_threshold = 2,
                                  delta_threshold = 4,
                                  skip_rho_sigma = T,
                                  verbose = F
  )

  # check final clustering
  print(plot_cell_clusters(monocle_cds_red, color_by = "as.factor(Cluster)"))
  #
  # for (i in colorval){
  #   print(plot_cell_clusters(monocle_cds_red, color_by = paste0('as.factor(', i, ')')))
  # }



  # perform differential expression -----------------------------------------

  # browser()

  # find expressed genes
  monocle_cds_expressed_genes <- rownames(subset(Biobase::featureData(monocle_cds), Biobase::featureData(monocle_cds)$num_cells_expressed >= 10))

  print("running differential expression test")
  tictoc::tic("finished differentiial expression with")

  diff_test_res <- differentialGeneTest(monocle_cds_red[monocle_cds_expressed_genes, ],
                                        fullModelFormulaStr = "~Cluster",
                                        cores = 6
  )
  toc()


  # select top 1000 signif. genes
  ordering_genes <- row.names(diff_test_res)[order(diff_test_res$qval)][1:1000]



  monocle_cds <- setOrderingFilter(monocle_cds, ordering_genes = ordering_genes)

  monocle_cds <- reduceDimension(monocle_cds, method = "DDRTree")

  monocle_cds <- orderCells(monocle_cds)

}

#' Preprocess a Monocle v2 object for heatmap based on provided pseudotime
#'
#' @param ptime
#' @param monocle_cds
#'
#' @return
#' @export
#'
#' @examples
process_monocle_child <- function(ptime, monocle_cds) {
  monocle_cds <- monocle_cds[, ptime$Sample_ID]

  old_ptime <- phenoData(monocle_cds)$Pseudotime

  ptime$ptime <- scales::rescale(ptime$ptime, range(old_ptime))

  phenoData(monocle_cds)$Pseudotime <- ptime$ptime

  monocle_cds_expressed_genes <- rownames(subset(Biobase::featureData(monocle_cds), Biobase::featureData(monocle_cds)$num_cells_expressed >= 10))

  print("running differential expression test")
  tictoc::tic("finished differentiial expression with")

  diff_test_res <- differentialGeneTest(monocle_cds,
                                        fullModelFormulaStr = "~sm.ns(Pseudotime)",
                                        cores = 6
  )

  toc()

  return(list(monocle_cds = monocle_cds, diff_test_res = diff_test_res))
}

#' Plot a heatmap of differentially expressed genes over a query and set of reference pseudotimes
#'
#' @param monocle_list
#' @param query_name
#'
#' @return
#' @export
#'
#' @examples
plot_monocle_all_ptimes <- function(monocle_list, query_name, ...) {

  monocle_list[[query_name]] <- monocle_list[[query_name]]

  sig_gene_names <- dplyr::filter(monocle_list[[query_name]]$diff_test_res, pval < 0.05) %>%
    dplyr::arrange(pval) %>%
    dplyr::pull(gene_short_name) %>%
    identity()

  sig_gene_names <- sig_gene_names[1:100]
  monocle_list[[query_name]]$monocle_cds <- monocle_list[[query_name]]$monocle_cds[sig_gene_names, ]

  message("creating main heatmap")
  monocle_list[[query_name]]$heatmap <- seuratTools::plot_pseudotime_heatmap(monocle_list[[query_name]]$monocle_cds,
                                                                cores = 6,
                                                                show_rownames = T,
                                                                return_heatmap = TRUE,
                                                                ...
  )

  # retrieve order of gene names in first heatmap
  gene_order <- monocle_list[[query_name]]$heatmap[["tree_row"]][["labels"]][monocle_list[[query_name]]$heatmap[["tree_row"]][["order"]]]



  monocle_list[[query_name]]$monocle_cds <- monocle_list[[query_name]]$monocle_cds[gene_order,]

  message("recreating main heatmap")
  monocle_list[[query_name]]$heatmap <- seuratTools::plot_pseudotime_heatmap(monocle_list[[query_name]]$monocle_cds,
                                                                cores = 6,
                                                                show_rownames = T,
                                                                return_heatmap = TRUE,
                                                                cluster_rows = F,
                                                                ...
  )

  # reference_cds <- monocle_list[names(monocle_list) != query_name]
  reference_names <- names(monocle_list)[names(monocle_list) != query_name]

  # browser()
  for (i in reference_names){
    monocle_list[[i]]$monocle_cds <- monocle_list[[i]]$monocle_cds[gene_order,]

    message(paste0("creating reference heatmap ", i))
    monocle_list[[i]]$heatmap  <- seuratTools::plot_pseudotime_heatmap(monocle_list[[i]]$monocle_cds,
                                                          cores = 6,
                                                          show_rownames = T,
                                                          return_heatmap = TRUE,
                                                          cluster_rows = F,
                                                          ...
    )
  }

  return(monocle_list)

}

#' Arrange processed monocle CDS objects containing heatmaps
#'
#' @param cds_list
#'
#' @return
#' @export
#'
#' @examples
arrange_ptime_heatmaps <- function(cds_list) {
  cds_heatmaps <- purrr::map(cds_list, ~purrr::pluck(.x, "heatmap", 4))

  message("arranging heatmaps")
  arranged_heatmaps <- cowplot::plot_grid(
    plotlist = cds_heatmaps,
    labels = names(cds_heatmaps),
    align = "h",
    nrow = 1) +
    theme(plot.margin=unit(c(1,1,1.5,1.2),"cm"))
}

#' Plot Expression of a Given Feature of a set of Pseudotimes
#'
#' @param cds_list
#' @param features
#'
#' @return
#' @export
#'
#' @examples
plot_feature_in_ref_query_ptime <- function(cds_list, features = c("RXRG"), color_by = "State", relative_expr = FALSE){
  sub_cds_list <- purrr::map(cds_list, ~.x$monocle_cds[features,])
  feature_plots_in_ptime <- purrr::map(sub_cds_list, ~monocle::plot_genes_in_pseudotime(.x, color_by = color_by, relative_expr = relative_expr))
  refquery_ptime_plot <- cowplot::plot_grid(plotlist = feature_plots_in_ptime, ncol = 1, labels = names(feature_plots_in_ptime))
  plot(refquery_ptime_plot)
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

#' Plot Pseudotime Heatmap
#'
#' @param cds_subset
#' @param cluster_rows
#' @param hclust_method
#' @param num_clusters
#' @param hmcols
#' @param add_annotation_row
#' @param add_annotation_col
#' @param show_rownames
#' @param use_gene_short_name
#' @param norm_method
#' @param scale_max
#' @param scale_min
#' @param trend_formula
#' @param return_heatmap
#' @param cores
#'
#' @return
#' @export
#'
#' @examples
plot_pseudotime_heatmap <- function(cds_subset, cluster_rows = TRUE, hclust_method = "ward.D2",
          num_clusters = 6, hmcols = NULL, add_annotation_row = NULL,
          add_annotation_col = NULL, show_rownames = FALSE, use_gene_short_name = TRUE,
          norm_method = c("log", "vstExprs"), scale_max = 3, scale_min = -3,
          trend_formula = "~sm.ns(Pseudotime, df=3)", return_heatmap = FALSE,
          cores = 1, ...)
{
  num_clusters <- min(num_clusters, nrow(cds_subset))
  pseudocount <- 1
  newdata <- data.frame(Pseudotime = seq(min(Biobase::pData(cds_subset)$Pseudotime),
                                         max(Biobase::pData(cds_subset)$Pseudotime), length.out = 100))
  m <- genSmoothCurves(cds_subset, cores = cores, trend_formula = trend_formula,
                       relative_expr = T, new_data = newdata)
  m = m[!apply(m, 1, sum) == 0, ]
  norm_method <- match.arg(norm_method)
  if (norm_method == "vstExprs" && is.null(cds_subset@dispFitInfo[["blind"]]$disp_func) ==
      FALSE) {
    m = vstExprs(cds_subset, expr_matrix = m)
  }
  else if (norm_method == "log") {
    m = log10(m + pseudocount)
  }
  # m = m[!apply(m, 1, sd) == 0, ]
  m = Matrix::t(scale(Matrix::t(m), center = TRUE))
  m = m[is.na(row.names(m)) == FALSE, ]
  m[is.nan(m)] = 0
  m[m > scale_max] = scale_max
  m[m < scale_min] = scale_min
  heatmap_matrix <- m
  row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
  row_dist[is.na(row_dist)] <- 1
  if (is.null(hmcols)) {
    bks <- seq(-3.1, 3.1, by = 0.1)
    hmcols <- colorRamps::blue2green2red(length(bks) - 1)
  }
  else {
    bks <- seq(-3.1, 3.1, length.out = length(hmcols))
  }
  ph <- pheatmap(heatmap_matrix, useRaster = T, cluster_cols = FALSE,
                 cluster_rows = cluster_rows, show_rownames = F, show_colnames = F,
                 clustering_distance_rows = row_dist, clustering_method = hclust_method,
                 cutree_rows = num_clusters, silent = TRUE, filename = NA,
                 breaks = bks, border_color = NA, color = hmcols, ...)
  if (cluster_rows) {
    annotation_row <- data.frame(Cluster = factor(cutree(ph$tree_row,
                                                         num_clusters)))
  }
  else {
    annotation_row <- NULL
  }
  if (!is.null(add_annotation_row)) {
    old_colnames_length <- ncol(annotation_row)
    annotation_row <- cbind(annotation_row, add_annotation_row[row.names(annotation_row),
                                                               ])
    colnames(annotation_row)[(old_colnames_length + 1):ncol(annotation_row)] <- colnames(add_annotation_row)
  }
  if (!is.null(add_annotation_col)) {
    annotation_col <- add_annotation_col
  }
  else {
    annotation_col <- NA
  }
  if (use_gene_short_name == TRUE) {
    if (is.null(Biobase::fData(cds_subset)$gene_short_name) == FALSE) {
      feature_label <- as.character(Biobase::fData(cds_subset)[row.names(heatmap_matrix),
                                                      "gene_short_name"])
      feature_label[is.na(feature_label)] <- row.names(heatmap_matrix)
      feature_label <- make.names(row.names(heatmap_matrix), unique = T)
      row_ann_labels <- as.character(Biobase::fData(cds_subset)[row.names(annotation_row),
                                                       "gene_short_name"])
      row_ann_labels[is.na(row_ann_labels)] <- row.names(annotation_row)
      row_ann_labels <- make.names(row.names(annotation_row), unique = T)
    }
    else {
      # feature_label <- row.names(heatmap_matrix)
      feature_label <- make.names(row.names(heatmap_matrix), unique=TRUE)
      row_ann_labels <- row.names(annotation_row)
    }
  }
  else {
    # feature_label <- row.names(heatmap_matrix)
    feature_label <- make.names(row.names(heatmap_matrix), unique=TRUE)
    if (!is.null(annotation_row))
      row_ann_labels <- row.names(annotation_row)
  }
  row.names(heatmap_matrix) <- feature_label
  if (!is.null(annotation_row))
    row.names(annotation_row) <- row_ann_labels
  colnames(heatmap_matrix) <- c(1:ncol(heatmap_matrix))
  ph_res <- pheatmap(heatmap_matrix[, ], useRaster = T, cluster_cols = FALSE,
                     cluster_rows = cluster_rows, show_rownames = show_rownames,
                     show_colnames = F, clustering_distance_rows = row_dist,
                     clustering_method = hclust_method, cutree_rows = num_clusters,
                     annotation_row = annotation_row, annotation_col = annotation_col,
                     treeheight_row = 20, breaks = bks, color = hmcols,
                     border_color = NA, silent = TRUE, filename = NA, ...)
  grid::grid.rect(gp = grid::gpar("fill", col = NA))
  grid::grid.draw(ph_res$gtable)
  if (return_heatmap) {
    return(ph_res)
  }
}

