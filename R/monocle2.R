#' Convert a Seurat V3 object to a Monocle v2 object
#'
#' @param seu
#'
#' @return
#' @export
#'
#' @examples
convert_seuv3_to_monoclev2 <- function(seu, return_census = FALSE, sig_slice = 1000) {

  # Load Seurat object
  # Extract data, phenotype data, and feature data from the SeuratObject
  data <- as(as.matrix(seu@assays$RNA@data), "sparseMatrix")

  data <- floor(data)

  pd <- new("AnnotatedDataFrame", data = seu@meta.data)

  fData <- data.frame(gene_short_name = row.names(data), row.names = row.names(data))
  fd <- new("AnnotatedDataFrame", data = fData)

  # Construct monocle cds
  monocle_cds <- monocle::newCellDataSet(data,
                                         phenoData = pd,
                                         featureData = fd,
                                         lowerDetectionLimit = 0.5,
                                         expressionFamily = VGAM::negbinomial.size()
  )

  if(return_census){
    rpc_matrix <- monocle::relative2abs(monocle_cds, method = "num_genes")

    monocle_cds <- monocle::newCellDataSet(as(rpc_matrix, "sparseMatrix"),
                                           phenoData = pd,
                                           featureData = fd,
                                           lowerDetectionLimit = 0.5,
                                           expressionFamily = VGAM::negbinomial.size()
    )

  }


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

  # find expressed genes
  monocle_cds_expressed_genes <- rownames(subset(Biobase::featureData(monocle_cds), Biobase::featureData(monocle_cds)$num_cells_expressed >= 10))

  print("running differential expression test")
  tictoc::tic("finished differentiial expression with")

  diff_test_res <- differentialGeneTest(monocle_cds_red[monocle_cds_expressed_genes, ],
                                        fullModelFormulaStr = "~Cluster",
                                        cores = 6
  )
  tictoc::toc()


  # select top 1000 signif. genes
  ordering_genes <- row.names(diff_test_res)[order(diff_test_res$qval)][1:sig_slice]

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
process_monocle_child <- function(ptime, monocle_cds, trend_formula = "~sm.ns(Pseudotime, df=3)") {
  monocle_cds <- monocle_cds[, ptime$sample_id]

  old_ptime <- phenoData(monocle_cds)$Pseudotime

  ptime$ptime <- scales::rescale(ptime$ptime, range(old_ptime))

  phenoData(monocle_cds)$Pseudotime <- ptime$ptime

  monocle_cds_expressed_genes <- rownames(subset(Biobase::featureData(monocle_cds), Biobase::featureData(monocle_cds)$num_cells_expressed >= 10))

  print("running differential expression test")
  tictoc::tic("finished differentiial expression with")

  diff_test_res <- differentialGeneTest(monocle_cds,
                                        fullModelFormulaStr = trend_formula,
                                        cores = 6
  )

  tictoc::toc()

  return(list(monocle_cds = monocle_cds, diff_test_res = diff_test_res))
}

#' Plot a Set of Heatmaps Based on a provided pseudotime
#'
#' @param monocle_list
#' @param query_name
#' @param sig_slice
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plot_all_ptimes <- function(monocle_list, query_name, sig_slice = 1000, ...) {

  sig_gene_names <- dplyr::filter(monocle_list[[query_name]]$diff_test_res, pval < 0.05) %>%
    dplyr::arrange(pval) %>%
    dplyr::slice(1:sig_slice) %>%
    dplyr::pull(gene_short_name) %>%
    identity()

  # sig_gene_names <- sig_gene_names[1:100]
  monocle_list[[query_name]]$monocle_cds <- monocle_list[[query_name]]$monocle_cds[sig_gene_names, ]

  message("creating main heatmap")
  monocle_list[[query_name]]$heatmap_matrix <- suppressWarnings(seuratTools::calc_pseudotime_heatmap(monocle_list[[query_name]]$monocle_cds,
                                                                                                      sig_gene_names = sig_gene_names,
                                                                                                      cores = 1,
                                                                                                      show_rownames = T,
                                                                                                      return_heatmap = TRUE,
                                                                                                      ...
  ))

  # # retrieve order of gene names in first heatmap
  heatmap_gene_order <- rownames(monocle_list[[query_name]]$heatmap_matrix)

  monocle_list[[query_name]]$monocle_cds <- monocle_list[[query_name]]$monocle_cds[heatmap_gene_order,]

  reference_names <- names(monocle_list)[names(monocle_list) != query_name]


  for (i in reference_names){
    monocle_list[[i]]$monocle_cds <- monocle_list[[i]]$monocle_cds[heatmap_gene_order,]

    message(paste0("creating reference heatmap ", i))
    monocle_list[[i]]$heatmap_matrix  <- suppressWarnings(seuratTools::calc_pseudotime_heatmap(monocle_list[[i]]$monocle_cds,
                                                                                        cores = 1,
                                                                                        show_rownames = F,
                                                                                        return_heatmap = TRUE,
                                                                                        cluster_rows = F,
                                                                                        trend_formula = "~sm.ns(Pseudotime, df=1)",
                                                                                        ...
    ))

  }


  monocle_list <- cross_check_heatmaps(monocle_list, query_name)

  return(monocle_list)

}

#' Title
#'
#' @param monocle_list
#' @param query_name
#' @param set_row_order
#' @param ...
#' @param cluster_rows
#'
#' @return
#' @export
#'
#' @examples
cross_check_heatmaps <- function(monocle_list, query_name, set_row_order = F, cluster_rows = T, ...){
  reference_names = names(monocle_list)[!names(monocle_list) == query_name]

  if (!set_row_order == F){
    monocle_list[[query_name]]$heatmap_matrix <- monocle_list[[query_name]]$heatmap_matrix[set_row_order,]
  }

  common_genes <- purrr::map(monocle_list, ~rownames(purrr::pluck(.x, "heatmap_matrix"))) %>%
    purrr::reduce(intersect)

  if (!set_row_order == F){
    common_genes <- common_genes[match(set_row_order, common_genes)]
    common_genes <- common_genes[!is.na(common_genes)]

    cluster_rows = F
  }


  monocle_list[[query_name]]$heatmap_matrix <- monocle_list[[query_name]]$heatmap_matrix[common_genes,]
  query_results <- seuratTools::plot_pseudotime_heatmap(monocle_list[[query_name]]$heatmap_matrix, query_name, cluster_rows = cluster_rows, heatmap_width = 12, ...)

  monocle_list[[query_name]]$heatmap <- query_results$heatmap

  for (i in reference_names){
    monocle_list[[i]]$heatmap_matrix <- monocle_list[[i]]$heatmap_matrix[common_genes,]
    monocle_list[[i]]$heatmap <- plot_pseudotime_heatmap(monocle_list[[i]]$heatmap_matrix, i, cluster_rows = F, query_set = F, ...)

  }

  heatmap_gene_order <- rownames(monocle_list[[query_name]]$heatmap_matrix)
  #
  mod_diff_test_res <- filter_rows_to_top(monocle_list[[query_name]]$diff_test_res, "gene_short_name", heatmap_gene_order)
  #
  annotation_row <- tibble::tibble(feature = heatmap_gene_order)
  #
  mod_diff_test_res <- dplyr::left_join(mod_diff_test_res, annotation_row, by = c("gene_short_name" = "feature"))
  #
  monocle_list[[query_name]]$diff_test_res <- mod_diff_test_res

  clusters <- dendextend::cutree(query_results$row_dend, k = 6) %>%
    tibble::enframe("gene_short_name", "cluster")

  labels <- labels(query_results$row_dend) %>%
    tibble::enframe("order", "gene_short_name") %>%
    dplyr::left_join(clusters, by = "gene_short_name")

  ordered_diff_test_res <- dplyr::left_join(labels, mod_diff_test_res, by = "gene_short_name")

  monocle_list[[query_name]]$ordered_diff_test_res <- ordered_diff_test_res

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
arrange_ptime_heatmaps <- function(cds_list, cds_name) {
  cds_heatmaps <- purrr::map(cds_list, ~purrr::pluck(.x, "heatmap"))

  query_heatmap <- cds_heatmaps[[cds_name]]

  reference_heatmaps <- cds_heatmaps[!names(cds_heatmaps) == cds_name]

  heatmaplist <- query_heatmap + purrr::reduce(reference_heatmaps, `+`)

  return(heatmaplist)
}

#' Plot Expression of a Given Feature of a set of Pseudotimes
#'
#' @param cds_list
#' @param features
#' @param color_by
#' @param relative_expr
#' @param ...
#' @param min_expr
#'
#' @return
#' @export
#'
#' @examples
plot_feature_in_ref_query_ptime <- function(cds_list, features = c("RXRG"), color_by = "State", relative_expr = FALSE, min_expr = 0.5, ...){
  sub_cds_list <- purrr::map(cds_list, ~.x$monocle_cds[features,])
  feature_plots_in_ptime <- purrr::map(sub_cds_list, monocle::plot_genes_in_pseudotime, trend_formula = "~sm.ns(Pseudotime, df=3)", color_by = color_by, relative_expr = relative_expr, min_expr = min_expr)
  refquery_ptime_plot <- cowplot::plot_grid(plotlist = feature_plots_in_ptime, ncol = 2, labels = names(feature_plots_in_ptime))
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

#' calculate pseudotime heatmaps
#'
#' @param cds_subset
#' @param cluster_rows
#' @param dend_k
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
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
calc_pseudotime_heatmap <- function(cds_subset, cluster_rows = TRUE, dend_k = 6, hclust_method = "ward.D2",
                                             num_clusters = 6, hmcols = NULL, add_annotation_row = NULL,
                                             add_annotation_col = NULL, show_rownames = FALSE, use_gene_short_name = TRUE,
                                             norm_method = c("log", "vstExprs"), scale_max = 3, scale_min = -3,
                                             trend_formula = "~sm.ns(Pseudotime, df=3)", return_heatmap = FALSE,
                                             cores = 1, ...)
{

  Biobase::fData(cds_subset) <- data.frame(row.names = rownames(cds_subset), gene_short_name = rownames(cds_subset))

  num_clusters <- min(num_clusters, nrow(cds_subset))
  pseudocount <- 1
  newdata <- data.frame(Pseudotime = seq(min(Biobase::pData(cds_subset)$Pseudotime),
                                         max(Biobase::pData(cds_subset)$Pseudotime), length.out = 100))
  m <- monocle::genSmoothCurves(cds_subset, cores = cores, trend_formula = trend_formula,
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

  # drop every row (gene) that has a standard deviation of zero
  # m = m[!apply(m, 1, sd) == 0, ]
  m = Matrix::t(scale(Matrix::t(m), center = TRUE))
  m = m[is.na(row.names(m)) == FALSE, ]
  m[is.nan(m)] = 0
  m[m > scale_max] = scale_max
  m[m < scale_min] = scale_min
  heatmap_matrix <- m

    return(heatmap_matrix)
}

#' plot pseudotime heatmap
#'
#' @param heatmap_matrix
#' @param dend_k
#' @param cluster_rows
#' @param hmcols
#' @param heatmap_title
#' @param seriation
#' @param row_font
#' @param heatmap_height
#' @param query_set
#' @param heatmap_width
#'
#' @return
#' @export
#'
#' @examples
plot_pseudotime_heatmap <- function(heatmap_matrix, heatmap_title, dend_k = 6, cluster_rows = T, query_set = T, hmcols = NULL, seriation = F, row_font = 4, heatmap_height = 96, heatmap_width = 8){

  row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
  row_dist[is.na(row_dist)] <- 1

  if (is.null(hmcols)) {
    bks <- seq(-3.1, 3.1, by = 0.1)
    hmcols <- colorRamps::blue2green2red(length(bks) - 1)
  }
  else {
    bks <- seq(-3.1, 3.1, length.out = length(hmcols))
  }

  if (seriation){
    o1 = seriation::seriate(row_dist, method = "GW")
    row_dend = as.dendrogram(o1[[1]])
  } else {
    row_dend = hclust(row_dist, method = "ward.D2")
    row_dend = dendextend::color_branches(row_dend, k = dend_k)
  }

  if (cluster_rows){
    row_split = dend_k
  } else {
    row_split = NULL
    row_dend = dendextend::rotate(row_dend, rownames(heatmap_matrix))
  }

  ph_res <- ComplexHeatmap::Heatmap(heatmap_matrix,
                                    # Remove name from fill legend
                                    name = "",
                                    # Keep original row/col order
                                    row_order = rownames(heatmap_matrix),
                                    column_order = colnames(heatmap_matrix),
                                    # Add left annotation (legend with tumor/normal)
                                    # right_annotation = gene_ann,
                                    # ACTUAL SPLIT by sample group
                                    row_split = row_split,
                                    show_row_names = TRUE,
                                    show_column_names = FALSE,
                                    show_column_dend = FALSE,
                                    row_title = NULL,
                                    cluster_rows = row_dend,
                                    show_row_dend = TRUE,
                                    row_names_gp = grid::gpar(fontsize = c(row_font)),
                                    heatmap_height = unit(heatmap_height, "cm"),
                                    heatmap_width = unit(heatmap_width, "cm"),
                                    column_title = heatmap_title,
                                    row_dend_width = unit(4, "cm"))

  if (query_set){
    query_results <- list(heatmap = ph_res, row_dend = row_dend)
    return(query_results)
  } else{
    return(ph_res)
  }

}

#' Create Single Cell Experiment from Tibbles
#'
#' @param counts
#' @param colData
#' @param experimentdata
#'
#' @return
#' @export
#'
#' @examples
sce_from_tibbles <- function(counts, colData, experimentdata = NULL){
  featuredata <- data.frame(rownames(counts), row.names = rownames(counts))

  counts <- data.frame(counts)
  counts <- as.matrix(counts)

  colData <- data.frame(colData)
  rownames(colData) <- gsub("-", ".", colData$sample_id)
  colData <- colData[colnames(counts),]
  sce <- SingleCellExperiment::SingleCellExperiment(assays=list(counts=counts), colData=colData, rowData=featuredata, metadata=experimentdata)

  return(sce)
}

#' Run Census
#'
#' @param sce a single cell experiment object
#' @param census_output_file name to be given to output file
#'
#' @return
#' @export
#'
#' @examples
run_census <- function(sce, census_output_file){

  # create new celldataset
  pd <- new("AnnotatedDataFrame", data=data.frame(colData(sce)))
  fd <- new("AnnotatedDataFrame", data=data.frame(rowData(sce)))

  HSMM <- monocle::newCellDataSet(counts(sce),
                         phenoData = pd,
                         featureData = fd,
                         lowerDetectionLimit=0.1,
                         expressionFamily=tobit(Lower=0.1))

  rpc_matrix <- monocle::relative2abs(HSMM, method = "num_genes")

  # Now, make a new CellDataSet using the RNA counts
  HSMM <- monocle::newCellDataSet(as(rpc_matrix, "sparseMatrix"),
                         phenoData = pd,
                         featureData = fd,
                         lowerDetectionLimit=1,
                         expressionFamily=negbinomial.size())

  return(HSMM)

  #print output of census to csv prior to monocle workflow

  write.csv(as.matrix(Biobase::exprs(HSMM)), census_output_file)

  census_meta_file <- gsub("census_matrix.csv", "census_meta.csv", census_output_file)

  write.csv(pData(HSMM), census_meta_file, row.names = FALSE)

  return(as.matrix(Biobase::exprs(HSMM)))

}

