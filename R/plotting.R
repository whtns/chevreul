
#' Plot pseudotime over multiple branches
#'
#' @param cds
#' @param branches
#' @param branches_name
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
plot_multiple_branches_heatmap <- function(cds,
                                           branches,
                                           branches_name = NULL,
                                           cluster_rows = TRUE,
                                           hclust_method = "ward.D2",
                                           num_clusters = 6,

                                           hmcols = NULL,

                                           add_annotation_row = NULL,
                                           add_annotation_col = NULL,
                                           show_rownames = FALSE,
                                           use_gene_short_name = TRUE,

                                           norm_method = c("vstExprs", "log"),
                                           scale_max=3,
                                           scale_min=-3,

                                           trend_formula = '~sm.ns(Pseudotime, df=3)',

                                           return_heatmap=FALSE,
                                           cores=1){
  pseudocount <- 1
  if(!(all(branches %in% Biobase::pData(cds)$State)) & length(branches) == 1){
    stop('This function only allows to make multiple branch plots where branches is included in the pData')
  }

  branch_label <- branches
  if(!is.null(branches_name)){
    if(length(branches) != length(branches_name)){
      stop('branches_name should have the same length as branches')
    }
    branch_label <- branches_name
  }

  #test whether or not the states passed to branches are true branches (not truncks) or there are terminal cells
  g <- cds@minSpanningTree
  m <- NULL
  # branche_cell_num <- c()
  for(branch_in in branches) {
    branches_cells <- row.names(subset(Biobase::pData(cds), State == branch_in))
    root_state <- subset(Biobase::pData(cds), Pseudotime == 0)[, 'State']
    root_state_cells <- row.names(subset(Biobase::pData(cds), State == root_state))

    if(cds@dim_reduce_type != 'ICA') {
      root_state_cells <- unique(paste('Y_', cds@auxOrderingData$DDRTree$pr_graph_cell_proj_closest_vertex[root_state_cells, ], sep = ''))
      branches_cells <- unique(paste('Y_', cds@auxOrderingData$DDRTree$pr_graph_cell_proj_closest_vertex[branches_cells, ], sep = ''))
    }
    root_cell <- root_state_cells[which(degree(g, v = root_state_cells) == 1)]
    tip_cell <- branches_cells[which(degree(g, v = branches_cells) == 1)]

    traverse_res <- traverseTree(g, root_cell, tip_cell)
    path_cells <- names(traverse_res$shortest_path[[1]])

    if(cds@dim_reduce_type != 'ICA') {
      pc_ind <- cds@auxOrderingData$DDRTree$pr_graph_cell_proj_closest_vertex
      path_cells <- row.names(pc_ind)[paste('Y_', pc_ind[, 1], sep = '') %in% path_cells]
    }

    cds_subset <- cds[, path_cells]

    newdata <- data.frame(Pseudotime = seq(0, max(Biobase::pData(cds_subset)$Pseudotime),length.out = 100))

    tmp <- genSmoothCurves(cds_subset, cores=cores, trend_formula = trend_formula,
                           relative_expr = T, new_data = newdata)
    if(is.null(m))
      m <- tmp
    else
      m <- cbind(m, tmp)
  }

  #remove genes with no expression in any condition
  m=m[!apply(m,1,sum)==0,]

  norm_method <- match.arg(norm_method)

  # FIXME: this needs to check that vst values can even be computed. (They can only be if we're using NB as the expressionFamily)
  if(norm_method == 'vstExprs' && is.null(cds@dispFitInfo[["blind"]]$disp_func) == FALSE) {
    m = vstExprs(cds, expr_matrix=m)
  }
  else if(norm_method == 'log') {
    m = log10(m+pseudocount)
  }

  # Row-center the data.
  m=m[!apply(m,1,sd)==0,]
  m=Matrix::t(scale(Matrix::t(m),center=TRUE))
  m=m[is.na(row.names(m)) == FALSE,]
  m[is.nan(m)] = 0
  m[m>scale_max] = scale_max
  m[m<scale_min] = scale_min

  heatmap_matrix <- m

  row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix)))/2)
  row_dist[is.na(row_dist)] <- 1

  if(is.null(hmcols)) {
    bks <- seq(-3.1,3.1, by = 0.1)
    hmcols <- colorRamps::blue2green2red(length(bks) - 1)
  }
  else {
    bks <- seq(-3.1,3.1, length.out = length(hmcols))
  }

  ph <- pheatmap(heatmap_matrix,
                 useRaster = T,
                 cluster_cols=FALSE,
                 cluster_rows=T,
                 show_rownames=F,
                 show_colnames=F,
                 clustering_distance_rows=row_dist,
                 clustering_method = hclust_method,
                 cutree_rows=num_clusters,
                 silent=TRUE,
                 filename=NA,
                 breaks=bks,
                 color=hmcols)

  annotation_col <- data.frame(Branch=factor(rep(rep(branch_label, each = 100))))
  annotation_row <- data.frame(Cluster=factor(cutree(ph$tree_row, num_clusters)))
  col_gaps_ind <- c(1:(length(branches) - 1)) * 100

  if(!is.null(add_annotation_row)) {
    old_colnames_length <- ncol(annotation_row)
    annotation_row <- cbind(annotation_row, add_annotation_row[row.names(annotation_row), ])
    colnames(annotation_row)[(old_colnames_length+1):ncol(annotation_row)] <- colnames(add_annotation_row)
    # annotation_row$bif_time <- add_annotation_row[as.character(Biobase::fData(absolute_cds[row.names(annotation_row), ])$gene_short_name), 1]
  }


  if (use_gene_short_name == TRUE) {
    if (is.null(Biobase::fData(cds)$gene_short_name) == FALSE) {
      feature_label <- as.character(Biobase::fData(cds)[row.names(heatmap_matrix), 'gene_short_name'])
      feature_label[is.na(feature_label)] <- row.names(heatmap_matrix)

      row_ann_labels <- as.character(Biobase::fData(cds)[row.names(annotation_row), 'gene_short_name'])
      row_ann_labels[is.na(row_ann_labels)] <- row.names(annotation_row)
    }
    else {
      feature_label <- row.names(heatmap_matrix)
      row_ann_labels <- row.names(annotation_row)
    }
  }
  else {
    feature_label <- row.names(heatmap_matrix)
    row_ann_labels <- row.names(annotation_row)
  }

  row.names(heatmap_matrix) <- feature_label
  row.names(annotation_row) <- row_ann_labels


  colnames(heatmap_matrix) <- c(1:ncol(heatmap_matrix))

  if(!(cluster_rows)) {
    annotation_row <- NA
  }

  ph_res <- pheatmap(heatmap_matrix[, ], #ph$tree_row$order
                     useRaster = T,
                     cluster_cols = FALSE,
                     cluster_rows = cluster_rows,
                     show_rownames=show_rownames,
                     show_colnames=F,
                     #scale="row",
                     clustering_distance_rows=row_dist, #row_dist
                     clustering_method = hclust_method, #ward.D2
                     cutree_rows=num_clusters,
                     # cutree_cols = 2,
                     annotation_row=annotation_row,
                     annotation_col=annotation_col,
                     gaps_col = col_gaps_ind,
                     treeheight_row = 20,
                     breaks=bks,
                     fontsize = 12,
                     color=hmcols,
                     silent=TRUE,
                     border_color = NA,
                     filename=NA
  )

  grid::grid.rect(gp=grid::gpar("fill", col=NA))
  grid::grid.draw(ph_res$gtable)
  if (return_heatmap){
    return(ph_res)
  }
}
