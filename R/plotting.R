
#' Cross plot vars
#'
#' @param seu
#' @param resolution
#' @param mycols
#'
#' @return
#' @export
#'
#' @examples
cross_plot_vars <- function(seu, resolution, mycols) {

  if ("integrated" %in% names(seu@assays)) {
    active_assay <- "integrated"
  }
  else {
    active_assay <- "RNA"
  }

  cluster_resolution = paste0(active_assay,
                              "_snn_res.", resolution)
  mycols <- gsub("^seurat$", cluster_resolution, mycols)
  newcolname = paste(mycols, collapse = "_by_")

  newdata <- seu[[mycols]] %>%
    tidyr::unite(!!newcolname, mycols)

  Idents(seu) <- newdata

  return(seu)

}

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


#' Plot RNA velocity Computed by Velocyto.R
#'
#' @param seu
#' @param reduction
#' @param arrow.scale
#' @param cell.colors
#' @param plot_format
#' @param velocity
#'
#' @return
#' @export
#'
#' @examples
plot_velocity_arrows <- function(seu, velocity, reduction = "umap", cell.colors, plot_format = "grid", arrow.scale = 3){

  velocity <- seu@misc$vel

  emb <- Embeddings(object = seu, reduction = reduction)

  cell.alpha=1.0; cell.cex=1; fig.height=4; fig.width=4.5;

  if (plot_format == "arrow") {
    velocyto.R::show.velocity.on.embedding.cor(emb, velocity, n=100, scale='sqrt',
                                               cell.colors=velocyto.R::ac(cell.colors, alpha=cell.alpha),
                                               cex=cell.cex, arrow.scale=arrow.scale, arrow.lwd=1)
  } else if (plot_format == "grid"){
    #Alternatively, the same function can be used to calculate a velocity vector field:
    velocyto.R::show.velocity.on.embedding.cor(emb, velocity, n=100, scale='sqrt',
                                               cell.colors=velocyto.R::ac(cell.colors, alpha=cell.alpha),
                                               cex=cell.cex, arrow.scale=arrow.scale,
                                               show.grid.flow=TRUE, min.grid.cell.mass=0.5,
                                               grid.n=20, arrow.lwd=2)
  }


}

#' Plot RNA velocity trajectory computed by Velocyto.R
#'
#' @param seu
#' @param reduction
#' @param format either an arrow or grid format
#' @param arrow.scale
#'
#' @return
#' @export
#'
#' @examples
plot_velocity_trajectory <- function(seu, reduction = "umap", format = "arrow", arrow.scale = 3){

  emb <- Embeddings(object = seu, reduction = reduction)
  ident.colors <- (scales::hue_pal())(n = length(x = levels(x = seu)))
  names(x = ident.colors) <- levels(x = seu)
  cell.colors <- ident.colors[Idents(object = seu)]
  names(x = cell.colors) <- colnames(x = seu)

  cell.alpha=1.0; cell.cex=1; fig.height=4; fig.width=4.5;

  velocyto.R::show.velocity.on.embedding.eu(emb, seu@misc$vel, n=40, scale='sqrt', cell.colors=ac(cell.colors,alpha=cell.alpha),
                                            cex=cell.cex, nPcs=30, sigma=2.5, show.trajectories=TRUE, diffusion.steps=400,
                                            n.trajectory.clusters=15, ntop.trajectories=1,
                                            embedding.knn=T, control.for.neighborhood.density=TRUE,n.cores=6)

}

#' Plot Metadata Variables
#'
#' @param seu
#' @param embedding
#' @param group
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plot_var <- function(seu, embedding = "umap", group = "batch", dims = c(1,2), highlight = NULL, ...){

  metadata <- tibble::as_tibble(seu[[]][Seurat::Cells(seu),], rownames = "sID")
  cellid <- metadata[["sID"]]
  key <- rownames(metadata)

  if (embedding == "umap"){
    dims = c(1,2)
  } else if (embedding == "tsne"){
    dims = c(1,2)
  }

  dims <- as.numeric(dims)

  d <- Seurat::DimPlot(object = seu, dims = dims, reduction = embedding, group.by = group, pt.size = 1.0, ...) +
    aes(key = key, cellid = cellid) +
    # gghighlight()
    # theme(legend.text=element_text(size=10)) +
    NULL

  plotly_plot <- plotly::ggplotly(d, tooltip = "cellid", height  = 500) %>%
    # htmlwidgets::onRender(javascript) %>%
    # plotly::highlight(on = "plotly_selected", off = "plotly_relayout") %>%
    plotly_settings() %>%
    plotly::toWebGL() %>%
    # plotly::partial_bundle() %>%
    identity()

}

#' Plotly settings
#'
#' @param plotly_plot
#'
#' @return
#' @export
#'
#' @examples
plotly_settings <- function(plotly_plot){
  plotly_plot %>%
    plotly::layout(dragmode = "lasso") %>%
    plotly::config(
      toImageButtonOptions = list(
        format = "png",
        filename = "myplot",
        width = 600,
        height = 700
      )) %>%
    identity()
}


#' plot Violin plot
#'
#' @param seu
#' @param plot_var
#' @param plot_vals
#' @param features
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plot_violin <- function(seu, plot_var = "batch", plot_vals = NULL, features = "RXRG",
                        ...){
  if (is.null(plot_vals)) {
    plot_vals = unique(seu[[]][[plot_var]])
    plot_vals <- plot_vals[!is.na(plot_vals)]
  }
  seu <- seu[, seu[[]][[plot_var]] %in% plot_vals]
  vln_plot <- Seurat::VlnPlot(seu, features = features, group.by = plot_var,
                              ...) + labs(title = "Expression Values for each cell are normalized by that cell's total expression then multiplied by 10,000 and natural-log transformed") +
    stat_summary(fun.y = mean, geom = "line", size = 4,
                 colour = "black")
  return(vln_plot)
}

#' plot heatmap
#'
#' @param seu
#' @param features
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plot_heatmap <- function(seu, features = "RXRG", ...){

  if ("integrated" %in% names(seu@assays)) {
    default_assay = "integrated"
  } else {
    default_assay = "RNA"
  }

  hm <- Seurat::DoHeatmap(seu, features = features, assay = default_assay, ...) +
    labs(title = "Expression Values for each cell are normalized by that cell's total expression then multiplied by 10,000 and natural-log transformed") +
    NULL

  return(hm)
}



#' Plot Features
#' This is a great function
#' @param seu
#' @param embedding
#' @param features
#' @param dims
#'
#' @return
#' @export
#'
#' @examples
plot_feature <- function(seu, embedding, features, dims = c(1,2), return_plotly = TRUE){

  metadata <- tibble::as_tibble(seu[[]][Seurat::Cells(seu),], rownames = "sID")

  cellid <- metadata[["sID"]]
  key <- rownames(metadata)

  if (embedding %in% c("tsne", "umap")){
    dims = c(1,2)
  }

  dims <- as.numeric(dims)

  fp <- Seurat::FeaturePlot(object = seu, features = features, dims = dims, reduction = embedding, pt.size = 1.0, blend = FALSE)	+
    aes(key = key, cellid = cellid, alpha = 0.7)

  if (return_plotly == FALSE) return(fp)

  plotly_plot <- plotly::ggplotly(fp, tooltip = "cellid", height = 500) %>%
    plotly_settings() %>%
    plotly::toWebGL() %>%
    # plotly::partial_bundle() %>%
    identity()

}

#' Plot Rides
#'
#' @param seu
#' @param features
#'
#' @return
#' @export
#'
#' @examples
plot_ridge <- function(seu, features){

  cc_genes_path <- "~/single_cell_projects/resources/regev_lab_cell_cycle_genes.txt"
  cc.genes <- readLines(con = cc_genes_path)
  s.genes <- cc.genes[1:43]
  g2m.genes <- cc.genes[44:97]

  seu <- CellCycleScoring(object = seu, s.genes, g2m.genes,
                          set.ident = TRUE)

  RidgePlot(object = seu, features = features)

  # plotly::ggplotly(r, height = 750)
  #
}


#' Plot Cluster Marker Genes
#'
#' @param seu
#' @param resolution
#' @param selected_clusters
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plot_markers <- function(seu, resolution, num_markers = 5, selected_clusters = NULL, return_plot = FALSE, ...){

  Idents(seu) <- seu[[resolution]]

  markers <- seu@misc$markers[[resolution]] %>%
    dplyr::top_n(n = num_markers, wt = logFC) %>%
    identity()

  if(!is.null(selected_clusters)){
    seu <- seu[,Idents(seu) %in% selected_clusters]
    markers <- markers %>%
      dplyr::filter(group %in% selected_clusters)
  }

  markers <- dplyr::pull(markers, feature)

  markerplot <- DotPlot(seu, features = unique(markers), group.by = resolution, dot.scale = 3) +
    theme(axis.text.x = element_text(size = 10, angle = 45, vjust = 1, hjust=1),
          axis.text.y = element_text(size = 10)) +
    scale_y_discrete(position = "right") +
    coord_flip() +
    NULL

  if (return_plot) return(markerplot)

  plot_height = (150*num_markers)
  plot_width = (100*length(levels(Idents(seu))))

  plotly_plot <- plotly::ggplotly(markerplot, height = plot_height, width = plot_width) %>%
    plotly_settings() %>%
    plotly::toWebGL() %>%
    # plotly::partial_bundle() %>%
    identity()

}

#' Plot Read Count
#'
#' @param seu
#' @param plot_type
#'
#' @return
#' @export
#'
#' @examples
plot_readcount <- function(seu, plot_type){

  seu_tbl <- tibble::rownames_to_column(seu[[]], "SID")

  rc_plot <- ggplot(seu_tbl, aes(x=reorder(SID, -nCount_RNA), y = nCount_RNA, fill = !!as.symbol(plot_type))) +
    # scale_y_continuous(breaks = seq(0, 8e7, by = 5e5)) +
    scale_y_log10() +
    geom_bar(position = "identity", stat = "identity") +
    # geom_text(data=subset(agg_qc_wo_na, Sample %in% thresholded_cells & align_type == "paired_total"),
    #   aes(Sample, count, label=Sample)) +
    # geom_text(data=subset(agg_qc_wo_na, Sample %in% low_read_count_cells & align_type == "paired_aligned_one"),
    #   aes(Sample, count, label=Sample)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    # scale_fill_manual(values = c( "low_read_count"="tomato", "keep"="gray" ), guide = FALSE ) +
    labs(title = "Paired Unique Reads", x = "Sample") +
    NULL

  rc_plot <- rc_plot %>%
    plotly::toWebGL() %>%
    # plotly::partial_bundle() %>%
    identity()

}

