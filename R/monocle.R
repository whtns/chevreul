#' Convert a Seurat Object to a Monocle Cell Data Set
#'
#' @param seu
#'
#' @return
#' @export
#'
#' @examples
convert_seu_to_cds <- function(seu) {

  ### Building the necessary parts for a basic cds

  # part one, gene annotations

  gene_annotation <- data.frame(gene_short_name = rownames(seu), row.names = rownames(seu))

  # part two, cell information

  cell_metadata <- seu[[]]

  # part three, counts sparse matrix

  expression_matrix <- GetAssayData(seu, slot = "counts")


  ### Construct the basic cds object

  cds_from_seurat <- monocle3::new_cell_data_set(expression_matrix,
                                       cell_metadata = cell_metadata,
                                       gene_metadata = gene_annotation)


  ### Construct and assign the made up partition

  recreate.partition <- c(rep(1, length(cds_from_seurat@colData@rownames)))
  names(recreate.partition) <- cds_from_seurat@colData@rownames
  recreate.partition <- as.factor(recreate.partition)

  cds_from_seurat@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition


  ### Assign the cluster info

  cds_from_seurat@clusters@listData[["UMAP"]][["clusters"]] <- seu[[]][[paste0(Seurat::DefaultAssay(seu), '_snn_res.1')]]


  ### Could be a space-holder, but essentially fills out louvain parameters

  cds_from_seurat@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"


  ### Assign UMAP coordinate

  cds_from_seurat@reducedDims@listData[["UMAP"]] <- Embeddings(seu, "umap")


  ### Assign feature loading for downstream module analysis

  cds_from_seurat@preprocess_aux$gene_loadings <- Loadings(seu, "pca")


  ### Learn graph, this step usually takes a significant period of time for larger samples

  print("Learning graph, which can take a while depends on the sample")

  cds_from_seurat <- monocle3::learn_graph(cds_from_seurat, use_partition = T)


  ### Plot cluster info with trajectory

  print("Plotting clusters")

  # pdf(sprintf("%s/clusters.with.trajectory.%s.pdf", output.dir, Dim), width = 10, height = 10)
  # clus <- monocle3::plot_cells(cds_from_seurat,
  #                    color_cells_by = paste0(DefaultAssay(seu), '_snn_res.1'),
  #                    label_cell_groups = FALSE,
  #                    label_groups_by_cluster=FALSE,
  #                    label_leaves=FALSE,
  #                    label_branch_points=FALSE)
  # clus
  return(cds_from_seurat)

}

#' Plot a Monocle Cell Data Set
#'
#' @param cds
#' @param resolution
#'
#' @return
#' @export
#'
#' @examples
plot_cds <- function(cds, resolution){

  key <- seq(1, length(colnames(cds)))
  cellid <- colnames(cds)

  cds[['key']] = key
  cds[['cellid']] = cellid

  if (grepl("integrated", colnames(colData(cds)))){
    default_assay = "integrated"
  } else {
    default_assay = "RNA"
  }

  cds_plot <- monocle3::plot_cells(cds,
                                   label_cell_groups = FALSE,
                                   label_groups_by_cluster = FALSE,
                                   label_leaves = FALSE,
                                   label_branch_points = FALSE,
                                   color_cells_by = paste0(default_assay, "_snn_res.", resolution)) +
    # aes(key = key, cellid = cellid) +
    NULL


  plotly::ggplotly(cds_plot, height = 750) %>%
    plotly::layout(dragmode = "lasso") %>%
    identity()


}

