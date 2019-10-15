#' Convert a Seurat Object to a Monocle Cell Data Set
#'
#' @param seu
#'
#' @return
#' @export
#'
#' @examples
convert_seu_to_cds <- function(seu, resolution = 1) {

  ### Building the necessary parts for a basic cds

  # part two, counts sparse matrix


  if (any(grepl("integrated", names(seu[[]])))){
    default_assay = "integrated"
  } else {
    default_assay = "RNA"
  }

  DefaultAssay(seu) <- default_assay

  expression_matrix <- Seurat::GetAssayData(seu, slot = "data", assay = default_assay)

  count_matrix <- Seurat::GetAssayData(seu, slot = "counts", assay = "RNA")

  count_matrix <- count_matrix[row.names(expression_matrix),]
  count_matrix <- count_matrix[,Matrix::colSums(count_matrix) != 0]

  # part three, gene annotations

  gene_annotation <- data.frame(gene_short_name = rownames(count_matrix),
                                row.names = rownames(count_matrix))

  # part one, cell information

  cell_metadata <- seu[[]][colnames(count_matrix),]

  seu <- seu[,colnames(count_matrix)]

  ### Construct the basic cds object
  cds_from_seurat <- monocle3::new_cell_data_set(expression_data = count_matrix,
                                       cell_metadata = cell_metadata,
                                       gene_metadata = gene_annotation)

  cds_from_seurat <- cds_from_seurat[,colnames(seu)]

  # estimate size factors
  cds_from_seurat <- cds_from_seurat[, colSums(as.matrix(monocle3::exprs(cds_from_seurat))) != 0]
  cds_from_seurat <- monocle3::estimate_size_factors(cds_from_seurat)


  ### Construct and assign the made up partition

  recreate.partition <- c(rep(1, length(cds_from_seurat@colData@rownames)))
  names(recreate.partition) <- cds_from_seurat@colData@rownames
  recreate.partition <- as.factor(recreate.partition)

  cds_from_seurat@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition

  ### Could be a space-holder, but essentially fills out louvain parameters
  cds_from_seurat@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"
  cds_from_seurat@reducedDims@listData[["UMAP"]] <- Embeddings(seu, "umap")
  # cds_from_seurat@reducedDims@listData[["PCA"]] <- Embeddings(seu, "pca")
  cds_from_seurat@preprocess_aux$gene_loadings <- Loadings(seu, "pca")

  cds_from_seurat <- learn_graph_by_resolution(cds_from_seurat, seu, resolution = resolution)

  return(cds_from_seurat)

}


#' Learn Monocle Graph by Resolution
#'
#' @param cds
#' @param resolution
#'
#' @return
#' @export
#'
#' @examples
learn_graph_by_resolution <- function(cds, seu, resolution = 1){

  ### Assign the cluster info
  if (any(grepl("integrated", names(cds@colData)))){
    default_assay = "integrated"
  } else {
    default_assay = "RNA"
  }

  cds <- monocle3::cluster_cells(cds)

  clusters <- seu[[paste0(default_assay, '_snn_res.', resolution)]]

  clusters <- purrr::set_names(clusters[[1]], rownames(clusters))

  cds <- assign_clusters_to_cds(cds, clusters = clusters)
  print("Learning graph, which can take a while depends on the sample")
  cds <- monocle3::learn_graph(cds, use_partition = T)

  return(cds)
}

#' Plot a Monocle Cell Data Set
#'
#' @param cds
#' @param resolution
#' @param color_cells_by
#'
#' @return
#' @export
#'
#' @examples
plot_cds <- function(cds, resolution, color_cells_by = "louvain_cluster"){

  key <- seq(1, length(colnames(cds)))
  cellid <- colnames(cds)

  cds[['key']] = key
  cds[['cellid']] = cellid

  if (any(grepl("integrated", names(cds@colData)))){
    default_assay = "integrated"
  } else {
    default_assay = "RNA"
  }

  if (color_cells_by == "louvain_cluster"){
    color_cells_by = paste0(default_assay, "_snn_res.", resolution)
  }


  cds_plot <- monocle3::plot_cells(cds,
                                   label_cell_groups = FALSE,
                                   label_groups_by_cluster = FALSE,
                                   label_leaves = FALSE,
                                   label_branch_points = FALSE,
                                   color_cells_by = color_cells_by) +
    # aes(key = key, cellid = cellid) +
    NULL


  plotly::ggplotly(cds_plot, height = 400) %>%
    # plotly::layout(dragmode = "lasso") %>%
    identity()


}

#' Plot a Monocle Cell Data Set
#'
#' @param cds
#' @param resolution
#' @param color_cells_by
#'
#' @return
#' @export
#'
#' @examples
plot_pseudotime <- function(cds, resolution, color_cells_by = "louvain_cluster"){

  key <- seq(1, length(colnames(cds)))
  cellid <- colnames(cds)

  cds[['key']] = key
  cds[['cellid']] = cellid

  if (any(grepl("integrated", colnames(cds@colData)))){
    default_assay = "integrated"
  } else {
    default_assay = "RNA"
  }

  if (color_cells_by == "louvain_cluster"){
    color_cells_by = paste0(default_assay, "_snn_res.", resolution)
  }

  cds_plot <- monocle3::plot_cells(cds,
                                   label_cell_groups = FALSE,
                                   label_groups_by_cluster = FALSE,
                                   label_leaves = FALSE,
                                   label_branch_points = FALSE,
                                   color_cells_by = color_cells_by) +
    # aes(key = key, cellid = cellid) +
    NULL

  print(cds_plot)

}


