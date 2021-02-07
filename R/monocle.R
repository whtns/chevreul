#' Convert a Seurat Object to a Monocle Cell Data Set
#'
#' @param seu
#'
#' @return
#' @export
#'
#' @examples
#' processed_seu <- clustering_workflow(human_gene_transcript_seu)
#' cds <- convert_seu_to_cds(processed_seu)
convert_seu_to_cds <- function(seu, resolution = 1) {
  print(resolution)

  # # drop sample_name metadata column as that is reserved by monocle3
  # seu$sample_name <- NULL
  #
  # cds <- SeuratWrappers::as.cell_data_set(seu) %>%
  #   monocle3::estimate_size_factors()
  #
  # rowData(cds)$gene_short_name <- rownames(cds)
  #
  # return(cds)

  # Building the necessary parts for a basic cds

  # part two, counts sparse matrix

  if ("integrated" %in% names(seu@assays)) {
    default_assay <- "integrated"
  } else {
    default_assay <- "gene"
  }

  DefaultAssay(seu) <- default_assay

  expression_matrix <- Seurat::GetAssayData(seu, slot = "data", assay = default_assay)

  count_matrix <- Seurat::GetAssayData(seu, slot = "counts", assay = "gene")

  count_matrix <- count_matrix[row.names(expression_matrix), ]
  count_matrix <- count_matrix[, Matrix::colSums(count_matrix) != 0]

  # part three, gene annotations

  gene_annotation <- data.frame(
    gene_short_name = rownames(count_matrix),
    row.names = rownames(count_matrix)
  )

  # part one, cell information
  cell_metadata <- seu[[]][colnames(count_matrix), ]

  # drop metadata column 'sample_name' for monocle plotting functions if present
  if (any(stringr::str_detect(colnames(cell_metadata), "sample_name"))) {
    cell_metadata <- subset(cell_metadata, select = -sample_name)
  }

  seu <- seu[, colnames(count_matrix)]

  ### Construct the basic cds object
  cds_from_seurat <- monocle3::new_cell_data_set(
    expression_data = count_matrix,
    cell_metadata = cell_metadata,
    gene_metadata = gene_annotation
  )

  cds_from_seurat <- cds_from_seurat[, colnames(seu)]

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


  cds_from_seurat <- monocle3::preprocess_cds(cds_from_seurat, method = "PCA", norm_method = "none")
  cds_from_seurat <- monocle3::reduce_dimension(cds_from_seurat, reduction_method = "UMAP")

  # reducedDim(cds_from_seurat, "PCA") <- Embeddings(seu, "pca")
  reducedDim(cds_from_seurat, "UMAP") <- Embeddings(seu, "umap")
  # cds_from_seurat@reducedDims@listData[["UMAP"]] <- Embeddings(seu, "umap")
  # cds_from_seurat@reducedDims@listData[["PCA"]] <- Embeddings(seu, "pca")
  cds_from_seurat@preprocess_aux$gene_loadings <- Loadings(seu, "pca")

  cds_from_seurat <- learn_graph_by_resolution(cds_from_seurat, seu, resolution = resolution)

  return(cds_from_seurat)
}

#' Assign Clusters to CDS
#'
#' @param cds
#' @param clusters
#'
#' @return
#' @export
#'
#' @examples
assign_clusters_to_cds <- function(cds, clusters) {
  clusters <- clusters[colnames(cds)]

  cds@clusters@listData[["UMAP"]][["clusters"]] <- clusters
  names(cds@clusters@listData[["UMAP"]][["clusters"]]) <- cds@colData@rownames

  return(cds)
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
learn_graph_by_resolution <- function(cds, seu, resolution = 1) {

  ### Assign the cluster info
  if (any(grepl("integrated", names(cds@colData)))) {
    default_assay <- "integrated"
  } else {
    default_assay <- "gene"
  }

  cds <- monocle3::cluster_cells(cds)

  clusters <- seu[[paste0(default_assay, "_snn_res.", resolution)]]

  clusters <- purrr::set_names(clusters[[1]], rownames(clusters))

  cds <- assign_clusters_to_cds(cds, clusters = clusters)
  print("Learning graph, which can take a while depends on the sample")
  cds <- monocle3::learn_graph(cds, use_partition = T)

  return(cds)
}

#' Plot a Monocle Cell Data Set
#'
#' @param cds
#' @param color_cells_by
#'
#' @return
#' @export
#'
#' @examples
plot_cds <- function(cds, color_cells_by = NULL, genes = NULL) {
  key <- colnames(cds)

  cds[["key"]] <- key

  if (any(grepl("integrated", names(cds@colData)))) {
    default_assay <- "integrated"
  } else {
    default_assay <- "gene"
  }

  # if (color_cells_by == "louvain_cluster"){
  #   color_cells_by = paste0(default_assay, "_snn_res.", resolution)
  # }


  cds_plot <- plot_cells(cds,
    genes = genes,
    label_cell_groups = FALSE,
    label_groups_by_cluster = FALSE,
    label_leaves = FALSE,
    label_branch_points = FALSE,
    color_cells_by = color_cells_by, key = key, cellid = key,
    cell_size = 0.75
  )
  NULL


  cds_plot <-
    cds_plot %>%
    plotly::ggplotly(height = 400) %>%
    plotly_settings() %>%
    plotly::toWebGL() %>%
    # plotly::partial_bundle() %>%
    identity()
}

#' Plot pseudotime on a Monocle Cell Data Set
#'
#' @param cds
#' @param resolution
#' @param color_cells_by
#'
#' @return
#' @export
#'
#' @examples
plot_pseudotime <- function(cds, resolution, color_cells_by = NULL, genes = NULL) {
  key <- seq(1, length(colnames(cds)))
  cellid <- colnames(cds)

  cds[["key"]] <- key
  cds[["cellid"]] <- cellid

  if (any(grepl("integrated", colnames(cds@colData)))) {
    default_assay <- "integrated"
  } else {
    default_assay <- "gene"
  }


  cds_plot <- monocle3::plot_cells(cds,
    show_trajectory_graph = TRUE,
    genes = genes,
    label_cell_groups = FALSE,
    label_groups_by_cluster = FALSE,
    label_leaves = FALSE,
    label_branch_points = FALSE,
    color_cells_by = color_cells_by,
    cell_size = 0.75
  ) +
    # aes(key = key, cellid = cellid) +
    NULL

  cds_plot <-
    cds_plot %>%
    plotly::ggplotly(height = 400) %>%
    plotly_settings() %>%
    plotly::toWebGL() %>%
    # plotly::partial_bundle() %>%
    identity()

  # print(cds_plot)
}

#' Plot feature expression on a Monocle Cell Data Set
#'
#' @param cds
#' @param resolution
#'
#' @return
#' @export
#'
#' @examples
plot_monocle_features <- function(cds, resolution, genes = NULL, ...) {
  key <- seq(1, length(colnames(cds)))
  cellid <- colnames(cds)

  cds[["key"]] <- key
  cds[["cellid"]] <- cellid

  if (any(grepl("integrated", colnames(cds@colData)))) {
    default_assay <- "integrated"
  } else {
    default_assay <- "gene"
  }

  cds_plot <- plot_cells(cds,
    genes = genes,
    label_cell_groups = FALSE,
    label_groups_by_cluster = FALSE,
    label_leaves = FALSE,
    label_branch_points = FALSE,
    cell_size = 0.75
  ) +
    # aes(key = key, cellid = cellid) +
    NULL

  cds_plot <-
    cds_plot %>%
    plotly::ggplotly(height = 400) %>%
    plotly::ggplotly(height = 400) %>%
    plotly_settings() %>%
    plotly::toWebGL() %>%
    # plotly::partial_bundle() %>%
    identity()

  # print(cds_plot)
}


#' Plot cells of a monocle cell data set
#'
#' @param cds
#' @param x
#' @param y
#' @param reduction_method
#' @param color_cells_by
#' @param group_cells_by
#' @param genes
#' @param show_trajectory_graph
#' @param trajectory_graph_color
#' @param trajectory_graph_segment_size
#' @param norm_method
#' @param label_cell_groups
#' @param label_groups_by_cluster
#' @param group_label_size
#' @param labels_per_group
#' @param label_branch_points
#' @param label_roots
#' @param label_leaves
#' @param graph_label_size
#' @param cell_size
#' @param cell_stroke
#' @param alpha
#' @param min_expr
#' @param rasterize
#' @param scale_to_range
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plot_cells <- function(cds, x = 1, y = 2, reduction_method = c(
                         "UMAP", "tSNE",
                         "PCA", "LSI", "Aligned"
                       ), color_cells_by = "cluster", group_cells_by = c(
                         "cluster",
                         "partition"
                       ), genes = NULL, show_trajectory_graph = TRUE,
                       trajectory_graph_color = "grey28", trajectory_graph_segment_size = 0.75,
                       norm_method = c("log", "size_only"), label_cell_groups = TRUE,
                       label_groups_by_cluster = TRUE, group_label_size = 2, labels_per_group = 1,
                       label_branch_points = TRUE, label_roots = TRUE, label_leaves = TRUE,
                       graph_label_size = 2, cell_size = 0.35, cell_stroke = I(cell_size / 2),
                       alpha = 1, min_expr = 0.1, rasterize = FALSE, scale_to_range = FALSE, ...) {
  reduction_method <- match.arg(reduction_method)
  assertthat::assert_that(methods::is(cds, "cell_data_set"))
  assertthat::assert_that(!is.null(reducedDims(cds)[[reduction_method]]),
    msg = paste(
      "No dimensionality reduction for", reduction_method,
      "calculated.", "Please run reduce_dimensions with",
      "reduction_method =", reduction_method, "before attempting to plot."
    )
  )
  low_dim_coords <- reducedDims(cds)[[reduction_method]]
  assertthat::assert_that(ncol(low_dim_coords) >= max(x, y),
    msg = paste(
      "x and/or y is too large. x and y must",
      "be dimensions in reduced dimension", "space."
    )
  )
  if (!is.null(color_cells_by)) {
    assertthat::assert_that(color_cells_by %in% c(
      "cluster",
      "partition", "pseudotime"
    ) | color_cells_by %in%
      names(colData(cds)), msg = paste(
      "color_cells_by must one of",
      "'cluster', 'partition', 'pseudotime,", "or a column in the colData table."
    ))
    if (color_cells_by == "pseudotime") {
      tryCatch(
        {
          pseudotime(cds, reduction_method = reduction_method)
        },
        error = function(x) {
          stop(paste(
            "No pseudotime for", reduction_method,
            "calculated. Please run order_cells with",
            "reduction_method =", reduction_method, "before attempting to color by pseudotime."
          ))
        }
      )
    }
  }
  assertthat::assert_that(!is.null(color_cells_by) || !is.null(markers),
    msg = paste(
      "Either color_cells_by or markers must",
      "be NULL, cannot color by both!"
    )
  )
  norm_method <- match.arg(norm_method)
  group_cells_by <- match.arg(group_cells_by)
  assertthat::assert_that(!is.null(color_cells_by) || !is.null(genes),
    msg = paste(
      "Either color_cells_by or genes must be",
      "NULL, cannot color by both!"
    )
  )
  if (show_trajectory_graph && is.null(monocle3::principal_graph(cds)[[reduction_method]])) {
    message("No trajectory to plot. Has learn_graph() been called yet?")
    show_trajectory_graph <- FALSE
  }
  gene_short_name <- NA
  sample_name <- NA
  data_dim_1 <- NA
  data_dim_2 <- NA
  if (rasterize) {
    plotting_func <- ggrastr::geom_point_rast
  }
  else {
    plotting_func <- ggplot2::geom_point
  }
  S_matrix <- reducedDims(cds)[[reduction_method]]
  data_df <- data.frame(S_matrix[, c(x, y)])
  colnames(data_df) <- c("data_dim_1", "data_dim_2")
  data_df$sample_name <- row.names(data_df)
  data_df <- as.data.frame(cbind(data_df, colData(cds)))
  if (group_cells_by == "cluster") {
    data_df$cell_group <- tryCatch(
      {
        clusters(cds, reduction_method = reduction_method)[data_df$sample_name]
      },
      error = function(e) {
        NULL
      }
    )
  }
  else if (group_cells_by == "partition") {
    data_df$cell_group <- tryCatch(
      {
        partitions(cds, reduction_method = reduction_method)[data_df$sample_name]
      },
      error = function(e) {
        NULL
      }
    )
  }
  else {
    stop("Error: unrecognized way of grouping cells.")
  }
  if (color_cells_by == "cluster") {
    data_df$cell_color <- tryCatch(
      {
        clusters(cds, reduction_method = reduction_method)[data_df$sample_name]
      },
      error = function(e) {
        NULL
      }
    )
  }
  else if (color_cells_by == "partition") {
    data_df$cell_color <- tryCatch(
      {
        partitions(cds, reduction_method = reduction_method)[data_df$sample_name]
      },
      error = function(e) {
        NULL
      }
    )
  }
  else if (color_cells_by == "pseudotime") {
    data_df$cell_color <- tryCatch(
      {
        pseudotime(cds, reduction_method = reduction_method)[data_df$sample_name]
      },
      error = function(e) {
        NULL
      }
    )
  }
  else {
    data_df$cell_color <- colData(cds)[
      data_df$sample_name,
      color_cells_by
    ]
  }
  if (show_trajectory_graph) {
    ica_space_df <- t(cds@principal_graph_aux[[reduction_method]]$dp_mst) %>%
      as.data.frame() %>%
      dplyr::select_(
        prin_graph_dim_1 = x,
        prin_graph_dim_2 = y
      ) %>%
      dplyr::mutate(
        sample_name = rownames(.),
        sample_state = rownames(.)
      )
    dp_mst <- cds@principal_graph[[reduction_method]]
    edge_df <- dp_mst %>%
      igraph::as_data_frame() %>%
      dplyr::select_(
        source = "from",
        target = "to"
      ) %>%
      dplyr::left_join(ica_space_df %>%
        dplyr::select_(
          source = "sample_name", source_prin_graph_dim_1 = "prin_graph_dim_1",
          source_prin_graph_dim_2 = "prin_graph_dim_2"
        ),
      by = "source"
      ) %>%
      dplyr::left_join(ica_space_df %>%
        dplyr::select_(
          target = "sample_name", target_prin_graph_dim_1 = "prin_graph_dim_1",
          target_prin_graph_dim_2 = "prin_graph_dim_2"
        ),
      by = "target"
      )
  }
  markers_exprs <- NULL
  expression_legend_label <- NULL
  if (!is.null(genes)) {
    if (!is.null(dim(genes)) && dim(genes) >= 2) {
      markers <- unlist(genes[, 1], use.names = FALSE)
    }
    else {
      markers <- genes
    }
    markers_rowData <- as.data.frame(subset(
      rowData(cds),
      gene_short_name %in% markers | row.names(rowData(cds)) %in%
        markers
    ))
    if (nrow(markers_rowData) == 0) {
      stop("None of the provided genes were found in the cds")
    }
    if (nrow(markers_rowData) >= 1) {
      cds_exprs <- SingleCellExperiment::counts(cds)[row.names(markers_rowData), ,
        drop = FALSE
      ]
      cds_exprs <- Matrix::t(Matrix::t(cds_exprs) / monocle3::size_factors(cds))
      if (!is.null(dim(genes)) && dim(genes) >= 2) {
        genes <- as.data.frame(genes)
        row.names(genes) <- genes[, 1]
        genes <- genes[row.names(cds_exprs), ]
        agg_mat <- as.matrix(monocle3::aggregate_gene_expression(cds,
          genes,
          norm_method = norm_method, scale_agg_values = FALSE
        ))
        if (dim(agg_mat)[2] == 1) agg_mat <- t(agg_mat)

        markers_exprs <- agg_mat
        markers_exprs <- reshape2::melt(markers_exprs)
        colnames(markers_exprs)[1:2] <- c(
          "feature_id",
          "cell_id"
        )
        if (is.factor(genes[, 2])) {
          markers_exprs$feature_id <- factor(markers_exprs$feature_id,
            levels = levels(genes[, 2])
          )
        }
        markers_exprs$feature_label <- markers_exprs$feature_id
        norm_method <- "size_only"
        expression_legend_label <- "Expression score"
      }
      else {
        cds_exprs@x <- round(10000 * cds_exprs@x) / 10000
        markers_exprs <- matrix(cds_exprs, nrow = nrow(markers_rowData))
        colnames(markers_exprs) <- colnames(SingleCellExperiment::counts(cds))
        row.names(markers_exprs) <- row.names(markers_rowData)
        markers_exprs <- reshape2::melt(markers_exprs)
        colnames(markers_exprs)[1:2] <- c(
          "feature_id",
          "cell_id"
        )
        markers_exprs <- merge(markers_exprs, markers_rowData,
          by.x = "feature_id", by.y = "row.names"
        )
        if (is.null(markers_exprs$gene_short_name)) {
          markers_exprs$feature_label <- as.character(markers_exprs$feature_id)
        }
        else {
          markers_exprs$feature_label <- as.character(markers_exprs$gene_short_name)
        }
        markers_exprs$feature_label <- ifelse(is.na(markers_exprs$feature_label) |
          !as.character(markers_exprs$feature_label) %in%
            markers, as.character(markers_exprs$feature_id),
        as.character(markers_exprs$feature_label)
        )
        markers_exprs$feature_label <- factor(markers_exprs$feature_label,
          levels = markers
        )
        if (norm_method == "size_only") {
          expression_legend_label <- "Expression"
        } else {
          expression_legend_label <- "log10(Expression)"
        }
      }
      if (scale_to_range) {
        markers_exprs <- dplyr::group_by(
          markers_exprs,
          feature_label
        ) %>%
          dplyr::mutate(
            max_val_for_feature = max(value),
            min_val_for_feature = min(value)
          ) %>%
          dplyr::mutate(value = 100 *
            (value - min_val_for_feature) / (max_val_for_feature -
              min_val_for_feature))
        expression_legend_label <- "% Max"
      }
    }
  }
  if (label_cell_groups && is.null(color_cells_by) == FALSE) {
    if (is.null(data_df$cell_color)) {
      if (is.null(genes)) {
        message(paste(
          color_cells_by, "not found in colData(cds), cells will",
          "not be colored"
        ))
      }
      text_df <- NULL
      label_cell_groups <- FALSE
    }
    else {
      if (is.character(data_df$cell_color) || is.factor(data_df$cell_color)) {
        if (label_groups_by_cluster && is.null(data_df$cell_group) ==
          FALSE) {
          text_df <- data_df %>%
            dplyr::group_by(cell_group) %>%
            dplyr::mutate(cells_in_cluster = dplyr::n()) %>%
            dplyr::group_by(cell_color, add = TRUE) %>%
            dplyr::mutate(per = dplyr::n() / cells_in_cluster)
          median_coord_df <- text_df %>% dplyr::summarize(
            fraction_of_group = dplyr::n(),
            text_x = stats::median(x = data_dim_1),
            text_y = stats::median(x = data_dim_2)
          )
          text_df <- suppressMessages(text_df %>% dplyr::select(per) %>%
            dplyr::distinct())
          text_df <- suppressMessages(dplyr::inner_join(
            text_df,
            median_coord_df
          ))
          text_df <- text_df %>%
            dplyr::group_by(cell_group) %>%
            dplyr::top_n(labels_per_group, per)
        }
        else {
          text_df <- data_df %>%
            dplyr::group_by(cell_color) %>%
            dplyr::mutate(per = 1)
          median_coord_df <- text_df %>% dplyr::summarize(
            fraction_of_group = dplyr::n(),
            text_x = stats::median(x = data_dim_1),
            text_y = stats::median(x = data_dim_2)
          )
          text_df <- suppressMessages(text_df %>% dplyr::select(per) %>%
            dplyr::distinct())
          text_df <- suppressMessages(dplyr::inner_join(
            text_df,
            median_coord_df
          ))
          text_df <- text_df %>%
            dplyr::group_by(cell_color) %>%
            dplyr::top_n(labels_per_group, per)
        }
        text_df$label <- as.character(text_df %>% dplyr::pull(cell_color))
      }
      else {
        message(paste(
          "Cells aren't colored in a way that allows them to",
          "be grouped."
        ))
        text_df <- NULL
        label_cell_groups <- FALSE
      }
    }
  }
  if (!is.null(markers_exprs) && nrow(markers_exprs) > 0) {
    data_df <- merge(data_df, markers_exprs,
      by.x = "sample_name",
      by.y = "cell_id"
    )
    data_df$value <- with(data_df, ifelse(value >= min_expr,
      value, NA
    ))
    na_sub <- data_df[is.na(data_df$value), ]
    if (norm_method == "size_only") {
      g <- ggplot(data = data_df, aes(
        x = data_dim_1,
        y = data_dim_2
      )) +
        plotting_func(aes(
          data_dim_1,
          data_dim_2
        ),
        size = I(cell_size), stroke = I(cell_stroke),
        color = "grey80", alpha = alpha, data = na_sub
        ) +
        plotting_func(aes(color = value),
          size = I(cell_size),
          stroke = I(cell_stroke), na.rm = TRUE
        ) +
        viridis::scale_color_viridis(
          option = "plasma",
          name = expression_legend_label, na.value = "grey80",
          end = 0.8, alpha = alpha
        ) +
        guides(alpha = FALSE) +
        facet_wrap(~feature_label)
    }
    else {
      g <- ggplot(data = data_df, aes(
        x = data_dim_1,
        y = data_dim_2
      )) +
        plotting_func(aes(
          data_dim_1,
          data_dim_2
        ),
        size = I(cell_size), stroke = I(cell_stroke),
        color = "grey80", data = na_sub, alpha = alpha
        ) +
        plotting_func(aes(color = log10(value + min_expr)),
          size = I(cell_size), stroke = I(cell_stroke),
          na.rm = TRUE, alpha = alpha
        ) +
        viridis::scale_color_viridis(
          option = "plasma",
          name = expression_legend_label, na.value = "grey80",
          end = 0.8, alpha = alpha
        ) +
        guides(alpha = FALSE) +
        facet_wrap(~feature_label)
    }
  }
  else {
    g <- ggplot(data = data_df, aes(x = data_dim_1, y = data_dim_2))
    if (color_cells_by %in% c("cluster", "partition")) {
      if (is.null(data_df$cell_color)) {
        g <- g + geom_point(
          color = I("gray"), size = I(cell_size),
          stroke = I(cell_stroke), na.rm = TRUE, alpha = I(alpha)
        )
        message(paste(
          "cluster_cells() has not been called yet, can't",
          "color cells by cluster"
        ))
      }
      else {
        g <- g + geom_point(aes(color = cell_color, ...),
          size = I(cell_size), stroke = I(cell_stroke),
          na.rm = TRUE, alpha = alpha
        )
      }
      g <- g + guides(color = guide_legend(
        title = color_cells_by,
        override.aes = list(size = 4)
      ))
    }
    else if (class(data_df$cell_color) == "numeric") {
      g <- g + geom_point(aes(color = cell_color, ...),
        size = I(cell_size),
        stroke = I(cell_stroke), na.rm = TRUE, alpha = alpha
      )
      g <- g + viridis::scale_color_viridis(
        name = color_cells_by,
        option = "C"
      )
    }
    else {
      g <- g + geom_point(aes(color = cell_color, ...),
        size = I(cell_size),
        stroke = I(cell_stroke), na.rm = TRUE, alpha = alpha
      )
      g <- g + guides(color = guide_legend(
        title = color_cells_by,
        override.aes = list(size = 4)
      ))
    }
  }
  if (show_trajectory_graph) {
    g <- g + geom_segment(aes_string(
      x = "source_prin_graph_dim_1",
      y = "source_prin_graph_dim_2", xend = "target_prin_graph_dim_1",
      yend = "target_prin_graph_dim_2"
    ),
    size = trajectory_graph_segment_size,
    color = I(trajectory_graph_color), linetype = "solid",
    na.rm = TRUE, data = edge_df
    )
    if (label_branch_points) {
      mst_branch_nodes <- branch_nodes(cds)
      branch_point_df <- ica_space_df %>%
        dplyr::slice(match(
          names(mst_branch_nodes),
          sample_name
        )) %>%
        dplyr::mutate(branch_point_idx = seq_len(dplyr::n()))
      g <- g + geom_point(aes_string(
        x = "prin_graph_dim_1",
        y = "prin_graph_dim_2"
      ),
      shape = 21, stroke = I(trajectory_graph_segment_size),
      color = "white", fill = "black", size = I(graph_label_size *
        1.5), na.rm = TRUE, branch_point_df
      ) + geom_text(aes_string(
        x = "prin_graph_dim_1",
        y = "prin_graph_dim_2", label = "branch_point_idx"
      ),
      size = I(graph_label_size), color = "white",
      na.rm = TRUE, branch_point_df
      )
    }
    if (label_leaves) {
      mst_leaf_nodes <- leaf_nodes(cds)
      leaf_df <- ica_space_df %>%
        dplyr::slice(match(
          names(mst_leaf_nodes),
          sample_name
        )) %>%
        dplyr::mutate(leaf_idx = seq_len(dplyr::n()))
      g <- g + geom_point(aes_string(
        x = "prin_graph_dim_1",
        y = "prin_graph_dim_2"
      ),
      shape = 21, stroke = I(trajectory_graph_segment_size),
      color = "black", fill = "lightgray", size = I(graph_label_size *
        1.5), na.rm = TRUE, leaf_df
      ) + geom_text(aes_string(
        x = "prin_graph_dim_1",
        y = "prin_graph_dim_2", label = "leaf_idx"
      ),
      size = I(graph_label_size), color = "black",
      na.rm = TRUE, leaf_df
      )
    }
    if (label_roots) {
      mst_root_nodes <- monocle3:::root_nodes(cds)
      root_df <- ica_space_df %>%
        dplyr::slice(match(
          names(mst_root_nodes),
          sample_name
        )) %>%
        dplyr::mutate(root_idx = seq_len(dplyr::n()))
      g <- g + geom_point(aes_string(
        x = "prin_graph_dim_1",
        y = "prin_graph_dim_2"
      ),
      shape = 21, stroke = I(trajectory_graph_segment_size),
      color = "black", fill = "white", size = I(graph_label_size *
        1.5), na.rm = TRUE, root_df
      ) + geom_text(aes_string(
        x = "prin_graph_dim_1",
        y = "prin_graph_dim_2", label = "root_idx"
      ),
      size = I(graph_label_size), color = "black",
      na.rm = TRUE, root_df
      )
    }
  }
  if (label_cell_groups) {
    g <- g + ggrepel::geom_text_repel(data = text_df, mapping = aes_string(
      x = "text_x",
      y = "text_y", label = "label"
    ), size = I(group_label_size))
    if (is.null(markers_exprs)) {
      g <- g + theme(legend.position = "none")
    }
  }
  g <- g + monocle3:::monocle_theme_opts() + xlab(paste(
    reduction_method,
    x
  )) + ylab(paste(reduction_method, y)) + theme(legend.key = element_blank()) +
    theme(panel.background = element_rect(fill = "white"))
  g
}


#' Title
#'
#' @param cds
#' @param cells
#' @param pr_deg_ids
#' @param seu_resolution
#' @param collapse_rows
#' @param collapse_cols
#' @param resolution
#' @param group.by
#' @param group.bar.height
#' @param cluster_columns
#' @param column_split
#' @param col_dendrogram
#' @param mm_col_dend
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
monocle_module_heatmap <- function(cds, pr_deg_ids, seu_resolution, cells = NULL, collapse_rows = TRUE,
                                   resolution = 10^seq(-6, -1), group.by = "batch", group.bar.height = 0.01,
                                   cluster_columns = FALSE, cluster_rows = TRUE,
                                   column_split = NULL, col_dendrogram = "ward.D2", mm_col_dend = 30, ...) {

  if (any(grepl("integrated", colnames(cds@colData)))) {
    default_assay <- "integrated"
  } else {
    default_assay <- "gene"
  }

  seu_resolution <- paste0(default_assay, "_snn_res.", seu_resolution)

  cds <- cds[pr_deg_ids, ]
  gene_module_df <- monocle3::find_gene_modules(cds, resolution = resolution) %>%
    dplyr::arrange(module)

  cell_group_df <- tibble::tibble(
    cell = row.names(colData(cds)),
    cell_group = colData(cds)[[seu_resolution]]
  ) %>%
    dplyr::mutate(cell_group = cell)

  if (collapse_rows != TRUE) {
    heatmap_row_df <-
      dplyr::select(gene_module_df, id) %>%
      dplyr::mutate(module = id)

    module_levels <- levels(gene_module_df$module)
    col <- scales::hue_pal()(length(module_levels))
    names(col) <- module_levels

    col <- list(module = col[gene_module_df$module])

    row_ha <- ComplexHeatmap::rowAnnotation(module = gene_module_df$module, col = col)
  } else {
    heatmap_row_df <- gene_module_df

    module_levels <- levels(gene_module_df$module)
    col <- scales::hue_pal()(length(module_levels))
    names(col) <- module_levels

    col <- list(module = col[gene_module_df$module])

    row_ha <- ComplexHeatmap::rowAnnotation(module = unique(gene_module_df$module), col = col)
  }

  agg_mat <- monocle3::aggregate_gene_expression(cds, heatmap_row_df)

  # reorder aggregation matrix by pseudotime
  col_order <- sort(monocle3::pseudotime(cds))
  agg_mat <- agg_mat[, names(col_order)]

  group.by <- group.by %||% "batch"
  cells <- cells %||% colnames(x = cds)

  pseudotime_tbl <-
    monocle3::pseudotime(cds) %>%
    tibble::enframe("sample_id", "pseudotime")

  groups.use <- colData(cds)[cells, group.by, drop = FALSE] %>%
    as.data.frame() %>%
    tibble::rownames_to_column("sample_id") %>%
    dplyr::mutate(across(where(is.character), as.factor)) %>%
    dplyr::left_join(pseudotime_tbl, by = "sample_id") %>%
    dplyr::arrange(pseudotime) %>%
    data.frame(row.names = 1) %>%
    identity()


  groups.use[is.na(groups.use)] <- min(groups.use$pseudotime, na.rm = TRUE)

  groups.use.factor <- groups.use[sapply(groups.use, is.factor)]
  ha_cols.factor <- NULL
  if (length(groups.use.factor) > 0) {
    ha_col_names.factor <- lapply(groups.use.factor, levels)
    ha_cols.factor <- purrr::map(ha_col_names.factor, ~ (scales::hue_pal())(length(.x))) %>%
      purrr::map2(ha_col_names.factor, set_names)
  }
  groups.use.numeric <- groups.use[sapply(groups.use, is.numeric)]
  ha_cols.numeric <- NULL
  if (length(groups.use.numeric) > 0) {
    numeric_col_fun <- function(myvec, color) {
      circlize::colorRamp2(range(myvec), c("white", color))
    }
    ha_col_names.numeric <- names(groups.use.numeric)
    ha_col_hues.numeric <- (scales::hue_pal())(length(ha_col_names.numeric))
    ha_cols.numeric <- purrr::map2(
      groups.use[ha_col_names.numeric],
      ha_col_hues.numeric, numeric_col_fun
    )
  }
  ha_cols <- c(ha_cols.factor, ha_cols.numeric)
  column_ha <- ComplexHeatmap::HeatmapAnnotation(
    df = groups.use,
    height = unit(group.bar.height, "points"), col = ha_cols
  )

  module_heatmap <- ComplexHeatmap::Heatmap(as.matrix(agg_mat),
    name = "log expression",
    top_annotation = column_ha,
    left_annotation = row_ha,
    cluster_columns = cluster_columns,
    cluster_rows = cluster_rows,
    show_column_names = FALSE,
    column_dend_height = unit(mm_col_dend, "mm"),
    # column_split = column_split,
    column_title = NULL,
    ...
  )

  # module_heatmap <- iheatmapr::iheatmap(as.matrix(agg_mat), col_labels = TRUE, row_labels = TRUE, cluster_rows = "hclust", cluster_cols = NULL)

  return(list(module_table = gene_module_df, module_heatmap = module_heatmap, agg_mat = agg_mat))
}

#' Flip Pseudotime
#'
#' @param cds
#'
#' @return
#' @export
#'
#' @examples
flip_pseudotime <- function(cds) {
  # pull original ptime
  orig_pseudotime <- monocle3::pseudotime(cds)
  orig_pseudotime[is.infinite(orig_pseudotime)] <- NA

  # sort ptime
  forward_pseudotime <- sort(orig_pseudotime[!is.na(orig_pseudotime)])

  # rev ptime and flip names
  rev_pseudotime <- rev(forward_pseudotime)
  names(rev_pseudotime) <- names(forward_pseudotime)

  rev_pseudotime <- c(rev_pseudotime, orig_pseudotime[is.na(orig_pseudotime)])

  # sort rev ptime by original order
  rev_pseudotime <- rev_pseudotime[names(orig_pseudotime)]

  # assign flipped ptime to cds
  cds@principal_graph_aux$UMAP$pseudotime <- rev_pseudotime

  return(cds)
}

export_pseudotime <- function(cds, root_cells) {
  root_cells <- root_cells %>%
    set_names(.) %>%
    tibble::enframe("sample_id", "root_cell") %>%
    dplyr::mutate(root_cell = 1)

  monocle_pt <- monocle3::pseudotime(cds) %>%
    tibble::enframe("sample_id", "pseudotime") %>%
    dplyr::arrange(pseudotime) %>%
    dplyr::left_join(root_cells, by = "sample_id")

  return(monocle_pt)
}
