#' Convert a Seurat Object to a Monocle Cell Data Set
#'
#' @param object
#'
#' @return
#'
#' @examples
#' processed_object <- clustering_workflow(human_gene_transcript_object)
#' cds <- convert_object_to_cds(processed_object)
convert_object_to_cds <- function(object, resolution = 1, min_expression = 0.05) {
    print(resolution)

    # Building the necessary parts for a basic cds

    # part two, counts sparse matrix

    if ("integrated" %in% names(object@assays)) {
        default_assay <- "integrated"
    } else {
        default_assay <- "gene"
    }

    DefaultAssay(object) <- default_assay

    expression_matrix <- Seurat::GetAssayData(object, slot = "data", assay = "gene")

    count_matrix <- Seurat::GetAssayData(object, slot = "counts", assay = "gene")

    count_matrix <- count_matrix[row.names(expression_matrix), ]
    count_matrix <- count_matrix[, Matrix::colSums(count_matrix) != 0]

    # part three, gene annotations

    gene_annotation <- data.frame(
        gene_short_name = rownames(count_matrix),
        row.names = rownames(count_matrix)
    )

    # part one, cell information
    cell_metadata <- get_cell_metadata(object)[colnames(count_matrix), ]

    # drop metadata column 'sample_name' for monocle plotting functions if present
    if (any(stringr::str_detect(colnames(cell_metadata), "sample_name"))) {
        cell_metadata <- subset(cell_metadata, select = -sample_name)
    }

    object <- object[, colnames(count_matrix)]

    ### Construct the basic cds object
    cds_from_object <- monocle3::new_cell_data_set(
        expression_data = count_matrix,
        cell_metadata = cell_metadata,
        gene_metadata = gene_annotation
    )

    cds_from_object <- cds_from_object[, colnames(object)]

    # estimate size factors
    cds_from_object <- cds_from_object[, colSums(as.matrix(monocle3::exprs(cds_from_object))) != 0]
    cds_from_object <- monocle3::estimate_size_factors(cds_from_object)


    ### Construct and assign the made up partition

    recreate.partition <- c(rep(1, length(cds_from_object@colData@rownames)))
    names(recreate.partition) <- cds_from_object@colData@rownames
    recreate.partition <- as.factor(recreate.partition)

    cds_from_object@clusters@listData[["UMAP"]][["partitions"]] <- recreate.partition

    ### Could be a space-holder, but essentially fills out louvain parameters
    cds_from_object@clusters@listData[["UMAP"]][["louvain_res"]] <- "NA"


    cds_from_object <- monocle3::preprocess_cds(cds_from_object, method = "PCA", norm_method = "none")
    cds_from_object <- monocle3::reduce_dimension(cds_from_object, reduction_method = "UMAP")

    # # reducedDim(cds_from_object, "PCA") <- Embeddings(object, "pca")
    reducedDim(cds_from_object, "UMAP") <- Embeddings(object, "umap")
    # # cds_from_object@reducedDims@listData[["UMAP"]] <- Embeddings(object, "umap")
    # # cds_from_object@reducedDims@listData[["PCA"]] <- Embeddings(object, "pca")
    cds_from_object@preprocess_aux$gene_loadings <- Loadings(object, "pca")

    # cds_from_object <- learn_graph_by_resolution(cds_from_object, object, resolution = resolution)

    return(cds_from_object)
}

#' Assign Clusters to CDS
#'
#' @param cds
#' @param clusters
#'
#' @return
#'
#' @examples
assign_clusters_to_cds <- function(cds, clusters) {
    clusters <- clusters[rownames(cds@colData)]

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
#'
#' @examples
learn_graph_by_resolution <- function(cds, object, resolution = 1) {
    ### Assign the cluster info
    if (any(grepl("integrated", names(cds@colData)))) {
        default_assay <- "integrated"
    } else {
        default_assay <- "gene"
    }

    cds <- monocle3::cluster_cells(cds)

    clusters <- object[[paste0(default_assay, "_snn_res.", resolution)]]

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
#'
#' @examples
plot_cds <- function(cds, color_cells_by = NULL, genes = NULL, return_plotly = TRUE) {
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

    if (return_plotly) {
        cds_plot <-
            cds_plot %>%
            ggplotly(height = 400) %>%
            plotly_settings() %>%
            toWebGL() %>%
            # partial_bundle() %>%
            identity()
    }

    return(cds_plot)
}

#' Plot pseudotime on a Monocle Cell Data Set
#'
#' @param cds
#' @param resolution
#' @param color_cells_by
#'
#' @return
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
        ggplotly(height = 400) %>%
        plotly_settings() %>%
        toWebGL() %>%
        # partial_bundle() %>%
        identity()

    # print(cds_plot)
}

#' Plot feature expression on a Monocle Cell Data Set
#'
#' @param cds
#' @param resolution
#'
#' @return
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
        ggplotly(height = 400) %>%
        ggplotly(height = 400) %>%
        plotly_settings() %>%
        toWebGL() %>%
        # partial_bundle() %>%
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
    } else {
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
    } else if (group_cells_by == "partition") {
        data_df$cell_group <- tryCatch(
            {
                partitions(cds, reduction_method = reduction_method)[data_df$sample_name]
            },
            error = function(e) {
                NULL
            }
        )
    } else {
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
    } else if (color_cells_by == "partition") {
        data_df$cell_color <- tryCatch(
            {
                partitions(cds, reduction_method = reduction_method)[data_df$sample_name]
            },
            error = function(e) {
                NULL
            }
        )
    } else if (color_cells_by == "pseudotime") {
        data_df$cell_color <- tryCatch(
            {
                pseudotime(cds, reduction_method = reduction_method)[data_df$sample_name]
            },
            error = function(e) {
                NULL
            }
        )
    } else {
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
            mutate(
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
            left_join(
                ica_space_df %>%
                    dplyr::select_(
                        source = "sample_name", source_prin_graph_dim_1 = "prin_graph_dim_1",
                        source_prin_graph_dim_2 = "prin_graph_dim_2"
                    ),
                by = "source"
            ) %>%
            left_join(
                ica_space_df %>%
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
        } else {
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
            } else {
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
                } else {
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
                markers_exprs <- group_by(
                    markers_exprs,
                    feature_label
                ) %>%
                    mutate(
                        max_val_for_feature = max(value),
                        min_val_for_feature = min(value)
                    ) %>%
                    mutate(value = 100 *
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
        } else {
            if (is.character(data_df$cell_color) || is.factor(data_df$cell_color)) {
                if (label_groups_by_cluster && is.null(data_df$cell_group) ==
                    FALSE) {
                    text_df <- data_df %>%
                        group_by(cell_group) %>%
                        mutate(cells_in_cluster = dplyr::n()) %>%
                        group_by(cell_color, add = TRUE) %>%
                        mutate(per = dplyr::n() / cells_in_cluster)
                    median_coord_df <- text_df %>% summarize(
                        fraction_of_group = dplyr::n(),
                        text_x = stats::median(x = data_dim_1),
                        text_y = stats::median(x = data_dim_2)
                    )
                    text_df <- suppressMessages(text_df %>% select(per) %>%
                        distinct())
                    text_df <- suppressMessages(dplyr::inner_join(
                        text_df,
                        median_coord_df
                    ))
                    text_df <- text_df %>%
                        group_by(cell_group) %>%
                        dplyr::top_n(labels_per_group, per)
                } else {
                    text_df <- data_df %>%
                        group_by(cell_color) %>%
                        mutate(per = 1)
                    median_coord_df <- text_df %>% summarize(
                        fraction_of_group = dplyr::n(),
                        text_x = stats::median(x = data_dim_1),
                        text_y = stats::median(x = data_dim_2)
                    )
                    text_df <- suppressMessages(text_df %>% select(per) %>%
                        distinct())
                    text_df <- suppressMessages(dplyr::inner_join(
                        text_df,
                        median_coord_df
                    ))
                    text_df <- text_df %>%
                        group_by(cell_color) %>%
                        dplyr::top_n(labels_per_group, per)
                }
                text_df$label <- as.character(text_df %>% pull(cell_color))
            } else {
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
                plotting_func(
                    aes(
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
        } else {
            g <- ggplot(data = data_df, aes(
                x = data_dim_1,
                y = data_dim_2
            )) +
                plotting_func(
                    aes(
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
    } else {
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
            } else {
                g <- g + geom_point(aes(color = cell_color, ...),
                    size = I(cell_size), stroke = I(cell_stroke),
                    na.rm = TRUE, alpha = alpha
                )
            }
            g <- g + guides(color = guide_legend(
                title = color_cells_by,
                override.aes = list(size = 4)
            ))
        } else if (class(data_df$cell_color) == "numeric") {
            g <- g + geom_point(aes(color = cell_color, ...),
                size = I(cell_size),
                stroke = I(cell_stroke), na.rm = TRUE, alpha = alpha
            )
            g <- g + viridis::scale_color_viridis(
                name = color_cells_by,
                option = "C"
            )
        } else {
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
        g <- g + geom_segment(
            aes_string(
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
                mutate(branch_point_idx = seq_len(dplyr::n()))
            g <- g + geom_point(
                aes_string(
                    x = "prin_graph_dim_1",
                    y = "prin_graph_dim_2"
                ),
                shape = 21, stroke = I(trajectory_graph_segment_size),
                color = "white", fill = "black", size = I(graph_label_size *
                    1.5), na.rm = TRUE, branch_point_df
            ) + geom_text(
                aes_string(
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
                mutate(leaf_idx = seq_len(dplyr::n()))
            g <- g + geom_point(
                aes_string(
                    x = "prin_graph_dim_1",
                    y = "prin_graph_dim_2"
                ),
                shape = 21, stroke = I(trajectory_graph_segment_size),
                color = "black", fill = "lightgray", size = I(graph_label_size *
                    1.5), na.rm = TRUE, leaf_df
            ) + geom_text(
                aes_string(
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
                mutate(root_idx = seq_len(dplyr::n()))
            g <- g + geom_point(
                aes_string(
                    x = "prin_graph_dim_1",
                    y = "prin_graph_dim_2"
                ),
                shape = 21, stroke = I(trajectory_graph_segment_size),
                color = "black", fill = "white", size = I(graph_label_size *
                    1.5), na.rm = TRUE, root_df
            ) + geom_text(
                aes_string(
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


#' Threshold monocle genes
#'
#' @param object
#' @param cds
#' @param min_expression
#'
#' @return
#'
#' @examples
threshold_monocle_genes <- function(object, cds, min_expression = 0.05) {
    agg_mat <- Seurat::GetAssayData(object, assay = "gene") %>%
        as.matrix()

    lgl_agg_mat <- agg_mat > min_expression

    percent_cells <- rowSums(lgl_agg_mat) / dim(lgl_agg_mat)[2]

    monocle3::fData(cds)$percent_cells <- percent_cells

    return(cds)
}

#' Find Modules from monocle
#'
#' Find modules in monocle cell data set
#'
#' @param cds
#' @param cells
#' @param pr_deg_ids
#' @param object_resolution
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
#' @importFrom circlize colorRamp2
#'
#' @return
monocle_module_heatmap <- function(cds, pr_deg_ids, object_resolution, cells = NULL, collapse_rows = TRUE,
    resolution = 10^seq(-6, -1), group.by = "batch", group.bar.height = 0.01,
    cluster_columns = FALSE, cluster_rows = TRUE,
    column_split = NULL, col_dendrogram = "ward.D2", mm_col_dend = 30, min_percent = 0.05, ...) {
    collapse_rows <- switch(collapse_rows,
        modules = TRUE,
        genes = FALSE
    )

    if (any(grepl("integrated", colnames(cds@colData)))) {
        default_assay <- "integrated"
    } else {
        default_assay <- "gene"
    }

    object_resolution <- paste0(default_assay, "_snn_res.", object_resolution)

    cds <- cds[pr_deg_ids, ]

    cds2 <- monocle3::preprocess_cds(cds) %>%
        monocle3::reduce_dimension(
            max_components = 2,
            reduction_method = "UMAP"
        )

    thresholded_genes <- monocle3::fData(cds) %>%
        tibble::as_tibble() %>%
        mutate(id = gene_short_name) %>%
        filter(percent_cells > 0.05)

    cds <- cds[rownames(cds) %in% thresholded_genes$id, ]

    gene_module_df <- monocle3::find_gene_modules(cds2, resolution = resolution) %>%
        arrange(module) %>%
        filter(id %in% rownames(cds)) %>%
        left_join(thresholded_genes, by = "id") %>%
        identity()

    cell_group_df <- tibble::tibble(
        cell = row.names(colData(cds)),
        cell_group = colData(cds)[[object_resolution]]
    ) %>%
        mutate(cell_group = cell)

    if (collapse_rows != TRUE) {
        heatmap_row_df <-
            select(gene_module_df, id) %>%
            mutate(module = id)

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
        enframe("sample_id", "pseudotime")

    groups.use <- colData(cds)[cells, group.by, drop = FALSE] %>%
        as.data.frame() %>%
        rownames_to_column("sample_id") %>%
        mutate(across(where(is.character), as.factor)) %>%
        left_join(pseudotime_tbl, by = "sample_id") %>%
        arrange(pseudotime) %>%
        data.frame(row.names = 1) %>%
        identity()


    groups.use[is.na(groups.use)] <- min(groups.use$pseudotime, na.rm = TRUE)

    groups.use.factor <- groups.use[sapply(groups.use, is.factor)]
    ha_cols.factor <- NULL
    if (length(groups.use.factor) > 0) {
        ha_col_names.factor <- lapply(groups.use.factor, levels)
        ha_cols.factor <- map(ha_col_names.factor, ~ (scales::hue_pal())(length(.x))) %>%
            purrr::map2(ha_col_names.factor, set_names)
    }
    groups.use.numeric <- groups.use[sapply(groups.use, is.numeric)]
    ha_cols.numeric <- NULL
    if (length(groups.use.numeric) > 0) {
        numeric_col_fun <- function(myvec, color) {
            colorRamp2(range(myvec), c("white", color))
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
#' @param cds a cell data set object
#'
#' @return
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
        enframe("sample_id", "root_cell") %>%
        mutate(root_cell = 1)

    monocle_pt <- monocle3::pseudotime(cds) %>%
        enframe("sample_id", "pseudotime") %>%
        arrange(pseudotime) %>%
        left_join(root_cells, by = "sample_id")

    return(monocle_pt)
}

#' Monocle UI Module
#'
#' @param id
#'
#' @noRd
monocleui <- function(id) {
  ns <- NS(id)
  tagList(
    fluidRow(
      chevreulBox(
        title = "Seurat Data",
        plotlyOutput(ns("objectdimplot"), height = 500),
        width = 6
        # plotDimRedui(ns("plotdimred")
      ),
      chevreulBox(
        title = "Pobjectdotime Settings",
        actionButton(ns("subsetSeurat"), "Subset Seurat before Pobjectdotime Calculation"),
        actionButton(ns("calcCDS"), "Calculate Pobjectdotime"),
        sliderInput(ns("cdsResolution"), "Resolution of clustering algorithm (affects number of clusters)",
                    min = 0.2, max = 2, step = 0.2, value = 0.6
        ),
        actionButton(ns("subsetCells"), "Subset Monocle Object After Pobjectdotime Calculation"),
        uiOutput(ns("rootCellsui")),
        actionButton(ns("plotPobjectdotime"), "Calculate Pobjectdotime With Root Cells"),
        downloadButton(ns("downloadPT"), "Export Pobjectdotime"),
        checkboxInput(ns("flipPtime"), "Invert Pobjectdotime", value = TRUE),
        width = 6
      )
    ),
    chevreulBox(
      title = "Embedding Plot",
      selectizeInput(ns("plottype1"), "Variable to Plot", choices = c(Louvain = "louvain"), selected = "Louvain", multiple = TRUE),
      selectizeInput(ns("customFeature1"), "Gene or transcript expression by which to color the plot",
                     choices = NULL, multiple = FALSE
      ),
      uiOutput(ns("moduleSelect1")),
      plotlyOutput(ns("monoclePlot1")),
      width = 6
    ),
    chevreulBox(
      title = "Embedding Plot",
      selectizeInput(ns("plottype2"), "Variable to Plot", choices = c(Louvain = "louvain"), selected = "Louvain", multiple = TRUE),
      selectizeInput(ns("customFeature2"), "gene or transcript on which to color the plot",
                     choices = NULL, multiple = FALSE
      ),
      uiOutput(ns("moduleSelect2")),
      plotlyOutput(ns("monoclePlot2")),
      width = 6
    ),
    fluidRow(
      chevreulBox(
        title = "calculate pseudotime",
        radioButtons(ns("diffexFeature"), "Feature for differential expression", choices = c("gene", "transcript")),
        actionButton(ns("calcPtimeGenes"), "Find Pobjectdotime Correlated Genes"),
        sliderInput(ns("qvalThreshold"), "Set q value threshold for module calculation", min = 0.01, 0.1, value = 0.05, step = 0.01),
        textOutput("pseudotimeMessages"),
        uiOutput(ns("partitionSelect")),
        uiOutput(ns("genePlotQuery2")),
        DT::DTOutput(ns("ptimeGenesDT")),
        downloadButton(ns("downloadGenesDT"), "Download data as csv"),
        # uiOutput(ns("ptimeGenes")),
        width = 6
      ),
      chevreulBox(
        title = "Plot Feature Expression over Pobjectdotime",
        plotlyOutput(ns("ptimeGenesLinePlot")),
        width = 6,
        height = 650
      )
    ),
    chevreulBox(
      title = "Heatmap",
      uiOutput(ns("colAnnoVarui")),
      radioButtons(ns("heatmapRows"), "annotate heatmap rows by genes or modules?", choices = c("modules", "genes")),
      downloadButton(ns("downloadPlot"), "Download Heatmap"),
      downloadButton(ns("downloadCds"), "Download celldataset"),
      plotOutput(ns("monocleHeatmap"), width = "800px", height = "1200px")
    ),
    chevreulBox(
      title = "Modules",
      plotOutput(ns("modulePlot")),
      div(DT::dataTableOutput(ns("moduleTable")), style = "font-size: 75%")
    )
  )
}

#' Monocle Server Module
#'
#' @param input
#' @param output
#' @param session
#' @param cds
#' @param object
#' @param plot_types
#' @param resolution
#'
#' @noRd
monocle <- function(input, output, session, object, plot_types, featureType,
                    organism_type, reductions) {
  ns <- session$ns

  # markermarker
  w <- waiter::Waiter$new(ns("monocleHeatmap"),
                          html = waiter::spin_loaders(id = 1, color = "black", style = "position:relative;margin:auto;"),
                          color = waiter::transparent(.5)
  )

  output$colAnnoVarui <- renderUI({
    req(object())

    selectizeInput(ns("colAnnoVar"), "Column Annotation(s)",
                   choices = colnames(get_cell_metadata(object())), selected = "batch", multiple = TRUE
    )
  })

  cds_rvs <- reactiveValues(selected = c(traj = TRUE, ptime = FALSE, diff_features = FALSE))
  cds_plot_types <- reactiveVal(c(Pobjectdotime = "pseudotime", Module = "module"))
  myplot_types <- reactive({
    c(purrr::flatten_chr(plot_types()), cds_plot_types())
  })

  # to be able to subset, create a new copy of the object

  object_monocle <- reactiveVal()

  observe({
    req(object())
    object_monocle(object())
  })

  louvain_resolution <- reactive({
    if (query_assay(object(), "integrated")) {
      assay <- "integrated"
    } else {
      assay <- "gene"
    }

    paste0(assay, "_snn_res.", input$cdsResolution)
  })

  objectdimplot <- reactive({
    req(object_monocle())

    plot_var(object_monocle(), embedding = "umap", group = louvain_resolution(), return_plotly = TRUE)
  })

  output$objectdimplot <- renderPlotly({
    objectdimplot()
  })

  # callModule(plotDimRed, "plotdimred", object, plot_types, featureType,
  #            organism_type, reductions)


  observeEvent(input$subsetSeurat, {
    req(object_monocle())

    d <- event_data("plotly_selected", priority = "event")
    if (is.null(d)) {
      msg <- "Click and drag events (i.e. select/lasso) appear here (double-click to clear)"
      print(d)
    } else {
      print(d$key)
      print(d)
      subset_monocle <- object_monocle()[, d$key]
      object_monocle(subset_monocle)
    }
  })

  observeEvent(input$calcCDS, {
    req(object_monocle())
    cds_rvs$selected <- c(traj = TRUE, ptime = FALSE, diff_features = FALSE)
    cds <- convert_object_to_cds(object(), resolution = input$cdsResolution)
    # cds <- convert_object_to_cds(object_monocle(), resolution = input$cdsResolution)
    cds <- cds[, colnames(cds) %in% colnames(object_monocle())]

    cds <- threshold_monocle_genes(object_monocle(), cds)

    cds <- learn_graph_by_resolution(cds, object_monocle(),
                                     resolution = input$cdsResolution
    )
    updateSelectizeInput(session, "plottype1", selected = "louvain", choices = myplot_types())
    updateSelectizeInput(session, "customFeature1", choices = rownames(cds), server = TRUE)
    updateSelectizeInput(session, "plottype2", selected = "louvain", choices = myplot_types())
    updateSelectizeInput(session, "customFeature2", choices = rownames(cds), server = TRUE)
    cds_rvs$traj <- cds
  })

  selected_plot <- reactiveVal()

  output$monoclePlot1 <- renderPlotly({
    req(input$plottype1)
    req(cds_rvs$traj)
    w$show()
    print(cds_rvs$selected)
    if (input$plottype1 == "louvain") {
      cluster_resolution <- reactive({
        if (any(stringr::str_detect(colnames(colData(cds_rvs$traj)), "integrated"))) {
          paste0("integrated", "_snn_res.", input$cdsResolution)
        } else {
          paste0("gene", "_snn_res.", input$cdsResolution)
        }
      })
      plot_cds(cds_rvs$traj, color_cells_by = cluster_resolution())
    } else if (input$plottype1 == "pseudotime") {
      plot_pseudotime(cds_rvs$traj, color_cells_by = "pseudotime", resolution = input$cdsResolution)
    } else if (input$plottype1 == "feature") {
      plot_monocle_features(cds_rvs$traj, genes = input$customFeature1, monocle_heatmap()$agg_mat)
    } else if (input$plottype1 == "module") {
      print(monocle_heatmap()$module_table)
      print(input$plotModule1)
      genes <- monocle_heatmap()$module_table %>%
        filter(module %in% input$plotModule1) %>%
        dplyr::mutate(module = factor(module))
      plot_monocle_features(cds_rvs$traj, genes = genes, monocle_heatmap()$agg_mat)
    } else {
      plot_cds(cds_rvs$traj, color_cells_by = input$plottype1)
    }
  })

  output$monoclePlot2 <- renderPlotly({
    req(input$plottype2)
    req(cds_rvs$traj)
    w$show()
    print(cds_rvs$selected)
    if (input$plottype2 == "louvain") {
      cluster_resolution <- reactive({
        if (any(stringr::str_detect(colnames(colData(cds_rvs$traj)), "integrated"))) {
          paste0("integrated", "_snn_res.", input$cdsResolution)
        } else {
          paste0("gene", "_snn_res.", input$cdsResolution)
        }
      })
      plot_cds(cds_rvs$traj, color_cells_by = cluster_resolution())
    } else if (input$plottype2 == "pseudotime") {
      plot_pseudotime(cds_rvs$traj, color_cells_by = "pseudotime", resolution = input$cdsResolution)
    } else if (input$plottype2 == "feature") {
      plot_monocle_features(cds_rvs$traj, genes = input$customFeature2, monocle_heatmap()$agg_mat)
    } else if (input$plottype2 == "module") {
      print(monocle_heatmap()$module_table)
      print(input$plotModule2)

      genes <- monocle_heatmap()$module_table %>%
        filter(module %in% input$plotModule2) %>%
        dplyr::mutate(module = factor(module))
      plot_monocle_features(cds_rvs$traj, genes = genes, monocle_heatmap()$agg_mat)
    } else {
      plot_cds(cds_rvs$traj, color_cells_by = input$plottype2)
    }
  })

  cdsbrush <- reactive({
    req(cds_rvs$traj)
    d <- event_data("plotly_selected")
    if (is.null(d)) {
      msg <- "Click and drag events (i.e. select/lasso) appear here (double-click to clear)"
      return(d)
    } else {
      # selected_cells <- colnames(cds_rvs$traj)[as.numeric(d$key)]
      d$key
    }
  })

  observeEvent(input$subsetCells, {
    req(cds_rvs$traj)
    print(cdsbrush())
    cds_rvs$traj <- cds_rvs$traj[, cdsbrush()]
  })

  output$rootCellsui <- renderUI({
    selectizeInput(ns("rootCells"), "Choose Root Cells", choices = c("Choose Root Cells" = "", colnames(cds_rvs$traj)), multiple = TRUE)
  })

  exported_pseudotime <- reactiveVal()

  observeEvent(input$plotPobjectdotime, {
    req(cds_rvs$traj)
    req(input$rootCells)
    cds_rvs$traj <- monocle3::order_cells(cds_rvs$traj, root_cells = input$rootCells)

    # # select only first partition
    # cds_rvs$traj <- cds_rvs$traj[, monocle3::partitions(cds_rvs$traj) == 1]


    if (input$flipPtime) {
      cds_rvs$traj <- flip_pseudotime(cds_rvs$traj)
    }
    updateSelectizeInput(session, "plottype1", selected = "pseudotime", choices = myplot_types())
    updateSelectizeInput(session, "plottype2", selected = "pseudotime", choices = myplot_types())
    cds_rvs$selected <- c(traj = FALSE, ptime = TRUE, diff_features = FALSE)
    # markermarker

    exported_pseudotime(export_pseudotime(cds_rvs$traj, input$rootCells))
  })


  output$downloadPT <- downloadHandler(
    filename = function() {
      paste("pseudotime-", Sys.Date(), ".csv", sep = "")
    },
    content = function(file) {
      readr::write_csv(exported_pseudotime(), file)
    }
  )

  # markermarker
  observeEvent(input$calcPtimeGenes, {
    req(input$diffexFeature)
    if (req(cds_rvs$selected["ptime"])) {
      # markermarker
      # cds_rvs$traj  <- swap_counts_from_feature(cds_rvs$traj, input$diffexFeature)

      showModal(modalDialog(
        title = "Calculating Pobjectdotime Correlated Features",
        "This may take a few minutes!"
      ))
      cds_rvs$traj@metadata[["diff_features"]] <- monocle3::graph_test(cds_rvs$traj, neighbor_graph = "principal_graph", cores = 4, expression_family = "negbinom")

      cds_rvs$selected <- c(traj = FALSE, ptime = FALSE, diff_features = TRUE)

      removeModal()
    }
  })

  cds_pr_test_res <- reactive({
    if (req(cds_rvs$selected["diff_features"])) {
      cds_rvs$traj@metadata$diff_features %>%
        # subset(q_value < 0.05) %>%
        arrange(q_value) %>%
        dplyr::select(-status) %>%
        # dplyr::filter %>%
        identity()
    }
  })

  observe({
    req(cds_pr_test_res())
    if (req(cds_rvs$selected["diff_features"])) {
      output$genePlotQuery2 <- renderUI({
        selectizeInput(ns("genePlotQuery1"), "Pick Gene to Plot on Pobjectdotime", choices = rownames(cds_pr_test_res()), multiple = TRUE, selected = rownames(cds_pr_test_res())[1])
      })

      output$partitionSelect <- renderUI({
        selectizeInput(ns("partitions"), "Select a Partition to Plot", choices = levels(monocle3::partitions(cds_rvs$traj)), multiple = FALSE)
      })
    }
  })

  observe({
    req(cds_pr_test_res())
    req(input$genePlotQuery1)
    if (req(cds_rvs$selected["diff_features"])) {
      output$ptimeGenesLinePlot <- renderPlotly({
        genes_in_pseudotime <- prep_plot_genes_in_pseudotime(cds_rvs$traj, input$genePlotQuery1, input$cdsResolution)
        genes_in_pseudotime <-
          genes_in_pseudotime %>%
          ggplotly(height = 600) %>%
          plotly_settings() %>%
          toWebGL() %>%
          # partial_bundle() %>%
          identity()
      })

      output$ptimeGenesDT <- DT::renderDT({
        DT::datatable(cds_pr_test_res(),
                      extensions = "Buttons",
                      options = list(dom = "Bftp", buttons = c("copy", "csv"), scrollX = "100px", scrollY = "400px", pageLength = 200, paging = TRUE)
        )
      })

      output$downloadGenesDT <- downloadHandler(
        filename = function() {
          paste("diffex_ptime-", Sys.Date(), ".csv", sep = "")
        },
        content = function(file) {
          write.csv(cds_pr_test_res(), file)
        }
      )
    }
  })

  monocle_heatmap <- reactive({
    req(cds_rvs$traj)
    req(input$colAnnoVar)

    heatmap_genes <- cds_pr_test_res() %>%
      filter(q_value < input$qvalThreshold)

    monocle_module_heatmap(cds_rvs$traj, rownames(heatmap_genes), input$cdsResolution, collapse_rows = input$heatmapRows, group.by = input$colAnnoVar)
  }) %>%
    bindCache(cds_rvs$traj, input$cdsResolution, input$heatmapRows)

  module_choices <- reactive({
    module_choices <- as.character(unique(monocle_heatmap()$module_table$module))
    # names(module_choices) <- paste("Module", module_choices)
  })

  output$moduleSelect1 <- renderUI({
    selectizeInput(ns("plotModule1"), "gene module to plot (if computed)", choices = module_choices(), multiple = TRUE)
  })
  output$moduleSelect2 <- renderUI({
    selectizeInput(ns("plotModule2"), "gene module to plot (if computed)", choices = module_choices(), multiple = TRUE)
  })

  observe({
    output$monocleHeatmap <- renderPlot({
      monocle_heatmap()$module_heatmap
    })

    output$moduleTable <- DT::renderDT({
      DT::datatable(monocle_heatmap()$module_table,
                    extensions = "Buttons",
                    options = list(dom = "Bft", buttons = c(
                      "copy",
                      "csv"
                    ), scrollX = "100px", scrollY = "400px")
      )
    })

    output$modulePlot <- renderPlot({
      ggplot(monocle_heatmap()$module_table, aes(dim_1, dim_2, color = module)) +
        geom_point()
    })
  })

  output$downloadPlot <- downloadHandler(
    filename = function() {
      paste("heatmap", ".pdf", sep = "")
    },
    content = function(file) {
      ggsave(file, ggplotify::as.ggplot(monocle_heatmap()$module_heatmap), width = 16, height = 12)
    }
  )

  output$downloadCds <- downloadHandler(
    filename = function() {
      paste("cds", ".rds", sep = "")
    },
    content = function(file) {
      saveRDS(cds_rvs$traj, file)
    }
  )
}

#' Prep plot genes in pseudotime
#'
#' @param cds
#' @param mygenes
#' @param resolution
#'
#' @return
#'
#' @examples
prep_plot_genes_in_pseudotime <- function(cds, mygenes, resolution, partition = FALSE) {
  if (partition) {
    partition_cells <- monocle3::partitions(cds)
    # partition_cells <-  split(names(partition_cells), partition_cells)[[input$partitions]]
    partition_cells <- split(names(partition_cells), partition_cells)[[1]]

    cds <- cds[, colnames(cds) %in% partition_cells]
  }

  cds <- cds[rownames(cds) %in% mygenes, ]

  if (any(grepl("integrated", colnames(colData(cds))))) {
    default_assay <- "integrated"
  } else {
    default_assay <- "gene"
  }

  color_cells_by <- paste0(default_assay, "_snn_res.", resolution)

  gene_ptime_plot <- monocle3::plot_genes_in_pseudotime(cds,
                                                        color_cells_by = color_cells_by,
                                                        min_expr = 0.5
  )

  return(gene_ptime_plot)
}

#' Swap counts from Feature
#'
#' @param cds
#' @param featureType
#'
#' @return
#'
#' @examples
swap_counts_from_feature <- function(cds, featureType) {
  print(featureType)
  #
  #   if (featureType == "transcript"){
  #     rowData(cds[[featureType]])$gene_short_name <- rownames(cds[[featureType]])
  #   }

  assay(cds$traj, withDimnames = FALSE) <- assay(cds[[featureType]])
  rowData(cds$traj) <- rowData(cds[[featureType]])
  rownames(cds$traj) <- rownames(cds[[featureType]])
  cds$traj@preprocess_aux$gene_loadings <- cds[[featureType]]@preprocess_aux$gene_loadings
  # counts(cds$traj) <- counts(cds[[featureType]])
  cds$traj
}
