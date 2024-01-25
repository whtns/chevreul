#' Unite metadata
#'
#'
#'
#' @param seu A seurat object
#' @param metavars A feature or variable to combine
#'
#' @return
#' @export
#'
#' @examples
unite_metadata <- function(seu, metavars) {
    newcolname <- paste(metavars, collapse = "_by_")

    newdata <- seu[[metavars]] %>%
        tidyr::unite(!!newcolname, metavars) %>%
        tibble::deframe()

    Idents(seu) <- newdata

    return(seu)
}

#' Plot monocle pseudotime over multiple branches
#'
#' Plots heatmap to de
#'
#' @param cds CellDataSet for the experiment
#' @param branches The terminal branches on the developmental tree to be investigated.
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
    scale_max = 3,
    scale_min = -3,
    trend_formula = "~sm.ns(Pseudotime, df=3)",
    return_heatmap = FALSE,
    cores = 1) {
    pseudocount <- 1
    if (!(all(branches %in% Biobase::pData(cds)$State)) & length(branches) == 1) {
        stop("This function only allows to make multiple branch plots where branches is included in the pData")
    }

    branch_label <- branches
    if (!is.null(branches_name)) {
        if (length(branches) != length(branches_name)) {
            stop("branches_name should have the same length as branches")
        }
        branch_label <- branches_name
    }

    # test whether or not the states passed to branches are true branches (not truncks) or there are terminal cells
    g <- cds@minSpanningTree
    m <- NULL
    # branche_cell_num <- c()
    for (branch_in in branches) {
        branches_cells <- row.names(subset(Biobase::pData(cds), State == branch_in))
        root_state <- subset(Biobase::pData(cds), Pseudotime == 0)[, "State"]
        root_state_cells <- row.names(subset(Biobase::pData(cds), State == root_state))

        if (cds@dim_reduce_type != "ICA") {
            root_state_cells <- unique(paste("Y_", cds@auxOrderingData$DDRTree$pr_graph_cell_proj_closest_vertex[root_state_cells, ], sep = ""))
            branches_cells <- unique(paste("Y_", cds@auxOrderingData$DDRTree$pr_graph_cell_proj_closest_vertex[branches_cells, ], sep = ""))
        }
        root_cell <- root_state_cells[which(degree(g, v = root_state_cells) == 1)]
        tip_cell <- branches_cells[which(degree(g, v = branches_cells) == 1)]

        traverse_res <- traverseTree(g, root_cell, tip_cell)
        path_cells <- names(traverse_res$shortest_path[[1]])

        if (cds@dim_reduce_type != "ICA") {
            pc_ind <- cds@auxOrderingData$DDRTree$pr_graph_cell_proj_closest_vertex
            path_cells <- row.names(pc_ind)[paste("Y_", pc_ind[, 1], sep = "") %in% path_cells]
        }

        cds_subset <- cds[, path_cells]

        newdata <- data.frame(Pseudotime = seq(0, max(Biobase::pData(cds_subset)$Pseudotime), length.out = 100))

        tmp <- genSmoothCurves(cds_subset,
            cores = cores, trend_formula = trend_formula,
            relative_expr = T, new_data = newdata
        )
        if (is.null(m)) {
            m <- tmp
        } else {
            m <- cbind(m, tmp)
        }
    }

    # remove genes with no expression in any condition
    m <- m[!apply(m, 1, sum) == 0, ]

    norm_method <- match.arg(norm_method)

    # FIXME: this needs to check that vst values can even be computed. (They can only be if we're using NB as the expressionFamily)
    if (norm_method == "vstExprs" && is.null(cds@dispFitInfo[["blind"]]$disp_func) == FALSE) {
        m <- vstExprs(cds, expr_matrix = m)
    } else if (norm_method == "log") {
        m <- log10(m + pseudocount)
    }

    # Row-center the data.
    m <- m[!apply(m, 1, sd) == 0, ]
    m <- Matrix::t(scale(Matrix::t(m), center = TRUE))
    m <- m[is.na(row.names(m)) == FALSE, ]
    m[is.nan(m)] <- 0
    m[m > scale_max] <- scale_max
    m[m < scale_min] <- scale_min

    heatmap_matrix <- m

    row_dist <- as.dist((1 - cor(Matrix::t(heatmap_matrix))) / 2)
    row_dist[is.na(row_dist)] <- 1

    if (is.null(hmcols)) {
        bks <- seq(-3.1, 3.1, by = 0.1)
        hmcols <- colorRamps::blue2green2red(length(bks) - 1)
    } else {
        bks <- seq(-3.1, 3.1, length.out = length(hmcols))
    }

    ph <- pheatmap(heatmap_matrix,
        useRaster = T,
        cluster_cols = FALSE,
        cluster_rows = T,
        show_rownames = F,
        show_colnames = F,
        clustering_distance_rows = row_dist,
        clustering_method = hclust_method,
        cutree_rows = num_clusters,
        silent = TRUE,
        filename = NA,
        breaks = bks,
        color = hmcols
    )

    annotation_col <- data.frame(Branch = factor(rep(rep(branch_label, each = 100))))
    annotation_row <- data.frame(Cluster = factor(cutree(ph$tree_row, num_clusters)))
    col_gaps_ind <- c(1:(length(branches) - 1)) * 100

    if (!is.null(add_annotation_row)) {
        old_colnames_length <- ncol(annotation_row)
        annotation_row <- cbind(annotation_row, add_annotation_row[row.names(annotation_row), ])
        colnames(annotation_row)[(old_colnames_length + 1):ncol(annotation_row)] <- colnames(add_annotation_row)
        # annotation_row$bif_time <- add_annotation_row[as.character(Biobase::fData(absolute_cds[row.names(annotation_row), ])$gene_short_name), 1]
    }


    if (use_gene_short_name == TRUE) {
        if (is.null(Biobase::fData(cds)$gene_short_name) == FALSE) {
            feature_label <- as.character(Biobase::fData(cds)[row.names(heatmap_matrix), "gene_short_name"])
            feature_label[is.na(feature_label)] <- row.names(heatmap_matrix)

            row_ann_labels <- as.character(Biobase::fData(cds)[row.names(annotation_row), "gene_short_name"])
            row_ann_labels[is.na(row_ann_labels)] <- row.names(annotation_row)
        } else {
            feature_label <- row.names(heatmap_matrix)
            row_ann_labels <- row.names(annotation_row)
        }
    } else {
        feature_label <- row.names(heatmap_matrix)
        row_ann_labels <- row.names(annotation_row)
    }

    row.names(heatmap_matrix) <- feature_label
    row.names(annotation_row) <- row_ann_labels


    colnames(heatmap_matrix) <- c(1:ncol(heatmap_matrix))

    if (!(cluster_rows)) {
        annotation_row <- NA
    }

    ph_res <- pheatmap(heatmap_matrix[, ], # ph$tree_row$order
        useRaster = T,
        cluster_cols = FALSE,
        cluster_rows = cluster_rows,
        show_rownames = show_rownames,
        show_colnames = F,
        # scale="row",
        clustering_distance_rows = row_dist, # row_dist
        clustering_method = hclust_method, # ward.D2
        cutree_rows = num_clusters,
        # cutree_cols = 2,
        annotation_row = annotation_row,
        annotation_col = annotation_col,
        gaps_col = col_gaps_ind,
        treeheight_row = 20,
        breaks = bks,
        fontsize = 12,
        color = hmcols,
        silent = TRUE,
        border_color = NA,
        filename = NA
    )

    grid::grid.rect(gp = grid::gpar("fill", col = NA))
    grid::grid.draw(ph_res$gtable)
    if (return_heatmap) {
        return(ph_res)
    }
}


#' Plot Metadata Variables
#'
#' Plots static or interactive plot where each point represents a cell metadata
#' variable whose position on the map depends on cell embeddings determined by the
#' reduction technique used
#'
#' @param seu A Seurat object
#' @param embedding The dimensional reduction technique to be used
#' @param group Name of one or more metadata columns to group (color) cells by.
#' @param dims Dimensions to plot, must be a two-length numeric vector
#' @param highlight A list of character or numeric vectors of cells to highlight
#' @param pt.size Adjust point size on the plot
#' @param return_plotly Convert plot to interactive web-based graph
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#'
#'
#'
#' # static mode
#' plot_var(human_gene_transcript_seu, group = "batch", return_plotly = FALSE)
#'
#' # interactive plotly plot
#' plotly_plot <- plot_var(human_gene_transcript_seu, group = "batch")
#' print(plotly_plot)
#'
plot_var <- function(seu, group = "batch", embedding = "umap", dims = c(1, 2), highlight = NULL, pt.size = 1.0, return_plotly = FALSE, ...) {
    Seurat::DefaultAssay(seu) <- "gene"

    metadata <- seu[[]][Seurat::Cells(seu), ]
    key <- rownames(metadata)

    if (embedding == "umap") {
        dims <- c(1, 2)
    } else if (embedding == "tsne") {
        dims <- c(1, 2)
    }

    dims <- as.numeric(dims)

    d <- Seurat::DimPlot(object = seu, dims = dims, reduction = embedding, group.by = group, pt.size = pt.size, ...) +
        aes(key = key, cellid = key) +
        # theme(legend.text=element_text(size=10)) +
        NULL

    if (return_plotly == FALSE) {
        return(d)
    }

    plotly_plot <- plotly::ggplotly(d, tooltip = "cellid", height = 500) %>%
        # htmlwidgets::onRender(javascript) %>%
        # plotly::highlight(on = "plotly_selected", off = "plotly_relayout") %>%
        plotly_settings() %>%
        plotly::toWebGL() %>%
        # plotly::partial_bundle() %>%
        identity()
}

#' Plotly settings
#'
#' Change settings of a plotly plot
#'
#' @param plotly_plot  A plotly plot
#' @param width Default set to '600'
#' @param height Default set to '700'
#'
#' @return
#'
#' @examples
plotly_settings <- function(plotly_plot, width = 600, height = 700) {
    plotly_plot %>%
        plotly::layout(dragmode = "lasso") %>%
        plotly::config(
            toImageButtonOptions = list(
                format = "svg",
                filename = "myplot",
                width = width,
                height = height
            )
        ) %>%
        identity()
}


#' Plot Violin plot
#'
#' Plots a Violin plot of a single data (gene expression, metrics, etc.)
#' grouped by a metadata variable
#'
#' @param seu A Seurat object
#' @param plot_var Variable to group (color) cells by
#' @param plot_vals
#' @param features Features to plot
#' @param assay Name of assay to use, defaults to the active assay
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#'
#' plot_violin(human_gene_transcript_seu, plot_var = "batch", features = c("NRL", "GNAT2"))
#'
plot_violin <- function(seu, plot_var = "batch", plot_vals = NULL, features = "RXRG", assay = "gene", ...) {
    if (is.null(plot_vals)) {
        plot_vals <- unique(seu[[]][[plot_var]])
        plot_vals <- plot_vals[!is.na(plot_vals)]
    }
    seu <- seu[, seu[[]][[plot_var]] %in% plot_vals]
    vln_plot <- Seurat::VlnPlot(seu, features = features, group.by = plot_var, assay = assay, pt.size = 1, ...) +
        geom_boxplot(width = 0.2) +
        # labs(title = "Expression Values for each cell are normalized by that cell's total expression then multiplied by 10,000 and natural-log transformed") +
        # stat_summary(fun.y = mean, geom = "line", size = 4, colour = "black") +
        NULL

    print(vln_plot)
}


#' Plot Feature
#'
#' Plots gene or transcript expression overlaid on a given embedding.
#' If multiple features are supplied the joint density of all features
#' will be plotted using [Nebulosa](https://www.bioconductor.org/packages/devel/bioc/html/Nebulosa.html)
#'
#' @param seu A Seurat object
#' @param embedding Dimensional reduction technique to be used
#' @param features Features to plot
#' @param dims Dimensions to plot, must be a two-length numeric vector
#'
#' @return
#' @export
#' @importFrom ggplot2 aes
#'
#' @examples
#'
#' # static, single feature
#' plot_feature(human_gene_transcript_seu, embedding = "umap", features = c("NRL"), return_plotly = FALSE)
#' # static, multi-feature
#' plot_feature(human_gene_transcript_seu, embedding = "umap", features = c("RXRG", "NRL"), return_plotly = FALSE)
#' # interactive, multi-feature
#' plotly_plot <- plot_feature(human_gene_transcript_seu, embedding = "umap", features = c("RXRG", "NRL"))
#' print(plotly_plot)
#'
plot_feature <- function(seu, embedding = c("umap", "pca", "tsne"), features, dims = c(1, 2), return_plotly = FALSE, pt.size = 1.0) {
    Seurat::DefaultAssay(seu) <- "gene"

    metadata <- seu[[]][Seurat::Cells(seu), ]
    key <- rownames(metadata)

    if (embedding %in% c("tsne", "umap")) {
        dims <- c(1, 2)
    }

    dims <- as.numeric(dims)

    if (length(features) == 1) {
        fp <- Seurat::FeaturePlot(object = seu, features = features, dims = dims, reduction = embedding, pt.size = pt.size, blend = FALSE) +
            ggplot2::aes(key = key, cellid = key, alpha = 0.7)
    } else if (length(features) > 1) {
        nebulosa_plots <- Nebulosa::plot_density(object = seu, features = features, dims = dims, reduction = embedding, size = pt.size, joint = TRUE, combine = FALSE)

        fp <- dplyr::last(nebulosa_plots) +
            ggplot2::aes(key = key, cellid = key, alpha = 0.7)
    }

    if (return_plotly == FALSE) {
        return(fp)
    }

    plotly_plot <- plotly::ggplotly(fp, tooltip = "cellid", height = 500) %>%
        plotly_settings() %>%
        plotly::toWebGL() %>%
        # plotly::partial_bundle() %>%
        identity()
}

#' Plot cell cycle distribution grouped by metadata
#'
#' Plot ridge plots of G1, S, and G2M phases grouped by provided metadata
#'
#' @param seu A seurat object
#' @param features  Features to plot (gene expression, metrics, PC scores, anything that can be retreived by Seurat::FetchData)
#'
#' @return
#' @export
#'
#' @examples
plot_cell_cycle_distribution <- function(seu, features) {
    cc_genes_path <- "~/single_cell_projects/resources/regev_lab_cell_cycle_genes.txt"
    cc.genes <- readLines(con = cc_genes_path)
    s.genes <- cc.genes[1:43]
    g2m.genes <- cc.genes[44:97]

    seu <- CellCycleScoring(
        object = seu, s.genes, g2m.genes,
        set.ident = TRUE
    )

    RidgePlot(object = seu, features = features)
}


#' Plot Cluster Marker Genes
#'
#' Plot a dot plot of n marker features grouped by cell metadata
#' available methods are wilcoxon rank-sum test implemented in
#' [presto](https://github.com/immunogenomics/presto) and specificity scores implemented in [genesorteR](https://github.com/mahmoudibrahim/genesorteR)
#'
#' @param seu a seurat object
#' @param marker_method either "presto" or "genesorteR"
#' @param metavar the metadata variable from which to pick clusters
#' @param num_markers default is 5
#' @param selected_values
#' @param return_plotly whether to return an interactive ploly plot
#' @param featureType
#' @param hide_technical whether to exclude mitochondrial or ribosomal genes
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#'
#' # interactive mode using "presto"
#' plot_markers(human_gene_transcript_seu, metavar = "tech", marker_method = "presto", return_plotly = TRUE)
#'
#' # static mode using "presto"
#' plot_markers(human_gene_transcript_seu, metavar = "tech", marker_method = "genesorteR", return_plotly = FALSE)
#'
plot_markers <- function(seu, metavar = "batch", num_markers = 5, selected_values = NULL, return_plotly = FALSE, marker_method = "presto", seurat_assay = "gene", hide_technical = NULL, unique_markers = FALSE, p_val_cutoff = 1, ...) {
    Idents(seu) <- seu[[]][[metavar]]

    # by default only resolution markers are calculated in pre-processing
    seu <- find_all_markers(seu, metavar, seurat_assay = seurat_assay, p_val_cutoff = p_val_cutoff)

    marker_table <- seu@misc$markers[[metavar]][[marker_method]]

    markers <-
        marker_table %>%
        enframe_markers() %>%
        dplyr::mutate(dplyr::across(.fns = as.character))

    if (!is.null(hide_technical)) {
        markers <- purrr::map(markers, c)

        if (hide_technical == "pseudo") {
            markers <- purrr::map(markers, ~ .x[!.x %in% pseudogenes[[seurat_assay]]])
        } else if (hide_technical == "mito_ribo") {
            markers <- purrr::map(markers, ~ .x[!str_detect(.x, "^MT-")])
            markers <- purrr::map(markers, ~ .x[!str_detect(.x, "^RPS")])
            markers <- purrr::map(markers, ~ .x[!str_detect(.x, "^RPL")])
        } else if (hide_technical == "all") {
            markers <- purrr::map(markers, ~ .x[!.x %in% pseudogenes[[seurat_assay]]])
            markers <- purrr::map(markers, ~ .x[!str_detect(.x, "^MT-")])
            markers <- purrr::map(markers, ~ .x[!str_detect(.x, "^RPS")])
            markers <- purrr::map(markers, ~ .x[!str_detect(.x, "^RPL")])
        }

        min_length <- min(purrr::map_int(markers, length))

        markers <- purrr::map(markers, head, min_length) %>%
            dplyr::bind_cols()
    }

    if (unique_markers) {
        markers <-
            markers %>%
            dplyr::mutate(precedence = row_number()) %>%
            pivot_longer(-precedence, names_to = "group", values_to = "markers") %>%
            dplyr::arrange(markers, precedence) %>%
            dplyr::group_by(markers) %>%
            dplyr::filter(row_number() == 1) %>%
            dplyr::arrange(group, precedence) %>%
            tidyr::drop_na() %>%
            dplyr::group_by(group) %>%
            dplyr::mutate(precedence = row_number()) %>%
            tidyr::pivot_wider(names_from = "group", values_from = "markers") %>%
            dplyr::select(-precedence)
    }

    sliced_markers <-
        markers %>%
        dplyr::slice_head(n = num_markers) %>%
        tidyr::pivot_longer(everything(), names_to = "group", values_to = "feature") %>%
        dplyr::arrange(group) %>%
        dplyr::distinct(feature, .keep_all = TRUE) %>%
        # dplyr::top_n(n = num_markers, wt = logFC) %>%
        identity()

    if (!is.null(selected_values)) {
        seu <- seu[, Idents(seu) %in% selected_values]
        sliced_markers <- sliced_markers %>%
            dplyr::filter(group %in% selected_values) %>%
            dplyr::distinct(feature, .keep_all = TRUE)
    }
    # browser()
    vline_coords <- head(cumsum(table(sliced_markers$group)) + 0.5, -1)

    sliced_markers <- dplyr::pull(sliced_markers, feature)

    seu[[metavar]][is.na(seu[[metavar]])] <- "NA"
    Idents(seu) <- metavar

    markerplot <- DotPlot(seu, assay = "gene", features = sliced_markers, group.by = metavar, dot.scale = 3) +
        ggplot2::theme(
            axis.text.x = ggplot2::element_text(size = 10, angle = 45, vjust = 1, hjust = 1),
            axis.text.y = ggplot2::element_text(size = 10)
        ) +
        ggplot2::scale_y_discrete(position = "left") +
        ggplot2::scale_x_discrete(limits = sliced_markers) +
        ggplot2::geom_vline(xintercept = vline_coords, linetype = 2) +
        ggplot2::coord_flip() +
        NULL


    if (return_plotly == FALSE) {
        return(markerplot)
    }

    plot_height <- (150 * num_markers)
    plot_width <- (100 * length(levels(Idents(seu))))

    markerplot <- plotly::ggplotly(markerplot, height = plot_height, width = plot_width) %>%
        plotly_settings() %>%
        plotly::toWebGL() %>%
        # plotly::partial_bundle() %>%
        identity()

    return(list(plot = markerplot, markers = marker_table))
}

#' Plot Read Count
#'
#' Draw a box plot for read count data of a metadata variable
#'
#' @param seu A seurat object
#' @param metavar Metadata variable to plot. Default set to "nCount_RNA"
#' @param color.by Variable to color bins by. Default set to "batch"
#' @param yscale Scale of y axis. Default set to "linear"
#' @param return_plotly whether to return an interactive ploly plot. Default set to FALSE
#'
#' @return
#' @export
#'
#' @examples
#' # interactive plotly
#' plot_readcount(human_gene_transcript_seu, return_plotly = TRUE)
#' # static plot
#' plot_readcount(human_gene_transcript_seu, return_plotly = FALSE)
#'
#' @importFrom ggplot2 ggplot aes geom_bar theme labs scale_y_log10
plot_readcount <- function(seu, metavar = "nCount_RNA", color.by = "batch", yscale = "linear", return_plotly = FALSE, ...) {
    seu_tbl <- tibble::rownames_to_column(seu[[]], "SID") %>%
        dplyr::select(SID, !!as.symbol(metavar), !!as.symbol(color.by))

    rc_plot <-
        ggplot(seu_tbl, aes(
            x = reorder(SID, -!!as.symbol(metavar)),
            y = !!as.symbol(metavar), fill = !!as.symbol(color.by)
        )) +
        geom_bar(position = "identity", stat = "identity") +
        theme(axis.text.x = element_blank()) +
        labs(
            title = metavar,
            x = "Sample"
        ) +
        NULL

    if (yscale == "log") {
        rc_plot <-
            rc_plot +
            scale_y_log10()
    }

    if (return_plotly == FALSE) {
        return(rc_plot)
    }

    rc_plot <- plotly::ggplotly(rc_plot, tooltip = "cellid", height = 500) %>%
        # htmlwidgets::onRender(javascript) %>%
        # plotly::highlight(on = "plotly_selected", off = "plotly_relayout") %>%
        plotly_settings() %>%
        plotly::toWebGL() %>%
        # plotly::partial_bundle() %>%
        identity()
}


#' Plot Annotated Complexheatmap from Seurat object
#'
#' @param object A Seurat object
#' @param features Vector of features to plot. Features can come
#' @param cells Cells to retain
#' @param group.by  Name of one or more metadata columns to annotate columns by (for example, orig.ident)
#' @param layer
#' @param assay
#' @param group.bar.height
#' @param col_arrangement how to arrange columns whether with a dendrogram (Ward.D2, average, etc.) or exclusively by metadata category
#' @param column_split whether to split columns by metadat value
#' @param mm_col_dend height of column dendrogram
#' @param ... additional arguments passed to ComplexHeatmap::Heatmap
#'
#' @return
#' @export
#'
#' @examples
#'
#' # plot top 50 variable genes
#' top_50_features <- VariableFeatures(human_gene_transcript_seu)[1:50]
#' seu_complex_heatmap(human_gene_transcript_seu, features = top_50_features)
#'
seu_complex_heatmap <- function(seu, features = NULL, group.by = "ident", cells = NULL,
    layer = "scale.data", assay = NULL, group.bar.height = 0.01,
    column_split = NULL, col_arrangement = "ward.D2", mm_col_dend = 30, ...) {
    if (length(GetAssayData(seu, layer = "scale.data")) == 0) {
        message("seurat object has not been scaled. Please run `Seurat::ScaleData` to view a scaled heatmap; showing unscaled expression data")
        layer <- "data"
    }

    cells <- cells %||% colnames(x = seu)
    if (is.numeric(x = cells)) {
        cells <- colnames(x = seu)[cells]
    }
    assay <- assay %||% Seurat::DefaultAssay(object = seu)
    Seurat::DefaultAssay(object = seu) <- assay
    features <- features %||% VariableFeatures(object = seu)
    features <- rev(x = unique(x = features))
    possible.features <- rownames(x = GetAssayData(
        object = seu,
        layer = layer
    ))
    if (any(!features %in% possible.features)) {
        bad.features <- features[!features %in% possible.features]
        features <- features[features %in% possible.features]
        if (length(x = features) == 0) {
            stop(
                "No requested features found in the ", layer,
                " layer for the ", assay, " assay."
            )
        }
        warning(
            "The following features were omitted as they were not found in the ",
            layer, " layer for the ", assay, " assay: ", paste(bad.features,
                collapse = ", "
            )
        )
    }
    data <- as.data.frame(x = t(x = as.matrix(x = GetAssayData(
        object = seu,
        layer = layer
    )[features, cells, drop = FALSE])))
    seu <- suppressMessages(expr = StashIdent(
        object = seu,
        save.name = "ident"
    ))

    if (any(col_arrangement %in% c(
        "ward.D", "single", "complete", "average", "mcquitty",
        "median", "centroid", "ward.D2"
    ))) {
        if ("pca" %in% Seurat::Reductions(seu)) {
            cluster_columns <-
                Seurat::Embeddings(seu, "pca") %>%
                dist() %>%
                hclust(col_arrangement)
        } else {
            message("pca not computed for this dataset; cells will be clustered by displayed features")
            cluster_columns <- function(m) as.dendrogram(cluster::agnes(m), method = col_arrangement)
        }
    } else {
        cells <-
            seu %>%
            Seurat::FetchData(vars = col_arrangement) %>%
            dplyr::arrange(across(all_of(col_arrangement))) %>%
            rownames()

        data <- data[cells, ]

        group.by <- base::union(group.by, col_arrangement)

        cluster_columns <- FALSE
    }

    group.by <- group.by %||% "ident"
    groups.use <- seu[[group.by]][cells, , drop = FALSE]

    groups.use <- groups.use %>%
        tibble::rownames_to_column("sample_id") %>%
        dplyr::mutate(across(where(is.character), ~ str_wrap(str_replace_all(.x, ",", " "), 10))) %>%
        dplyr::mutate(across(where(is.character), as.factor)) %>%
        data.frame(row.names = 1) %>%
        identity()

    # factor colors
    groups.use.factor <- groups.use[sapply(groups.use, is.factor)]
    ha_cols.factor <- NULL
    if (length(groups.use.factor) > 0) {
        ha_col_names.factor <- lapply(groups.use.factor, levels)

        ha_cols.factor <- purrr::map(ha_col_names.factor, ~ scales::hue_pal()(length(.x))) %>%
            purrr::map2(ha_col_names.factor, purrr::set_names)
    }

    # numeric colors
    groups.use.numeric <- groups.use[sapply(groups.use, is.numeric)]
    ha_cols.numeric <- NULL
    if (length(groups.use.numeric) > 0) {
        numeric_col_fun <- function(myvec, color) {
            circlize::colorRamp2(range(myvec), c("white", color))
        }

        ha_col_names.numeric <- names(groups.use.numeric)
        ha_col_hues.numeric <- scales::hue_pal()(length(ha_col_names.numeric))

        ha_cols.numeric <- purrr::map2(groups.use[ha_col_names.numeric], ha_col_hues.numeric, numeric_col_fun)
    }

    ha_cols <- c(ha_cols.factor, ha_cols.numeric)

    column_ha <- ComplexHeatmap::HeatmapAnnotation(df = groups.use, height = grid::unit(group.bar.height, "points"), col = ha_cols)

    hm <- ComplexHeatmap::Heatmap(t(data),
        name = "log expression", top_annotation = column_ha,
        cluster_columns = cluster_columns,
        show_column_names = FALSE,
        column_dend_height = grid::unit(mm_col_dend, "mm"),
        column_split = column_split,
        column_title = NULL,
        ...
    )

    return(hm)
}



#' Plot Transcript Composition
#'
#' plot the proportion of reads of a given gene map to each transcript
#'
#' @param seu A seurat object
#' @param gene_symbol Gene symbol of gene of intrest
#' @param group.by Name of one or more metadata columns to annotate columns by
#' (for example, orig.ident)
#' @param standardize
#' @param drop_zero Drop zero values
#'
#' @return
#' @export
#'
#' @examples
#' plot_transcript_composition(human_gene_transcript_seu, "RXRG", group.by = "gene_snn_res.0.6")
#'
plot_transcript_composition <- function(seu, gene_symbol, group.by = "batch", standardize = FALSE, drop_zero = FALSE) {
    transcripts <- annotables::grch38 %>%
        dplyr::filter(symbol == gene_symbol) %>%
        dplyr::left_join(annotables::grch38_tx2gene, by = "ensgene") %>%
        dplyr::pull(enstxp)

    metadata <- seu@meta.data
    metadata$sample_id <- NULL
    metadata <-
        metadata %>%
        tibble::rownames_to_column("sample_id") %>%
        dplyr::select(sample_id, group.by = {{ group.by }})

    data <- FetchData(seu$transcript, vars = transcripts)

    data <- expm1(as.matrix(data))

    data <-
        data %>%
        as.data.frame() %>%
        tibble::rownames_to_column("sample_id") %>%
        tidyr::pivot_longer(
            cols = starts_with("ENST"),
            names_to = "transcript",
            values_to = "expression"
        ) %>%
        dplyr::left_join(metadata, by = "sample_id") %>%
        dplyr::mutate(
            group.by = as.factor(group.by),
            transcript = as.factor(transcript)
        )

    data <- dplyr::group_by(data, group.by, transcript)

    # drop zero values

    if (drop_zero) {
        data <- dplyr::filter(data, expression != 0)
    }

    data <- dplyr::summarize(data, expression = mean(expression))

    position <- ifelse(standardize, "fill", "stack")

    p <- ggplot(
        data = data,
        aes(x = group.by, y = expression, fill = transcript)
    ) +
        # stat_summary(fun = "mean", geom = "col") +
        geom_col(stat = "identity", position = position) +
        theme_minimal() +
        theme(
            axis.title.x = element_blank(),
            axis.text.x = element_text(
                angle = 45, hjust = 1, vjust = 1, size = 12
            )
        ) +
        labs(title = paste("Mean expression by", group.by, "-", gene_symbol), subtitle = "data scaled by library size then ln transformed") +
        NULL

    return(list(plot = p, data = data))
}

#' Plot All Transcripts
#'
#' plot expression all transcripts for an input gene superimposed on an embedding
#'
#' @param seu A seurat object
#' @param features gene or vector of transcripts
#' @param embedding umap
#' @param from_gene whether to look up transcripts for an input gene
#'
#' @return
#' @export
#'
#' @examples
#'
#' processed_seu <- clustering_workflow(human_gene_transcript_seu)
#' transcripts_to_plot <- genes_to_transcripts("RXRG")
#' plot_all_transcripts(processed_seu, features = transcripts_to_plot)
#'
plot_all_transcripts <- function(seu, features, embedding = "umap", from_gene = TRUE, combine = TRUE) {
    if (from_gene) {
        features <- genes_to_transcripts(features)
    }

    features <- features[features %in% rownames(seu[["transcript"]])]

    # transcript_cols <- as.data.frame(t(as.matrix(seu[["transcript"]][features,])))

    transcript_cols <- FetchData(seu, features)

    seu <- AddMetaData(seu, transcript_cols)

    plot_out <- purrr::map(paste0("transcript_", features), ~ plot_feature(seu,
        embedding = embedding,
        features = .x, return_plotly = FALSE
    )) %>%
        purrr::set_names(features)

    if (combine) {
        plot_out <- wrap_plots(plot_out)
    }

    return(plot_out)
}


#' Plot CC
#'
#' @param object
#' @param group.by
#' @param assay
#' @param label
#' @param phase_levels
#' @param kept_phases
#'
#' @return
#' @export
#'
#' @examples
plot_cc <- function(object = NULL, group.by = "gene_snn_res.0.6", assay = "gene", label = "_filtered_", color.by = "batch", mytitle = "Title", faceted = TRUE, pt.size = 1, phase_levels = c("g1", "g1_s", "s", "s_g2", "g2", "g2_m", "post_mitotic", "hsp", "g1_stress", "other", "s_2"), kept_phases = NULL){

  kept_phases <- kept_phases %||% phase_levels

    heatmap_features  <-
      object@misc$markers[[group.by]][["presto"]]

    cluster_order <- levels(object@meta.data[[group.by]]) %>%
      set_names(.)

    object@meta.data[[group.by]] <-
      factor(object@meta.data[[group.by]], levels = cluster_order)

    group_by_clusters <- object@meta.data[[group.by]]

    object@meta.data$clusters <- names(cluster_order[group_by_clusters])

    object@meta.data$clusters <- factor(object@meta.data$clusters, levels = unique(setNames(names(cluster_order), cluster_order)[levels(object@meta.data[[group.by]])]))

    heatmap_features <-
      heatmap_features %>%
      dplyr::arrange(Cluster) %>%
      group_by(Cluster) %>%
      slice_head(n=6) %>%
      dplyr::filter(Gene.Name %in% rownames(GetAssayData(object, layer = "gene", slot = "scale.data"))) %>%
      dplyr::filter(Gene.Name %in% VariableFeatures(object)) %>%
      identity()

  heatmap_features <-
    heatmap_features %>%
    dplyr::ungroup() %>%
    dplyr::distinct(Gene.Name, .keep_all = TRUE)

  # row_ha = ComplexHeatmap::rowAnnotation(term = rev(heatmap_features$term))

  object_heatmap <- ggplotify::as.ggplot(
    seu_complex_heatmap(object,
                        features = heatmap_features$Gene.Name,
                        group.by = c("G2M.Score", "S.Score", "clusters"),
                        col_arrangement = c("clusters"),
                        cluster_rows = FALSE,
                        column_split =  sort(object@meta.data$clusters),
                        row_split = rev(heatmap_features$Cluster),
                        row_title_rot = 0,
                        # row_split = sort(object@meta.data$clusters)
    )) +
    labs(title = mytitle) +
    theme()


  # browser()
  labels <- data.frame(clusters=unique(object[[]][["clusters"]]), label =unique(object[[]][["clusters"]])) %>%
    # dplyr::rename({{group.by}} := cluster) %>%
    identity()

  cc_data <- FetchData(object, c("clusters", "G2M.Score", "S.Score", "Phase", color.by))

  centroid_data <-
    cc_data %>%
    dplyr::group_by(clusters) %>%
    dplyr::summarise(mean_x = mean(S.Score), mean_y = mean(G2M.Score)) %>%
    dplyr::mutate(clusters = factor(clusters, levels = levels(cc_data$clusters))) %>%
    dplyr::mutate(centroid = "centroids") %>%
    identity()

  centroid_plot <-
    cc_data %>%
    ggplot(aes(x = `S.Score`, y = `G2M.Score`, group = .data[["clusters"]], color = .data[[color.by]])) +
    geom_point(size = pt.size) +
    theme_light() +
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank()
    ) +
    geom_point(data = centroid_data, aes(x = mean_x, y = mean_y, fill = .data[["clusters"]]), size = 6, alpha = 0.7, shape = 23,  colour="black") +
    NULL

  facet_cell_cycle_plot <-
    cc_data %>%
    ggplot(aes(x = `S.Score`, y = `G2M.Score`, group = .data[["clusters"]], color = .data[[color.by]])) +
    geom_point(size = pt.size) +
    geom_point(data = centroid_data, aes(x = mean_x, y = mean_y, fill = .data[["clusters"]]), size = 6, alpha = 0.7, shape = 23,  colour="black") +
    facet_wrap(~.data[["clusters"]], ncol = 2) +
    theme_light() +
    geom_label(data = labels,
               aes(label = label),
               # x = Inf,
               # y = -Inf,
               x = max(cc_data$S.Score)+0.05,
               y = max(cc_data$G2M.Score)-0.1,
               hjust=1,
               vjust=1,
               inherit.aes = FALSE) +
    theme(
      strip.background = element_blank(),
      strip.text.x = element_blank()
    ) +
    # guides(color = "none") +
    NULL

  appender <- function(string) str_wrap(string, width = 40)

  if(faceted){
    # collage_plots <- list(object_heatmap, facet_cell_cycle_plot)
    return(facet_cell_cycle_plot)
  } else {
    # collage_plots <- list(object_heatmap, centroid_plot)
    return(centroid_plot)
  }
#
#     mypatch <- patchwork::wrap_plots(collage_plots) +
#       # plot_layout(widths = c(16, 4)) +
#       # plot_layout(design = layout) +
#       NULL

}
