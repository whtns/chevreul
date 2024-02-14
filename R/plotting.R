#' Unite metadata
#'
#' @param object A object
#' @param group_bys A feature or variable to combine
#'
#' @return an object with Idents formed from concatenation of group_bys
#' @export
setGeneric("unite_metadata", function(object, group_bys) standardGeneric("unite_metadata"))

setMethod(
    "unite_metadata", "Seurat",
    function(object, group_bys) {
        newcolname <- paste(group_bys, collapse = "_by_")
        newdata <- object[[group_bys]] %>%
            tidyr::unite(!!newcolname, group_bys) %>%
            tibble::deframe()
        Idents(object) <- newdata
        return(object)
    }
)

setMethod(
    "unite_metadata", "SingleCellExperiment",
    function(object, group_bys) {
        newcolname <- paste(group_bys, collapse = "_by_")
        newdata <- colData(object)[group_bys] %>%
            as.data.frame() %>%
            tidyr::unite(!!newcolname, group_bys) %>%
            tibble::deframe()
        # Idents(object) <- newdata
        return(object)
    }
)

#' Plot Metadata Variables
#'
#' Plots static or interactive plot where each point represents a cell metadata
#' variable whose position on the map depends on cell embeddings determined by the
#' reduction technique used
#'
#' @param object A Seurat object
#' @param embedding The dimensional reduction technique to be used
#' @param group Name of one or more metadata columns to group (color) cells by.
#' @param dims Dimensions to plot, must be a two-length numeric vector
#' @param highlight A list of character or numeric vectors of cells to highlight
#' @param pt.size Adjust point size on the plot
#' @param return_plotly Convert plot to interactive web-based graph
#' @param ... extra parameters passed to ggplot
#'
#' @return a ggplot
#' @export
#'
#' @examples
#' # static mode
#' \dontrun{
#' plot_var(human_gene_transcript_object, group = "batch", return_plotly = FALSE)
#' }
#'
#' # interactive plotly plot
#' \dontrun{
#' plotly_plot <- plot_var(human_gene_transcript_object, group = "batch")
#' }
#' \dontrun{
#' print(plotly_plot)
#' }
setGeneric("plot_var", function(object, group = "batch", embedding = "umap", dims = c(1, 2), highlight = NULL, pt.size = 1.0, return_plotly = FALSE, ...) {
    standardGeneric("plot_var")
})

setMethod(
    "plot_var", "Seurat",
    function(object, group = "batch", embedding = "umap", dims = c(1, 2), highlight = NULL, pt.size = 1.0, return_plotly = FALSE, ...) {
        Seurat::DefaultAssay(object) <- "gene"

        # metadata <- tibble::as_tibble(pull_metadata(object)[Seurat::Cells(object),], rownames = "sID")
        # cellid <- metadata[["sID"]]
        # key <- rownames(metadata)

        metadata <- pull_metadata(object)[Seurat::Cells(object), ]
        key <- rownames(metadata)

        if (embedding == "umap") {
            dims <- c(1, 2)
        } else if (embedding == "tsne") {
            dims <- c(1, 2)
        }

        dims <- as.numeric(dims)

        d <- Seurat::DimPlot(object = object, dims = dims, reduction = embedding, group.by = group, pt.size = pt.size, ...) +
            aes(key = key, cellid = key) +
            # theme(legend.text=element_text(size=10)) +
            NULL

        if (return_plotly == FALSE) {
            return(d)
        }

        plotly_plot <- ggplotly(d, tooltip = "cellid", height = 500) %>%
            # htmlwidgets::onRender(javascript) %>%
            # highlight(on = "plotly_selected", off = "plotly_relayout") %>%
            plotly_settings() %>%
            toWebGL() %>%
            # partial_bundle() %>%
            identity()
    }
)

setMethod(
    "plot_var", "SingleCellExperiment",
    function(object, group = "batch", embedding = "UMAP", dims = c(1, 2), highlight = NULL, pt.size = 1.0, return_plotly = FALSE, ...) {
        metadata <- pull_metadata(object)
        key <- rownames(metadata)

        if (embedding == "UMAP") {
            dims <- c(1, 2)
        } else if (embedding == "TSNE") {
            dims <- c(1, 2)
        }

        dims <- as.numeric(dims)

        d <- scater::plotReducedDim(object = object, dimred = embedding, ncomponents = 2, color_by = group, ...) +
            aes(key = key, cellid = key) +
            # theme(legend.text=element_text(size=10)) +
            NULL

        if (return_plotly == FALSE) {
            return(d)
        }

        plotly_plot <- ggplotly(d, tooltip = "cellid", height = 500) %>%
            # htmlwidgets::onRender(javascript) %>%
            # highlight(on = "plotly_selected", off = "plotly_relayout") %>%
            plotly_settings() %>%
            toWebGL() %>%
            # partial_bundle() %>%
            identity()
    }
)

#' Plotly settings
#'
#' Change settings of a plotly plot
#'
#' @param plotly_plot  A plotly plot
#' @param width Default set to '600'
#' @param height Default set to '700'
#'
#'@noRd
plotly_settings <- function(plotly_plot, width = 600, height = 700) {
    plotly_plot %>%
        layout(dragmode = "lasso") %>%
        config(toImageButtonOptions = list(format = "svg", filename = "myplot", width = width, height = height)) %>%
        identity()
}

#' Plot Violin plot
#'
#' Plots a Violin plot of a single data (gene expression, metrics, etc.)
#' grouped by a metadata variable
#'
#' @param object A Seurat object
#' @param plot_var Variable to group (color) cells by
#' @param plot_vals plot values
#' @param features Features to plot
#' @param assay Name of assay to use, defaults to the active assay
#' @param ... extra parameters passed to ggplot2
#'
#' @return a violin plot
#' @export
#'
#' @examples
#'
#' plot_violin(human_gene_transcript_object, plot_var = "batch", features = c("NRL", "GNAT2"))
#'
setGeneric("plot_violin", function(object, plot_var = "batch", plot_vals = NULL, features = "RXRG", assay = "gene", ...) {
    standardGeneric("plot_violin")
})

setMethod(
    "plot_violin", "Seurat",
    function(object, plot_var = "batch", plot_vals = NULL, features = "RXRG", assay = "gene", ...) {
        if (is.null(plot_vals)) {
            plot_vals <- unique(object[[]][[plot_var]])
            plot_vals <- plot_vals[!is.na(plot_vals)]
        }
        object <- object[, object[[]][[plot_var]] %in% plot_vals]
        vln_plot <- Seurat::VlnPlot(object, features = features, group.by = plot_var, assay = assay, pt.size = 1, ...) + geom_boxplot(width = 0.2) + NULL
        print(vln_plot)
    }
)

setMethod(
    "plot_violin", "SingleCellExperiment",
    function(object, plot_var = "batch", plot_vals = NULL, features = "RXRG", assay = "gene", ...) {
        if (is.null(plot_vals)) {
            plot_vals <- unique(pull_metadata(object)[[plot_var]])
            plot_vals <- plot_vals[!is.na(plot_vals)]
        }
        object <- object[, pull_metadata(object)[[plot_var]] %in% plot_vals]
        vln_plot <- scater::plotExpression(object, features = features, x = plot_var, color_by = plot_var) + geom_boxplot(width = 0.2) + NULL
        print(vln_plot)
    }
)



#' Plot Feature
#'
#' Plots gene or transcript expression overlaid on a given embedding.
#' If multiple features are supplied the joint density of all features
#' will be plotted using [Nebulosa](https://www.bioconductor.org/packages/devel/bioc/html/Nebulosa.html)
#'
#' @param object A Seurat object
#' @param embedding Dimensional reduction technique to be used
#' @param features Features to plot
#' @param dims Dimensions to plot, must be a two-length numeric vector
#'
#' @return an embedding colored by a feature of interest
#' @export
#' @importFrom ggplot2 aes
setGeneric("plot_feature", function(object, embedding = c("umap", "pca", "tsne"), features, dims = c(1, 2), return_plotly = FALSE, pt.size = 1.0) {
    standardGeneric("plot_feature")
})

setMethod(
    "plot_feature", "Seurat",
    function(object, embedding = c("umap", "pca", "tsne"), features, dims = c(1, 2), return_plotly = FALSE, pt.size = 1.0) {
        embedding <- tolower(embedding)

        Seurat::DefaultAssay(object) <- "gene"

        metadata <- pull_metadata(object)[Seurat::Cells(object), ]
        key <- rownames(metadata)

        if (embedding %in% c("tsne", "umap")) {
            dims <- c(1, 2)
        }

        dims <- as.numeric(dims)

        if (length(features) == 1) {
            fp <- Seurat::FeaturePlot(object = object, features = features, dims = dims, reduction = embedding, pt.size = pt.size, blend = FALSE) +
                ggplot2::aes(key = key, cellid = key, alpha = 0.7)
        } else if (length(features) > 1) {
            nebulosa_plots <- Nebulosa::plot_density(object = object, features = features, dims = dims, reduction = embedding, size = pt.size, joint = TRUE, combine = FALSE)

            fp <- dplyr::last(nebulosa_plots) +
                ggplot2::aes(key = key, cellid = key, alpha = 0.7)
        }

        if (return_plotly == FALSE) {
            return(fp)
        }

        plotly_plot <- ggplotly(fp, tooltip = "cellid", height = 500) %>%
            plotly_settings() %>%
            toWebGL() %>%
            # partial_bundle() %>%
            identity()
    }
)

setMethod(
    "plot_feature", "SingleCellExperiment",
    function(object, embedding = c("umap", "pca", "tsne"), features, dims = c(1, 2), return_plotly = FALSE, pt.size = 1.0) {
        embedding <- toupper(embedding)

        metadata <- pull_metadata(object)
        key <- rownames(metadata)

        if (embedding %in% c("TSNE", "UMAP")) {
            dims <- c(1, 2)
        }

        dims <- as.numeric(dims)

        if (length(features) == 1) {
            fp <- scater::plotReducedDim(object = object, color_by = features, dimred = embedding) +
                ggplot2::aes(key = key, cellid = key, alpha = 0.7)
        } else if (length(features) > 1) {
            nebulosa_plots <- Nebulosa::plot_density(object = object, features = features, dims = dims, reduction = embedding, size = pt.size, joint = TRUE, combine = FALSE)

            fp <- dplyr::last(nebulosa_plots) +
                ggplot2::aes(key = key, cellid = key, alpha = 0.7)
        }

        if (return_plotly == FALSE) {
            return(fp)
        }

        plotly_plot <- ggplotly(fp, tooltip = "cellid", height = 500) %>%
            plotly_settings() %>%
            toWebGL() %>%
            # partial_bundle() %>%
            identity()
    }
)

#' Plot cell cycle distribution grouped by metadata
#'
#' Plot ridge plots of G1, S, and G2M phases grouped by provided metadata
#'
#' @param object A single cell object
#' @return a ggplot of cell cycle scores
#' @export
setGeneric("plot_cell_cycle_distribution", function(object) standardGeneric("plot_cell_cycle_distribution"))

setMethod(
    "plot_cell_cycle_distribution", "Seurat",
    function(object) {
        s.genes <- cc.genes[["s.genes"]]
        g2m.genes <- cc.genes[["g2m.genes"]]
        object <- CellCycleScoring(object = object, s.genes, g2m.genes, set.ident = TRUE)
        return(object)
    }
)

setMethod(
    "plot_cell_cycle_distribution", "SingleCellExperiment",
    function(object) {
        # s.genes <- cc.genes[["s.genes"]]
        # g2m.genes <- cc.genes[["g2m.genes"]]
        hs_pairs0 <- cc.genes.cyclone
        assignments <- cyclone(object, hs_pairs0, gene.names = rownames(object))
        colData(object)[colnames(assignments$scores)] <- assignments$scores
        colData(object)["Phase"] <- assignments$phases
        return(object)
    }
)


#' Plot Cluster Marker Genes
#'
#' Plot a dot plot of n marker features grouped by cell metadata
#' available methods are wilcoxon rank-sum test
#'
#' @param object a object
#' @param marker_method "wilcox"
#' @param group_by the metadata variable from which to pick clusters
#' @param num_markers default is 5
#' @param selected_values selected values to display
#' @param return_plotly whether to return an interactive ploly plot
#' @param featureType gene or transcript
#' @param hide_technical whether to exclude mitochondrial or ribosomal genes
#' @param ... extra parameters passed to ggplot2
#'
#' @return a ggplot with marker genes from group_by
#' @export
#'
#' @examples
#'
#' interactive mode using "wilcox"
#' \dontrun{plot_markers(human_gene_transcript_object, group_by = "tech", marker_method = "wilcox", return_plotly = TRUE)}
#'
setGeneric("plot_markers", function(object, group_by = "batch", num_markers = 5, selected_values = NULL, return_plotly = FALSE, marker_method = "wilcox", object_assay = "gene", hide_technical = NULL, unique_markers = FALSE, p_val_cutoff = 1, ...) standardGeneric("plot_markers"))

setMethod(
    "plot_markers", "Seurat",
    function(object, group_by = "batch", num_markers = 5, selected_values = NULL, return_plotly = FALSE, marker_method = "wilcox", object_assay = "gene", hide_technical = NULL, unique_markers = FALSE, p_val_cutoff = 1, ...) {
        Idents(object) <- pull_metadata(object)[[group_by]]
        object <- find_all_markers(object, group_by, object_assay = object_assay, p_val_cutoff = p_val_cutoff)
        marker_table <- Misc(object)$markers[[group_by]][[marker_method]]
        markers <- marker_table %>%
            enframe_markers() %>%
            mutate(dplyr::across(.fns = as.character))
        if (!is.null(hide_technical)) {
            markers <- map(markers, c)
            if (hide_technical == "pobjectdo") {
                markers <- map(markers, ~ .x[!.x %in% pobjectdogenes[[object_assay]]])
            } else if (hide_technical == "mito_ribo") {
                markers <- map(markers, ~ .x[!str_detect(.x, "^MT-")])
                markers <- map(markers, ~ .x[!str_detect(.x, "^RPS")])
                markers <- map(markers, ~ .x[!str_detect(.x, "^RPL")])
            } else if (hide_technical == "all") {
                markers <- map(markers, ~ .x[!.x %in% pobjectdogenes[[object_assay]]])
                markers <- map(markers, ~ .x[!str_detect(.x, "^MT-")])
                markers <- map(markers, ~ .x[!str_detect(.x, "^RPS")])
                markers <- map(markers, ~ .x[!str_detect(.x, "^RPL")])
            }
            min_length <- min(purrr::map_int(markers, length))
            markers <- map(markers, head, min_length) %>% dplyr::bind_cols()
        }
        if (unique_markers) {
            markers <- markers %>%
                mutate(precedence = row_number()) %>%
                pivot_longer(-precedence, names_to = "group", values_to = "markers") %>%
                arrange(markers, precedence) %>%
                group_by(markers) %>%
                filter(row_number() == 1) %>%
                arrange(group, precedence) %>%
                drop_na() %>%
                group_by(group) %>%
                mutate(precedence = row_number()) %>%
                tidyr::pivot_wider(names_from = "group", values_from = "markers") %>%
                select(-precedence)
        }
        sliced_markers <- markers %>%
            dplyr::slice_head(n = num_markers) %>%
            tidyr::pivot_longer(everything(), names_to = "group", values_to = "feature") %>%
            arrange(group) %>%
            distinct(feature, .keep_all = TRUE) %>%
            identity()
        if (!is.null(selected_values)) {
            object <- object[, Idents(object) %in% selected_values]
            sliced_markers <- sliced_markers %>%
                filter(group %in% selected_values) %>%
                distinct(feature, .keep_all = TRUE)
        }
        vline_coords <- head(cumsum(table(sliced_markers$group)) + 0.5, -1)
        sliced_markers <- pull(sliced_markers, feature)
        object[[group_by]][is.na(object[[group_by]])] <- "NA"
        Idents(object) <- group_by

        markerplot <- DotPlot(object, assay = "gene", features = sliced_markers, group.by = group_by, dot.scale = 3) + ggplot2::theme(axis.text.x = ggplot2::element_text(size = 10, angle = 45, vjust = 1, hjust = 1), axis.text.y = ggplot2::element_text(size = 10)) + ggplot2::scale_y_discrete(position = "left") + ggplot2::scale_x_discrete(limits = sliced_markers) + ggplot2::geom_vline(xintercept = vline_coords, linetype = 2) + ggplot2::coord_flip() + NULL

        if (return_plotly == FALSE) {
            return(markerplot)
        }
        plot_height <- (150 * num_markers)
        plot_width <- (100 * length(levels(Idents(object))))
        markerplot <- ggplotly(markerplot, height = plot_height, width = plot_width) %>%
            plotly_settings() %>%
            toWebGL() %>%
            identity()
        return(list(plot = markerplot, markers = marker_table))
    }
)

setMethod(
    "plot_markers", "SingleCellExperiment",
    function(object, group_by = "batch", num_markers = 5, selected_values = NULL, return_plotly = FALSE, marker_method = "wilcox", object_assay = "gene", hide_technical = NULL, unique_markers = FALSE, p_val_cutoff = 1, ...) {
        # Idents(object) <- pull_metadata(object)[[group_by]]
        object <- find_all_markers(object, group_by, object_assay = object_assay, p_val_cutoff = p_val_cutoff)
        marker_table <- metadata(object)$markers[[group_by]][[marker_method]]
        markers <- marker_table %>%
            enframe_markers() %>%
            mutate(dplyr::across(.fns = as.character))
        if (!is.null(hide_technical)) {
            markers <- map(markers, c)
            if (hide_technical == "pobjectdo") {
                markers <- map(markers, ~ .x[!.x %in% pobjectdogenes[[object_assay]]])
            } else if (hide_technical == "mito_ribo") {
                markers <- map(markers, ~ .x[!str_detect(.x, "^MT-")])
                markers <- map(markers, ~ .x[!str_detect(.x, "^RPS")])
                markers <- map(markers, ~ .x[!str_detect(.x, "^RPL")])
            } else if (hide_technical == "all") {
                markers <- map(markers, ~ .x[!.x %in% pobjectdogenes[[object_assay]]])
                markers <- map(markers, ~ .x[!str_detect(.x, "^MT-")])
                markers <- map(markers, ~ .x[!str_detect(.x, "^RPS")])
                markers <- map(markers, ~ .x[!str_detect(.x, "^RPL")])
            }
            min_length <- min(purrr::map_int(markers, length))
            markers <- map(markers, head, min_length) %>% dplyr::bind_cols()
        }
        if (unique_markers) {
            markers <- markers %>%
                mutate(precedence = row_number()) %>%
                pivot_longer(-precedence, names_to = "group", values_to = "markers") %>%
                arrange(markers, precedence) %>%
                group_by(markers) %>%
                filter(row_number() == 1) %>%
                arrange(group, precedence) %>%
                drop_na() %>%
                group_by(group) %>%
                mutate(precedence = row_number()) %>%
                tidyr::pivot_wider(names_from = "group", values_from = "markers") %>%
                select(-precedence)
        }
        sliced_markers <- markers %>%
            dplyr::slice_head(n = num_markers) %>%
            tidyr::pivot_longer(everything(), names_to = "group", values_to = "feature") %>%
            arrange(group) %>%
            distinct(feature, .keep_all = TRUE) %>%
            identity()
        if (!is.null(selected_values)) {
            object <- object[, pull_metadata(object)[[group_by]] %in% selected_values]
            sliced_markers <- sliced_markers %>%
                filter(group %in% selected_values) %>%
                distinct(feature, .keep_all = TRUE)
        }
        vline_coords <- head(cumsum(table(sliced_markers$group)) + 0.5, -1)
        sliced_markers <- pull(sliced_markers, feature)
        object[[group_by]][is.na(object[[group_by]])] <- "NA"
        markerplot <- scater::plotDots(object, features = sliced_markers, group = group_by) +
            ggplot2::theme(axis.text.x = ggplot2::element_text(size = 10, angle = 45, vjust = 1, hjust = 1), axis.text.y = ggplot2::element_text(size = 10)) +
            # ggplot2::scale_y_discrete(position = "left") +
            # ggplot2::scale_x_discrete(limits = sliced_markers) +
            ggplot2::geom_hline(yintercept = vline_coords, linetype = 2) +
            NULL
        if (return_plotly == FALSE) {
            return(markerplot)
        }
        plot_height <- (150 * num_markers)
        plot_width <- (100 * length(levels(as.factor(pull_metadata(object)[[group_by]]))))
        markerplot <- ggplotly(markerplot, height = plot_height, width = plot_width) %>%
            plotly_settings() %>%
            toWebGL() %>%
            identity()
        return(list(plot = markerplot, markers = marker_table))
    }
)


#' Plot Read Count
#'
#' Draw a box plot for read count data of a metadata variable
#'
#' @param object A object
#' @param group_by Metadata variable to plot. Default set to "nCount_RNA"
#' @param color.by Variable to color bins by. Default set to "batch"
#' @param yscale Scale of y axis. Default set to "linear"
#' @param return_plotly whether to return an interactive ploly plot. Default set to FALSE
#' @param ... extra args passed to ggplot2
#'
#' @return a histogram of read counts
#' @export
#'
#' @examples
#' interactive plotly
#' \dontrun{plot_readcount(human_gene_transcript_object, return_plotly = TRUE)}
#'
#' static plot
#' \dontrun{plot_readcount(human_gene_transcript_object, return_plotly = FALSE)}
#'
#' @importFrom ggplot2 ggplot aes geom_bar theme labs scale_y_log10
setGeneric("plot_readcount", function(object, group_by = "nCount_RNA", color.by = "batch", yscale = "linear", return_plotly = FALSE, ...) standardGeneric("plot_readcount"))

setMethod(
    "plot_readcount", "Seurat",
    function(object, group_by = "nCount_RNA", color.by = "batch", yscale = "linear", return_plotly = FALSE, ...) {
        object_tbl <- rownames_to_column(pull_metadata(object), "SID") %>% select(SID, !!as.symbol(group_by), !!as.symbol(color.by))
        rc_plot <- ggplot(object_tbl, aes(x = reorder(SID, -!!as.symbol(group_by)), y = !!as.symbol(group_by), fill = !!as.symbol(color.by))) +
            geom_bar(position = "identity", stat = "identity") +
            theme(axis.text.x = element_blank()) +
            labs(title = group_by, x = "Sample") +
            NULL
        if (yscale == "log") {
            rc_plot <- rc_plot + scale_y_log10()
        }
        if (return_plotly == FALSE) {
            return(rc_plot)
        }
        rc_plot <- ggplotly(rc_plot, tooltip = "cellid", height = 500) %>%
            plotly_settings() %>%
            toWebGL() %>%
            identity()
    }
)

setMethod(
    "plot_readcount", "SingleCellExperiment",
    function(object, group_by = "nCount_RNA", color.by = "batch", yscale = "linear", return_plotly = FALSE, ...) {
        object_tbl <- rownames_to_column(pull_metadata(object), "SID") %>% select(SID, !!as.symbol(group_by), !!as.symbol(color.by))
        rc_plot <- ggplot(object_tbl, aes(x = reorder(SID, -!!as.symbol(group_by)), y = !!as.symbol(group_by), fill = !!as.symbol(color.by))) +
            geom_bar(position = "identity", stat = "identity") +
            theme(axis.text.x = element_blank()) +
            labs(title = group_by, x = "Sample") +
            NULL
        if (yscale == "log") {
            rc_plot <- rc_plot + scale_y_log10()
        }
        if (return_plotly == FALSE) {
            return(rc_plot)
        }
        rc_plot <- ggplotly(rc_plot, tooltip = "cellid", height = 500) %>%
            plotly_settings() %>%
            toWebGL() %>%
            identity()
    }
)

#' Plot Annotated Complexheatmap from Seurat object
#'
#' @param object A Seurat object
#' @param features Vector of features to plot. Features can come
#' @param cells Cells to retain
#' @param group.by  Name of one or more metadata columns to annotate columns by (for example, orig.ident)
#' @param layer "counts" for raw data "scale.data" for log-normalized data
#' @param assay assay to display
#' @param group.bar.height height for group bars
#' @param col_arrangement how to arrange columns whether with a dendrogram (Ward.D2, average, etc.) or exclusively by metadata category
#' @param column_split whether to split columns by metadat value
#' @param mm_col_dend height of column dendrogram
#' @param ... additional arguments passed to ComplexHeatmap::Heatmap
#'
#' @return a complexheatmap
#' @export
#'
#' @examples
#' \dontrun{top_50_features <- get_variable_features(human_gene_transcript_object)[1:50]}
#' \dontrun{make_complex_heatmap(human_gene_transcript_object, features = top_50_features)}
#'
setGeneric("make_complex_heatmap", function(object, features = NULL, group.by = "ident", cells = NULL, layer = "scale.data", assay = NULL, group.bar.height = 0.01, column_split = NULL, col_arrangement = "ward.D2", mm_col_dend = 30, ...) standardGeneric("make_complex_heatmap"))

setMethod(
    "make_complex_heatmap", "Seurat",
    function(object, features = NULL, group.by = "ident", cells = NULL, layer = "scale.data", assay = NULL, group.bar.height = 0.01, column_split = NULL, col_arrangement = "ward.D2", mm_col_dend = 30, ...) {
        if (length(GetAssayData(object, layer = "scale.data")) == 0) {
            message("object has not been scaled. Please run `Seurat::ScaleData` to view a scaled heatmap; showing unscaled expression data")
            layer <- "data"
        }
        cells <- cells %||% colnames(x = object)
        if (is.numeric(x = cells)) {
            cells <- colnames(x = object)[cells]
        }
        assay <- assay %||% Seurat::DefaultAssay(object = object)
        Seurat::DefaultAssay(object = object) <- assay
        features <- features %||% get_variable_features(object = object)
        features <- rev(x = unique(x = features))
        possible.features <- rownames(x = GetAssayData(object = object, layer = layer))
        if (any(!features %in% possible.features)) {
            bad.features <- features[!features %in% possible.features]
            features <- features[features %in% possible.features]
            if (length(x = features) == 0) {
                stop("No requested features found in the ", layer, " layer for the ", assay, " assay.")
            }
            warning("The following features were omitted as they were not found in the ", layer, " layer for the ", assay, " assay: ", paste(bad.features, collapse = ", "))
        }
        data <- as.data.frame(x = t(x = as.matrix(x = GetAssayData(object = object, layer = layer)[features, cells, drop = FALSE])))
        object <- suppressMessages(expr = StashIdent(object = object, save.name = "ident"))
        if (any(col_arrangement %in% c("ward.D", "single", "complete", "average", "mcquitty", "median", "centroid", "ward.D2"))) {
            if ("pca" %in% Seurat::Reductions(object)) {
                cluster_columns <- Seurat::Embeddings(object, "pca") %>%
                    dist() %>%
                    hclust(col_arrangement)
            } else {
                message("pca not computed for this dataset; cells will be clustered by displayed features")
                cluster_columns <- function(m) as.dendrogram(cluster::agnes(m), method = col_arrangement)
            }
        } else {
            cells <- object %>%
                Seurat::FetchData(vars = col_arrangement) %>%
                arrange(across(all_of(col_arrangement))) %>%
                rownames()
            data <- data[cells, ]
            group.by <- base::union(group.by, col_arrangement)
            cluster_columns <- FALSE
        }
        group.by <- group.by %||% "ident"
        groups.use <- object[[group.by]][cells, , drop = FALSE]
        groups.use <- groups.use %>%
            rownames_to_column("sample_id") %>%
            mutate(across(where(is.character), ~ stringr::str_wrap(stringr::str_replace_all(.x, ",", " "), 10))) %>%
            mutate(across(where(is.character), as.factor)) %>%
            data.frame(row.names = 1) %>%
            identity()
        groups.use.factor <- groups.use[sapply(groups.use, is.factor)]
        ha_cols.factor <- NULL
        if (length(groups.use.factor) > 0) {
            ha_col_names.factor <- lapply(groups.use.factor, levels)
            ha_cols.factor <- map(ha_col_names.factor, ~ (scales::hue_pal())(length(.x))) %>% purrr::map2(ha_col_names.factor, purrr::set_names)
        }
        groups.use.numeric <- groups.use[sapply(groups.use, is.numeric)]
        ha_cols.numeric <- NULL
        if (length(groups.use.numeric) > 0) {
            numeric_col_fun <- function(myvec, color) {
                circlize::colorRamp2(range(myvec), c("white", color))
            }
            ha_col_names.numeric <- names(groups.use.numeric)
            ha_col_hues.numeric <- (scales::hue_pal())(length(ha_col_names.numeric))
            ha_cols.numeric <- purrr::map2(groups.use[ha_col_names.numeric], ha_col_hues.numeric, numeric_col_fun)
        }
        ha_cols <- c(ha_cols.factor, ha_cols.numeric)
        column_ha <- ComplexHeatmap::HeatmapAnnotation(df = groups.use, height = grid::unit(group.bar.height, "points"), col = ha_cols)
        hm <- ComplexHeatmap::Heatmap(t(data), name = "log expression", top_annotation = column_ha, cluster_columns = cluster_columns, show_column_names = FALSE, column_dend_height = grid::unit(mm_col_dend, "mm"), column_split = column_split, column_title = NULL, ...)
        return(hm)
    }
)

setMethod(
    "make_complex_heatmap", "SingleCellExperiment",
    function(object, features = NULL, group.by = "ident", cells = NULL, layer = "scale.data", assay = NULL, group.bar.height = 0.01, column_split = NULL, col_arrangement = "ward.D2", mm_col_dend = 30, ...) {
        assay_method <- switch(layer,
            counts = "counts",
            scale.data = "logcounts"
        )


        cells <- cells %||% colnames(x = object)
        if (is.numeric(x = cells)) {
            cells <- colnames(x = object)[cells]
        }
        assay <- assay %||% mainExpName(object)
        if (!assay == mainExpName(object)) {
            object <- swapAltExp(object, name = assay)
        }
        features <- features %||% scran::getTopHVGs(object)
        features <- rev(unique(features))
        possible.features <- rownames(x = SummarizedExperiment::assay(object, assay_method))
        if (any(!features %in% possible.features)) {
            bad.features <- features[!features %in% possible.features]
            features <- features[features %in% possible.features]
            if (length(x = features) == 0) {
                stop("No requested features found in the ", layer, " layer for the ", assay, " assay.")
            }
            warning("The following features were omitted as they were not found in the ", assay_method, " layer for the ", assay, " assay: ", paste(bad.features, collapse = ", "))
        }
        data <- as.data.frame(x = t(x = as.matrix(x = assay(object, assay_method)[features, cells, drop = FALSE])))

        # object <- suppressMessages(expr = StashIdent(object = object, save.name = "ident"))

        if (any(col_arrangement %in% c("ward.D", "single", "complete", "average", "mcquitty", "median", "centroid", "ward.D2"))) {
            if ("PCA" %in% reducedDimNames(object)) {
                cluster_columns <- reducedDim(human_gene_transcript_sce, "PCA") %>%
                    dist() %>%
                    hclust(col_arrangement)
            } else {
                message("pca not computed for this dataset; cells will be clustered by displayed features")
                cluster_columns <- function(m) as.dendrogram(cluster::agnes(m), method = col_arrangement)
            }
        } else {
            cells <- colData(object)[col_arrangement] %>%
              as.data.frame() %>%
                arrange(across(all_of(col_arrangement))) %>%
                rownames()
            data <- data[cells, ]
            group.by <- base::union(group.by, col_arrangement)
            cluster_columns <- FALSE
        }
        group.by <- group.by %||% "ident"
        groups.use <- colData(object)[group.by] %>% as.data.frame()
        groups.use <- groups.use %>%
            rownames_to_column("sample_id") %>%
            mutate(across(where(is.character), ~ stringr::str_wrap(stringr::str_replace_all(.x, ",", " "), 10))) %>%
            mutate(across(where(is.character), as.factor)) %>%
            data.frame(row.names = 1) %>%
            identity()
        groups.use.factor <- groups.use[sapply(groups.use, is.factor)]
        ha_cols.factor <- NULL
        if (length(groups.use.factor) > 0) {
            ha_col_names.factor <- lapply(groups.use.factor, levels)
            ha_cols.factor <- map(ha_col_names.factor, ~ (scales::hue_pal())(length(.x))) %>% purrr::map2(ha_col_names.factor, purrr::set_names)
        }
        groups.use.numeric <- groups.use[sapply(groups.use, is.numeric)]
        ha_cols.numeric <- NULL
        if (length(groups.use.numeric) > 0) {
            numeric_col_fun <- function(myvec, color) {
                circlize::colorRamp2(range(myvec), c("white", color))
            }
            ha_col_names.numeric <- names(groups.use.numeric)
            ha_col_hues.numeric <- (scales::hue_pal())(length(ha_col_names.numeric))
            ha_cols.numeric <- purrr::map2(groups.use[ha_col_names.numeric], ha_col_hues.numeric, numeric_col_fun)
        }
        ha_cols <- c(ha_cols.factor, ha_cols.numeric)
        column_ha <- ComplexHeatmap::HeatmapAnnotation(df = groups.use, height = grid::unit(group.bar.height, "points"), col = ha_cols)
        hm <- ComplexHeatmap::Heatmap(t(data), name = "log expression", top_annotation = column_ha, cluster_columns = cluster_columns, show_column_names = FALSE, column_dend_height = grid::unit(mm_col_dend, "mm"), column_split = column_split, column_title = NULL, ...)
        return(hm)
    }
)

#' Plot Transcript Composition
#'
#' plot the proportion of reads of a given gene map to each transcript
#'
#' @param object A object
#' @param gene_symbol Gene symbol of gene of intrest
#' @param group.by Name of one or more metadata columns to annotate columns by
#' (for example, orig.ident)
#' @param standardize whether to standardize values
#' @param drop_zero Drop zero values
#'
#' @return a stacked barplot of transcript counts
#' @export
#'
#' @examples
#' \dontrun{plot_transcript_composition(human_gene_transcript_object, "RXRG",
#' group.by = "gene_snn_res.0.6")}
#'
setGeneric("plot_transcript_composition", function(object, gene_symbol, group.by = "batch", standardize = FALSE, drop_zero = FALSE) standardGeneric("plot_transcript_composition"))

setMethod(
    "plot_transcript_composition", "Seurat",
    function(object, gene_symbol, group.by = "batch", standardize = FALSE, drop_zero = FALSE) {
        transcripts <- annotables::grch38 %>%
            filter(symbol == gene_symbol) %>%
            left_join(annotables::grch38_tx2gene, by = "ensgene") %>%
            pull(enstxp)
        metadata <- object@meta.data
        metadata$sample_id <- NULL
        metadata <- metadata %>%
            rownames_to_column("sample_id") %>%
            select(sample_id, group.by = {{ group.by }})
        data <- FetchData(object$transcript, vars = transcripts)
        data <- expm1(as.matrix(data))
        data <- data %>%
            as.data.frame() %>%
            rownames_to_column("sample_id") %>%
            tidyr::pivot_longer(cols = starts_with("ENST"), names_to = "transcript", values_to = "expression") %>%
            left_join(metadata, by = "sample_id") %>%
            mutate(group.by = as.factor(group.by), transcript = as.factor(transcript))
        data <- group_by(data, group.by, transcript)
        if (drop_zero) {
            data <- filter(data, expression != 0)
        }
        data <- summarize(data, expression = mean(expression))
        position <- ifelse(standardize, "fill", "stack")
        p <- ggplot(data = data, aes(x = group.by, y = expression, fill = transcript)) +
            geom_col(stat = "identity", position = position) +
            theme_minimal() +
            theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12)) +
            labs(title = paste("Mean expression by", group.by, "-", gene_symbol), subtitle = "data scaled by library size then ln transformed") +
            NULL
        return(list(plot = p, data = data))
    }
)

setMethod(
    "plot_transcript_composition", "SingleCellExperiment",
    function(object, gene_symbol, group.by = "batch", standardize = FALSE, drop_zero = FALSE) {
        transcripts <- annotables::grch38 %>%
            filter(symbol == gene_symbol) %>%
            left_join(annotables::grch38_tx2gene, by = "ensgene") %>%
            pull(enstxp)
        metadata <- pull_metadata(object)
        metadata$sample_id <- NULL
        metadata <- metadata %>%
            rownames_to_column("sample_id") %>%
            select(sample_id, group.by = {{ group.by }})

        transcripts <- transcripts[transcripts %in% rownames(altExp(object, "transcript"))]

        data <- counts(altExp(object, "transcript"))[transcripts, ] %>%
            as.matrix() %>%
            t()

        data <- data %>%
            as.data.frame() %>%
            rownames_to_column("sample_id") %>%
            tidyr::pivot_longer(cols = starts_with("ENST"), names_to = "transcript", values_to = "expression") %>%
            left_join(metadata, by = "sample_id") %>%
            mutate(group.by = as.factor(group.by), transcript = as.factor(transcript))
        data <- group_by(data, group.by, transcript)
        if (drop_zero) {
            data <- filter(data, expression != 0)
        }
        data <- summarize(data, expression = mean(expression))
        position <- ifelse(standardize, "fill", "stack")
        p <- ggplot(data = data, aes(x = group.by, y = expression, fill = transcript)) +
            geom_col(stat = "identity", position = position) +
            theme_minimal() +
            theme(axis.title.x = element_blank(), axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1, size = 12)) +
            labs(title = paste("Mean expression by", group.by, "-", gene_symbol)) +
            NULL
        return(list(plot = p, data = data))
    }
)

#' Plot All Transcripts
#'
#' plot expression all transcripts for an input gene superimposed on an embedding
#'
#' @param object A object
#' @param features gene or vector of transcripts
#' @param embedding umap
#' @param from_gene whether to look up transcripts for an input gene
#' @param combine TRUE
#'
#' @return a list of embedding plots colored by a feature of interest
#' @export
#'
#' @examples
#'
#' \dontrun{processed_object <- clustering_workflow(human_gene_transcript_object)
#' transcripts_to_plot <- genes_to_transcripts("RXRG")
#' plot_all_transcripts(processed_object, features = transcripts_to_plot)}
#'
setGeneric("plot_all_transcripts", function(object, features, embedding = "umap", from_gene = TRUE, combine = TRUE) standardGeneric("plot_all_transcripts"))

setMethod(
    "plot_all_transcripts", "Seurat",
    function(object, features, embedding = "umap", from_gene = TRUE, combine = TRUE) {
        if (from_gene) {
            features <- genes_to_transcripts(features)
        }
        features <- features[features %in% rownames(object[["transcript"]])]
        transcript_cols <- FetchData(object, features)
        object <- AddMetaData(object, transcript_cols)
        plot_out <- map(paste0("transcript_", features), ~ plot_feature(object, embedding = embedding, features = .x, return_plotly = FALSE)) %>% purrr::set_names(features)
        if (combine) {
            plot_out <- wrap_plots(plot_out)
        }
        return(plot_out)
    }
)

setMethod(
    "plot_all_transcripts", "SingleCellExperiment",
    function(object, features, embedding = "UMAP", from_gene = TRUE, combine = TRUE) {
        if (from_gene) {
            features <- genes_to_transcripts(features)
        }
        features <- features[features %in% rownames(altExp(object, "transcript"))]
        transcript_cols <- assay(altExp(object, "transcript"))[features, ]
        colData(object)[features] <- t(as.matrix(transcript_cols))
        # plot_out <- scater::plotReducedDim(altExp(object), features = features, dimred = embedding)
        plot_out <- map(paste0(features), ~ plot_feature(object, embedding = embedding, features = .x, return_plotly = FALSE)) %>% purrr::set_names(features)
        if (combine) {
            plot_out <- wrap_plots(plot_out)
        }
        return(plot_out)
    }
)
