#' Preprocess Seurat Object
#'
#' Performs standard pre-processing workflow for scRNA-seq data
#'
#' @param assay Assay to use
#' @param scale Perform linear transformation 'Scaling'
#' @param normalize Perform normalization
#' @param features Identify highly variable features
#' @param legacy_settings Use legacy settings
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#'
#' panc8[["gene"]] <- seurat_preprocess(panc8[["gene"]])
#'
seurat_preprocess <- function(assay, scale = TRUE, normalize = TRUE, features = NULL, legacy_settings = FALSE, ...) {
    # Normalize data

    if (legacy_settings) {
        message("using legacy settings")

        logtransform_exp <- as.matrix(log1p(Seurat::GetAssayData(assay)))

        assay <- Seurat::SetAssayData(assay, slot = "data", logtransform_exp) %>%
            Seurat::ScaleData(features = rownames(.))

        return(assay)
    }

    if (normalize) {
        assay <- Seurat::NormalizeData(assay, verbose = FALSE, ...)
    }

    # Filter out only variable genes
    assay <- Seurat::FindVariableFeatures(assay, selection.method = "vst", verbose = FALSE, ...)

    # Regress out unwanted sources of variation
    if (scale) {
        assay <- Seurat::ScaleData(assay, features = rownames(assay), ...)
    }

    return(assay)
}

#' Find All Markers
#'
#' Find all markers at a range of resolutions
#'
#' @param seu A seurat object.
#' @param metavar A metadata variable to group by.
#' @param seurat_assay Assay to use, Default "gene".
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' markers_stashed_seu <- find_all_markers(panc8)
#' marker_genes <- Misc(markers_stashed_seu, "markers")
#' str(marker_genes)
find_all_markers <- function(seu, metavar = NULL, seurat_assay = "gene", ...) {
    if (is.null(metavar)) {
        resolutions <- colnames(seu[[]])[grepl(paste0(seurat_assay, "_snn_res."), colnames(seu[[]]))]

        cluster_index <- grepl(paste0(seurat_assay, "_snn_res."), colnames(seu[[]]))

        if (!any(cluster_index)) {
            warning("no clusters found in metadata. runnings seurat_cluster")
            seu <- seurat_cluster(seu, resolution = seq(0.2, 2.0, by = 0.2))
        }

        clusters <- seu[[]][, cluster_index]

        cluster_levels <- purrr::map_int(clusters, ~ length(unique(.x)))
        cluster_levels <- cluster_levels[cluster_levels > 1]

        clusters <- dplyr::select(clusters, dplyr::one_of(names(cluster_levels)))
        metavar <- names(clusters)
    }

    new_markers <- purrr::map(metavar, stash_marker_features, seu, seurat_assay = seurat_assay, ...)
    names(new_markers) <- metavar

    old_markers <- seu@misc$markers[!names(seu@misc$markers) %in% names(new_markers)]

    seu@misc$markers <- c(old_markers, new_markers)

    return(seu)
}

enframe_markers <- function(marker_table) {
    marker_table %>%
        dplyr::select(Gene.Name, Cluster) %>%
        dplyr::mutate(rn = row_number()) %>%
        tidyr::pivot_wider(names_from = Cluster, values_from = Gene.Name) %>%
        dplyr::select(-rn)
}

#' Stash Marker Genes in a Seurat Object
#'
#' Marker Genes will be stored in slot `@misc$markers`
#'
#' @param metavar A metadata variable to group by
#' @param seu A seurat object
#' @param seurat_assay An assay to use
#' @param top_n Use top n genes, Default "200"
#' @param p_val_cutoff p value cut-off, Default value is "0.5"
#'
#' @return
#'
#' @examples
#'
#' seu <- stash_marker_features(metavar = "batch", seu, seurat_assay = "gene")
#'
stash_marker_features <- function(metavar, seu, seurat_assay, top_n = 200, p_val_cutoff = 0.5) {
    message(paste0("stashing presto markers for ", metavar))

    markers <- list()
    markers$presto <-
        presto::wilcoxauc(seu, metavar, seurat_assay = seurat_assay) %>%
        dplyr::group_by(group) %>%
        dplyr::filter(padj < p_val_cutoff) %>%
        dplyr::top_n(n = top_n, wt = logFC) %>%
        dplyr::arrange(group, desc(logFC)) %>%
        dplyr::select(Gene.Name = feature, Average.Log.Fold.Change = logFC, Adjusted.pvalue = padj, avgExpr, Cluster = group)

    # message(paste0("stashing genesorteR markers for ", metavar))
    #
    # markers$genesorteR <- tryCatch(
    #   {
    #     gs <- genesorteR::sortGenes(
    #       Seurat::GetAssayData(seu, assay = seurat_assay, slot = "data"),
    #       tidyr::replace_na(seu[[]][[metavar]], "NA")
    #     )
    #
    #     pp <- genesorteR::getPValues(gs)
    #
    #     genesorter_table <- genesorteR::getTable(gs, pp, adjpval_cutoff = p_val_cutoff)
    #
    #   },
    #   error = function(e) {
    #     message(sprintf("Error in %s: %s", deparse(e[["call"]]), e[["message"]]))
    #     NULL
    #   },
    #   finally = {
    #   }
    # )

    return(markers)
}
