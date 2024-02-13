#' Preprocess Single Cell Object
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
#' @return a preprocessed single cell object
#' @export
#'
#' @examples
#'
#' panc8[["gene"]] <- object_preprocess(panc8[["gene"]])
#'
setGeneric("object_preprocess", function (assay, scale = TRUE, normalize = TRUE, features = NULL, legacy_settings = FALSE, ...)  standardGeneric("object_preprocess"))

setMethod("object_preprocess", "Seurat",
          function (assay, scale = TRUE, normalize = TRUE, features = NULL, legacy_settings = FALSE, ...)
          {
            if (legacy_settings) {
              message("using legacy settings")
              logtransform_exp <- as.matrix(log1p(Seurat::GetAssayData(assay)))
              assay <- Seurat::SetAssayData(assay, slot = "data", logtransform_exp) %>% Seurat::ScaleData(features = rownames(.))
              return(assay)
            }
            if (normalize) {
              assay <- Seurat::NormalizeData(assay, verbose = FALSE, ...)
            }
            assay <- Seurat::FindVariableFeatures(assay, selection.method = "vst", verbose = FALSE, ...)
            if (scale) {
              assay <- Seurat::ScaleData(assay, features = rownames(assay), ...)
            }
            return(assay)
          }
)

setMethod("object_preprocess", "SingleCellExperiment",
          function (assay, scale = TRUE, normalize = TRUE, features = NULL, legacy_settings = FALSE, ...)
          {
            # assay <- computeLibraryFactors(assay)
            # if (normalize) {
            #   assay <- logNormCounts(assay)
            # }
            #
            # assay <- modelGeneVar(assay)
            #
            # if (scale) {
            # }

            clusters <- quickCluster(assay)
            assay <- computeSumFactors(assay, clusters=clusters)
            # summary(sizeFactors(assay))

            assay <- logNormCounts(assay)

            dec <- modelGeneVar(assay)

            # Get the top 10% of genes.
            top.hvgs <- getTopHVGs(dec, prop=0.1)
            assay <- runPCA(assay, subset_row=top.hvgs)
            # ncol(reducedDim(assayd, "PCA"))
            output <- getClusteredPCs(reducedDim(assay))
            # reducedDim(assay, "PCAsub") <- reducedDim(assay, "PCA")[,1:npcs,drop=FALSE]
            # npcs
            # In this case, using the PCs that we chose from getClusteredPCs().
            g1 <- buildSNNGraph(assay, use.dimred="PCA")


            return(assay)
          }


)

#' Find All Markers
#'
#' Find all markers at a range of resolutions
#'
#' @param object An object.
#' @param group_by A metadata variable to group by.
#' @param object_assay Assay to use, Default "gene".
#' @param ...
#'
#' @return a single cell object containing marker genes
#' @export
#'
#' @examples
#' markers_stashed_object <- find_all_markers(panc8)
#' marker_genes <- Misc(markers_stashed_object, "markers")
#' str(marker_genes)
setGeneric("find_all_markers", function(object, group_by = NULL, object_assay = "gene", ...) standardGeneric("find_all_markers"))

setMethod(
    "find_all_markers", "Seurat",
    function(object, group_by = NULL, object_assay = "gene", ...) {
        if (is.null(group_by)) {
            resolutions <- colnames(pull_metadata(object))[grepl(paste0(object_assay, "_snn_res."), colnames(pull_metadata(object)))]
            cluster_index <- grepl(paste0(object_assay, "_snn_res."), colnames(pull_metadata(object)))
            if (!any(cluster_index)) {
                warning("no clusters found in metadata. runnings object_cluster")
                object <- object_cluster(object, resolution = seq(0.2, 2, by = 0.2))
            }
            clusters <- pull_metadata(object)[, cluster_index]
            cluster_levels <- purrr::map_int(clusters, ~ length(unique(.x)))
            cluster_levels <- cluster_levels[cluster_levels > 1]
            clusters <- select(clusters, dplyr::one_of(names(cluster_levels)))
            group_by <- names(clusters)
        }
        new_markers <- map(group_by, ~ stash_marker_features(object, .x, object_assay = object_assay, ...))
        names(new_markers) <- group_by
        old_markers <- Misc(object)$markers[!names(Misc(object)$markers) %in% names(new_markers)]
        Misc(object)$markers <- c(old_markers, new_markers)
        return(object)
    }
)

setMethod(
    "find_all_markers", "SingleCellExperiment",
    function(object, group_by = NULL, object_assay = "gene", ...) {
        if (is.null(group_by)) {
            resolutions <- colnames(pull_metadata(object))[grepl(paste0(object_assay, "_snn_res."), colnames(pull_metadata(object)))]
            cluster_index <- grepl(paste0(object_assay, "_snn_res."), colnames(pull_metadata(object)))
            if (!any(cluster_index)) {
                warning("no clusters found in metadata. runnings object_cluster")
                object <- object_cluster(object, resolution = seq(0.2, 2, by = 0.2))
            }
            clusters <- pull_metadata(object)[, cluster_index]
            cluster_levels <- purrr::map_int(clusters, ~ length(unique(.x)))
            cluster_levels <- cluster_levels[cluster_levels > 1]
            clusters <- select(clusters, dplyr::one_of(names(cluster_levels)))
            group_by <- names(clusters)
        }
        new_markers <- map(group_by, ~ stash_marker_features(object, .x, object_assay = object_assay, ...))
        names(new_markers) <- group_by
        old_markers <- metadata(object)$markers[!names(metadata(object)$markers) %in% names(new_markers)]
        metadata(object)$markers <- c(old_markers, new_markers)
        return(object)
    }
)


#' Title
#'
#' @param marker_table a table of marker genes
#'
#' @return a table of marker genes
#'
enframe_markers <- function(marker_table) {
    marker_table %>%
        select(Gene.Name, Cluster) %>%
        mutate(rn = row_number()) %>%
        tidyr::pivot_wider(names_from = Cluster, values_from = Gene.Name) %>%
        select(-rn)
}

#' Stash Marker Genes in a Seurat Object
#'
#' Marker Genes will be stored in slot `@misc$markers`
#'
#' @param group_by A metadata variable to group by
#' @param object A object
#' @param object_assay An assay to use
#' @param top_n Use top n genes, Default "200"
#' @param p_val_cutoff p value cut-off, Default value is "0.5"
#'
#' @return a single cell object with marker genes
#'
#' @examples
#'
#' object <- stash_marker_features(group_by = "batch", object, object_assay = "gene")
#'
setGeneric("stash_marker_features", function(object, group_by, object_assay = "gene", top_n = 200, p_val_cutoff = 0.5) standardGeneric("stash_marker_features"))

setMethod(
    "stash_marker_features", "Seurat",
    function(object, group_by, object_assay = "gene", top_n = 200, p_val_cutoff = 0.5) {
        message(paste0("stashing presto markers for ", group_by))
        markers <- list()
        markers$presto <- presto::wilcoxauc(object, group_by, seurat_assay = object_assay) %>%
            group_by(group) %>%
            filter(padj < p_val_cutoff) %>%
            dplyr::top_n(n = top_n, wt = logFC) %>%
            arrange(group, desc(logFC)) %>%
            select(Gene.Name = feature, Average.Log.Fold.Change = logFC, Adjusted.pvalue = padj, avgExpr, Cluster = group)
        return(markers)
    }
)

setMethod(
    "stash_marker_features", "SingleCellExperiment",
    function(object, group_by, object_assay = "gene", top_n = 200, p_val_cutoff = 0.5) {
        message(paste0("stashing presto markers for ", group_by))
        markers <- list()
        markers$presto <- presto::wilcoxauc(object, group_by, seurat_assay = object_assay) %>%
            group_by(group) %>%
            filter(padj < p_val_cutoff) %>%
            dplyr::top_n(n = top_n, wt = logFC) %>%
            arrange(group, desc(logFC)) %>%
            select(Gene.Name = feature, Average.Log.Fold.Change = logFC, Adjusted.pvalue = padj, avgExpr, Cluster = group)
        return(markers)
    }
)
