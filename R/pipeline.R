#' Run SingleCellExperiment Integration
#'
#' Run batch correction, followed by:
#' 1) stashing of batches in metadata 'batch'
#' 2) clustering with resolution 0.2 to 2.0 in increments of 0.2
#' 3) saving to <proj_dir>/output/sce/<feature>_object_<suffix>.rds
#'
#' @param suffix a suffix to be appended to a file save in output dir
#' @param object_list List of objects to be integrated
#' @param resolution Range of resolution
#' @param organism Default "human"
#' @param annotate_cell_cycle whether to score cell cycle phases
#' @param reduction pca, umap, or tsne
#' @param ... extra args passed to object_integrate
#' @param annotate_percent_mito logical scalar whether to annotate mitochondrial percentage
#'
#' @return an integrated SingleCellExperiment object
object_integration_pipeline <- function(object_list, resolution = seq(0.2, 2.0, by = 0.2), suffix = "", organism = "human", annotate_cell_cycle = FALSE, annotate_percent_mito = FALSE, reduction = "PCA", ...) {
    experiment_names <- names(object_list)

    organisms <- case_when(
        grepl("Hs", experiment_names) ~ "human",
        grepl("Mm", experiment_names) ~ "mouse"
    )

    names(organisms) <- experiment_names

    organisms[is.na(organisms)] <- organism

    integrated_object <- object_integrate(object_list, organism = organism, ...)

    integrated_object <- object_reduce_dimensions(integrated_object, ...)

    # cluster merged objects
    integrated_object <- object_cluster(integrated_object, resolution = resolution, algorithm = algorithm, reduction = reduction, ...)

    experiment <- "gene"
    integrated_object <- find_all_markers(integrated_object, experiment = experiment)

    #   enriched_object <- tryCatch(getEnrichedPathways(integrated_object), error = function(e) e)
    #   enrichr_available <- !any(class(enriched_object) == "error")
    #   if(enrichr_available){
    #     integrated_object <- enriched_object
    #   }

    # annotate cell cycle scoring to objects
    if (annotate_cell_cycle) {
        integrated_object <- annotate_cell_cycle(integrated_object, ...)
    }

    # annotate mitochondrial percentage in object metadata
    if (annotate_percent_mito) {
        integrated_object <- add_percent_mito(integrated_object, ...)
    }

    # annotate excluded cells
    # integrated_object <- annotate_excluded(integrated_object, excluded_cells)

    return(integrated_object)
}

# object_pipeline ------------------------------

#' Run SingleCellExperiment Pipeline
#'
#' This functions allows you to Preprocess, Cluster and Reduce Dimensions for a single object.
#'
#' @param object A SingleCellExperiment object
#' @param experiment Assay of interest in SingleCellExperiment object
#' @param resolution Resolution for clustering cells. Default set to 0.6.
#' @param reduction Dimensional reduction object
#' @param organism Organism
#' @param ... extra parameters passed to internal functions
#'
#' @return a processed SingleCellExperiment object
object_pipeline <- function(object, experiment = "gene", resolution = 0.6, reduction = "PCA", organism = "human", ...) {
    object <- object_preprocess(object, scale = TRUE, ...)
    for (experiment in altExpNames(object)) {
        altExp(object, experiment) <- object_preprocess(altExp(object, experiment), scale = TRUE, ...)
    }

    # PCA
    object <- object_reduce_dimensions(object, ...)

    object <- object_cluster(object = object, resolution = resolution, reduction = reduction, ...)

    object <- find_all_markers(object, experiment = "gene")

    # if (feature == "gene"){
    #   enriched_object <- tryCatch(getEnrichedPathways(object), error = function(e) e)
    #   enrichr_available <- !any(class(enriched_object) == "error")
    #   if(enrichr_available){
    #     object <- enriched_object
    #   }
    # }


    # annotate cell cycle scoring to objects
    object <- annotate_cell_cycle(object, ...)

    # annotate mitochondrial percentage in object metadata
    object <- add_percent_mito(object, ...)

    return(object)
}

#' Run Louvain Clustering at Multiple Resolutions
#'
#' @param object A SingleCellExperiment objects
#' @param resolution Clustering resolution
#' @param custom_clust custom cluster
#' @param reduction Set dimensional reduction object
#' @param algorithm 1
#' @param ... extra args passed to single cell packages
#'
#' @return a SingleCellExperiment object with louvain clusters
object_cluster <- function(object = object, resolution = 0.6, custom_clust = NULL, reduction = "PCA", algorithm = 1, ...) {
    message(glue("[{format(Sys.time(), '%H:%M:%S')}] Clustering Cells..."))
    if (length(resolution) > 1) {
        for (i in resolution) {
            message(glue("clustering at {i} resolution"))
            cluster_labels <- clusterCells(object,
                                                                         use.dimred = reduction,
                                                                         BLUSPARAM = NNGraphParam(cluster.fun = "louvain", cluster.args = list(resolution = i))
            )
            colData(object)[[glue("gene_snn_res.{i}")]] <- cluster_labels
        }
    } else if (length(resolution) == 1) {
        message(glue("clustering at {resolution} resolution"))
        cluster_labels <- clusterCells(object,
                                                                     use.dimred = reduction,
                                                                     BLUSPARAM = NNGraphParam(cluster.fun = "louvain", cluster.args = list(resolution = resolution))
        )
        
        
        colData(object)[[glue("gene_snn_res.{resolution}")]] <- cluster_labels
    }
    
    return(object)
}

#' Dimensional Reduction
#'
#' Run PCA, TSNE and UMAP on a singlecell objects
#' perplexity should not be bigger than 3 * perplexity < nrow(X) - 1, see details for interpretation
#'
#' @param object A SingleCellExperiment object
#' @param experiment Experiment of interest to be processed
#' @param ... Extra parameters passed to object_reduce_dimensions
#'
#' @return a SingleCellExperiment object with embeddings
object_reduce_dimensions <- function(object, experiment = "gene", ...) {
    num_samples <- dim(object)[[2]]
    if (num_samples < 50) {
        npcs <- num_samples - 1
    } else {
        npcs <- 50
    }
    if ("gene" == experiment) {
        object <- runPCA(x = object, subset_row = getTopHVGs(stats = object), ncomponents = npcs, ...)
    } else {
        object <- runPCA(x = object, altexp = experiment, subset_row = getTopHVGs(stats = object), ncomponents = npcs, ...)
    }
    
    if ((ncol(object) - 1) > 3 * 30) {
        if ("gene" == experiment) {
            object <- runTSNE(x = object, dimred = "PCA", n_dimred = seq(30))
        } else {
            object <- runTSNE(x = object, altexp = experiment, dimred = "PCA", n_dimred = seq(30))
        }
        if ("gene" == experiment) {
            object <- runUMAP(x = object, dimred = "PCA", n_dimred = seq(30))
        } else {
            object <- runUMAP(x = object, altexp = experiment, dimred = "PCA", n_dimred = seq(30))
        }
    }
    return(object)
}

#' Give a new project name to a SingleCellExperiment object
#'
#' @param object A SingleCellExperiment object
#' @param new_name New name to assign
#'
#' @return a renamed SingleCellExperiment object
#' @export
#' @examples
#' 
#' 
#' data(small_example_dataset)
#' rename_object(small_example_dataset, "new_name")
rename_object <- function(object, new_name) {
    metadata(object)["project.name"] <- new_name
    return(object)
}
