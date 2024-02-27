#' Split SingleCellExperiment by colData variable
#'
#' @param x single cell object
#' @param f colData variable as a string
#'
#' @return a list of singlecellexperiments name by colData value
#' @export
#'
#' @examples
#' chevreul_sce <- chevreuldata::human_gene_transcript_sce()
#' splitByCol(chevreul_sce, "batch")
splitByCol <- function(x, f = "batch") {

  f <- colData(x)[[f]]

  i <- split(seq_along(f), f)

  v <- vector(mode = "list", length = length(i))

  names(v) <- names(i)

  for (n in names(i)) { v[[n]] <- x[, i[[n]]] }

  return(v)

}


#' Merge Small SingleCellExeriment Objects
#'
#' @param object_list List of two or more singlecell objects
#' @param k.filter minimum cell number for integration
#'
#' @return a single cell object
#' @export
#' @examples
#' chevreul_sce <- chevreuldata::human_gene_transcript_sce()
#' merge_small_objects(
#' "small_batch1" = chevreul_sce[,1:40],
#' "small_batch2" = chevreul_sce[,41:80],
#' "large_batch" = chevreul_sce[,81:300])
#'
merge_small_objects <- function(..., k.filter = 50) {
  object_list <- list(...)

    # check if any singlecell objects are too small and if so merge with the first singlecell objects
    object_dims <- map(object_list, dim) %>%
        map_lgl(~ .x[[2]] < k.filter)

    small_objects <- object_list[object_dims]

    object_list <- object_list[!object_dims]

    object_list[[1]] <- reduce(c(small_objects, object_list[[1]]), correctExperiments, PARAM=NoCorrectParam())

    return(object_list)
}


#' Batch Correct Multiple Single Cell Objects
#'
#' @param object_list List of two or more single cell objects
#' @param method Default "cca"
#' @param organism human or mouse
#' @param ... extra args passed to object_reduce_dimensions
#'
#' @return an integrated single cell object
#' @importFrom batchelor correctExperiments
#' @export
#' @examples
#' chevreul_sce <- chevreuldata::human_gene_transcript_sce()
#' batches <- splitByCol(chevreul_sce, "batch")
#' object_integrate(batches)
object_integrate <- function(object_list, organism = "human", ...) {

    geneCorrected <- correctExperiments(object_list)
    mainExpName(geneCorrected) <- "integrated"

    geneMerged <- correctExperiments(object_list)

    transcriptBatches <- map(object_list, swapAltExp, "transcript")

    transcriptMerged <- correctExperiments(transcriptBatches, PARAM=NoCorrectParam())

    altExp(geneCorrected, "gene") <- geneMerged
    altExp(geneCorrected, "transcript") <- transcriptMerged

    geneCorrected <- record_experiment_data(geneCorrected, experiment_name = "integrated", organism = organism)

    return(geneCorrected)
}

#' Run Louvain Clustering at Multiple Resolutions
#'
#' @param object A single cell objects
#' @param resolution Clustering resolution
#' @param custom_clust custom cluster
#' @param reduction Set dimensional reduction object
#' @param algorithm 1
#' @param ... extra args passed to single cell packages
#'
#' @return a single cell object with louvain clusters
#' @export
#' @importFrom bluster NNGraphParam
#' @importFrom scran clusterCells
#' @importFrom glue glue
#' @examples
#' chevreul_sce <- chevreuldata::human_gene_transcript_sce()
#' object_cluster(chevreul_sce)
object_cluster <- function(object = object, resolution = 0.6, custom_clust = NULL, reduction = "PCA", algorithm = 1, ...) {
        message(paste0("[", format(Sys.time(), "%H:%M:%S"), "] Clustering Cells..."))
        if (length(resolution) > 1) {
            for (i in resolution) {
                message(paste0("clustering at ", i, " resolution"))
                cluster_labels <- clusterCells(object,
                    use.dimred = reduction,
                    BLUSPARAM = NNGraphParam(cluster.fun = "louvain", cluster.args = list(resolution = i))
                )
                colData(object)[[glue("gene_snn_res.{i}")]] <- cluster_labels
            }
        } else if (length(resolution) == 1) {
            message(paste0("clustering at ", resolution, " resolution"))
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
#' @param experiment Assay of interest to be run on the singlecell objects
#' @param reduction Set dimensional reduction object
#' @param ... Extra parameters passed to object_reduce_dimensions
#'
#' @return a single cell object with embeddings
#' @export
#' @examples
#' chevreul_sce <- chevreuldata::human_gene_transcript_sce()
#' object_reduce_dimensions(chevreul_sce)
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
                object <- runTSNE(x = object, dimred = "PCA", n_dimred = 1:30)
            } else {
                object <- runTSNE(x = object, altexp = experiment, dimred = "PCA", n_dimred = 1:30)
            }
            if ("gene" == experiment) {
                object <- runUMAP(x = object, dimred = "PCA", n_dimred = 1:30)
            } else {
                object <- runUMAP(x = object, altexp = experiment, dimred = "PCA", n_dimred = 1:30)
            }
        }
        return(object)
    }

#' Give a new project name to a single cell object
#'
#' @param object A SingleCellExperiment object
#' @param new_name New name to assign
#'
#' @return a renamed single cell object
#' @export
#' @examples
#' chevreul_sce <- chevreuldata::human_gene_transcript_sce()
#' rename_object(chevreul_sce, "new_name")
rename_object <- function (object, new_name)
          {
            metadata(object)["project.name"] <- new_name
            return(object)
          }

#' Reintegrate (filtered) singlecell objectss
#'
#' This function takes a SCE object and perfroms teh below steps
#' 1) split by batch
#' 2) integrate
#' 3) run integration pipeline and save
#'
#' @param object A singlecell objects
#' @param suffix to be appended to file saved in output dir
#' @param reduction to use default is pca
#' @param ... extra args passed to object_integration_pipeline
#'
#' @return a single cell object
#' @export
#' @examples
#' chevreul_sce <- chevreuldata::human_gene_transcript_sce()
#' reintegrate_object(chevreul_sce)
reintegrate_object <- function (object, suffix = "", reduction = "PCA", ...)
          {
            organism <- metadata(object)$experiment$organism
            experiment_name <- metadata(object)$experiment$experiment_name
            objects <- splitByCol(object, "batch")
            object <- object_integration_pipeline(objects, suffix = suffix, ...)
            object <- record_experiment_data(object, experiment_name, organism)
          }
