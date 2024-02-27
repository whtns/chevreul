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
#'
#' @return an integrated single cell object
#' @export
#'
#' @examples
#' chevreul_sce <- chevreuldata::human_gene_transcript_sce()
#' batches <- splitByCol(chevreul_sce, "batch")
#' integrated_object <- object_integration_pipeline(batches)
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
#' @return a processed single cell object
#' @export
#'
#' @examples
#' chevreul_sce <- chevreuldata::human_gene_transcript_sce()
#'
#' processed_object <- object_pipeline(chevreul_sce)
#'
object_pipeline <- function(object, experiment = "gene", resolution = 0.6, reduction = "PCA", organism = "human", ...) {

  object <- object_preprocess(object, scale = TRUE, ...)
        for (experiment in altExpNames(object)[!altExpNames(object) == "velocity"]) {
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
