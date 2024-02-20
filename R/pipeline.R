# seurat integration ------------------------------

#' Run Seurat Integration
#'
#' Run batch correction, followed by:
#' 1) stashing of batches in metadata 'batch'
#' 2) clustering with resolution 0.2 to 2.0 in increments of 0.2
#' 3) saving to <proj_dir>/output/sce/<feature>_object_<suffix>.rds
#'
#' @param suffix a suffix to be appended to a file save in output dir
#' @param object_list List of objects to be integrated
#' @param resolution Range of resolution
#' @param algorithm Algorithm for modularity optimization. Default 1:original Louvain algorithm
#' @param organism Default "human"
#' @param annotate_cell_cycle whether to score cell cycle phases
#' @param reduction pca, umap, or tsne
#' @param ... extra args passed to object_integrate
#'
#' @return an integrated single cell object
#' @export
#'
#' @examples
#'
#' batches <- panc8 %>%
#'     Seurat::SplitObject(split.by = "tech")
#'
#' integrated_object <- object_integration_pipeline(batches)
object_integration_pipeline <- function(object_list, resolution = seq(0.2, 2.0, by = 0.2), suffix = "", algorithm = 1, organism = "human", annotate_cell_cycle = FALSE, annotate_percent_mito = FALSE, reduction = "pca", ...) {
        experiment_names <- names(object_list)

        organisms <- case_when(
            grepl("Hs", experiment_names) ~ "human",
            grepl("Mm", experiment_names) ~ "mouse"
        )

        names(organisms) <- experiment_names

        organisms[is.na(organisms)] <- organism

        integrated_object <- object_integrate(object_list, organism = organism, ...)

        # cluster merged objects
        if ("harmony" %in% names(integrated_object@reductions)) reduction <- "harmony"
        integrated_object <- object_cluster(integrated_object, resolution = resolution, algorithm = algorithm, reduction = reduction, ...)


        object_assay <- "gene"
        if ("harmony" %in% names(integrated_object@reductions)) object_assay <- "integrated"
        integrated_object <- find_all_markers(integrated_object, object_assay = object_assay)

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

#' Run Seurat Pipeline
#'
#' This functions allows you to Preprocess, Cluster and Reduce Dimensions for a single object.
#'
#' @param object A Seurat object
#' @param assay Assay of interest in Seurat object
#' @param resolution Resolution for clustering cells. Default set to 0.6.
#' @param reduction Dimensional reduction object
#' @param organism Organism
#' @param ... extra parameters passed to internal functions
#'
#' @return a processed single cell object
#' @export
#'
#' @examples
#'
#' processed_object <- object_pipeline(panc8)
#'
object_pipeline <- function(object, assay = "gene", resolution = 0.6, reduction = "PCA", organism = "human", ...) {
        assays <- names(object@assays)

        assays <- assays[assays %in% c("gene", "transcript")]

        for (assay in assays) {
            object[[assay]] <- object_preprocess(object[[assay]], scale = TRUE, ...)
        }

        # PCA
        object <- object_reduce_dimensions(object, reduction = reduction, ...)

        object <- object_cluster(object = object, resolution = resolution, reduction = reduction, ...)

        object <- find_all_markers(object, object_assay = "gene")

        # if (feature == "gene"){
        #   enriched_object <- tryCatch(getEnrichedPathways(object), error = function(e) e)
        #   enrichr_available <- !any(class(enriched_object) == "error")
        #   if(enrichr_available){
        #     object <- enriched_object
        #   }
        # }



        # annotate cell cycle scoring to objects
        object <- annotate_cell_cycle(object, organism = organism, ...)

        # annotate mitochondrial percentage in object metadata
        object <- add_percent_mito(object, organism = organism)

        return(object)
    }
