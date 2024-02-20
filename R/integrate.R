#' Integrate small datasets with harmony
#'
#' @param object_list List of two or more singlecell objectss to integrate
#'
#' @return an integrated single cell object
#' @export
#'
harmony_integrate <- function(object_list) {
    object_list.integrated <- purrr::reduce(object_list, merge)
    object_list.integrated@assays[["integrated"]] <- object_list.integrated@assays[["gene"]]
    DefaultAssay(object_list.integrated) <- "integrated"
    object_list.integrated <- object_preprocess(object_list.integrated)
    object_list.integrated <- object_reduce_dimensions(object_list.integrated)
    object_list.integrated <- RunHarmony(object_list.integrated, group.by.vars = "batch", assay.use = "integrated")
    object_list.integrated <- RunUMAP(object_list.integrated, assay = "integrated", reduction = "harmony", dims = 1:30)
    object_list.integrated <- FindNeighbors(object_list.integrated, assay = "integrated", reduction = "harmony", dims = 1:30)
    object_list.integrated
}

#' Merge Small SingleCellExeriment Objects
#'
#' @param object_list List of two or more singlecell objects
#' @param k.filter minimum cell number for integration
#'
#' @return a single cell object
#' @export
#'
merge_small_objects <- function(object_list, k.filter = 50) {
    # check if any singlecell objects are too small and if so merge with the first singlecell objects
    object_dims <- map(object_list, dim) %>%
        purrr::map_lgl(~ .x[[2]] < k.filter)

    small_objects <- object_list[object_dims]

    object_list <- object_list[!object_dims]

    object_list[[1]] <- purrr::reduce(c(small_objects, object_list[[1]]), merge)

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
#'
object_integrate <- function(object_list, method = "cca", organism = "human", ...) {

  universe <-
    map(object_list, rownames) %>%
    purrr::reduce(intersect) %>%
    identity()

  # object_list <- merge_small_objects(object_list)
  object_list.integrated <- batchelor::correctExperiments(object_list)

      # see https://github.com/satijalab/seurat/issues/6341------------------------------
    cells_per_batch <- sapply(object_list, ncol)
    min_k_weight <- min(cells_per_batch) - 1
    min_k_weight <- ifelse(min_k_weight < 100, min_k_weight, 100)

    object_list.integrated <- record_experiment_data(object_list.integrated, experiment_name = "integrated", organism = organism)

    return(object_list.integrated)
}

#' Batch Correct Multiple SingleCellExperiments
#'
#' @param object_list List of two or more singlecell objects
#' @param method Default "cca"
#' @param organism human or mouse
#' @param ... extra args passed to object_reduce_dimensions
#'
#' @return an integrated single cell object
#' @export
#'
sce_integrate <- function(object_list, method = "cca", organism = "human", ...) {
  # To construct a reference we will identify ‘anchors’ between the individual datasets. First, we split the combined object into a list, with each dataset as an element.

  # Prior to finding anchors, we perform standard preprocessing (log-normalization), and identify variable features individually for each. Note that Seurat v3 implements an improved method for variable feature selection based on a variance stabilizing transformation ("vst")

  for (i in 1:length(x = object_list)) {
    object_list[[i]][["gene"]] <- object_preprocess(object_list[[i]][["gene"]], scale = TRUE)

    object_list[[i]]$batch <- names(object_list)[[i]]
  }

  object_list <- merge_small_objects(object_list)

  if (method == "rpca") {
    # scale and run pca for each separate batch in order to use reciprocal pca instead of cca
    features <- SelectIntegrationFeatures(object.list = object_list)
    object_list <- map(object_list, Seurat::ScaleData, features = features)
    object_list <- map(object_list, Seurat::RunPCA, features = features)
    object_list.anchors <- FindIntegrationAnchors(object.list = object_list, reduction = "rpca", dims = 1:30)
  } else if (method == "cca") {
    # Next, we identify anchors using the FindIntegrationAnchors function, which takes a list of Seurat objects as input.
    object_list.anchors <- Seurat::FindIntegrationAnchors(object.list = object_list, dims = 1:30, k.filter = 50)
  }

  # proceed with integration
  # object_list.integrated <- IntegrateData(anchorset = object_list.anchors, dims = 1:30)

  # see https://github.com/satijalab/seurat/issues/6341------------------------------
  cells_per_batch <- sapply(object_list, ncol)
  min_k_weight <- min(cells_per_batch) - 1
  min_k_weight <- ifelse(min_k_weight < 100, min_k_weight, 100)

  object_list.integrated <- tryCatch(IntegrateData(anchorset = object_list.anchors, dims = 1:30, k.weight = min_k_weight), error = function(e) e)
  run_harmony <- any(class(object_list.integrated) == "error")
  if (run_harmony) {
    object_list.integrated <- harmony_integrate(object_list)
  }

  #   enriched_object <- tryCatch(getEnrichedPathways(object), error = function(e) e)
  #   enrichr_available <- !any(class(enriched_object) == "error")
  #   if(enrichr_available){
  #     object <- enriched_object
  #   }

  # Next, we identify anchors using the FindIntegrationAnchors function, which takes a list of Seurat objects as input.

  # #stash batches
  Idents(object_list.integrated) <- "batch"
  object_list.integrated[["batch"]] <- Idents(object_list.integrated)

  # switch to integrated assay. The variable features of this assay are
  # automatically set during IntegrateData
  Seurat::DefaultAssay(object = object_list.integrated) <- "integrated"

  # if not integrated with harmony run the standard workflow for visualization and clustering
  if (!"harmony" %in% names(object_list.integrated@reductions)) {
    object_list.integrated <- Seurat::ScaleData(object = object_list.integrated, verbose = FALSE)
    object_list.integrated <- object_reduce_dimensions(object_list.integrated, ...)
  }

  object_list.integrated <- record_experiment_data(object_list.integrated, experiment_name = "integrated", organism = organism)

  return(object_list.integrated)
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
object_cluster <- function(object = object, resolution = 0.6, custom_clust = NULL, reduction = "PCA", algorithm = 1, ...) {
        message(paste0("[", format(Sys.time(), "%H:%M:%S"), "] Clustering Cells..."))
        # return list of graph object with KNN (SNN?)
        # object <- scran::buildKNNGraph(x = object, use.dimred = reduction)
        if (length(resolution) > 1) {
            for (i in resolution) {
                message(paste0("clustering at ", i, " resolution"))
                # object <- Seurat::FindClusters(object = object, resolution = i, algorithm = algorithm, ...)
                cluster_labels <- scran::clusterCells(object,
                    use.dimred = reduction,
                    BLUSPARAM = NNGraphParam(cluster.fun = "louvain", cluster.args = list(resolution = i))
                )
                # colLabels(object) <- cluster_labels
                colData(object)[[glue::glue("gene_snn_res.{i}")]] <- cluster_labels
            }
        } else if (length(resolution) == 1) {
            message(paste0("clustering at ", resolution, " resolution"))
            # object <- Seurat::FindClusters(object = object, resolution = resolution, algorithm = algorithm, ...)
            cluster_labels <- scran::clusterCells(object,
                use.dimred = reduction,
                BLUSPARAM = NNGraphParam(cluster.fun = "louvain", cluster.args = list(resolution = resolution))
            )


            colData(object)[[glue("gene_snn_res.{resolution}")]] <- cluster_labels
        }
        # if (!is.null(custom_clust)) {
        #
        #   clusters <- tibble::tibble(sample_id = rownames(get_cell_metadata(object))) %>% rownames_to_column("order") %>% dplyr::inner_join(custom_clust, by = "sample_id") %>% pull(cluster) %>% identity()
        #   Idents(object = object) <- clusters
        #   return(object)
        # }
        return(object)
    }

#' Read in Gene and Transcript Seurat Objects
#'
#' @param proj_dir path to project directory
#' @param prefix default "unfiltered"
#'
#' @return a single cell object
#' @export
load_object_path <- function(proj_dir = getwd(), prefix = "unfiltered") {
    object_regex <- paste0(paste0(".*/", prefix, "_object.rds"))

    object_path <- path(proj_dir, "output", "seurat") %>%
        dir_ls(regexp = object_regex)

    if (!rlang::is_empty(object_path)) {
        return(object_path)
    }

    stop("'", object_path, "' does not exist",
        paste0(" in current working directory ('", getwd(), "')"),
        ".",
        call. = FALSE
    )
}


#' Load Seurat Files from a signle project path
#'
#' @param proj_dir project directory
#' @param ... extra args passed to load_object_path
#'
#' @return a single cell object
load_object_from_proj <- function(proj_dir, ...) {
        object_file <- load_object_path(proj_dir, ...)
        object_file <- readRDS(object_file)
    }

#' Dimensional Reduction
#'
#' Run PCA, TSNE and UMAP on a singlecell objects
#' perplexity should not be bigger than 3 * perplexity < nrow(X) - 1, see details for interpretation
#'
#' @param object A Seurat object
#' @param assay Assay of interest to be run on the singlecell objects
#' @param reduction Set dimensional reduction object
#' @param legacy_settings Use legacy settings
#' @param ... Extra parameters passed to object_reduce_dimensions
#'
#' @return a single cell object with embeddings
#' @export
object_reduce_dimensions <- function(object, assay = "gene", reduction = "PCA", legacy_settings = FALSE, ...) {
        if ("integrated" %in% names(object@assays)) {
            assay <- "integrated"
        } else {
            assay <- "gene"
        }
        num_samples <- dim(object)[[2]]
        if (num_samples < 50) {
            npcs <- num_samples - 1
        } else {
            npcs <- 50
        }
        if (legacy_settings) {
            message("using legacy settings")

            if ("gene" == assay) {
                object <- scater::runPCA(x = object, subset_row = rownames(object))
            } else {
                object <- scater::runPCA(x = object, altexp = assay, subset_row = rownames(object))
            }
        } else {
            if ("gene" == assay) {
                object <- scater::runPCA(x = object, subset_row = scran::getTopHVGs(stats = object), ncomponents = npcs, ...)
            } else {
                object <- scater::runPCA(x = object, altexp = assay, subset_row = scran::getTopHVGs(stats = object), ncomponents = npcs, ...)
            }
        }
        if (reduction == "harmony") {
            object <- RunHarmony(object, "batch")
        }
        if ((ncol(object) - 1) > 3 * 30) {
            if ("gene" == assay) {
                object <- scater::runTSNE(x = object, dimred = reduction, n_dimred = 1:30)
            } else {
                object <- scater::runTSNE(x = object, altexp = assay, dimred = reduction, n_dimred = 1:30)
            }
            if ("gene" == assay) {
                object <- scater::runUMAP(x = object, dimred = reduction, n_dimred = 1:30)
            } else {
                object <- scater::runUMAP(x = object, altexp = assay, dimred = reduction, n_dimred = 1:30)
            }
        }
        return(object)
    }

#'
#' Give a new project name to a single cell object
#'
#' @param object A Seurat object
#' @param new_name New name to assign
#'
#' @return a renamed single cell object
#' @export
rename_object <- function (object, new_name)
          {
            metadata(object)["project.name"] <- new_name
            return(object)
          }

#' Filter a List of Seurat Objects
#'
#' Filter Seurat Objects by custom variable and reset assay to uncorrected "gene"
#'
#' @param objects single cell projects
#' @param filter_var filter variable
#' @param filter_val filter values
#' @param .drop whether to drop from single cell object
#'
#' @return a list of single cell objects
#' @export
filter_merged_objects <- function(objects, filter_var, filter_val, .drop = F) {
    objects <- map(objects, ~ filter_merged_object(object = .x, filter_var = filter_var, filter_val = filter_val, .drop = .drop))
}


#' Filter a Single Seurat Object
#'
#' @param object A singlecell objects
#' @param filter_var filter variable
#' @param filter_val filter values
#' @param .drop whether to drop values
#'
#' @return a single cell object
#' @export
filter_merged_object <- function(object, filter_var, filter_val, .drop = .drop) {
    if (.drop) {
        mycells <- get_cell_metadata(object)[[filter_var]] == filter_val
    } else {
        mycells <- get_cell_metadata(object)[[filter_var]] == filter_val | is.na(get_cell_metadata(object)[[filter_var]])
    }
    mycells <- colnames(object)[mycells]
    object <- object[, mycells]
    return(object)
}


#' Reintegrate (filtered) singlecell objectss
#'
#' 1) split by batch
#' 2) integrate
#' 3) run integration pipeline and save
#'
#' @param object A singlecell objects
#' @param feature gene or transcript
#' @param suffix to be appended to file saved in output dir
#' @param reduction to use default is pca
#' @param algorithm 1
#' @param ... extra args passed to object_integration_pipeline
#'
#' @return a single cell object
#' @export
reintegrate_object <- function (object, feature = "gene", suffix = "", reduction = "PCA", algorithm = 1, ...)
          {
            # Seurat::DefaultAssay(object) <- "gene"
            organism <- metadata(object)$experiment$organism
            experiment_name <- metadata(object)$experiment$experiment_name
            # object <- Seurat::DietSeurat(object, counts = TRUE, data = TRUE, scale.data = FALSE)
            objects <- batchelor::divideIntoBatches(object, object$batch, byrow = FALSE, restrict = NULL)[["batches"]]
            object <- object_integration_pipeline(objects, feature = feature, suffix = suffix, algorithm = algorithm, ...)
            object <- record_experiment_data(object, experiment_name, organism)
          }
