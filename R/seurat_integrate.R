old_harmony_integrate <- function (seu_list)
  {
    seu_list.integrated <- purrr::reduce(seu_list, merge)
    seu_list.integrated <- seurat_preprocess(seu_list.integrated)
    seu_list.integrated <- seurat_reduce_dimensions(seu_list.integrated)
    seu_list.integrated <- RenameAssays(seu_list.integrated,
                                        gene = "RNA")
    seu_list.integrated@assays[["integrated"]] <- seu_list.integrated@assays[["RNA"]]
    seu_list.integrated <- harmony::RunHarmony(seu_list.integrated,
                                               group.by.vars = "batch")
    seu_list.integrated <- RunUMAP(seu_list.integrated, reduction = "harmony",
                                   dims = 1:30)
    seu_list.integrated <- FindNeighbors(seu_list.integrated,
                                         reduction = "harmony", dims = 1:30) %>% FindClusters()
    seu_list.integrated <- RenameAssays(seu_list.integrated,
                                        RNA = "gene")
    seu_list.integrated
  }

#' Integrate small datasets with harmony
#'
#' @param seu_list List of two or more seurat objects to integrate
#'
#' @return
#' @export
#'
#' @examples
harmony_integrate <- function(seu_list){
  seu_list.integrated <- purrr::reduce(seu_list, merge)
  seu_list.integrated@assays[["integrated"]] <- seu_list.integrated@assays[["gene"]]
  DefaultAssay(seu_list.integrated) <- "integrated"
  seu_list.integrated <- seurat_preprocess(seu_list.integrated)
  seu_list.integrated <- seurat_reduce_dimensions(seu_list.integrated)
  seu_list.integrated <- harmony::RunHarmony(seu_list.integrated, group.by.vars = "batch", assay.use = "integrated")
  seu_list.integrated <- RunUMAP(seu_list.integrated, assay = "integrated", reduction = "harmony", dims = 1:30)
  seu_list.integrated <- FindNeighbors(seu_list.integrated, assay = "integrated", reduction = "harmony", dims = 1:30)
  seu_list.integrated
}

#' Merge Small Seurat Objects
#'
#' @param seu_list List of two or more seurat objects
#'
#' @return
#' @export
#'
#' @examples
merge_small_seus <- function(seu_list, k.filter = 50) {
  # check if any seurat objects are too small and if so merge with the first seurat object
  seu_dims <- purrr::map(seu_list, dim) %>%
    purrr::map_lgl(~ .x[[2]] < k.filter)

  small_seus <- seu_list[seu_dims]

  seu_list <- seu_list[!seu_dims]

  seu_list[[1]] <- purrr::reduce(c(small_seus, seu_list[[1]]), merge)

  return(seu_list)
}


#' Batch Correct Multiple Seurat Objects
#'
#' @param seu_list List of two or more seurat objects
#' @param method Default "cca"
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
seurat_integrate <- function(seu_list, method = "cca", organism = "human", ...) {

  # To construct a reference we will identify ‘anchors’ between the individual datasets. First, we split the combined object into a list, with each dataset as an element.

  # Prior to finding anchors, we perform standard preprocessing (log-normalization), and identify variable features individually for each. Note that Seurat v3 implements an improved method for variable feature selection based on a variance stabilizing transformation ("vst")

  for (i in 1:length(x = seu_list)) {
    seu_list[[i]][["gene"]] <- seurat_preprocess(seu_list[[i]][["gene"]], scale = TRUE)

    seu_list[[i]]$batch <- names(seu_list)[[i]]
  }

  seu_list <- merge_small_seus(seu_list)

  if (method == "rpca") {
    # scale and run pca for each separate batch in order to use reciprocal pca instead of cca
    features <- SelectIntegrationFeatures(object.list = seu_list)
    seu_list <- purrr::map(seu_list, Seurat::ScaleData, features = features)
    seu_list <- purrr::map(seu_list, Seurat::RunPCA, features = features)
    seu_list.anchors <- FindIntegrationAnchors(object.list = seu_list, reduction = "rpca", dims = 1:30)
  } else if (method == "cca") {
    # Next, we identify anchors using the FindIntegrationAnchors function, which takes a list of Seurat objects as input.
    seu_list.anchors <- Seurat::FindIntegrationAnchors(object.list = seu_list, dims = 1:30, k.filter = 50)
  }

  # proceed with integration
  # seu_list.integrated <- IntegrateData(anchorset = seu_list.anchors, dims = 1:30)

  # see https://github.com/satijalab/seurat/issues/6341------------------------------
  cells_per_batch <- sapply(seu_list, ncol)
  min_k_weight = min(cells_per_batch) - 1
  min_k_weight <- ifelse(min_k_weight < 100, min_k_weight, 100)

  seu_list.integrated <- tryCatch(IntegrateData(anchorset = seu_list.anchors, dims = 1:30, k.weight = min_k_weight), error = function(e) e)
  run_harmony <- any(class(seu_list.integrated) == "error")
  if(run_harmony){
    seu_list.integrated <- harmony_integrate(seu_list)
  }

  #   enriched_seu <- tryCatch(getEnrichedPathways(seu), error = function(e) e)
  #   enrichr_available <- !any(class(enriched_seu) == "error")
  #   if(enrichr_available){
  #     seu <- enriched_seu
  #   }

  # Next, we identify anchors using the FindIntegrationAnchors function, which takes a list of Seurat objects as input.

  # #stash batches
  Idents(seu_list.integrated) <- "batch"
  seu_list.integrated[["batch"]] <- Idents(seu_list.integrated)

  # switch to integrated assay. The variable features of this assay are
  # automatically set during IntegrateData
  Seurat::DefaultAssay(object = seu_list.integrated) <- "integrated"

  # if not integrated with harmony run the standard workflow for visualization and clustering
  if (!"harmony" %in% names(seu_list.integrated@reductions)){
    seu_list.integrated <- Seurat::ScaleData(object = seu_list.integrated, verbose = FALSE)
    seu_list.integrated <- seurat_reduce_dimensions(seu_list.integrated, ...)
  }

  seu_list.integrated <- record_experiment_data(seu_list.integrated, experiment_name = "integrated", organism = organism)

  return(seu_list.integrated)
}


#' Run Louvain Clustering at Multiple Resolutions
#'
#' @param seu A seurat object
#' @param resolution Clustering resolution
#' @param custom_clust
#' @param reduction Set dimensional reduction object
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
seurat_cluster <- function(seu = seu, resolution = 0.6, custom_clust = NULL, reduction = "pca", algorithm = 1, ...) {
  message(paste0("[", format(Sys.time(), "%H:%M:%S"), "] Clustering Cells..."))
  seu <- FindNeighbors(object = seu, dims = 1:30, reduction = reduction)

  if (length(resolution) > 1) {
    for (i in resolution) {
      message(paste0("clustering at ", i, " resolution"))
      seu <- Seurat::FindClusters(object = seu, resolution = i, algorithm = algorithm, ...)
    }
  } else if (length(resolution) == 1) {
    message(paste0("clustering at ", resolution, " resolution"))
    seu <- Seurat::FindClusters(object = seu, resolution = resolution, algorithm = algorithm, ...)
  }

  if (!is.null(custom_clust)) {
    seu <- Seurat::StashIdent(object = seu, save.name = "old.ident")
    clusters <- tibble::tibble("sample_id" = rownames(seu[[]])) %>%
      tibble::rownames_to_column("order") %>%
      dplyr::inner_join(custom_clust, by = "sample_id") %>%
      dplyr::pull(cluster) %>%
      identity()

    Idents(object = seu) <- clusters


    return(seu)
  }

  return(seu)
}

#' Read in Gene and Transcript Seurat Objects
#'
#' @param proj_dir path to project directory
#' @param prefix default "unfiltered"
#'
#' @return
#' @export
#'
#'
#' @examples
load_seurat_path <- function(proj_dir = getwd(), prefix = "unfiltered") {
  seu_regex <- paste0(paste0(".*/", prefix, "_seu.rds"))

  seu_path <- fs::path(proj_dir, "output", "seurat") %>%
    fs::dir_ls(regexp = seu_regex)

  if (!rlang::is_empty(seu_path)) {
    return(seu_path)
  }

  stop("'", seu_path, "' does not exist",
    paste0(" in current working directory ('", getwd(), "')"),
    ".",
    call. = FALSE
  )
}


#' Load Seurat Files from a signle project path
#'
#' @param proj_dir
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
load_seurat_from_proj <- function(proj_dir, ...) {
  seu_file <- load_seurat_path(proj_dir, ...)

  seu_file <- readRDS(seu_file)
}

#' Dimensional Reduction
#'
#' Run PCA, TSNE and UMAP on a seurat object
#' perplexity should not be bigger than 3 * perplexity < nrow(X) - 1, see details for interpretation
#'
#' @param seu A Seurat object
#' @param assay Assay of interest to be run on the seurat object
#' @param reduction Set dimensional reduction object
#' @param legacy_settings Use legacy settings
#' @param ... Extra parameters passed to seurat_reduce_dimensions
#'
#' @return
#' @export
#'
#' @examples
seurat_reduce_dimensions <- function(seu, assay = "gene", reduction = "pca", legacy_settings = FALSE, ...) {
  if ("integrated" %in% names(seu@assays)) {
    assay <- "integrated"
  } else {
    assay <- "gene"
  }

  num_samples <- dim(seu)[[2]]

  if (num_samples < 50) {
    npcs <- num_samples - 1
  } else {
    npcs <- 50
  }

  if (legacy_settings) {
    message("using legacy settings")
    seu <- Seurat::RunPCA(seu, assay = assay, features = rownames(seu))
  } else {
    # seu <- Seurat::RunPCA(object = seu, do.print = FALSE, npcs = npcs, ...)
    seu <- Seurat::RunPCA(object = seu, assay = assay, features = Seurat::VariableFeatures(object = seu), do.print = FALSE, npcs = npcs, ...)
  }

  if (reduction == "harmony") {
    seu <- harmony::RunHarmony(seu, "batch")
  }

  if ((ncol(seu) - 1) > 3 * 30) {
    seu <- Seurat::RunTSNE(object = seu, assay = assay, reduction = reduction, dims = 1:30)
    seu <- Seurat::RunUMAP(object = seu, assay = assay, reduction = reduction, dims = 1:30)
  }

  return(seu)
}

#'
#' Give a new project name to a seurat object
#'
#' @param seu A Seurat object
#' @param new_name New name to assign
#'
#' @return
#' @export
#'
#' @examples
#'
rename_seurat <- function(seu, new_name) {
  seu@project.name <- new_name
  return(seu)
}

#' Reset default assay
#'
#' Function to reset the default assay of a seurat object
#'
#' @param seu A seurat object
#' @param new_assay Assay to set as default assay
#'
#' @return
#' @export
#'
#' @examples
SetDefaultAssay <- function(seu, new_assay) {
  Seurat::DefaultAssay(seu) <- new_assay
  return(seu)
}



#' Filter a List of Seurat Objects
#'
#' Filter Seurat Objects by custom variable and reset assay to uncorrected "gene"
#'
#' @param seus
#' @param filter_var
#' @param filter_val
#' @param .drop
#'
#' @return
#' @export
#'
#' @examples
filter_merged_seus <- function(seus, filter_var, filter_val, .drop = F) {
  seus <- purrr::map(seus, ~ filter_merged_seu(seu = .x, filter_var = filter_var, filter_val = filter_val, .drop = .drop))
}


#' Filter a Single Seurat Object
#'
#' @param seu A seurat object
#' @param filter_var
#' @param filter_val
#' @param .drop
#'
#' @return
#' @export
#'
#' @examples
filter_merged_seu <- function(seu, filter_var, filter_val, .drop = .drop) {
  if (.drop) {
    mycells <- seu[[]][[filter_var]] == filter_val
  } else {
    mycells <- seu[[]][[filter_var]] == filter_val | is.na(seu[[]][[filter_var]])
  }
  mycells <- colnames(seu)[mycells]
  seu <- seu[, mycells]
  return(seu)
}


#' Reintegrate (filtered) seurat objects
#'
#' 1) split by batch
#' 2) integrate
#' 3) run integration pipeline and save
#'
#' @param seu A seurat object
#' @param suffix to be appended to file saved in output dir
#' @param reduction to use default is pca
#'
#' @return
#' @export
#'
#' @examples
#' panc8$batch <- panc8$gene$tech
#' reintegrate_seu(panc8)
reintegrate_seu <- function(seu, feature = "gene", suffix = "", reduction = "pca", algorithm = 1, ...) {
  Seurat::DefaultAssay(seu) <- "gene"

  organism <- Misc(seu)$experiment$organism
  experiment_name <- Misc(seu)$experiment$experiment_name

  seu <- Seurat::DietSeurat(seu, counts = TRUE, data = TRUE, scale.data = FALSE)
  seus <- Seurat::SplitObject(seu, split.by = "batch")
  seu <- seurat_integration_pipeline(seus, feature = feature, suffix = suffix, algorithm = algorithm, ...)

  seu <- record_experiment_data(seu, experiment_name, organism)

  # integration_workflow <- function(batches, excluded_cells = NULL, resolution = seq(0.2, 2.0, by = 0.2), experiment_name = "default_experiment", organism = "human", ...) {
}
