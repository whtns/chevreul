#' Clustering Workflow
#'
#' Cluster and Reduce Dimensions of a seurat object
#'
#' @param feature_seus list of seurat objedevelcts named according to feature of interest ("gene" or "transcript")
#' @param excluded_cells named list of cells to exclude
#' @param resolution resolution(s) to use for clustering cells
#' @param organism
#' @param experiment_name
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' clustered_seu <- clustering_workflow(seurat_pancrease_reduced)
mod_clustering_workflow <- function(seu = NULL, excluded_cells, resolution = seq(0.2, 2.0, by = 0.2), organism = "human", experiment_name = "default_experiment", ...){

  seu <- mod_seurat_pipeline(seu, resolution = resolution, organism = organism, ...)

  seu <- record_experiment_data(seu, experiment_name, organism)

  # save_seurat(feature_seus, proj_dir = proj_dir, ...)

}

#' Run Seurat Pipeline
#'
#' Preprocess, Cluster and Reduce Dimensions for a single seurat object
#'
#' @param seu
#' @param resolution
#'
#' @return
#' @export
#'
#' @examples
mod_seurat_pipeline <- function(seu, feature = "none", resolution=0.6, reduction = "pca", annotate_cell_cycle = TRUE, annotate_percent_mito = TRUE, organism = "human", ...){

  seu <- mod_seurat_preprocess(seu, scale = T, ...)

  # PCA
  seu <- mod_seurat_reduce_dimensions(seu, check_duplicates = FALSE, ...)

  seu <- mod_seurat_cluster(seu = seu, resolution = resolution, reduction = reduction, ...)

  seu_assays = Seurat::Assays(seu)
  names(seu_assays) = seu_assays

  seu@misc$markers <- purrr::map(seu_assays, ~mod_find_all_markers(seu, assay = .x, resolution = resolution))

  if (feature == "gene"){
    enriched_seu <- tryCatch(getEnrichedPathways(seu), error = function(e) e)
    enrichr_available <- !any(class(enriched_seu) == "error")
    if(enrichr_available){
      seu <- enriched_seu
    }
  }

  # annotate low read count category in seurat metadata
  # seu <- seuratTools::add_read_count_col(seu)

  # annotate cell cycle scoring to seurat objects
  if (annotate_cell_cycle){
    seu <- mod_annotate_cell_cycle(seu, assay = "gene", ...)
  }

  # annotate mitochondrial percentage in seurat metadata
  if (annotate_percent_mito){
    seu <- mod_add_percent_mito(seu, assay = "gene", organism)
  }

  return(seu)
}

#' Dimensional Reduction
#'
#' Run PCA, TSNE and UMAP on a seurat object
#' perplexity should not be bigger than 3 * perplexity < nrow(X) - 1, see details for interpretation
#'
#' @param seu
#' @param assays
#' @param legacy_settings
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
mod_seurat_reduce_dimensions <- function(seu, assays = NULL, legacy_settings = FALSE, ...) {

  assays %||% Seurat::Assays(seu)

  num_samples <- dim(seu)[[2]]

  if (num_samples < 50){
    npcs = num_samples - 1
  } else {
    npcs = 50
  }

  if (legacy_settings){
    seu <- Seurat::RunPCA(seu, features = rownames(seu))
  } else {
    for (i in assays){
      reduction.name = paste0(i,".pca")
      seu <- Seurat::RunPCA(object = seu, assay = i, do.print = FALSE, npcs = npcs, reduction.name = reduction.name, ...)
    }
  }

  if ((ncol(seu) -1) > 3*30){
    for (i in assays){
      reduction.name = paste0(i,".tsne")
      seu <- Seurat::RunTSNE(object = seu, reduction = paste0(i,".pca"), dims = 1:30, reduction.name = reduction.name, ...)
    }

    for (i in assays){
      reduction.name = paste0(i,".umap")
      seu <- Seurat::RunUMAP(object = seu, reduction = paste0(i,".pca"), dims = 1:30, reduction.name = reduction.name, ...)
    }
  }

  return(seu)

}

#' Run Louvain Clustering at Multiple Resolutions
#'
#' @param seu
#' @param resolution
#' @param custom_clust
#' @param reduction
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
mod_seurat_cluster <- function(seu = seu, assays = NULL, resolution = 0.6, custom_clust = NULL, reduction = "pca", algorithm = 1, ...){

  assays %||% Seurat::Assays(seu)

  message(paste0("[", format(Sys.time(), "%H:%M:%S"), "] Clustering Cells..."))
  for (i in assays){
    reduction.name = paste0(i, ".pca")
    seu <- FindNeighbors(object = seu, assay = i, reduction = reduction.name, dims = 1:10)
  }

  if (length(resolution) > 1){
    for (i in assays){
      graph.name = paste0(i,"_snn")
      for (j in resolution){
        message(paste0("clustering ", i, " at ", j, " resolution"))
        seu <- Seurat::FindClusters(object = seu, graph.name = graph.name, resolution = j, algorithm = algorithm, ...)
      }
    }

  } else if (length(resolution) == 1){
    message(paste0("clustering at ", resolution, " resolution"))
    for (i in assays){
      graph.name = paste0(i,"_nn")
      seu <- Seurat::FindClusters(object = seu, graph.name = graph.name, resolution = resolution, algorithm = algorithm, ...)
    }
  }

  return(seu)
}

#' Find All Markers at a range of resolutions
#'
#' @param seu
#' @param metavar
#'
#' @return
#' @export
#'
#' @examples
mod_find_all_markers <- function(seu, assay, metavar = NULL, ...){
  if (is.null(metavar)){
    if ("integrated" %in% names(seu@assays)) {
      default_assay = "integrated"
    } else {
      default_assay = assay
    }

    resolutions <- colnames(seu[[]])[grepl(paste0(default_assay, "_snn_res."), colnames(seu[[]]))]

    cluster_index <- grepl(paste0(default_assay, "_snn_res."), colnames(seu[[]]))

    if(!any(cluster_index)) {
      stop("no clusters found in metadata. Please run seurat_cluster")
    }

    clusters <- seu[[]][,cluster_index]

    cluster_levels <- purrr::map_int(clusters, ~length(unique(.x)))
    cluster_levels <- cluster_levels[cluster_levels > 1]

    clusters <- dplyr::select(clusters, dplyr::one_of(names(cluster_levels)))
    metavar = names(clusters)
  }

  new_markers <- purrr::map(metavar, mod_stash_marker_features, seu, assay)
  names(new_markers) <- metavar

  old_markers <- seu@misc$markers[!names(seu@misc$markers) %in% names(new_markers)]

  markers <- c(old_markers, new_markers)

  return(markers)

}

#' Preprocess Seurat Object
#'
#' @param seu
#' @param scale
#' @param normalize
#' @param features
#' @param legacy_settings
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
mod_seurat_preprocess <- function(seu, scale=TRUE, normalize = TRUE, features = NULL, legacy_settings = FALSE, ...){
  # Normalize data

  if (legacy_settings){

    logtransform_exp <- as.matrix(log1p(Seurat::GetAssayData(seu)))

    seu <- Seurat::SetAssayData(seu, slot = "data", logtransform_exp) %>%
      Seurat::ScaleData(features = rownames(.))

    return(seu)
  }

  if (normalize){
    for (i in Seurat::Assays(seu)){
      seu <- Seurat::NormalizeData(object = seu, assay = i, verbose = FALSE, ...)
    }
  }

  # Filter out only variable genes
  for (i in Seurat::Assays(seu)){
    seu <- Seurat::FindVariableFeatures(object = seu, assay = i, selection.method = "vst", verbose = FALSE, ...)
  }

  # Regress out unwanted sources of variation
  if (scale){
    for (i in Seurat::Assays(seu)){
      seu <- Seurat::ScaleData(object = seu, assay = i, ...)
    }

  }

  return(seu)
}

#' Stash Marker Genes in a Seurat Object
#'
#' Marker Genes will be stored in slot `@misc$markers`
#'
#' @param metavar
#' @param seu
#' @param top_n
#'
#' @return
#' @export
#'
#' @examples
mod_stash_marker_features <- function(metavar, seu, assay, top_n = 200){
  markers <- list()
  markers$presto <- presto::wilcoxauc(seu, metavar, seurat_assay = assay) %>%
    dplyr::group_by(group) %>%
    dplyr::filter(padj < 0.05) %>%
    dplyr::top_n(n = top_n, wt = logFC) %>%
    dplyr::arrange(group, desc(logFC)) %>%
    dplyr::select(feature, group) %>%
    dplyr::mutate(rn = row_number()) %>%
    tidyr::pivot_wider(names_from = group, values_from = feature) %>%
    dplyr::select(-rn)

  markers$genesorteR <- genesorteR::sortGenes(
    GetAssayData(seu, assay = assay, slot = "data"),
    seu[[]][[metavar]]
  )

  markers$genesorteR <-
    apply(markers$genesorteR$specScore, 2, function(x) names(head(sort(x, decreasing = TRUE), n = top_n))) %>%
    as.data.frame()

  return(markers)

}

#' Annotate Cell Cycle
#'
#' Annotate Cell Cycle for Gene and Transcript Seurat Objects
#'
#' @param seu_list A list of seurat objects
#'
#' @return
#' @export
#'
#' @examples
mod_annotate_cell_cycle <- function(seu, assay, organism = "human", ...){


  # setdefaultassay to "RNA"
  Seurat::DefaultAssay(seu) <- assay

  s_genes <- cc.genes$s.genes
  g2m_genes <- cc.genes$g2m.genes

  s_transcripts <- cc.transcripts$s.genes
  g2m_transcripts <- cc.transcripts$g2m.genes

  if(organism == "mouse"){
    s_genes <- stringr::str_to_title(s_genes)
    g2m_genes <- stringr::str_to_title(g2m_genes)
    s_transcripts <- genes_to_transcripts(s_genes, organism = organism)
    g2m_transcripts <- genes_to_transcripts(g2m_genes, organism = organism)
  }

  if (assay == "gene"){
    seu <- CellCycleScoring(seu, s.features = s_genes, g2m.features = g2m_genes, set.ident = FALSE)

  } else if (assay == "transcript"){
    seu <- CellCycleScoring(seu, s.features = s_transcripts, g2m.features = g2m_transcripts, set.ident = FALSE)
  }

}

#' Annotate Low Read Count Category
#'
#'  Add a Read Count Categorical Variable to Seurat Object (based on nCount_RNA)
#'
#' @param seu A seurat object
#' @param feature
#' @param organism
#'
#' @return
#' @export
#'
#' @examples
mod_add_percent_mito <- function(seu, assay, organism){

  DefaultAssay(seu) <- assay

  mito_features <- mito_features[[organism]][[assay]]

  mito_features <- mito_features[mito_features %in% rownames(seu)]

  seu[["percent.mt"]] <- PercentageFeatureSet(seu, features = mito_features)

  return(seu)

}

#' Integration Workflow
#'
#' Integrate multiple seurat objects and save to file
#'
#' @param batches seurat objects for each all batches provided as a list. If named, the resulting integrated object will be identified with corresponding values in 'batch' metadata
#' @param excluded_cells named list of cells to exclude
#' @param resolution value(s) to control the clustering resolution via `Seurat::FindMarkers`
#' @param experiment_name arbitrary name to identify experiment
#' @param organism either "human" or "mouse"
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
#' batches <- seurat_pancreas_reduced %>%
#'   purrr::map(Seurat::SplitObject, split.by = "dataset") %>%
#'   purrr::transpose()
#'
#' inegrated_seu <- integration_workflow(batches)
mod_integration_workflow <- function(batches, excluded_cells, resolution = seq(0.2, 2.0, by = 0.2), experiment_name = "default_experiment", organism = "human", ...) {

  # names(child_proj_dirs) <- gsub("_proj", "", fs::path_file(child_proj_dirs))

  # load seurat objects from 'child' projects
  # seus <- purrr::map(child_proj_dirs, load_seurat_from_proj)

  organisms <- purrr::map(batches, list("meta.data", "organism", 1))

  if (any(map_lgl(organisms, is.null))){
    organisms <- case_when(
      grepl("Hs", names(batches)) ~ "human",
      grepl("Mm", names(batches)) ~ "mouse"
    )
    names(organisms) <- names(batches)
  }

  experiment_names <- names(batches)

  batches <- purrr::pmap(list(batches, experiment_names, organisms), record_experiment_data)

  #
  # batches <- purrr::transpose(batches)
  # for (i in names(batches)){
  #   batches[[i]] <- purrr::pmap(list(batches[[i]], experiment_names, organisms), record_experiment_data)
  # }
  # batches <- purrr::transpose(batches)

  if (all(purrr::map(batches, list("misc", "experiment", "organism")) == "human")){
    integrated_seu <- mod_seurat_integration_pipeline(batches, resolution = resolution, organism = "human", ...)

    integrated_seu@misc$batches <- names(batches)

  } else if (all(purrr::map(batches, list("misc", "experiment", "organism")) == "mouse")){
    integrated_seu <- mod_seurat_integration_pipeline(batches, resolution = resolution, organism = "mouse", ...)

    integrated_seu@misc$batches <- names(batches)

  }  else {

    # mouse_seu_list <- batches[grepl("Mm", names(batches))]
    mouse_seu_list <- batches[names(organisms[organisms == "mouse"])]
    # human_seu_list <- batches[grepl("Hs", names(batches))]
    human_seu_list <- batches[names(organisms[organisms == "human"])]
    merged_batches <- mod_cross_species_integrate(mouse_seu_list = mouse_seu_list, human_seu_list = human_seu_list)
    for (i in names(merged_batches)){
      merged_batches[[i]]@misc$batches <- names(batches)
    }
  }

  integrated_seu <- record_experiment_data(integrated_seu, experiment_name, organism)

  return(integrated_seu)

}

#' Run Seurat Integration
#'
#' run batch correction, followed by:
#' 1) stashing of batches in metadata 'batch'
#' 2) clustering with resolution 0.2 to 2.0 in increments of 0.2
#' 3) saving to <proj_dir>/output/sce/<feature>_seu_<suffix>.rds
#'
#' @param suffix a suffix to be appended to a file save in output dir
#' @param seu_list
#' @param resolution
#' @param feature
#' @param algorithm
#' @param organism
#' @param ...
#'
#' @return
#' @export
#'
#'
#' @examples
mod_seurat_integration_pipeline <- function(seu_list, feature = "none", resolution = seq(0.2, 2.0, by = 0.2), suffix = '', algorithm = 1, organism, annotate_cell_cycle = FALSE, annotate_percent_mito = FALSE, ...) {

  integrated_seu <- mod_seurat_integrate(seu_list, ...)

  integrated_assays <- stringr::str_subset(Seurat::Assays(integrated_seu), "integrated")
  names(integrated_assays) = integrated_assays

  # cluster merged seurat objects
  integrated_seu <- mod_seurat_cluster(integrated_seu, assays = integrated_assays, resolution = resolution, algorithm = algorithm, ...)

  integrated_seu@misc$markers <- purrr::map(integrated_assays, ~mod_find_all_markers(integrated_seu, assay = .x, resolution = resolution))

  if (feature == "gene"){
    enriched_seu <- tryCatch(getEnrichedPathways(integrated_seu), error = function(e) e)
    enrichr_available <- !any(class(enriched_seu) == "error")
    if(enrichr_available){
      integrated_seu <- enriched_seu
    }
  }

  # add read count column
  # integrated_seu <- mod_add_read_count_col(integrated_seu)

  # annotate cell cycle scoring to seurat objects
  if (annotate_cell_cycle){
    integrated_seu <- mod_annotate_cell_cycle(integrated_seu, assay = "gene.integrated", ...)
  }

  # annotate mitochondrial percentage in seurat metadata
  if (annotate_percent_mito){
    integrated_seu <- mod_add_percent_mito(integrated_seu, assay = "gene.integrated", ...)
  }

  #annotate excluded cells
  # integrated_seu <- annotate_excluded(integrated_seu, excluded_cells)

  return(integrated_seu)

}

#' Batch Correct Multiple Seurat Objects
#'
#' @param seu_list
#' @param method
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
mod_seurat_integrate <- function(seu_list, method = "cca", ...) {
  #
  # To construct a reference we will identify ‘anchors’ between the individual datasets. First, we split the combined object into a list, with each dataset as an element.

  # Prior to finding anchors, we perform standard preprocessing (log-normalization), and identify variable features individually for each. Note that Seurat v3 implements an improved method for variable feature selection based on a variance stabilizing transformation ("vst")
  seu_list <- purrr::map(seu_list, mod_seurat_preprocess, scale = TRUE, ...)


  # for (i in 1:length(x = seu_list)) {
  #   seu_list[[i]] <- seurat_preprocess(seu_list[[i]], scale = TRUE, ...)
  #   seu_list[[i]]$batch <- names(seu_list)[[i]]
  # }

  seu_list <- merge_small_seus(seu_list)


  assays_in_common <- reduce(map(seu_list, Assays), intersect)

  seu_list.integrated = list()
  for (i in assays_in_common){
    # Next, we identify anchors using the FindIntegrationAnchors function, which takes a list of Seurat objects as input.
    seu_list.anchors <- Seurat::FindIntegrationAnchors(object.list = seu_list, assay = rep(i, length(seu_list)), dims = 1:30, k.filter = 50)
    # proceed with integration
    seu_list.integrated[[i]]  <- IntegrateData(anchorset = seu_list.anchors, dims = 1:30, new.assay.name = paste0(i,".integrated"))
  }

  integrated_seu <- seu_list.integrated$gene
  integrated_seu[["transcript.integrated"]] <- seu_list.integrated$transcript[["transcript.integrated"]]

  # Next, we identify anchors using the FindIntegrationAnchors function, which takes a list of Seurat objects as input.

  # #stash batches
  Idents(integrated_seu) <- "batch"
  integrated_seu[["batch"]] <- Idents(integrated_seu)

  # switch to integrated assay. The variable features of this assay are
  # automatically set during IntegrateData
  # Seurat::DefaultAssay(object = seu_list.integrated) <- "integrated"


  integrated_assays <- stringr::str_subset(Seurat::Assays(integrated_seu), "integrated")

  for (i in integrated_assays){
    # Run the standard workflow for visualization and clustering
    integrated_seu <- Seurat::ScaleData(object = integrated_seu, assay = i, verbose = FALSE)
  }

  integrated_seu <- mod_seurat_reduce_dimensions(integrated_seu, assays = integrated_assays, ...)

  return(integrated_seu)
}


#' Title
#'
#' @param seu_list
#' @param feature
#' @param resolution
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
mod_update_seuratTools_object <- function(seu_list, feature, resolution = seq(0.2, 2.0, by = 0.2), ...){


  seu_list <- purrr::imap(seu_list, conform_assay_name)

  integrated_seu <- seu_list$gene
  integrated_seu@reductions$tsne <- NULL
  integrated_seu@reductions$umap <- NULL
  integrated_seu[["transcript"]] <- seu_list$transcript[["transcript"]]
  integrated_seu[["transcript.integrated"]] <- seu_list$transcript[["transcript.integrated"]]

  # #stash batches
  Idents(integrated_seu) <- "batch"
  integrated_seu[["batch"]] <- Idents(integrated_seu)

  # switch to integrated assay. The variable features of this assay are
  # automatically set during IntegrateData
  # Seurat::DefaultAssay(object = seu_list.integrated) <- "integrated"

  integrated_assays <- stringr::str_subset(Seurat::Assays(integrated_seu), "integrated")

  integrated_seu <- mod_seurat_reduce_dimensions(integrated_seu, assays = integrated_assays, ...)

  seuratTools_version <- integrated_seu@misc$experiment$technical_info$seuratTools_version

  seuratTools_version <- ifelse(is.null(seuratTools_version), 0.1, seuratTools_version)

  if(seuratTools_version < getNamespaceVersion("seuratTools")){

    if (!any(grepl("_snn_res", colnames(seu@meta.data)))){
      integrated_seu <- seurat_cluster(integrated_seu, resolution = resolution, reduction = "pca", ...)

    }
    integrated_seu <- find_all_markers(integrated_seu, resolution = resolution)
    integrated_seu <- record_experiment_data(integrated_seu)
    integrated_seu <- seu_calcn(integrated_seu)
  }

  DefaultAssay(integrated_seu) <- "gene.integrated"

  return(integrated_seu)

}

conform_assay_name <- function(seu, assay){

  seu[[assay]] <- seu[["RNA"]]
  integrated_assay = paste0(assay, ".integrated")
  seu[[integrated_assay]] <- seu[["integrated"]]

  DefaultAssay(seu) <- assay
  seu[["RNA"]] <- NULL
  seu[["integrated"]] <- NULL

  # names(seu@reductions) <- paste0(assay, ".", names(seu@reductions))

  return(seu)

}
