#' scvelo_assay
#'
#' run scvelo on a gene or transcript level object
#'
#' @param object a object
#' @param loom_path path to matching loom file
#' @param fit.quantile
#' @param ...
#'
#' @return a single cell object with RNA velocity calculated
#' @export
#'
#' @examples
run_scvelo <- function(object, loom_path, assay = "gene", fit.quantile = 0.05, check_loom = FALSE, ...) {
    # if(DefaultAssay(object) == "SCT"){
    #   object <-
    #     object %>%
    #     Seurat::FindVariableFeatures(nfeatures= 3000)
    # }

    ldat <- SeuratWrappers::ReadVelocity(file = loom_path)

    for (assay in names(ldat)) {
        colnames(ldat[[assay]]) <- str_replace(colnames(ldat[[assay]]), ".*:", "")
        colnames(ldat[[assay]]) <- str_replace(colnames(ldat[[assay]]), "x$", "-1")
    }

    bm <- Seurat::as.Seurat(x = ldat)
    bm <- bm[rownames(bm) %in% rownames(object), ]

    bm[[assay]] <- bm[["spliced"]]

    # subset bm by object.size
    bm <- bm[, colnames(bm) %in% colnames(object)]

    # subset object by ldat
    sub_object <- object[, colnames(object) %in% colnames(bm)]

    sub_object@assays[names(bm@assays)] <- bm@assays
    DefaultAssay(sub_object) <- assay
    Misc(sub_object)$vel <- NULL
    Misc(sub_object)[names(Misc(sub_object)) == "experiment"] <- NULL

    # sub_object <- SeuratObject::RenameAssays(sub_object, gene = "RNA")

    h5ad_path <- str_replace(loom_path, ".loom", ".h5ad")

    # sceasy::convertFormat(sub_object, from="seurat", to="anndata",
    #                       outFile=fs::path_expand(h5ad_path))

    convert_to_h5ad(sub_object, file_path = loom_path)

    return(sub_object)
}

#' convert a object to an on-disk anndata object
#'
#' @param object A object
#' @param file_path Path to file
#'
#' @return a path to an h5ad file
#' @export
#'
#' @examples
#'
#' convert_to_h5ad(human_gene_transcript_object, "inst/extdata/object.rds")
#'
convert_to_h5ad <- function(object, file_path) {
    h5object_path <- fs::path_ext_set(file_path, ".h5Seurat")
    message(h5object_path)
    SeuratDisk::SaveH5Seurat(object, filename = h5object_path, overwrite = TRUE)

    h5ad_path <- fs::path_ext_set(file_path, ".h5ad")

    message(h5ad_path)
    SeuratDisk::Convert(h5object_path, dest = h5ad_path, overwrite = TRUE)
}

#' scvelo_assay
#'
#' run scvelo on a gene or transcript level object
#'
#' @param object a object
#' @param loom_path path to matching loom file
#' @param group.by metadata to color plot
#' @param plot_method plotting method to use from scvelo
#' @param ...
#'
#' @return a single cell object with velocity calculated
#' @export
#'
#' @examples
prep_scvelo <- function(object, loom_path, velocity_mode = c("deterministic", "stochastic", "dynamical"), ...) {
    h5ad_path <- fs::path_ext_set(loom_path, ".h5ad")
    message(h5ad_path)
    adata_matches_object <- function(object, adata) {
        identical(sort(adata$obs_names$values), sort(colnames(object)))
    }

    if (fs::file_exists(h5ad_path)) {
        adata <- scvelo$read(fs::path_expand(h5ad_path))

        if (!adata_matches_object(object, adata)) {
            object <- run_scvelo(object, loom_path, ...)
        }
    } else {
        object <- run_scvelo(object, loom_path, ...)
    }

    adata <- scvelo$read(fs::path_expand(h5ad_path))
    # reticulate::source_python("scripts/rename_raw.py")
    # adata$raw$var$rename(columns = list('_index' = 'symbol'), inplace = True)

    scvelo$pp$filter_and_normalize(adata, min_shared_counts = 20L, n_top_genes = 2000L)

    scvelo$pp$moments(adata, n_pcs = 30L, n_neighbors = 30L)

    if (velocity_mode == "dynamical") {
        if (!"recover_dynamics" %in% adata$uns_keys()) {
            scvelo$tl$recover_dynamics(adata)
            reticulate::py_del_attr(adata, "raw")
            adata$write_h5ad(h5ad_path)
        }
        scvelo$tl$latent_time(adata)
    }

    scvelo$tl$velocity(adata, mode = velocity_mode)

    scvelo$tl$velocity_graph(adata)


    return(adata)
}

#' Plot scvelo on embedding plot
#'
#' @param adata
#' @param group.by
#' @param plot_method
#'
#' @return a matplotlib of RNA velocity
#' @export
#'
#' @examples
plot_scvelo <- function(adata, group.by = "batch", basis = "umap", plot_method = c("stream", "arrow", "dynamics"), ...) {
    num_cols <- length(unique(adata$obs[[group.by]]))

    mycols <- scales::hue_pal()(num_cols)

    if (plot_method == "stream") {
        scvelo$pl$velocity_embedding_stream(adata, basis = basis, palette = mycols, color = group.by, dpi = 200, figsize = c(20, 12), ...)
    } else if (plot_method == "arrow") {
        scvelo$pl$velocity_embedding(adata, basis = basis, palette = mycols, color = group.by, arrow_length = 3, arrow_size = 2, dpi = 200, figsize = c(20, 12), ...)
    } else if (plot_method == "dynamics") {
        scvelo$pl$scatter(adata, color = "latent_time", color_map = "gnuplot", figsize = c(20, 12), dpi = 200, ...)
    }

    # pyplot$show()
}

scvelo_expression <- function(adata, features = c("RXRG")) {
    scvelo$pl$velocity(adata, var_names = features, figsize = c(10, 10), dpi = 200)

    # pyplot$show()
}
