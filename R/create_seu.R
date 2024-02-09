#' Set Column Names from `tximport`
#'
#' @param txi
#' @param colnames
#'
#' @return
#' @export
#'
#' @examples
set_colnames_txi <- function(txi, colnames) {
    colnames(txi$counts) <- colnames
    colnames(txi$abundance) <- colnames
    colnames(txi$length) <- colnames
    return(txi)
}

#' Run \href{http://bioconductor.org/packages/release/bioc/html/tximport.html}{tximport} on a set of cells
#'
#' cells can be quantified using:
#' \itemize{
#'   \item stringtie
#'   \item salmon
#' }
#' @param proj_dir project directory
#' @param type stringtie or salmon
#' @param countsFromAbundance argument provided to tximport
#' @param edb
#'
#' @return
#' @export
#'
#' @examples
load_counts_by_tximport <- function(proj_dir, type = "salmon", countsFromAbundance = "scaledTPM", edb = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86) {
    sample_glob <- switch(type,
        kallisto = "*abundance.h5",
        salmon = "*quant.sf",
        stringtie = "*t_data.ctab"
    )

    sample_paths <- rlang::with_handlers(error = ~ rlang::abort("Can't find input files",
        parent = .
    ), sample_files <- fs::path(
        proj_dir, "output",
        type
    ) %>% fs::dir_ls(recurse = T, glob = sample_glob) %>%
        identity())

    tx2gene <- ensembldb::transcripts(edb, return.type = "data.frame")[, c("tx_id", "gene_id")] %>%
        dplyr::left_join(annotables::grch38, by = c("gene_id" = "ensgene")) %>%
        dplyr::select(tx_id, symbol) %>%
        tidyr::drop_na()

    txi_transcripts <- tximport::tximport(sample_files, type = type, tx2gene = tx2gene, txOut = T, countsFromAbundance = countsFromAbundance, ignoreTxVersion = TRUE)

    # sanitize transcript ids with trailing (.1, .2, etc)
    txi_transcripts <- purrr::map_if(
        txi_transcripts, is.matrix,
        ~ `rownames<-`(.x, stringr::str_remove(rownames(.x), "\\.[0-9]$"))
    )

    txi_genes <- tximport::summarizeToGene(txi_transcripts, tx2gene = tx2gene, ignoreTxVersion = TRUE)

    txi_transcripts$tx2gene <- tx2gene

    sample_names <- fs::path_file(fs::path_dir(sample_paths))

    txi_features <- map(list(gene = txi_genes, transcript = txi_transcripts), ~ set_colnames_txi(.x, sample_names))
}


#' Load Sample Metadata for a given project
#'
#'
#'
#' @param proj_dir
#'
#' @return
#' @export
#'
#' @examples
load_meta <- function(proj_dir) {
    # load metadata
    meta_file <- gsub("_proj", "_metadata.csv", path_file(proj_dir))
    meta_file <- fs::path(proj_dir, "data", meta_file)

    tpm_meta <- read_csv(meta_file)
}


#' Create a object from output of  \href{http://bioconductor.org/packages/release/bioc/html/tximport.html}{tximport} and a table of cell metadata
#'
#' @param txi output from load_counts_by_tximport
#' @param meta_tbl a tibble of cell metadata with cell ids as the first column
#' @param feature the feature level on which to summarize counts gene or transcript
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
object_from_tximport <- function(txi, meta_tbl, ...) {
    gene_expression <- as.matrix(txi$gene$counts)
    expid <- gsub("-[0-9]*", "", colnames(gene_expression))

    featuredata <- data.frame(
        feature = rownames(gene_expression),
        row.names = rownames(gene_expression)
    )

    meta_tbl <- data.frame(meta_tbl,
        row.names = meta_tbl[["sample_id"]]
    )

    meta_tbl <- meta_tbl[colnames(gene_expression), ]

    # create gene assay
    object <- Seurat::CreateSeuratObject(counts = gene_expression, project = expid, assay = "gene", meta.data = meta_tbl)
    object@assays[["gene"]] <- AddMetaData(object@assays[["gene"]], featuredata)

    if ("transcript" %in% names(txi)) {
        # create transcript assay
        transcript_expression <- as.matrix(txi$transcript$counts)
        object[["transcript"]] <- CreateAssayObject(transcript_expression)
    }

    # add default batch if missing
    object$batch <- object@project.name

    return(object)
}


#' Create a Seurat Object from a set of tibbles
#'
#' @param exp_tbl
#' @param feature
#' @param meta_tbl
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
object_from_tibbles <- function(exp_tbl, feature, meta_tbl, ...) {
    expid <- gsub("-.*", "", colnames(exp_tbl))

    featuredata <- data.frame(rownames(exp_tbl))
    rownames(featuredata) <- featuredata[, 1]
    if (feature == "transcript") {
        # gene_id <- tx2gene
        # featuredata$gene_symbol =
    }

    meta_tbl <- data.frame(meta_tbl)
    rownames(meta_tbl) <- meta_tbl[, "sample_id"]

    meta_tbl <- meta_tbl[colnames(exp_tbl), ]

    object <- Seurat::CreateSeuratObject(counts = exp_tbl, project = expid, assay = "gene", meta.data = meta_tbl)

    # add default batch if missing
    object$batch <- object@project.name

    return(object)
}


## ------------------------------------------------------------------------
# filter out low read count cells (threshold 1e5)

#' Filter our Cells from Seurat below read count threshold
#'
#' @param object A object
#' @param read_thresh
#'
#' @return
#' @export
#'
#' @examples
filter_low_rc_cells <- function(object, read_thresh = 1e5) {
    counts <- as.matrix(object@assays[["gene"]]@counts)

    counts <- colSums(counts)

    keep_cells <- counts[counts > read_thresh]

    removed_cells <- counts[counts <= read_thresh]
    print(removed_cells)

    object <- subset(object, cells = names(keep_cells))
}

#' Save object to <project>/output/sce/<feature>_object.rds
#'
#' @param ... named arguments specifying objects list of objects; default "gene" and "transcript"
#' @param prefix
#' @param proj_dir
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' save_object(gene = feature_objects$gene, transcript = feature_objects$transcript, proj_dir = proj_dir)
#'
#' save_object(gene = feature_objects$gene, transcript = feature_objects$transcript, prefix = "remove_nonPRs", proj_dir = proj_dir)
#' }
save_object <- function(object, prefix = "unfiltered", proj_dir = getwd()) {
    object_dir <- fs::path(proj_dir, "output", "seurat")

    fs::dir_create(object_dir)

    object_path <- fs::path(object_dir, paste0(prefix, "_object.rds"))

    # if (interactive()) {
    #   message(paste0("Do you want to save to ", fs::path_file(object_path)))
    #   confirm_save <- (menu(c("Yes", "No")) == 1)
    # } else {
    #   confirm_save <- TRUE
    # }
    #
    # if (!confirm_save){
    #   stop("aborting project save")
    # }

    message(paste0("saving to ", object_path))
    saveRDS(object, object_path)
    # if(prefix == "unfiltered"){
    #   Sys.chmod(object_path, "775")
    # }

    return(object)
}
