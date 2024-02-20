#' Set Column Names from `tximport`
#'
#' @param txi txi
#' @param colnames column names
#'
#' @return a list of matrices
#' @export
set_colnames_txi <- function(txi, colnames) {
    colnames(txi$counts) <- colnames
    colnames(txi$abundance) <- colnames
    colnames(txi$length) <- colnames
    return(txi)
}

#' Load Counts by Tximport
#'
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
#' @param edb ensembldb reference
#'
#' @return a list of feature count matrices for gene symbols and transcripts
#' @export
#'
#' @importFrom ensembldb transcripts
#' @importFrom EnsDb.Hsapiens.v86 EnsDb.Hsapiens.v86
#' @importFrom rlang with_handlers abort
#' @importFrom tximport tximport summarizeToGene
load_counts_by_tximport <- function(proj_dir, type = "salmon", countsFromAbundance = "scaledTPM", edb = EnsDb.Hsapiens.v86) {
    sample_glob <- switch(type,
        kallisto = "*abundance.h5",
        salmon = "*quant.sf",
        stringtie = "*t_data.ctab"
    )

    sample_paths <- with_handlers(error = ~ abort("Can't find input files",
        parent = .
    ), sample_files <- path(
        proj_dir, "output",
        type
    ) %>% dir_ls(recurse = T, glob = sample_glob) %>%
        identity())

    tx2gene <- transcripts(edb, return.type = "data.frame")[, c("tx_id", "gene_id")] %>%
        left_join(annotables::grch38, by = c("gene_id" = "ensgene")) %>%
        select(tx_id, symbol) %>%
        drop_na()

    txi_transcripts <- tximport(sample_files, type = type, tx2gene = tx2gene, txOut = T, countsFromAbundance = countsFromAbundance, ignoreTxVersion = TRUE)

    message("sanitize transcript ids with trailing (.1, .2, etc)")
    txi_transcripts <- map_if(
        txi_transcripts, is.matrix,
        ~ `rownames<-`(.x, str_remove(rownames(.x), "\\.[0-9]$"))
    )

    txi_genes <- summarizeToGene(txi_transcripts, tx2gene = tx2gene, ignoreTxVersion = TRUE)

    txi_transcripts$tx2gene <- tx2gene

    sample_names <- path_file(path_dir(sample_paths))

    txi_features <- map(list(gene = txi_genes, transcript = txi_transcripts), ~ set_colnames_txi(.x, sample_names))
}


#' Load Sample Metadata for a given project
#'
#'
#'
#' @param proj_dir a project directory
#'
#' @return a cell level metadata tibble
#' @export
#'
load_meta <- function(proj_dir) {
    meta_file <- gsub("_proj", "_metadata.csv", path_file(proj_dir))
    meta_file <- path(proj_dir, "data", meta_file)

    tpm_meta <- read_csv(meta_file)
}


#' Create a object from output of  \href{http://bioconductor.org/packages/release/bioc/html/tximport.html}{tximport} and a table of cell metadata
#'
#' @param txi output from load_counts_by_tximport
#' @param meta_tbl a tibble of cell metadata with cell ids as the first column
#'
#' @return a single cell object
#' @export
#'
object_from_tximport <- function(txi, meta_tbl) {
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
#' @param exp_tbl experiment tibble
#' @param feature feature type
#' @param meta_tbl metadata tibble
#'
#' @return a single cell object
#' @export
#'
object_from_tibbles <- function(exp_tbl, feature, meta_tbl) {
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

#' Filter our Cells from Seurat below read count threshold
#'
#' @param object A object
#' @param read_thresh read threshold
#'
#' @return a single cell object
#' @export
#'
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
#' @param object a single cell object
#' @param prefix a prefix for saving
#' @param proj_dir path to a project directory
#'
#' @return a path to an rds file containing a single cell object
#' @export
#'
#' @examples
#' \dontrun{
#' save_object(gene = feature_objects$gene, transcript = feature_objects$transcript, proj_dir = proj_dir)
#'
#' save_object(gene = feature_objects$gene, transcript = feature_objects$transcript, prefix = "remove_nonPRs", proj_dir = proj_dir)
#' }
save_object <- function(object, prefix = "unfiltered", proj_dir = getwd()) {
    object_dir <- path(proj_dir, "output", "seurat")

    dir.create(object_dir)

    object_path <- path(object_dir, paste0(prefix, "_object.rds"))

    message(paste0("saving to ", object_path))
    saveRDS(object, object_path)

    return(object)
}
