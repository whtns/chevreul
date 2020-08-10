

#' Set Column Names from `tximport`
#'
#' @param txi
#' @param colnames
#'
#' @return
#' @export
#'
#' @examples
set_colnames_txi <- function(txi, colnames){
  colnames(txi$counts) <- colnames
  colnames(txi$abundance) <- colnames
  colnames(txi$length) <- colnames
  return(txi)
}

#' Load Sample Counts then run Tximport
#'
#' @param countsFromAbundance
#' @param proj_dir
#' @param type
#' @param edb
#'
#' @return
#' @export
#'
#' @examples
load_counts_by_tximport <- function(proj_dir, type = "salmon", countsFromAbundance = "scaledTPM", edb = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86){

  sample_glob <- switch(type,
                        salmon = "*quant.sf",
                        stringtie = "*t_data.ctab")

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
  # for (i in seq_along(txi_transcripts)){
  #   if (class(txi_transcripts[[i]]) == "matrix"){
  #     rownames(txi_transcripts[[i]]) <- stringr::str_remove(rownames(txi_transcripts[[i]]), "\\.[0-9]$")
  #   }
  # }

  txi_genes <- tximport::summarizeToGene(txi_transcripts, tx2gene = tx2gene, ignoreTxVersion = TRUE)

  txi_transcripts$tx2gene <- tx2gene

  sample_names <- fs::path_file(fs::path_dir(sample_paths))

  txi_features <- purrr::map(list(gene = txi_genes, transcript = txi_transcripts), ~set_colnames_txi(.x, sample_names))


}


#' Load Project Sample Data
#'
#' @param meta_path
#'
#' @return
#' @export
#'
#' @examples
load_meta <- function(proj_dir){
  # load metadata
  meta_file <- gsub("_proj", "_metadata.csv", path_file(proj_dir))
  meta_file <- fs::path(proj_dir, "data", meta_file)

  tpm_meta <- read_csv(meta_file)
}


## ------------------------------------------------------------------------


#' Create a Seurat Object from return of tximport
#'
#' @param txi output from load_counts_from_stringtie
#' @param meta_tbl a tibble of cell metadata with cell ids as the first column
#' @param feature gene or transcript
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
seu_from_tximport <- function(txi, feature, meta_tbl, ...){
  # seu <- seu_from_tibbles(exp_tbl, feature, meta_tbl)

  exp_tbl <- as.matrix(txi$counts)

  expid <- gsub("-.*", "", colnames(exp_tbl))

  featuredata <- data.frame(feature = rownames(exp_tbl))
  rownames(featuredata) <- featuredata[,1]
  if (feature == "transcript"){
    featuredata <-
      featuredata %>%
      tibble::rownames_to_column("t_name") %>%
      dplyr::left_join(txi$tx2gene, by = "t_name") %>%
      tibble::column_to_rownames("t_name")
  }

  meta_tbl <- data.frame(meta_tbl)
  rownames(meta_tbl) <- meta_tbl[,"sample_id"]

  meta_tbl <- meta_tbl[colnames(exp_tbl),]

  seu <- Seurat::CreateSeuratObject(counts = exp_tbl, project = expid, assay = "RNA", meta.data = meta_tbl)

  seu@assays$RNA <- AddMetaData(seu@assays$RNA, featuredata)

  # add default batch if missing
  seu$batch <- seu@project.name

  return(seu)

}

#' create a multimodal seurat object from genes and transcript counts from tximport
#'
#' @param txi
#' @param meta_tbl
#' @param ...
#'
#' @return
#'
#'
#' @examples
multimodal_seu_from_tximport <- function (txi, meta_tbl, ...){
  gene_tbl <- as.matrix(txi$gene$counts)
  expid <- gsub("-.*", "", colnames(gene_tbl))
  genedata <- data.frame(feature = rownames(gene_tbl))
  rownames(genedata) <- genedata[, 1]
  meta_tbl <- data.frame(meta_tbl)
  rownames(meta_tbl) <- meta_tbl[, "sample_id"]
  meta_tbl <- meta_tbl[colnames(gene_tbl), ]
  seu <- Seurat::CreateSeuratObject(counts = gene_tbl, project = expid,
                                    assay = "gene", meta.data = meta_tbl)
  seu@assays$gene <- AddMetaData(seu@assays$gene, genedata)

  transcript_tbl <- as.matrix(txi$transcript$counts)
  transcriptdata <- data.frame(feature = rownames(transcript_tbl))
  rownames(transcriptdata) <- transcriptdata[, 1]

  seu[["transcript"]] <- CreateAssayObject(counts = transcript_tbl)
  seu@assays$transcript <- AddMetaData(seu@assays$transcript, transcriptdata)

  seu$batch <- seu@project.name
  return(seu)
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
seu_from_tibbles <- function(exp_tbl, feature, meta_tbl, ...){

  expid <- gsub("-.*", "", colnames(exp_tbl))

  featuredata <- data.frame(rownames(exp_tbl))
  rownames(featuredata) <- featuredata[,1]
  if (feature == "transcript"){
    # gene_id <- tx2gene
    # featuredata$gene_symbol =
  }

  meta_tbl <- data.frame(meta_tbl)
  rownames(meta_tbl) <- meta_tbl[,"sample_id"]

  meta_tbl <- meta_tbl[colnames(exp_tbl),]

  seu <- Seurat::CreateSeuratObject(counts = exp_tbl, project = expid, assay = "RNA", meta.data = meta_tbl)

  # add default batch if missing
  seu$batch <- seu@project.name

  return(seu)
}


## ------------------------------------------------------------------------
# filter out low read count cells (threshold 1e5)

#' Filter our Cells from Seurat below read count threshold
#'
#' @param seu A seurat object
#' @param read_thresh
#'
#' @return
#' @export
#'
#' @examples
filter_low_rc_cells <- function(seu, read_thresh = 1e5){

	counts <- as.matrix(seu@assays$RNA@counts)

	counts <- colSums(counts)

	keep_cells <- counts[counts > read_thresh]

	removed_cells <- counts[counts <= read_thresh]
	print(removed_cells)

	seu <- subset(seu, cells = names(keep_cells))
}

#' Save seurat object to <project>/output/sce/<feature>_seu.rds
#'
#' @param ... named arguments specifying seurat objects list of seurat objects; default "gene" and "transcript"
#' @param prefix
#' @param proj_dir
#'
#' @return
#' @export
#'
#' @examples
#' \dontrun{
#' save_seurat(gene = feature_seus$gene, transcript = feature_seus$transcript, proj_dir = proj_dir)
#'
#' save_seurat(gene = feature_seus$gene, transcript = feature_seus$transcript, prefix = "remove_nonPRs", proj_dir = proj_dir)
#'
#' }
save_seurat <- function(..., prefix = "unfiltered", proj_dir = getwd()){

  seu_list = list(...)
  if(is.null(names(seu_list))){
    seu_list = purrr::flatten(seu_list)
  }

  seurat_dir <- fs::path(proj_dir, "output", "seurat")

  fs::dir_create(seurat_dir)

  seu_path <- fs::path(seurat_dir, paste0(prefix, "_seu.rds"))


  if (interactive()) {
    message(paste0("Do you want to save to ", fs::path_file(seu_path)))
    confirm_save <- (menu(c("Yes", "No")) == 1)
  } else {
    confirm_save <- TRUE
  }

  if (!confirm_save){
    stop("aborting project save")
  }

  message(paste0("saving to ", fs::path_file(seu_path)))
  saveRDS(seu_list, seu_path)
  if(prefix == "unfiltered"){
    Sys.chmod(seu_path, "775")
  }

  return(seu_list)

}

