

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

#' Load Counts from Stringtie Files
#'
#' @param stringtie_paths
#' @param txOut
#' @param countsFromAbundance
#'
#' @return
#' @export
#'
#' @examples
load_counts_from_stringtie <- function(proj_dir, txOut, countsFromAbundance = "scaledTPM"){
  stringtie_paths <- rlang::with_handlers(
    error = ~ rlang::abort("Can't find input stringtie files (stringtie output with extension t._ctab)", parent  = .),
    stringtie_files <- fs::path(proj_dir, "output", "stringtie") %>%
      dir_ls(recursive = T) %>%
      path_filter("*t_data.ctab") %>%
      identity()
  )

  tmp <- read_tsv(stringtie_paths[1])
  tx2gene <- tmp[, c("t_name", "gene_name")]

  txi <- tximport::tximport(stringtie_files, type = "stringtie", tx2gene = tx2gene, txOut = txOut, countsFromAbundance = countsFromAbundance)

  sample_names <- path_file(path_dir(stringtie_paths))

  txi <- set_colnames_txi(txi, sample_names)

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

#Make SCE from census counts-------------------------------------------

#' Create a Seurat Object from return of tximport
#'
#' @param txi
#' @param meta_tbl
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
seu_from_tximport <- function(txi, meta_tbl, ...){
  exp_tbl <- as.matrix(txi$counts)

  seu <- seu_from_tibbles(exp_tbl, meta_tbl)

  return(seu)

}

#' Create a Seurat Object from a set of tibbles
#'
#' @param exp_tbl
#' @param meta_tbl
#' @param census_counts
#'
#' @return
#' @export
#'
#' @examples
seu_from_tibbles <- function(exp_tbl, meta_tbl, census_counts=NULL){
  # browser()

  expid <- gsub("-.*", "", colnames(exp_tbl))

  featuredata <- data.frame(rownames(exp_tbl))
  rownames(featuredata) <- featuredata[,1]

  meta_tbl <- data.frame(meta_tbl)
  rownames(meta_tbl) <- meta_tbl[,"sample_id"]

  # exp_tbl[1:3] <- purrr::map(exp_tbl[1:3], ~ .x[,(colnames(.x) %in% rownames(meta_tbl))])

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
	# browser()
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
<<<<<<< HEAD
#' \dontrun{
#' save_seurat(gene = feature_seus$gene, transcript = feature_seus$transcript, proj_dir = proj_dir)
#'
#' save_seurat(gene = feature_seus$gene, transcript = feature_seus$transcript, prefix = "remove_nonPRs", proj_dir = proj_dir)
#'
#' }
save_seurat <- function(..., prefix = "unfiltered", proj_dir = getwd()){

  seu_list = list(...)
  if(is.null(names(seu_list))){
    seu_list = flatten(seu_list)
  }
||||||| merged common ancestors
save_seurat <- function(seu, feature, suffix = "", proj_dir = getwd(), temp = F){

  if (temp == TRUE) return(seu)

  if(suffix != ""){
    suffix = paste0("_", suffix)
  }
=======
save_seurat <- function(seu, feature, prefix = "unfiltered", proj_dir = getwd()){

  prefix <- paste0(prefix, "_")
>>>>>>> 2d43e2725a02764750212f5c9acf866ac6e6a483

  seurat_dir <- fs::path(proj_dir, "output", "seurat")

  fs::dir_create(seurat_dir)

<<<<<<< HEAD
  seu_path <- fs::path(seurat_dir, paste0(prefix, "_seu.rds"))

  message(paste0("saving seurat objects to ", seu_path))
  saveRDS(seu_list, seu_path)
||||||| merged common ancestors
  seu_path <- fs::path(seurat_dir, paste0(feature, "_seu", suffix, ".rds"))
=======
  seu_path <- fs::path(seurat_dir, paste0(prefix, "_seu.rds"))
>>>>>>> 2d43e2725a02764750212f5c9acf866ac6e6a483

<<<<<<< HEAD
  return(seu_list)
||||||| merged common ancestors
  saveRDS(seu, seu_path)
  seu
=======
  message(paste0("saving seurat object to ", seu_path))
  saveRDS(seu, seu_path)
  seu
>>>>>>> 2d43e2725a02764750212f5c9acf866ac6e6a483
}

