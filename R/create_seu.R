

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
#'
#' @return
#' @export
#'
#' @examples
load_counts_from_stringtie <- function(proj_dir, txOut){
  stringtie_paths <- rlang::with_handlers(
    error = ~ abort("Can't find input stringtie files (stringtie output with extension t._ctab)", parent  = .),
    stringtie_files <- fs::path(proj_dir, "output", "stringtie") %>%
      dir_ls(recursive = T) %>%
      path_filter("*t_data.ctab") %>%
      identity()
  )

  tmp <- read_tsv(stringtie_paths[1])
  tx2gene <- tmp[, c("t_name", "gene_name")]

  txi <- tximport::tximport(stringtie_files, type = "stringtie", tx2gene = tx2gene, txOut = txOut, countsFromAbundance = "scaledTPM")

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

#' Create a Seurat Object from tximport and sample data
#'
#' @param txi
#' @param colData
#' @param census_counts
#'
#' @return
#' @export
#'
#' @examples
seu_from_tibbles <- function(txi, colData, census_counts=NULL){
	# browser()

	expid <- gsub("_.*", "", colData[1,1])

  featuredata <- data.frame(rownames(txi$counts))
  rownames(featuredata) <- featuredata[,1]

  txi$counts <- as.matrix(txi$counts)


  colData <- data.frame(colData)
  rownames(colData) <- colData[,1]

  txi[1:3] <- purrr::map(txi[1:3], ~ .x[,(colnames(.x) %in% rownames(colData))])

  colData <- colData[colnames(txi$counts),]

  seu <- CreateSeuratObject(counts = txi$counts, project = expid, assay = "RNA", meta.data = colData)

  # add default batch if missing
  seu$batch <- seu@project.name

  return(seu)
}


## ------------------------------------------------------------------------
# filter out low read count cells (threshold 1e5)

#' Filter our Cells from Seurat below read count threshold
#'
#' @param seu
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

	seu <- subset(seu, cells = names(keep_cells))
}


#' Save seurat object to <project>/output/sce/<feature>_seu.rds
#'
#' @param seu
#' @param feature
#'
#' @return
#' @export
#'
#' @examples
save_seurat <- function(seu, feature){
  sce_dir <- fs::path(proj_dir, "output", "sce")

  dir_create(sce_dir)

  seu_path <- fs::path(sce_dir, paste0(feature, "_seu.rds"))

  saveRDS(feature_seu, seu_path)
}


