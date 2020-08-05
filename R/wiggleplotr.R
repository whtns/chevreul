#' Create a tibble of big file paths by sample
#'
#' @param seu
#' @param proj_dir
#'
#' @return
#' @export
#'
#' @examples
load_bigwigs <- function(seu, proj_dir){

  bigwig_dir <- fs::path(proj_dir, "output", "HISAT2bigwig")

  if(!dir.exists(bigwig_dir)) stop("Sample coverage files (.bw) do not exist in the selected project")

  bigwig_tbl <-
    bigwig_dir %>%
    fs::dir_ls(glob = "*.bw", recurse = TRUE) %>%
    purrr::set_names(stringr::str_remove(fs::path_file(.), "_Aligned.sortedByCoord.out.bw")) %>%
    tibble::enframe("sample_id", "bigWig")

  if(!all(colnames(seu) %in% bigwig_tbl$sample_id)) stop("Sample coverage files (.bw) do not match samples in seurat object (check file names)")

  return(bigwig_tbl)

}

#' Plot BigWig Coverage for Genes of Interest Colored by a Given Variable
#'
#' @param genes_of_interest
#' @param cell_metadata
#' @param bigwig_tbl
#' @param var_of_interest
#' @param values_of_interest
#' @param edb
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plot_gene_coverage_by_var <- function(genes_of_interest = "RXRG",
                                      cell_metadata,
                                      bigwig_tbl,
                                      var_of_interest = NULL,
                                      values_of_interest = NULL,
                                      edb = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86,
                                      ...) {

  cell_metadata["sample_id"] <- NULL

  new_track_data <-
    cell_metadata %>%
    tibble::rownames_to_column("sample_id") %>%
    dplyr::select(sample_id,
                  condition = {{var_of_interest}},
                  track_id = {{var_of_interest}},
                  colour_group = {{var_of_interest}},
                  everything()) %>%
    dplyr::mutate(scaling_factor = 1, condition = as.factor(condition), colour_group = as.factor(colour_group)) %>%
    dplyr::left_join(bigwig_tbl, by = "sample_id") %>%
    identity()

  if (!is.null(values_of_interest)){
    new_track_data <-
      new_track_data %>%
      dplyr::filter(condition %in% values_of_interest)
  }

  coverage_plot <- wiggleplotr::plotCoverageFromEnsembldb(ensembldb = edb,
                            gene_names = genes_of_interest,
                            track_data = new_track_data,
                            heights = c(2,1),
                            alpha = 0.5,
                            fill_palette = scales::hue_pal()(length(levels(new_track_data$colour_group))),
                            ...
                            )

  return(coverage_plot)

}

#' Retrieve Metadata from Batch
#'
#' @param batch
#' @param projects_dir
#' @param db_path
#'
#' @return
#'
#' @examples
metadata_from_batch <- function(batch, projects_dir = "/dataVolume/storage/single_cell_projects",
                                db_path = "single-cell-projects.db"){

    mydb <- DBI::dbConnect(RSQLite::SQLite(), fs::path(projects_dir, db_path))

    projects_tbl <- DBI::dbReadTable(mydb, "projects_tbl") %>%
      dplyr::filter(!project_type %in% c("integrated_projects", "resources"))

    DBI::dbDisconnect(mydb)

    metadata <-
      projects_tbl %>%
      dplyr::filter(project_slug == batch) %>%
      dplyr::pull(project_path) %>%
      fs::path("data") %>%
      fs::dir_ls(glob = "*.csv") %>%
      identity()

}

