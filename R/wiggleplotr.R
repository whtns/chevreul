#' Prep Cell Metadata for Wiggleplotr
#'
#' @param cell_metadata_file
#' @param var_of_interest
#' @param bigwig_dir
#'
#' @return
#' @export
#'
#' @examples
prep_cell_metadata <- function(cell_metadata_file, var_of_interest, bigwig_dir){

  # browser()
  bigwigs <-
    bigwig_dir %>%
    fs::dir_ls(glob = "*.bw", recurse = TRUE) %>%
    purrr::set_names(stringr::str_extract(fs::path_file(.), "ds.*[0-9](?=.)")) %>%
    tibble::enframe("sample_id", "bigWig")

  metadata <- read_csv(cell_metadata_file)
  metadata <-
    metadata %>%
    dplyr::select(sample_id,
                  condition = {{var_of_interest}},
                  track_id = {{var_of_interest}},
                  colour_group = {{var_of_interest}},
                  everything()) %>%
    dplyr::mutate(scaling_factor = 1, condition = as.factor(condition), colour_group = as.factor(colour_group)) %>%
    dplyr::left_join(bigwigs, by = "sample_id") %>%
    identity()

}

#' Plot BigWig Coverage for Genes of Interest Colored by a Given Variable
#'
#' @param genes_of_interest
#' @param cell_metadata_file
#' @param var_of_interest
#' @param edb
#'
#' @return
#' @export
#'
#' @examples
plot_gene_coverage_by_var <- function(genes_of_interest = "RXRG",
                                      cell_metadata_file, var_of_interest,
                                      edb = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86) {

  myannofilterlist <- AnnotationFilterList(GeneNameFilter(genes_of_interest))

  exonstoplot <- exonsBy(edb, by = "tx", filter = myannofilterlist)
  cdstoplot <- cdsBy(edb, by = "tx", filter = myannofilterlist)

  transcript_annotations <-
    annotables::grch38 %>%
    dplyr::filter(symbol %in% genes_of_interest) %>%
    dplyr::left_join(annotables::grch38_tx2gene, by = "ensgene") %>%
    dplyr::select(transcript_id = enstxp, gene_name = "symbol", strand) %>%
    dplyr::distinct()

  bigwig_dir = cell_metadata_file %>%
    fs::path_dir() %>%
    fs::path_dir() %>%
    fs::path("output", "HISAT2bigwig")

  new_track_data <- prep_cell_metadata(cell_metadata_file, var_of_interest = var_of_interest, bigwig_dir = bigwig_dir)

  coverage_plot <- wiggleplotr::plotCoverage(exonstoplot, cdstoplot,
                                             transcript_annotations = transcript_annotations,
                                             track_data = new_track_data, heights = c(2, 1),
                                             mean_only = FALSE, alpha = 0.5,
                                             fill_palette = scales::hue_pal()(length(levels(new_track_data$colour_group)))
  )

}

#' Retrieve Metadata from Batch
#'
#' @param batch
#' @param projects_dir
#' @param db_path
#'
#' @return
#' @export
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

