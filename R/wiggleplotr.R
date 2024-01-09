#' Create a database of bigwigfiles
#'
#' Create a sqlite database of bigwig files matching cell ids in objects
#'
#' @param bam_files
#' @param bigwig_db
#'
#' @return
#' @export
#'
#' @examples
build_bigwig_db <- function(bam_files, bigwig_db = "~/.cache/chevreul/bw-files.db"){

  bam_files <- fs::path_abs(bam_files)

  bigwigfiles <- purrr::map_chr(bam_files, ~megadepth::bam_to_bigwig(.x, prefix = fs::path_ext_remove(.x), overwrite = TRUE)) %>%
    purrr::set_names(fs::path_file) %>%
    tibble::enframe("name", "bigWig") %>%
    dplyr::mutate(sample_id = stringr::str_remove(name, "_Aligned.sortedByCoord.out.bw")) %>%
    identity()

  con <- DBI::dbConnect(RSQLite::SQLite(), dbname = bigwig_db)

  DBI::dbWriteTable(con, "bigwigfiles", bigwigfiles, append = TRUE)

  DBI::dbDisconnect(con)

}

#' Load Bigwigs
#'
#' Load a tibble of bigwig file paths by cell id
#'
#' @param object A object
#' @param bigwig_db Sqlite database of bigwig files
#'
#' @return
#' @export
#'
#' @examples
load_bigwigs <- function(object, bigwig_db = "~/.cache/chevreul/bw-files.db") {
  con <- DBI::dbConnect(RSQLite::SQLite(), dbname = bigwig_db)

  bigwigfiles <- DBI::dbReadTable(con, "bigwigfiles") %>%
    dplyr::filter(sample_id %in% colnames(object)) %>%
    identity()

  missing_bigwigs <- colnames(object)[!(colnames(object) %in% bigwigfiles$sample_id)] %>%
    paste(collapse = ", ")

  warning(paste0("Sample coverage files ", missing_bigwigs, "(.bw) do not match samples in object (check file names)"))

  DBI::dbDisconnect(con)

  bigwigfiles <-
    bigwigfiles %>%
    dplyr::filter(sample_id %in% colnames(object))

  return(bigwigfiles)
}

#' Plot BigWig Coverage for Genes of Interest by a Given Variable
#'
#' Plot BigWig coverage for genes of interest colored by a given variable
#'
#' @param genes_of_interest Gene of interest
#' @param cell_metadata a dataframe with cell metadata from object
#' @param bigwig_tbl a tibble with colnames "name", "bigWig", and "sample_id" matching the filename, absolute path, and sample name of each cell in the cell_metadata
#' @param var_of_interest Variable to color by
#' @param values_of_interest
#' @param organism Organism
#' @param edb ensembldb object
#' @param heights The heights of each row in the grid of plot
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
plot_gene_coverage_by_var <- function(genes_of_interest = "RXRG",
                                      cell_metadata,
                                      bigwig_tbl,
                                      var_of_interest = "batch",
                                      values_of_interest = NULL,
                                      organism = "human",
                                      edb = NULL,
                                      heights = c(3, 1),
                                      scale_y = "log10",
                                      reverse_x = FALSE,
                                      start = NULL,
                                      end = NULL,
                                      summarize_transcripts = FALSE,
                                      ...) {
  if (organism == "mouse") {
    edb <- EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79
  } else {
    edb <- EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86
  }

  cell_metadata["sample_id"] <- NULL


  new_track_data <-
    cell_metadata %>%
    tibble::rownames_to_column("sample_id") %>%
    dplyr::select(sample_id,
      condition = {{ var_of_interest }},
      track_id = {{ var_of_interest }},
      colour_group = {{ var_of_interest }},
      everything()
    ) %>%
    dplyr::mutate(scaling_factor = 1) %>% # scales::rescale(nCount_RNA)
    dplyr::mutate(condition = as.factor(condition), colour_group = as.factor(colour_group)) %>%
    dplyr::left_join(bigwig_tbl, by = "sample_id") %>%
    dplyr::filter(!is.na(bigWig)) %>%
    identity()

  if (!is.null(values_of_interest)) {
    new_track_data <-
      new_track_data %>%
      dplyr::filter(condition %in% values_of_interest)
  }

  if (is.null(start) | is.null(end) | list(...)$rescale_introns == TRUE) {
    region_coords <- NULL
  } else {
    region_coords <- c(start, end)
  }

  coverage_plot_list <- wiggleplotr::plotCoverageFromEnsembldb(
    ensembldb = edb,
    gene_names = genes_of_interest,
    track_data = new_track_data,
    heights = heights,
    alpha = 0.5,
    fill_palette = scales::hue_pal()(length(levels(new_track_data$colour_group))),
    return_subplots_list = TRUE,
    region_coords = region_coords,
    ...
  )

  if (scale_y == "log10") {
    coverage_plot_list$coverage_plot <-
      coverage_plot_list$coverage_plot +
      scale_y_continuous(trans = scales::pobjectdo_log_trans(base = 10), breaks = 10^(0:4)) +
      NULL
  }

  if (reverse_x) {
    transformed_x_lim <- ggplot2::ggplot_build(coverage_plot_list$coverage_plot)$layout$panel_params[[1]]$x.range
    # transformed_x_lim[1] <- transformed_x_lim[1] - diff(transformed_x_lim)*0.2

    coverage_plot_list$coverage_plot <-
      coverage_plot_list$coverage_plot +
      coord_cartesian(
        xlim = rev(transformed_x_lim),
        expand = FALSE
      ) +
      NULL

    transformed_x_lim <- ggplot2::ggplot_build(coverage_plot_list$tx_structure)$layout$panel_params[[1]]$x.range
    # transformed_x_lim[1] <- transformed_x_lim[1] - diff(transformed_x_lim)*0.2

    coverage_plot_list$tx_structure <-
      coverage_plot_list$tx_structure +
      scale_x_reverse() +
      coord_cartesian(
        xlim = rev(transformed_x_lim),
        expand = FALSE
      ) +
      NULL
  }

  x_lim <- ggplot2::ggplot_build(coverage_plot_list$coverage_plot)$layout$panel_params[[1]]$x.range

  base_coverage <- coverage_plot_list$coverage_plot$data %>%
    dplyr::filter(!is.na(coverage)) %>%
    # tidyr::drop_na() %>%
    dplyr::group_by(sample_id, colour_group) %>%
    dplyr::summarize(sum = signif(sum(coverage), 4)) %>%
    dplyr::mutate(sum = sum / diff(x_lim)) %>%
    identity()

  coverage_plot <- patchwork::wrap_plots(coverage_plot_list, heights = heights, ncol = 1)

  return(list(plot = coverage_plot, table = base_coverage))
}
