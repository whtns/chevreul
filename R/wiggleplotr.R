#' Create a database of bigwigfiles
#'
#' Create a sqlite database of bigwig files matching cell ids in objects
#'
#' @param bam_files vector of paths to bam files
#' @param bigwig_db bigwig database
#'
#' @return a path to a bigwig file sqlite database
#' @export
build_bigwig_db <- function(bam_files, 
                            bigwig_db = "~/.cache/chevreul/bw-files.db") {
    bam_files <- normalizePath(bam_files)

    bigwigfiles <- map_chr(bam_files,
                           ~ bam_to_bigwig(.x,
                                           prefix = path_ext_remove(.x),
                                           overwrite = TRUE)) %>%
        set_names(path_file) %>%
        enframe("name", "bigWig") %>%
        mutate(sample_id = 
                   str_remove(name, "_Aligned.sortedByCoord.out.bw")) %>%
        identity()

    con <- dbConnect(SQLite(), dbname = bigwig_db)

    dbWriteTable(con, "bigwigfiles", bigwigfiles, append = TRUE)

    dbDisconnect(con)
}

#' Load Bigwigs
#'
#' Load a tibble of bigwig file paths by cell id
#'
#' @param object A object
#' @param bigwig_db Sqlite database of bigwig files
#'
#' @return a vector of bigwigs file paths
load_bigwigs <- function(object, bigwig_db = "~/.cache/chevreul/bw-files.db") {
    con <- dbConnect(SQLite(), dbname = bigwig_db)

    bigwigfiles <- dbReadTable(con, "bigwigfiles") %>%
        filter(sample_id %in% colnames(object)) %>%
        identity()

    missing_bigwigs <- colnames(object)[!(colnames(object) %in% 
                                              bigwigfiles$sample_id)] %>%
        paste(collapse = ", ")

    warning(paste0("Sample coverage files ", missing_bigwigs, 
                   "(.bw) do not match samples in object (check file names)"))

    dbDisconnect(con)

    bigwigfiles <-
        bigwigfiles %>%
        filter(sample_id %in% colnames(object))

    return(bigwigfiles)
}

#' Plot BigWig Coverage for Genes of Interest by a Given Variable
#'
#' Plot BigWig coverage for genes of interest colored by a given variable
#'
#' @param genes_of_interest Gene of interest
#' @param cell_metadata a dataframe with cell metadata from object
#' @param bigwig_tbl a tibble with colnames "name", "bigWig", and "sample_id"
#' matching the filename, absolute path, and sample name of each cell in the
#' cell_metadata
#' @param group_by Variable to color by
#' @param values_of_interest values of interest
#' @param organism Organism
#' @param edb ensembldb object
#' @param heights The heights of each row in the grid of plot
#' @param scale_y whether to scale coverage
#' @param reverse_x whether to reverse x axis
#' @param start start coordinates
#' @param end end coordinates
#' @param summarize_transcripts whether to summarize transcript counts
#' @param ... extra arguments passed to plotCoverageFromEnsembldb
#'
#' @return a ggplot with coverage faceted by group_by
plot_gene_coverage_by_var <- function(
        genes_of_interest = "NRL",
        cell_metadata,
        bigwig_tbl,
        group_by = "batch",
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
        edb <- EnsDb.Mmusculus.v79
    } else {
        edb <- EnsDb.Hsapiens.v86
    }

    cell_metadata["sample_id"] <- NULL


    new_track_data <-
        cell_metadata %>%
        rownames_to_column("sample_id") %>%
        select(sample_id,
            condition = {{ group_by }},
            track_id = {{ group_by }},
            colour_group = {{ group_by }},
            everything()
        ) %>%
        mutate(scaling_factor = 1) %>% # rescale(nCount_RNA)
        mutate(condition = as.factor(condition),
               colour_group = as.factor(colour_group)) %>%
        left_join(bigwig_tbl, by = "sample_id") %>%
        filter(!is.na(bigWig)) %>%
        identity()

    if (!is.null(values_of_interest)) {
        new_track_data <-
            new_track_data %>%
            filter(condition %in% values_of_interest)
    }

    if (is.null(start) | is.null(end) | list(...)$rescale_introns == TRUE) {
        region_coords <- NULL
    } else {
        region_coords <- c(start, end)
    }

    coverage_plot_list <- plotCoverageFromEnsembldb(
        ensembldb = edb,
        gene_names = genes_of_interest,
        track_data = new_track_data,
        heights = heights,
        alpha = 0.5,
        fill_palette = hue_pal()(length(levels(new_track_data$colour_group))),
        return_subplots_list = TRUE,
        region_coords = region_coords,
        ...
    )

    if (scale_y == "log10") {
        coverage_plot_list$coverage_plot <-
            coverage_plot_list$coverage_plot +
            scale_y_continuous(trans = pseudo_log_trans(base = 10), 
                               breaks = 10^(0:4)) +
            NULL
    }

    if (reverse_x) {
        transformed_x_lim <- ggplot_build(
            coverage_plot_list$coverage_plot)$layout$panel_params[[1]]$x.range

        coverage_plot_list$coverage_plot <-
            coverage_plot_list$coverage_plot +
            coord_cartesian(
                xlim = rev(transformed_x_lim),
                expand = FALSE
            ) +
            NULL

        transformed_x_lim <- ggplot_build(
            coverage_plot_list$tx_structure)$layout$panel_params[[1]]$x.range
        coverage_plot_list$tx_structure <-
            coverage_plot_list$tx_structure +
            scale_x_reverse() +
            coord_cartesian(
                xlim = rev(transformed_x_lim),
                expand = FALSE
            ) +
            NULL
    }

    x_lim <- ggplot_build(
        coverage_plot_list$coverage_plot)$layout$panel_params[[1]]$x.range

    base_coverage <- coverage_plot_list$coverage_plot$data %>%
        filter(!is.na(coverage)) %>%
        group_by(sample_id, colour_group) %>%
        summarize(sum = signif(sum(coverage), 4)) %>%
        mutate(sum = sum / diff(x_lim)) %>%
        identity()

    coverage_plot <- 
        wrap_plots(coverage_plot_list, heights = heights, ncol = 1)

    return(list(plot = coverage_plot, table = base_coverage))
}
