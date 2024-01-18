#' Gene enrichment using Enrichr.
#' @title Gene enrichment using Enrichr.
#' @description Gene enrichment using Enrichr.
#' @param genes Gene names or dataframe of gene names in first column and a
#' score between 0 and 1 in the other.
#' @param databases Databases to search.
#' @param URL_API URL to send requests to (Enrichr API).
#' See http://amp.pharm.mssm.edu/Enrichr/ for available databases.
#' @export
#' @return Returns a data frame of enrichment terms, p-values, ...
#' @author Wajid Jawaid, modified by Roman Hillje
.send_enrichr_query <- function(genes,
    databases = NULL,
    URL_API = NULL) {
    if (
        is.vector(genes) &
            !all(genes == "") &
            length(genes) != 0
    ) {
        temp <- httr::POST(
            url = URL_API,
            body = list(list = paste(genes, collapse = "\n"))
        )
    } else if (is.data.frame(genes)) {
        temp <- httr::POST(
            url = URL_API,
            body = list(list = paste(paste(genes[, 1], genes[, 2], sep = ","), collapse = "\n"))
        )
    } else {
        warning(
            paste0(
                "genes must be a non-empty vector of gene names or a dataframe ",
                "with genes and score."
            )
        )
    }
    httr::GET(url = "http://amp.pharm.mssm.edu/Enrichr/share")
    dfSAF <- options()$stringsAsFactors
    options(stringsAsFactors = FALSE)
    result <- future.apply::future_sapply(databases,
        USE.NAMES = TRUE,
        simplify = FALSE, function(x) {
            r <- httr::GET(
                url = "http://amp.pharm.mssm.edu/Enrichr/export",
                query = list(file = "API", backgroundType = x)
            )
            r <- gsub("&#39;", "'", intToUtf8(r$content))
            tc <- textConnection(r)
            r <- utils::read.table(tc,
                sep = "\t", header = TRUE, quote = "",
                comment.char = ""
            )
            close(tc)
            return(r)
        }
    )
    return(result)
}


#' Get enriched pathways based on marker genes from EnrichR.
#' @title Get enriched pathways based on marker genes from EnrichR.
#' @description This function uses the enrichR API to look for enriched pathways
#' in marker gene sets of samples and clusters.
#' @keywords Cerebro scRNAseq Seurat Enrichr
#' @param object Seurat object.
#' about sample; defaults to 'sample'.
#' @param column_cluster Column in object@meta.data that contains information
#' about cluster; defaults to 'cluster'.
#' @param databases Which databases to query. Use enrichR::listEnrichrDbs() to
#' check what databases are available.
#' @param adj_p_cutoff Cut-off for adjusted p-value of enriched pathways;
#' defaults to 0.05,
#' @param max_terms Save only first n entries of each database; defaults to 100.
#' @param URL_API URL to send requests to (Enrichr API). Allows to overwrite
#' default URL with an alternative taken from the Enrichr website in case the
#' original is out-of-service; defaults to
#' 'http://amp.pharm.mssm.edu/Enrichr/enrich'.
#' @export
#' @return Seurat object with Enrichr results for samples and clusters stored in
#' object@misc$enriched_pathways$enrichr
#' @import dplyr
#' @importFrom rlang .data
#' @author Roman Hillje, modified by Kevin Stachelek
#' @examples
#' pbmc <- readRDS(system.file("extdata/v1.2/seurat_pbmc.rds",
#'     package = "cerebroApp"
#' ))
#' pbmc <- getEnrichedPathways(
#'     object = pbmc,
#'     column_sample = "sample",
#'     column_cluster = "seurat_clusters",
#'     databases = c("GO_Biological_Process_2018", "GO_Cellular_Component_2018"),
#'     adj_p_cutoff = 0.01,
#'     max_terms = 100,
#'     URL_API = "http://amp.pharm.mssm.edu/Enrichr/enrich"
#' )
getEnrichedPathways <- function(object,
    column_cluster = "group",
    databases = c(
        "GO_Biological_Process_2018",
        "GO_Cellular_Component_2018",
        "GO_Molecular_Function_2018",
        "KEGG_2016",
        "WikiPathways_2016",
        "Reactome_2016",
        "Panther_2016",
        "Human_Gene_Atlas",
        "Mouse_Gene_Atlas"
    ),
    adj_p_cutoff = 0.05,
    max_terms = 100,
    URL_API = "http://amp.pharm.mssm.edu/Enrichr/enrich") {
    ## check if Seurat is installed
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        stop(
            "Package 'Seurat' needed for this function to work. Please install it.",
            call. = FALSE
        )
    }
    ## --------------------------------------------------------------------------##
    ## check if marker genes are present and stop if they aren't
    ## --------------------------------------------------------------------------##
    if (is.null(object@misc$markers)) {
        stop(
            "No marker genes found. Please run 'getMarkerGenes()' first.",
            call. = FALSE
        )
    }

    ## --------------------------------------------------------------------------##
    ## clusters
    ## - check if marker genes by cluster are available
    ## - extract marker genes by cluster
    ## - get cluster names and remove those for which no marker genes are available
    ## - create slot for annotation if doesn't already exist
    ## - annotate marker genes for each cluster in parallel
    ## - try up to three times to run enrichR annotation (fails sometimes)
    ## - filter results
    ## --------------------------------------------------------------------------##
    if (!is.null(object@misc$markers[[1]]$presto)) {
        if (is.data.frame(object@misc$markers[[1]]$presto)) {
            message(
                paste0(
                    "[", format(Sys.time(), "%H:%M:%S"),
                    "] Get enriched pathway for clusters..."
                )
            )

            ## remove clusters for which no marker genes were found
            markers_by_cluster <- object@misc$markers[[1]]$presto %>%
                # dplyr::filter(padj < 0.05) %>%
                identity()

            cluster_names <- names(markers_by_cluster)

            results_by_cluster <- future.apply::future_sapply(
                cluster_names,
                USE.NAMES = TRUE, simplify = FALSE,
                future.globals = FALSE, function(x) {
                    temp <- list()
                    attempt <- 1
                    while (
                        "Adjusted.P.value" %in% names(temp) == FALSE &&
                            attempt <= 3
                    ) {
                        attempt <- attempt + 1
                        try(
                            temp <- markers_by_cluster %>%
                                dplyr::pull(x) %>%
                                .send_enrichr_query(databases = databases, URL_API = URL_API)
                        )
                    }
                    #
                    results_2 <- sapply(names(temp),
                        USE.NAMES = TRUE,
                        simplify = FALSE, function(y) {
                            ## apply cut-off of adj. p-value and add database info as column
                            out <- temp[[y]] %>%
                                dplyr::filter(.data$Adjusted.P.value <= adj_p_cutoff) %>%
                                dplyr::mutate(db = y)
                            ## if there are more than max_terms entries...
                            if (nrow(out) > max_terms) {
                                out <- out %>% dplyr::top_n(-max_terms, .data$Adjusted.P.value)
                                ## if there are no entries left
                            } else if (nrow(out) == 0) {
                                out <- NULL
                            }
                            return(out)
                        }
                    )
                    ## remove dbs without any enriched entries
                    for (i in names(results_2))
                    {
                        if (is.null(results_2[[i]])) results_2[[i]] <- NULL
                    }
                    ## merge databases within each cluster
                    results_2 <- do.call(rbind, results_2)
                    return(results_2)
                }
            )

            ## remove clusters without any enriched entry in any database
            for (i in names(results_by_cluster))
            {
                if (is.null(results_by_cluster[[i]])) results_by_cluster[[i]] <- NULL
            }
            ## add cluster info as column
            for (i in names(results_by_cluster))
            {
                results_by_cluster[[i]] <- results_by_cluster[[i]] %>%
                    dplyr::mutate(group = i)
            }
            ## merge clusters into single table
            results_by_cluster <- do.call(rbind, results_by_cluster) %>%
                dplyr::select(.data$group, .data$db, dplyr::everything()) %>%
                dplyr::mutate(
                    cluster = factor(.data$group, levels = intersect(
                        cluster_names,
                        .data$group
                    )),
                    db = factor(.data$db, databases)
                )
            message(
                paste0(
                    "[", format(Sys.time(), "%H:%M:%S"), "] ", nrow(results_by_cluster),
                    " pathways passed the thresholds across all clusters and databases."
                )
            )
        } else if (object@misc$markers == "no_markers_found") {
            message(
                paste0(
                    "[", format(Sys.time(), "%H:%M:%S"),
                    "] Skipping pathway enrichment for cluster because no marker genes ",
                    "were identified for any cluster."
                )
            )
            results_by_cluster <- "no_markers_found"
        } else {
            warning(
                paste0(
                    "Unexpected data format of marker genes for clusters. Please submit ",
                    "an issue on GitHub: https://github.com/romanhaa/cerebroApp."
                )
            )
        }
    } else {
        message(
            paste0(
                "[", format(Sys.time(), "%H:%M:%S"),
                "] No marker genes for clusters available."
            )
        )
    }
    ## ---------------------------------------------------------------------------#
    ## merge results, add to Seurat object and return Seurat object
    ## ---------------------------------------------------------------------------#
    results <- list(
        by_cluster = results_by_cluster,
        parameters = list(
            databases = databases,
            adj_p_cutoff = adj_p_cutoff,
            max_terms = max_terms
        )
    )
    object@misc$enriched_pathways$enrichr <- results

    ## --------------------------------------------------------------------------##
    ## return Seurat object
    ## --------------------------------------------------------------------------##
    return(object)
}


#' Title
#'
#' @param enrich_by_cluster
#'
#' @return
#' @export
#'
#' @examples
format_pathway_table <- function(enrich_by_cluster, cluster, db) {
    enrich_by_cluster <-
        enrich_by_cluster %>%
        dplyr::filter(
            cluster == cluster,
            db == db
        ) %>%
        dplyr::select(3, 4, 5, 6, 10, 11) %>%
        dplyr::arrange(-Combined.Score) %>%
        dplyr::mutate(
            P.value = formatC(P.value, format = "e", digits = 2),
            Adjusted.P.value = formatC(Adjusted.P.value, format = "e", digits = 2),
            Combined.Score = formatC(Combined.Score, format = "f", digits = 2)
        ) %>%
        dplyr::rename(
            "p-value" = P.value,
            "adj. p-value" = Adjusted.P.value,
            "combined score" = Combined.Score,
        )

    enrich_by_cluster <-
        enrich_by_cluster %>%
        formattable::formattable(
            list("combined score" = formattable::color_bar("pink"))
        ) %>%
        formattable::as.datatable(
            filter = "top",
            selection = "none",
            escape = FALSE,
            autoHideNavigation = TRUE,
            rownames = FALSE,
            extensions = c("Buttons"),
            class = "cell-border stripe",
            options = list(
                columnDefs = list(list(visible = FALSE, targets = c(2, 5))),
                scrollX = TRUE,
                dom = "Bfrtip",
                lengthMenu = c(15, 30, 50, 100),
                pageLength = 15,
                buttons = list(
                    "colvis",
                    list(
                        extend = "collection",
                        text = "Download",
                        buttons = list(
                            list(
                                extend = "csv",
                                filename = "enriched_pathways_from_enrichr_by_cluster",
                                title = "Enriched pathways from Enrichr by cluster"
                            ),
                            list(
                                extend = "excel",
                                filename = "enriched_pathways_from_enrichr_by_cluster",
                                title = "Enriched pathways from Enrichr by cluster"
                            )
                        )
                    )
                )
            )
        )

    enrich_by_cluster <-
        enrich_by_cluster %>%
        DT::formatStyle(columns = c("combined score"), textAlign = "right")

    return(enrich_by_cluster)
}
