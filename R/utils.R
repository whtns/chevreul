#' Add Seurat Object Metadata
#'
#' Adds data to the seurat object to produce a seurat object with metadata added.
#'
#' @param seu A seurat object
#' @param datapath Path to file containing metadata
#'
#' @return
#' @export
#'
#' @examples
format_new_metadata <- function(seu, datapath) {
    new_meta <- read_csv(datapath) %>%
        dplyr::mutate(across(contains("snn"), as.factor))

    rowname_col <- colnames(new_meta)[1]

    new_meta <- tibble::column_to_rownames(new_meta, rowname_col)

    seu <- Seurat::AddMetaData(seu, new_meta)
    DefaultAssay(seu) <- "gene"
    ncalc <- Seurat:::CalcN(seu)
    seu$nFeature_RNA <- ncalc$nFeature
    seu$nCount_RNA <- ncalc$nCount

    return(seu)
}


#' Reformat Seurat Object Metadata
#'
#' Reformat Seurat Object Metadata by Coalesced Columns
#'
#' @param seu A Seurat object
#' @param cols Columns
#' @param new_col New columns
#'
#' @return
#' @export
#'
#' @examples
combine_cols <- function(seu, cols, new_col) {
    new_col <- janitor::make_clean_names(new_col)

    # ensure that new colname will not be dropped
    drop_cols <- cols[!cols == new_col]

    # make sure that none of the columns to be coalesced are entirely NA
    na_cols <- purrr::map_lgl(cols, ~ all(is.na(seu[[.x]])))
    cols <- cols[!na_cols]

    # check class of cols to be coalesced


    meta <- tibble::rownames_to_column(seu[[]]) %>%
        dplyr::mutate_at(vars(one_of(cols)), as.character) %>%
        dplyr::mutate(!!new_col := dplyr::coalesce(!!!syms(cols))) %>%
        dplyr::select(-drop_cols) %>%
        tibble::column_to_rownames(var = "rowname") %>%
        identity()
}

#' Filter Rows to Top
#'
#' @param df
#' @param column
#' @param values
#'
#' @return
#' @export
#'
#' @examples
filter_rows_to_top <- function(df, column, values) {
    matched_df <- df[df[[column]] %in% values, ]

    matched_df <- matched_df[match(values, matched_df[[column]]), ]

    unmatched_df <- df[!(df[[column]] %in% values), ]

    total_df <- list(matched_df = matched_df, unmatched_df = unmatched_df)
    total_df <- dplyr::bind_rows(total_df)

    return(total_df)
}

#' Collate list of variables to be plotted
#'
#' @param seu
#'
#' @return plot_types
#' @export
#' @examples
list_plot_types <- function(seu) {
    meta_types <- tibble::tibble(
        vars = colnames(seu[[]]),
        var_type = purrr::map_chr(purrr::map(seu[[]], pillar::new_pillar_type), 1),
        num_levels = unlist(purrr::map(seu[[]], ~ length(unique(.x))))
    )

    meta_types <- meta_types %>%
        # dplyr::filter(!grepl("_snn_res", vars)) %>%
        dplyr::mutate(meta_type = dplyr::case_when(
            var_type %in% c("int", "dbl") ~ "continuous",
            var_type %in% c("chr", "fct", "ord", "lgl", "glue") ~ "category"
        )) %>%
        dplyr::mutate(meta_type = ifelse(meta_type == "continuous" & num_levels < 30, "category", meta_type)) %>%
        dplyr::filter(num_levels > 1) %>%
        identity()

    continuous_vars <- meta_types %>%
        dplyr::filter(meta_type == "continuous") %>%
        dplyr::pull(vars)

    continuous_vars <- c("feature", continuous_vars) %>%
        purrr::set_names(stringr::str_to_title(stringr::str_replace_all(., "[[:punct:]]", " ")))


    category_vars <- meta_types %>%
        dplyr::filter(meta_type == "category") %>%
        dplyr::pull(vars) %>%
        purrr::set_names(stringr::str_to_title(stringr::str_replace_all(., "[^[:alnum:][:space:]\\.]", " ")))

    plot_types <- list(category_vars = category_vars, continuous_vars = continuous_vars)



    return(plot_types)
}

#' Get Transcripts in Seurat Object
#'
#' Get transcript ids in Seurat objects for one or more gene of interest
#'
#' @param seu A seurat object
#' @param gene Gene of intrest
#' @param organism Organism
#'
#' @return
#' @export
#'
#' @examples
#' RXRG_transcripts <- get_transcripts_from_seu(human_gene_transcript_seu, "RXRG")
#'
get_transcripts_from_seu <- function(seu, gene, organism = "human") {
    transcripts <- genes_to_transcripts(gene, organism)

    transcripts <- transcripts[transcripts %in% rownames(GetAssay(seu, "transcript"))]
}

#' Title
#'
#' @param cds
#' @param mygenes
#' @param resolution
#'
#' @return
#' @export
#'
#' @examples
prep_plot_genes_in_pseudotime <- function(cds, mygenes, resolution, partition = FALSE) {
    if (partition) {
        partition_cells <- monocle3::partitions(cds)
        # partition_cells <-  split(names(partition_cells), partition_cells)[[input$partitions]]
        partition_cells <- split(names(partition_cells), partition_cells)[[1]]

        cds <- cds[, colnames(cds) %in% partition_cells]
    }

    cds <- cds[rownames(cds) %in% mygenes, ]

    if (any(grepl("integrated", colnames(colData(cds))))) {
        default_assay <- "integrated"
    } else {
        default_assay <- "gene"
    }

    color_cells_by <- paste0(default_assay, "_snn_res.", resolution)

    gene_ptime_plot <- monocle3::plot_genes_in_pseudotime(cds,
        color_cells_by = color_cells_by,
        min_expr = 0.5
    )

    return(gene_ptime_plot)
}


#' Record Experiment Metadata
#'
#' @param object A seurat objet
#' @param experiment_name
#' @param organism
#'
#' @return
#' @export
#'
#' @examples
#' logged_seu <- record_experiment_data(human_gene_transcript_seu, experiment_name = "human_gene_transcript", organism = "mouse")
#' Misc(logged_seu, "experiment")
#' @importFrom purrr %||%
record_experiment_data <- function(object, experiment_name = "default_experiment", organism = "human") {
    if (!requireNamespace("Seurat", quietly = TRUE)) {
        stop("Package 'Seurat' needed for this function to work. Please install it.",
            call. = FALSE
        )
    }

    organism <- Seurat::Misc(object, "experiment")[["organism"]] %||% organism

    experiment_name <- Seurat::Misc(object, "experiment")[["experiment_name"]] %||% experiment_name

    message(paste0("[", format(Sys.time(), "%H:%M:%S"), "] Logging Technical Details..."))
    experiment <- list(
        experiment_name = experiment_name,
        organism = organism
    )
    experiment$date_of_export <- Sys.Date()
    experiment$date_of_analysis <- Sys.Date()

    experiment$parameters <- list(
        gene_nomenclature = "gene_symbol",
        discard_genes_expressed_in_fewer_cells_than = 10,
        keep_mitochondrial_genes = TRUE,
        variables_to_regress_out = "nCount_RNA",
        number_PCs = 30,
        tSNE_perplexity = 30,
        cluster_resolution = seq(0.2, 2.0, by = 0.2)
    )
    experiment$filtering <- list(
        UMI_min = 50,
        genes_min = 10
    )
    experiment$session_info <- list(
        capture.output(sessioninfo::session_info())
    )

    if (!is.null(object@version)) {
        experiment$seurat_version <- object@version
    }

    experiment$chevreul_version <- utils::packageVersion("chevreul")

    object@misc[["experiment"]] <- NULL
    object@misc[["experiment"]] <- experiment

    return(object)
}


#' Update a chevreul Object
#'
#' @param seu_path Path to a seurat object
#' @param feature
#' @param resolution
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
update_chevreul_object <- function(seu_path, feature, resolution = seq(0.2, 2.0, by = 0.2), return_seu = TRUE, ...) {
    message(seu_path)
    seu <- readRDS(seu_path)

    if (is.list(seu)) {
        seu <- convert_seu_list_to_multimodal(seu)
        # seu <- Seurat::UpdateSeuratObject(seu)
    } else if (all(names(seu@assays) == "RNA")) {
        seu <- RenameAssays(seu, RNA = "gene")
    } else if (identical(names(seu@assays), c("RNA", "integrated"))) {
        seu <- RenameAssays(seu, RNA = "gene")
    }

    seurat_version <- seu@misc$experiment$seurat_version

    if (packageVersion("Seurat") == "5.0.0" & (seurat_version < 5 || is.null(seurat_version))) {
        seu <- convert_v3_to_v5(seu)
    }

    seu <- propagate_spreadsheet_changes(seu@meta.data, seu)

    # set appropriate assay
    if ("integrated" %in% names(seu@assays)) {
        default_assay <- "integrated"
    } else {
        default_assay <- "gene"
    }

    DefaultAssay(seu) <- default_assay

    cluster_tag <- glue::glue("{DefaultAssay(seu)}_snn_res\\.")

    cluster_names <- str_subset(names(seu@meta.data), cluster_tag)
    new_cluster_names <- str_replace(cluster_names, cluster_tag, "cluster_resolution_")

    new_cluster_cols <- seu@meta.data[cluster_names]
    names(new_cluster_cols) <- new_cluster_names

    new_meta <- cbind(seu@meta.data, new_cluster_cols)

    seu@meta.data <- new_meta

    # names(seu@meta.data) <- new_cluster_names


    chevreul_version <- seu@misc$experiment$chevreul_version

    chevreul_version <- ifelse(is.null(chevreul_version), 0.1, chevreul_version)

    # update human gene symbols to grch38
    old_symbol <- "CTC-378H22.2"
    if (old_symbol %in% rownames(seu[["gene"]])) {
        for (i in names(seu@assays)[names(seu@assays) %in% c("gene", "integrated")]) {
            # seu <- update_human_gene_symbols(seu, assay = i)
        }
    }

    if (chevreul_version < getNamespaceVersion("chevreul")) {
        message(paste0(seu_path, " is out of date! updating..."))
        if (!any(grepl("_snn_res", colnames(seu@meta.data)))) {
            seu <- seurat_cluster(seu = seu, resolution = resolution, reduction = "pca", ...)
        }

        for (i in names(seu@assays)[names(seu@assays) %in% c("gene", "integrated")]) {
            seu <- find_all_markers(seu, seurat_assay = i)
        }

        seu <- record_experiment_data(seu, ...)
        seu <- seu_calcn(seu)
    }


    if (return_seu) {
        return(seu)
    } else {
        message(paste0("saving ", seu_path))
        # saveRDS(seu, gsub(".rds", "_multimodal.rds", seu_path))
        saveRDS(seu, seu_path)
    }
}


#' Calculate Read Count Metrics for a Seurat object
#'
#' Recalculate counts/features per cell for a seurat object
#'
#' @param seu A seurat object
#' @param assay Assay to use, Default = "gene"
#' @param slot
#'
#' @return
#' @export
#'
#' @examples
seu_calcn <- function(seu, assay = "gene", slot = "counts") {
    n.calc <- Seurat:::CalcN(object = GetAssay(seu, assay))
    if (!is.null(x = n.calc)) {
        names(x = n.calc) <- paste(names(x = n.calc), assay, sep = "_")
        seu[[names(x = n.calc)]] <- n.calc
    }

    return(seu)
}

#' Propagate Metadata Changes
#'
#' @param updated_table
#' @param seu
#'
#' @return
#' @export
#'
#' @examples
propagate_spreadsheet_changes <- function(updated_table, seu) {
    meta <- updated_table

    sample_ids <- rownames(meta)

    meta <- meta %>%
        dplyr::mutate(meta, across(contains("snn"), as.factor)) %>%
        mutate(across(where(is.ordered), ~ as.factor(as.character(.x))))

    rownames(meta) <- sample_ids

    seu@meta.data <- meta

    return(seu)
}

#' Create a database of chevreul projects
#'
#' Create a database containing chevreul projects
#'
#' @param cache_location Path to cache "~/.cache/chevreul"
#' @param sqlite_db Database to be created
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
create_project_db <- function(cache_location = "~/.cache/chevreul",
    sqlite_db = "single-cell-projects.db", verbose = TRUE) {
    if (!dir.exists(cache_location)) {
        dir.create(cache_location)
    }

    con <- DBI::dbConnect(RSQLite::SQLite(), fs::path(cache_location, sqlite_db))

    projects_tbl <- tibble::tibble(
        project_name = character(),
        project_path = character(),
        project_slug = character(),
        project_type = character(),
    )

    message(paste0("building table of chevreul projects at ", fs::path(cache_location, sqlite_db)))
    # DBI::dbWriteTable(con, "projects_tbl", projects_tbl)

    tryCatch({
        DBI::dbWriteTable(con, "projects_tbl", projects_tbl)
    }, warning = function(w) {
        message(sprintf("Warning in %s: %s", deparse(w[["call"]]), w[["message"]]))
    }, error = function(e) {
        message("projects db already exists!")
    }, finally = {
    })

    DBI::dbDisconnect(con)
}

#' Update a database of chevreul projects
#'
#' Add new/update existing projects to the database by recursing fully
#'
#' @param projects_dir The project directory to be updated
#' @param cache_location Path to cache "~/.cache/chevreul"
#' @param sqlite_db sqlite db
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
update_project_db <- function(projects_dir = NULL,
    cache_location = "~/.cache/chevreul",
    sqlite_db = "single-cell-projects.db",
    verbose = TRUE) {
    if (!dir.exists(cache_location)) {
        dir.create(cache_location)
    }

    con <- DBI::dbConnect(RSQLite::SQLite(), fs::path(cache_location, sqlite_db))

    projects_tbl <-
        fs::dir_ls(projects_dir, glob = "*.here", recurse = TRUE, fail = FALSE, all = TRUE) %>%
        fs::path_dir(.) %>%
        purrr::set_names(fs::path_file(.)) %>%
        tibble::enframe("project_name", "project_path") %>%
        dplyr::mutate(project_slug = stringr::str_remove(project_name, "_proj$")) %>%
        dplyr::mutate(project_type = fs::path_file(fs::path_dir(project_path))) %>%
        identity()

    current_projects_tbl <-
        DBI::dbReadTable(con, "projects_tbl") %>%
        dplyr::filter(fs::file_exists(project_path)) %>%
        dplyr::filter(!project_path %in% projects_tbl$project_path) %>%
        dplyr::bind_rows(projects_tbl) %>%
        dplyr::distinct(project_path, .keep_all = TRUE)

    DBI::dbWriteTable(con, "projects_tbl", projects_tbl, overwrite = TRUE)

    DBI::dbDisconnect(con)
}

#' Update a database of chevreul projects
#'
#' Append projects to datatbase
#'
#' @param new_project_path
#' @param projects_dir
#' @param cache_location Path to cache "~/.cache/chevreul"
#' @param sqlite_db sqlite db
#' @param verbose
#'
#'
#' @return
#' @export
#'
#' @examples
append_to_project_db <- function(new_project_path, projects_dir = NULL,
    cache_location = "~/.cache/chevreul",
    sqlite_db = "single-cell-projects.db",
    verbose = TRUE) {
    if (!dir.exists(cache_location)) {
        dir.create(cache_location)
    }

    con <- DBI::dbConnect(RSQLite::SQLite(), fs::path(cache_location, sqlite_db))

    projects_tbl <-
        new_project_path %>%
        purrr::set_names(fs::path_file(.)) %>%
        tibble::enframe("project_name", "project_path") %>%
        dplyr::mutate(project_slug = stringr::str_remove(project_name, "_proj$")) %>%
        dplyr::mutate(project_type = fs::path_file(fs::path_dir(project_path))) %>%
        identity()

    current_projects_tbl <-
        DBI::dbReadTable(con, "projects_tbl") %>%
        dplyr::filter(fs::file_exists(project_path)) %>%
        dplyr::filter(!project_path %in% projects_tbl$project_path) %>%
        dplyr::bind_rows(projects_tbl) %>%
        dplyr::distinct(project_path, .keep_all = TRUE)

    DBI::dbWriteTable(con, "projects_tbl", current_projects_tbl, overwrite = TRUE)

    DBI::dbDisconnect(con)
}

#' Read a database of chevreul projects
#'
#' Reads database of chevreul projects to a data frame
#'
#' @param projects_dir
#' @param cache_location Path to cache "~/.cache/chevreul"
#' @param sqlite_db sqlite db
#' @param verbose
#'
#' @return
#' @export
#'
#' @examples
read_project_db <- function(projects_dir = NULL,
    cache_location = "~/.cache/chevreul",
    sqlite_db = "single-cell-projects.db",
    verbose = TRUE) {
    if (!dir.exists(cache_location)) {
        dir.create(cache_location)
    }

    con <- DBI::dbConnect(RSQLite::SQLite(), fs::path(cache_location, sqlite_db))

    current_projects_tbl <-
        DBI::dbReadTable(con, "projects_tbl")

    DBI::dbDisconnect(con)

    return(current_projects_tbl)
}

#' Make Bigwig Database
#'
#'
#' @param new_project Project directory
#' @param cache_location Path to cache "~/.cache/chevreul"
#' @param sqlite_db sqlite db containing bw files
#'
#' @return
#' @export
#'
#' @examples
make_bigwig_db <- function(new_project = NULL, cache_location = "~/.cache/chevreul/", sqlite_db = "bw-files.db") {
    new_bigwigfiles <- fs::dir_ls(new_project, glob = "*.bw", recurse = TRUE) %>%
        purrr::set_names(fs::path_file(.)) %>%
        tibble::enframe("name", "bigWig") %>%
        dplyr::mutate(sample_id = stringr::str_remove(name, "_Aligned.sortedByCoord.out.*bw$")) %>%
        dplyr::filter(!stringr::str_detect(name, "integrated")) %>%
        dplyr::distinct(sample_id, .keep_all = TRUE) %>%
        identity()

    con <- DBI::dbConnect(RSQLite::SQLite(), dbname = fs::path(cache_location, sqlite_db))

    all_bigwigfiles <-
        dbReadTable(con, "bigwigfiles") %>%
        dplyr::bind_rows(new_bigwigfiles)

    DBI::dbWriteTable(con, "bigwigfiles", all_bigwigfiles, overwrite = TRUE)

    return(all_bigwigfiles)
}

#' Retrieve Metadata from Batch
#'
#' @param batch
#' @param projects_dir path to project dir
#' @param db_path path to .db file
#'
#' @return
#'
#' @examples
metadata_from_batch <- function(batch, projects_dir = "/dataVolume/storage/single_cell_projects",
    db_path = "single-cell-projects.db") {
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

#' Swap counts from Feature
#'
#' @param cds
#' @param featureType
#'
#' @return
#' @export
#'
#' @examples
swap_counts_from_feature <- function(cds, featureType) {
    print(featureType)
    #
    #   if (featureType == "transcript"){
    #     rowData(cds[[featureType]])$gene_short_name <- rownames(cds[[featureType]])
    #   }

    assay(cds$traj, withDimnames = FALSE) <- assay(cds[[featureType]])
    rowData(cds$traj) <- rowData(cds[[featureType]])
    rownames(cds$traj) <- rownames(cds[[featureType]])
    cds$traj@preprocess_aux$gene_loadings <- cds[[featureType]]@preprocess_aux$gene_loadings
    # counts(cds$traj) <- counts(cds[[featureType]])
    cds$traj
}

#' convert seurat list to multimodal object
#'
#' @param seu_list
#'
#' @return
#' @export
#'
#' @examples
convert_seu_list_to_multimodal <- function(seu_list) {
    colnames(seu_list[["gene"]]@meta.data) <- gsub("RNA_", "gene_", colnames(seu_list[["gene"]]@meta.data))

    multimodal_seu <- seu_list$gene
    multimodal_seu <- RenameAssays(multimodal_seu, RNA = "gene")

    if ("transcript" %in% names(seu_list)) {
        if (identical(length(Cells(seu_list$gene)), length(Cells(seu_list$transcript)))) {
            colnames(seu_list[["transcript"]]@meta.data) <- gsub("RNA_", "transcript_", colnames(seu_list[["transcript"]]@meta.data))
            multimodal_seu[["transcript"]] <- seu_list$transcript$RNA
            transcript_markers <- grepl("transcript_", names(seu_list$transcript@meta.data))
            transcript_cluster_cols <- seu_list[["transcript"]]@meta.data[transcript_markers]
            if (length(transcript_cluster_cols) > 0) {
                multimodal_seu <- AddMetaData(multimodal_seu, transcript_cluster_cols)
            }
        }
    }

    marker_names <- names(Misc(multimodal_seu)[["markers"]])

    if (!is.null(multimodal_seu@misc$markers)) {
        names(multimodal_seu@misc$markers) <- gsub("RNA", "gene", marker_names)
    }


    return(multimodal_seu)
}

#' Clean Vector of Chevreul Names
#'
#' Cleans names of seurat objects provided in a vector form
#'
#' @param myvec A vector of seurat object names
#'
#' @return
#' @export
#'
#' @examples
make_chevreul_clean_names <- function(myvec) {
    myvec %>%
        purrr::set_names(stringr::str_to_title(stringr::str_replace_all(., "[^[:alnum:][:space:]\\.]", " ")))
}


#' Convert seurat object to seurat V5 format
#'
#' Convert seurat object from v3 to v5 format
#'
#' @param seu_v3 a version 3 seurat object
#'
#' @return
#' @export
#'
#' @examples
#' convert_v3_to_v5(human_gene_transcript_seu)
convert_v3_to_v5 <- function(seu_v3) {
    seurat_version <- seu_v3@misc$experiment$seurat_version

    if (seurat_version < 5 || is.null(seurat_version)) {
        meta <- seu_v3@meta.data

        seu_v5 <- CreateSeuratObject(counts = seu_v3$gene@counts, data = seu_v3$gene@data, assay = "gene", meta.data = meta)

        transcript_assay.v5 <- CreateAssay5Object(counts = seu_v3$transcript@counts, data = seu_v3$transcript@data)
        seu_v5$transcript <- transcript_assay.v5

        seu_v5$gene <- seurat_preprocess(seu_v5$gene, normalize = FALSE)
        seu_v5$transcript <- seurat_preprocess(seu_v5$transcript, normalize = FALSE)

        # seu_v5 <- clustering_workflow(seu_v5)
        seu_v5@reductions <- seu_v3@reductions
        seu_v5@graphs <- seu_v3@graphs
        seu_v5@neighbors <- seu_v3@neighbors

        seu_v5@misc <- seu_v3@misc

        Idents(seu_v5) <- Idents(seu_v3)

        seu_v5@misc$experiment$seurat_version <- packageVersion("Seurat")
    } else {
        seu_v5 <- seu_v3
    }

    return(seu_v5)
}
