#' Run Differential Expression
#'
#'
#'
#' @param object
#' @param cluster1
#' @param cluster2
#' @param resolution
#' @param diffex_scheme
#'
#' @return
#' @export
#'
#' @examples
setGeneric("run_object_de", function(object, cluster1, cluster2, resolution, diffex_scheme = "louvain", featureType, tests = c("t", "wilcox", "bimod")) standardGeneric("run_object_de"))

setMethod(
    "run_object_de", "Seurat",
    function(object, cluster1, cluster2, resolution, diffex_scheme = "louvain", featureType, tests = c("t", "wilcox", "bimod")) {
        match.arg(tests)
        if (diffex_scheme == "louvain") {
            if ("integrated" %in% names(object@assays)) {
                active_assay <- "integrated"
            } else {
                active_assay <- "gene"
            }
            Idents(object) <- paste0(active_assay, "_snn_res.", resolution)
            object <- subset(object, idents = c(cluster1, cluster2))
        } else if (diffex_scheme == "feature") {
            object <- object[, c(cluster1, cluster2)]
            keep_cells <- c(cluster1, cluster2)
            new_idents <- c(rep(1, length(cluster1)), rep(2, length(cluster2)))
            names(new_idents) <- keep_cells
            new_idents <- new_idents[colnames(object)]
            Idents(object) <- new_idents
            cluster1 <- 1
            cluster2 <- 2
        }
        test_list <- vector("list", length(tests))
        for (test in tests) {
            print(test)
            de <- FindMarkers(object, assay = featureType, ident.1 = cluster1, ident.2 = cluster2, test.use = test)
            if (featureType == "transcript") {
                de_cols <- c("enstxp", "ensgene", "symbol", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")
                de <- de %>%
                    tibble::rownames_to_column("enstxp") %>%
                    dplyr::left_join(annotables::grch38_tx2gene, by = "enstxp") %>%
                    dplyr::left_join(annotables::grch38, by = "ensgene")
                if ("avg_logFC" %in% colnames(de)) {
                    de <- dplyr::mutate(de, avg_log2FC = log(exp(avg_logFC), 2))
                }
                de <- dplyr::select(de, any_of(de_cols))
            } else if (featureType == "gene") {
                de_cols <- c("ensgene", "symbol", "p_val", "avg_log2FC", "pct.1", "pct.2", "p_val_adj")
                de <- de %>%
                    tibble::rownames_to_column("symbol") %>%
                    dplyr::left_join(annotables::grch38, by = "symbol")
                if ("avg_logFC" %in% colnames(de)) {
                    de <- dplyr::mutate(de, avg_log2FC = log(exp(avg_logFC), 2))
                }
                de <- dplyr::select(de, any_of(de_cols))
            }
            test_list[[match(test, tests)]] <- de
        }
        names(test_list) <- tests
        return(test_list)
    }
)

setMethod(
    "run_object_de", "SingleCellExperiment",
    function(object, cluster1, cluster2, resolution, diffex_scheme = "louvain", featureType, tests = c("t", "wilcox", "bimod")) {
        match.arg(tests)

        if (featureType == "transcript") object <- altExp(object, "transcript")

        if (diffex_scheme == "louvain") {
            if ("integrated" %in% names(object@assays)) {
                active_assay <- "integrated"
            } else {
                active_assay <- "gene"
            }
            colLabels(object) <- colData(object)[[paste0(active_assay, "_snn_res.", resolution)]]
            object <- object[, colLabels(object) %in% c(cluster1, cluster2)]
            colLabels(object) <- factor(colLabels(object))
        } else if (diffex_scheme == "feature") {
            object <- object[, c(cluster1, cluster2)]
            keep_cells <- c(cluster1, cluster2)
            new_idents <- c(rep(1, length(cluster1)), rep(2, length(cluster2)))
            names(new_idents) <- keep_cells
            new_idents <- new_idents[colnames(object)]
            colLabels(object) <- new_idents
            cluster1 <- 1
            cluster2 <- 2
        }
        test_list <- vector("list", length(tests))
        for (test in tests) {
            print(test)
            de <- scran::findMarkers(object, test.type = test)
            if (featureType == "transcript") {
                de_cols <- c("enstxp", "ensgene", "symbol", "p_val" = "p.value", "avg_log2FC", "pct.1", "pct.2", "p_val_adj" = "FDR")
                de <- de[[1]] %>%
                    as.data.frame() %>%
                    tibble::rownames_to_column("enstxp") %>%
                    dplyr::left_join(annotables::grch38_tx2gene, by = "enstxp") %>%
                    dplyr::left_join(annotables::grch38, by = "ensgene")
                if ("summary.logFC" %in% colnames(de)) {
                    de <- dplyr::mutate(de, avg_log2FC = log(exp(summary.logFC), 2))
                }
                de <- dplyr::select(de, any_of(de_cols))
            } else if (featureType == "gene") {
                de_cols <- c("ensgene", "symbol", "p_val" = "p.value", "avg_log2FC", "pct.1", "pct.2", "p_val_adj" = "FDR")
                de <- de[[1]] %>%
                    as.data.frame() %>%
                    tibble::rownames_to_column("symbol") %>%
                    dplyr::left_join(annotables::grch38, by = "symbol")
                if ("summary.logFC" %in% colnames(de)) {
                    de <- dplyr::mutate(de, avg_log2FC = log(exp(summary.logFC), 2))
                }
                de <- dplyr::select(de, any_of(de_cols))
            }
            test_list[[match(test, tests)]] <- de
        }
        names(test_list) <- tests
        return(test_list)
    }
)




#' Run Enrichment Browser on Differentially Expressed Genes
#'
#' @param object
#' @param enrichment_method
#' @param ...
#' @param cluster_list
#' @param de_results
#'
#' @return
#' @export
#'
#' @examples
run_enrichmentbrowser <- function(object, cluster_list, de_results, enrichment_method = c("ora"), ...) {
    cluster1_cells <- cluster_list$cluster1
    cluster2_cells <- cluster_list$cluster2

    test_diffex_results <- de_results$t %>%
        dplyr::mutate(FC = log2(exp(avg_log2FC))) %>%
        dplyr::mutate(ADJ.PVAL = p_val_adj) %>%
        dplyr::distinct(symbol, .keep_all = TRUE) %>%
        tibble::column_to_rownames("symbol") %>%
        identity()

    # subset by supplied cell ids
    #
    object <- object[, c(cluster1_cells, cluster2_cells)]
    object <- object[rownames(object) %in% de_results$t$symbol, ]

    object[["gene"]]@meta.features <- test_diffex_results

    keep_cells <- c(cluster1_cells, cluster2_cells)
    new_idents <- c(rep(0, length(cluster1_cells)), rep(1, length(cluster2_cells)))
    names(new_idents) <- keep_cells
    new_idents <- new_idents[colnames(object)]
    Idents(object) <- new_idents

    counts <- GetAssayData(object, assay = "gene", slot = "counts")
    counts <- as.matrix(counts)
    mode(counts) <- "integer"

    rowData <- data.frame(
        FC = object[["gene"]][[]]$FC,
        ADJ.PVAL = object[["gene"]][[]]$ADJ.PVAL,
        row.names = rownames(object@assays[["gene"]])
    )

    colData <- as.data.frame(pull_metadata(object))

    se <- SummarizedExperiment::SummarizedExperiment(
        assays = list(counts = counts),
        rowData = rowData, colData = colData
    )

    se$GROUP <- forcats::fct_inseq(Idents(object))

    # se <- EnrichmentBrowser::deAna(se, grp = se$GROUP, de.method = "edgeR")
    se <- EnrichmentBrowser::idMap(se, org = "hsa", from = "SYMBOL", to = "ENTREZID")

    outdir <- fs::path("www", "enrichmentbrowser")

    # hsa.grn <- EnrichmentBrowser::compileGRN(org="hsa", db="kegg")

    # EnrichmentBrowser::ebrowser( meth=c("ora", "ggea"), perm=0, comb=TRUE,
    #           exprs=se, gs=go.gs, grn=hsa.grn, org="hsa", nr.show=3,
    #           out.dir=outdir, report.name=report.name, browse = FALSE)


    enrichment.res <- list()

    if ("ora" %in% enrichment_method) {
        enrichment.res$ora <- EnrichmentBrowser::sbea(
            method = "ora", se = se, gs = go.gs, perm = 0,
            alpha = 0.1
        )
        results <- enrichment.res$ora
        report.name <- "ora.html"
        EnrichmentBrowser::eaBrowse(results,
            html.only = TRUE, out.dir = outdir, graph.view = hsa.grn,
            report.name = report.name
        )
    }

    if ("gsea" %in% enrichment_method) {
        enrichment.res$gsea <- EnrichmentBrowser::sbea(
            method = "gsea", se = se, gs = msigdb.gs, perm = 100,
            alpha = 0.1
        )
        results <- enrichment.res$gsea
        report.name <- "gsea.html"
        EnrichmentBrowser::eaBrowse(results,
            html.only = TRUE, out.dir = outdir, graph.view = hsa.grn,
            report.name = report.name
        )
    }

    if ("nbea" %in% enrichment_method) {
        enrichment.res$nbea <- EnrichmentBrowser::nbea(method = "ggea", se = se, gs = go.gs, grn = hsa.grn)

        results <- enrichment.res$nbea
        report.name <- "nbea.html"
        EnrichmentBrowser::eaBrowse(results,
            html.only = TRUE, out.dir = outdir, graph.view = hsa.grn,
            report.name = report.name
        )
    }

    return(list(report = fs::path("enrichmentbrowser", report.name), results = results$res.tbl))
}

#' Prep Slider Values
#'
#' @param default_val
#'
#' @return
#' @export
#'
#' @examples
prep_slider_values <- function(default_val) {
    min <- round(default_val * 0.25, digits = 1)
    max <- round(default_val * 2.0, digits = 1)
    step <- 10^((ceiling(log10(default_val))) - 1)
    value <- default_val
    return(list(min = min, max = max, value = value, step = step))
}


#' Create a shiny app for a project on disk
#'
#' @param preset_project A preloaded project to start the app with
#' @param appTitle A title of the App
#' @param futureMb amount of Mb allocated to future package
#' @param preset_project
#' @param organism_type
#'
#' @return
#' @export
#'
#' @examples
objectApp <- function(preset_project, appTitle = "chevreul", organism_type = "human", db_path = "~/.cache/chevreul/single-cell-projects.db", futureMb = 13000) {
    print(packageVersion("chevreul"))
    future::plan(strategy = "multicore", workers = 6)
    future_size <- futureMb * 1024^2
    options(future.globals.maxSize = future_size)
    options(shiny.maxRequestSize = 40 * 1024^2)
    options(DT.options = list(
        pageLength = 2000, paging = FALSE,
        info = TRUE, searching = TRUE, autoWidth = F, ordering = TRUE,
        scrollX = TRUE, language = list(search = "Filter:")
    ))
    header <- shinydashboard::dashboardHeader(title = appTitle)
    sidebar <- shinydashboard::dashboardSidebar(
        uiOutput("projInput"),
        actionButton("loadProject", "Load Selected Project") %>%
            default_helper(type = "markdown", content = "overview"),
        textOutput("appTitle"),
        bookmarkButton(),
        shinyWidgets::prettyRadioButtons("organism_type",
            inline = TRUE,
            "Organism", choices = c("human", "mouse"), selected = organism_type
        ),
        shinyFiles::shinyFilesButton("objectUpload", "Load a Seurat Dataset",
            "Please select a .rds file",
            multiple = FALSE
        ),
        shinyFiles::shinySaveButton("saveSeurat", "Save Current Dataset",
            "Save file as...",
            filetype = list(rds = "rds")
        ),
        verbatimTextOutput("savefile"),
        shinydashboard::sidebarMenu(
            shinydashboard::menuItem("Integrate Projects",
                tabName = "integrateProjects", icon = icon("object-group")
            ), shinydashboard::menuItem("Reformat Metadata",
                tabName = "reformatMetadataDR", icon = icon("columns")
            ), shinydashboard::menuItem("Plot Data",
                tabName = "comparePlots", icon = icon("chart-bar"), selected = TRUE
            ), shinydashboard::menuItem("Heatmap/Violin Plots",
                tabName = "violinPlots", icon = icon("sort")
            ), shinydashboard::menuItem("Coverage Plots",
                tabName = "coveragePlots", icon = icon("mountain")
            ), shinydashboard::menuItem("Differential Expression",
                tabName = "diffex", icon = icon("magnet")
                # ), shinydashboard::menuItem("Pathway Enrichment Analysis",
                #   tabName = "pathwayEnrichment", icon = icon("sitemap")
            ), shinydashboard::menuItem("Find Markers",
                tabName = "findMarkers", icon = icon("bullhorn")
            ), shinydashboard::menuItem("Subset Seurat Input",
                tabName = "subsetSeurat", icon = icon("filter")
            ), shinydashboard::menuItem("All Transcripts",
                tabName = "allTranscripts", icon = icon("sliders-h")
            ), shinydashboard::menuItem("RNA Velocity",
                tabName = "plotVelocity", icon = icon("tachometer-alt")
            ), shinydashboard::menuItem("Monocle",
                tabName = "monocle", icon = icon("bullseye")
            ), shinydashboard::menuItem("Regress Features",
                tabName = "regressFeatures", icon = icon("eraser")
            ), shinydashboard::menuItem("Technical Information",
                tabName = "techInfo", icon = icon("cogs")
            )
        ),
        actionButton("changeEmbedAction",
            label = "Change Embedding Parameters"
        ),
        changeEmbedParamsui("changeembed"),
        width = 250
    )
    body <- shinydashboard::dashboardBody(
        waiter::use_waiter(),
        shinydashboard::tabItems(
            shinydashboard::tabItem(
                tabName = "comparePlots",
                h2("Compare Plots") %>%
                    default_helper(type = "markdown", content = "comparePlots"),
                plotDimRedui("plotdimred1"),
                plotDimRedui("plotdimred2"),
                plotReadCountui("plotreadcount1"),
                plotReadCountui("plotreadcount2"),
                chevreulBox(
                    title = "Selected Cells",
                    tableSelectedui("tableselected"),
                    width = 6
                ),
                plotClustree_UI("clustreePlot")
            ),
            shinydashboard::tabItem(
                tabName = "violinPlots",
                fluidRow(
                    plotHeatmapui("heatMap")
                ),
                fluidRow(
                    plotViolinui("violinPlot")
                )
            ),
            shinydashboard::tabItem(
                tabName = "coveragePlots",
                fluidRow(
                    plotCoverage_UI("coverageplots")
                )
            ),
            shinydashboard::tabItem(
                tabName = "integrateProjects",
                fluidRow(
                    integrateProjui("integrateproj"),
                    width = 12
                )
            ),
            shinydashboard::tabItem(
                tabName = "reformatMetadataDR",
                fluidRow(
                    reformatMetadataDRui("reformatMetadataDR")
                )
            ),
            shinydashboard::tabItem(
                tabName = "subsetSeurat",
                h2("Subset Seurat Input") %>%
                    default_helper(type = "markdown", content = "subsetSeurat"),
                plotDimRedui("subset"),
                chevreulBox(
                    title = "Subset Settings",
                    checkboxInput("legacySettingsSubset", "Use Legacy Settings", value = FALSE),
                    actionButton("subsetAction", "subset object by selected cells"),
                    actionButton("subsetCsv", "subset object by uploaded csv"),
                    fileInput("uploadCsv",
                        "Upload .csv file with cells to include",
                        accept = c(".csv")
                    ),
                    shinyjs::useShinyjs(),
                    # shinyjs::runcodeUI(code = "shinyjs::alert('Hello!')"),
                    textOutput("subsetMessages"),
                    width = 6
                ),
                chevreulBox(
                    title = "Selected Cells", tableSelectedui("subset"),
                    width = 6
                )
            ), shinydashboard::tabItem(
                tabName = "findMarkers",
                h2("Find Markers"),
                findMarkersui("findmarkers"),
                plotDimRedui("markerScatter")
            ), shinydashboard::tabItem(
                tabName = "allTranscripts",
                h2("All Transcripts"),
                plotDimRedui("alltranscripts2"),
                allTranscriptsui("alltranscripts1")
            ),
            shinydashboard::tabItem(
                tabName = "plotVelocity",
                h2("RNA Velocity"),
                fluidRow(
                    plotVelocityui("plotvelocity"),
                )
            ),
            shinydashboard::tabItem(
                tabName = "diffex",
                h2("Differential Expression") %>%
                    default_helper(type = "markdown", content = "diffex"),
                plotDimRedui("diffex"),
                diffexui("diffex")
            ),
            # shinydashboard::tabItem(
            #   tabName = "pathwayEnrichment",
            #   h2("Pathway Enrichment"),
            #   fluidRow(
            #     pathwayEnrichmentui("pathwayEnrichment")
            #   )
            # ),
            shinydashboard::tabItem(
                tabName = "regressFeatures",
                fluidRow(
                    chevreulBox(
                        title = "Regress Features",
                        actionButton(
                            "regressAction",
                            "Regress Seurat Objects By Genes"
                        ),
                        checkboxInput("runRegression",
                            "Run Regression?",
                            value = FALSE
                        ), radioButtons("priorGeneSet",
                            "Choose a marker gene set:",
                            choices = c(
                                "Apoptosis",
                                "Cell Cycle",
                                "Mitochondrial",
                                "Ribosomal"
                            )
                        ), selectizeInput("geneSet",
                            "List of genes",
                            choices = NULL, multiple = TRUE
                        ),
                        textInput("geneSetName", "Name for Gene Set"),
                        width = 12
                    ) %>%
                        default_helper(type = "markdown", content = "regressFeatures")
                )
            ), shinydashboard::tabItem(
                tabName = "monocle",
                h2("Monocle"),
                monocleui("monocle")
            ), shinydashboard::tabItem(
                tabName = "techInfo",
                h2("Technical Information"),
                h3(paste0("chevreul version: ", packageVersion("chevreul"))),
                techInfoui("techInfo")
            )
        )
    )

    ui <- function(request) {
        ui <- dashboardPage(
            header = header, sidebar = sidebar,
            body = body
        )
    }
    server <- function(input, output, session) {
        # shinyjs::runcodeServer()
        # observe({
        #   list_of_inputs <- reactiveValuesToList(input)
        #   list_of_inputs <<- reactiveValuesToList(input)
        #   print(list_of_inputs)
        #
        #   retained_inputs <- c("setProject")
        #
        #   list_of_inputs[!list_of_inputs %in% retained_inputs]
        # })
        # setBookmarkExclude(names(list_of_inputs))

        onBookmark(function(state) {
            state$values$uploadSeuratPath <- uploadSeuratPath()
        })

        onRestore(function(state) {
            uploadSeuratPath(state$values$uploadSeuratPath[[1]])
        })

        w <- waiter::Waiter$new()

        # lib.loc = "/dataVolume/storage/rpkgs/devel_install/"
        shinyhelper::observe_helpers(help_dir = system.file("helpers", package = "chevreul"))
        options(warn = -1)
        # shinylogs::track_usage(storage_mode = shinylogs::store_json(path = "logs/"))
        # projects_db <- "/dataVolume/storage/single_cell_projects/single_cell_projects.db"

        con <- reactive({
            DBI::dbConnect(
                RSQLite::SQLite(),
                db_path
            )
        })

        projList <- reactivePoll(4000, session, checkFunc = function() {
            if (file.exists(db_path)) {
                DBI::dbReadTable(con(), "projects_tbl") %>%
                    tibble::deframe()
            } else {
                ""
            }
        }, valueFunc = function() {
            DBI::dbReadTable(con(), "projects_tbl") %>%
                tibble::deframe()
        })

        output$projInput <- renderUI({
            selectizeInput("setProject", "Select Project to Load",
                choices = projList(), selected = preset_project,
                multiple = F
            )
        })

        proj_matrices <- reactive({
            create_proj_matrix(projList())
        })

        object <- reactiveVal()
        proj_dir <- reactiveVal()
        uploadSeuratPath <- reactiveVal()

        if (!is.null(preset_project)) {
            proj_dir(preset_project)
        }
        organism_type <- reactive({
            input$organism_type
        })
        plot_types <- reactive({
            list_plot_types(object())
        })
        observeEvent(input$loadProject, {
            proj_dir(input$setProject)
        })
        output$appTitle <- renderText({
            req(proj_dir())
            paste0("Loaded Project: ", fs::path_file(proj_dir()))
        })
        dataset_volumes <- reactive({
            print(proj_dir())
            dataset_volumes <- c(
                Home = fs::path(
                    proj_dir(),
                    "output", "seurat"
                ), "R Installation" = R.home(),
                shinyFiles::getVolumes()()
            )
        })
        observe({
            req(dataset_volumes())
            shinyFiles::shinyFileChoose(input, "objectUpload",
                roots = dataset_volumes(), session = session
            )
        })
        observeEvent(input$objectUpload, {
            req(dataset_volumes())

            file <- shinyFiles::parseFilePaths(
                dataset_volumes(),
                input$objectUpload
            )
            print(file)
            uploadSeuratPath(file$datapath)
        })

        observe({
            req(uploadSeuratPath())
            print("uploaded")
            shiny::withProgress(
                message = paste0("Uploading Data"),
                value = 0,
                {
                    shiny::incProgress(2 / 10)
                    print(uploadSeuratPath())
                    updated_object <- update_chevreul_object(object_path = uploadSeuratPath(), organism = organism)
                    object(updated_object)
                    shiny::incProgress(6 / 10)

                    organism <- case_when(
                        stringr::str_detect(uploadSeuratPath(), "Hs") ~ "human",
                        stringr::str_detect(uploadSeuratPath(), "Mm") ~ "mouse"
                    )

                    print(uploadSeuratPath())
                    print(names(object))
                    shiny::incProgress(8 / 10)
                }
            )
        })


        featureType <- reactive({
            "gene"
        })

        observe({
            req(dataset_volumes())
            shinyFiles::shinyFileSave(input, "saveSeurat",
                # roots = dataset_volumes(),
                roots = c(Home = fs::path(
                    proj_dir(),
                    "output", "seurat"
                )),
                session = session, restrictions = system.file(package = "base")
            )
        })
        subSeuratPath <- eventReactive(input$saveSeurat, {
            req(object())
            savefile <- shinyFiles::parseSavePath(
                dataset_volumes(),
                input$saveSeurat
            )
            return(savefile$datapath)
        })
        observeEvent(input$saveSeurat, {
            req(object())
            req(subSeuratPath())

            shiny::withProgress(
                message = paste0("Saving Data"),
                value = 0,
                {
                    Sys.sleep(6)
                    shiny::incProgress(2 / 10)
                    saveRDS(
                        object(),
                        subSeuratPath()
                    )
                    shiny::incProgress(10 / 10)
                }
            )
        })
        integrationResults <- callModule(
            integrateProj, "integrateproj",
            proj_matrices, object, proj_dir, con()
        )
        newprojList <- reactive({
            req(integrationResults())
            integration_path <- paste0(integrationResults())
            proj_dir(integration_path)
            newintegrated_project <- purrr::set_names(
                integration_path,
                fs::path_file(integration_path)
            )
            newprojList <- c(projList(), newintegrated_project)
        })
        observe({
            # print(newprojList())
            updateSelectizeInput(session, "setProject",
                label = "Select input label",
                choices = newprojList(),
                server = TRUE
            )
        })

        observe({
            reformatted_object <- callModule(reformatMetadataDR, "reformatMetadataDR", object, featureType)
            object(reformatted_object())
        })

        reductions <- reactive({
            req(object())
            names(object()@reductions)
            # c("pca", "tsne", "umap")
        })

        observe({
            req(object())

            callModule(plotDimRed, "plotdimred1", object, plot_types, featureType,
                organism_type = organism_type, reductions
            )
            callModule(plotDimRed, "plotdimred2", object, plot_types, featureType,
                organism_type = organism_type, reductions
            )
            callModule(plotDimRed, "diffex", object, plot_types, featureType,
                organism_type = organism_type, reductions
            )
            callModule(plotDimRed, "subset", object, plot_types, featureType,
                organism_type = organism_type, reductions
            )
            callModule(plotDimRed, "markerScatter", object, plot_types, featureType,
                organism_type = organism_type, reductions
            )
        })

        callModule(plotReadCount, "plotreadcount1", object, plot_types)
        callModule(plotReadCount, "plotreadcount2", object, plot_types)
        callModule(
            plotViolin, "violinPlot", object, featureType,
            organism_type
        )
        callModule(
            plotHeatmap, "heatMap", object, featureType,
            organism_type
        )

        callModule(
            plotCoverage, "coverageplots", object, plot_types, proj_dir, organism_type
        )
        callModule(plotClustree, "clustreePlot", object)
        callModule(tableSelected, "tableselected", object)
        diffex_selected_cells <- callModule(
            tableSelected, "diffex",
            object
        )

        callModule(pathwayEnrichment, "pathwayEnrichment", object, featureType)

        subset_selected_cells <- callModule(
            tableSelected, "subset",
            object
        )
        observeEvent(input$subsetAction, {
            req(subset_selected_cells())
            withCallingHandlers(
                {
                    shinyjs::html("subsetMessages", "")
                    message("Beginning")

                    subset_object <- object()[, colnames(object()) %in% subset_selected_cells()]
                    object(subset_object)
                    if (length(unique(object()[["batch"]])) > 1) {
                        message(paste0("reintegrating gene expression"))
                        reintegrated_object <- reintegrate_object(object(),
                            resolution = seq(0.2, 2, by = 0.2),
                            legacy_settings = input$legacySettingsSubset,
                            organism = Misc(object())$experiment$organism
                        )
                        object(reintegrated_object)
                    } else {
                        subset_object <- object_pipeline(object(), resolution = seq(0.2, 2, by = 0.2), legacy_settings = input$legacySettingsSubset) # markermarker
                        object(subset_object)
                    }
                    message("Complete!")
                },
                message = function(m) {
                    shinyjs::html(id = "subsetMessages", html = paste0(
                        "Subsetting Seurat Object: ",
                        m$message
                    ), add = FALSE)
                }
            )
        })
        observeEvent(input$subsetCsv, {
            req(input$subsetCsv)
            req(input$uploadCsv)
            withCallingHandlers(
                {
                    shinyjs::html("subsetMessages", "")
                    message("Beginning")
                    subset_object <- subset_by_meta(
                        input$uploadCsv$datapath,
                        object()
                    )
                    object(subset_object)
                    if (length(unique(object()[["batch"]])) > 1) {
                        message(paste0("reintegrating gene expression"))
                        reintegrated_object <- reintegrate_object(object(),
                            resolution = seq(0.2, 2, by = 0.2),
                            legacy_settings = input$legacySettingsSubset,
                            organism = Misc(object())$experiment$organism
                        )
                        object(reintegrated_object)
                    } else {
                        subset_object <- object_pipeline(object(),
                            resolution = seq(0.2,
                                2,
                                by = 0.2
                            ),
                            legacy_settings = input$legacySettingsSubset
                        )
                        object(subset_object)
                    }
                    message("Complete!")
                },
                message = function(m) {
                    shinyjs::html(id = "subsetMessages", html = paste0(
                        "Subsetting Seurat Object: ",
                        m$message
                    ), add = FALSE)
                }
            )
        })
        observeEvent(input$changeEmbedAction, {
            showModal(modalDialog(
                title = "Recalculating Embedding",
                "This process may take a minute or two!"
            ))
            object <- callModule(
                changeEmbedParams, "changeembed",
                object
            )
            removeModal()
        })
        callModule(findMarkers, "findmarkers", object, plot_types, featureType)
        diffex_results <- callModule(
            diffex, "diffex", object, featureType,
            diffex_selected_cells
        )

        observe({
            req(featureType())
            req(object())
            callModule(
                allTranscripts, "alltranscripts1", object, featureType,
                organism_type
            )

            callModule(plotDimRed, "alltranscripts2", object, plot_types, featureType,
                organism_type = organism_type, reductions
            )
        })

        prior_gene_set <- reactive({
            # req(input$priorGeneSet)
            req(object())

            if (is.null(input$priorGeneSet)) {
                ""
            } else if (input$priorGeneSet == "Apoptosis") {
                marker_genes <- c(
                    "CASP3", "CASP7", "BAX", "BAK1", "BID", "BBC3",
                    "BCL2", "MCL1"
                )

                marker_genes[marker_genes %in% rownames(object())]
            } else if (input$priorGeneSet == "Cell Cycle") {
                marker_genes <- c(
                    "MCM5", "PCNA", "TYMS", "FEN1", "MCM2", "MCM4",
                    "RRM1", "UNG", "GINS2", "MCM6", "CDCA7", "DTL",
                    "PRIM1", "UHRF1", "MLF1IP", "HELLS", "RFC2",
                    "RPA2", "NASP", "RAD51AP1", "GMNN", "WDR76",
                    "SLBP", "CCNE2", "UBR7", "POLD3", "MSH2",
                    "ATAD2", "RAD51", "RRM2", "CDC45", "CDC6",
                    "EXO1", "TIPIN", "DSCC1", "BLM", "CASP8AP2",
                    "USP1", "CLSPN", "POLA1", "CHAF1B", "BRIP1",
                    "E2F8"
                )

                marker_genes[marker_genes %in% rownames(object())]
            } else if (input$priorGeneSet == "Mitochondrial") {
                marker_genes <- mito_features[[organism_type()]][["gene"]]

                marker_genes[marker_genes %in% rownames(object())]
            } else if (input$priorGeneSet == "Ribosomal") {
                marker_genes <- ribo_features[[organism_type()]][["gene"]]

                marker_genes[marker_genes %in% rownames(object())]
            }
        })

        observe({
            updateSelectizeInput(session, "geneSet",
                choices = rownames(object()),
                selected = prior_gene_set(),
                server = TRUE
            )
        })
        observeEvent(input$regressAction, {
            req(object())
            showModal(modalDialog(
                title = "Regressing out provided list of features",
                "This process may take a minute or two!"
            ))
            regressed_object <- chevreul::regress_by_features(object(),
                feature_set = list(input$geneSet), set_name = janitor::make_clean_names(input$geneSetName),
                regress = input$runRegression
            )
            object(regressed_object)
            removeModal()
        })

        observe({
            req(reductions())
            callModule(
                monocle, "monocle", object, plot_types, featureType,
                organism_type, reductions
            )
        })


        observe({
            req(uploadSeuratPath())
            req(object())

            proj_path <- stringr::str_replace(uploadSeuratPath(), "output.*", "")

            proj_name <- fs::path_file(proj_path)
            print(proj_name)

            loom_path <- fs::path(proj_path, "output", "velocyto", paste0(proj_name, ".loom"))
            print(loom_path)
            # need to check if this file exists

            callModule(plotVelocity, "plotvelocity", object, loom_path)
        })

        callModule(techInfo, "techInfo", object)

        sessionId <- as.integer(runif(1, 1, 100000))
        output$sessionId <- renderText(paste0("Session id: ", sessionId))
        session$onSessionEnded(function() {
            cat(paste0("Ended: ", sessionId))
            observe(DBI::dbConnect(con()))
        })
    }

    # onStop(function() DBI::dbDisconnect(con))
    shinyApp(ui, server, enableBookmarking = "url")
}
