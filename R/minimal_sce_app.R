#' Create a minimal chevreul app using SingleCellExperiment input
#'
#' @param single_cell_object a singlecell object
#' @param loom_path path to a loom file
#' @param appTitle a title for the app
#' @param organism_type human or mouse
#' @param futureMb the megabytes available for the future package
#'
#' @return a minimal chevreul app
#' @export
#'
#' @examples
#' \donttest{
#' minimalSceApp(human_gene_transcript_sce)
#' }
#'
minimalSceApp <- function(single_cell_object = human_gene_transcript_sce, loom_path = NULL, appTitle = NULL, organism_type = "human", futureMb = 13000, bigwig_db = "~/.cache/chevreul/bw-files.db") {
    future::plan(strategy = "multicore", workers = 6)
    future_size <- futureMb * 1024^2
    options(future.globals.maxSize = future_size)
    options(shiny.maxRequestSize = 40 * 1024^2)
    options(DT.options = list(
        pageLength = 2000, paging = FALSE,
        info = TRUE, searching = TRUE, autoWidth = F, ordering = TRUE,
        scrollX = TRUE, language = list(search = "Filter:")
    ))
    header <- dashboardHeader(title = appTitle)
    sidebar <- dashboardSidebar(
        textOutput("appTitle"),
        uiOutput("featureType"),
        sidebarMenu(
            menuItem("Reformat Metadata",
                tabName = "reformatMetadata", icon = icon("columns")
            ), menuItem("Plot Data",
                tabName = "comparePlots", icon = icon("chart-bar"), selected = TRUE
            ), menuItem("Heatmap/Violin Plots",
                tabName = "violinPlots", icon = icon("sort")
                # ), menuItem("Coverage Plots",
                #   tabName = "coveragePlots", icon = icon("mountain")
            ), menuItem("Differential Expression",
                tabName = "diffex", icon = icon("magnet")
                # ), menuItem("Pathway Enrichment Analysis",
                #   tabName = "pathwayEnrichment", icon = icon("sitemap")
            ), menuItem("Find Markers",
                tabName = "findMarkers", icon = icon("bullhorn")
            ), menuItem("Subset Object Input",
                tabName = "subsetObject", icon = icon("filter")
            ), menuItem("All Transcripts",
                tabName = "allTranscripts", icon = icon("sliders-h")
            ), menuItem("Regress Features",
                tabName = "regressFeatures", icon = icon("eraser")
            ), menuItem("Technical Information",
                tabName = "techInfo", icon = icon("cogs")
            )
        ),
        actionButton("changeEmbedAction",
            label = "Change Embedding Parameters"
        ),
        changeEmbedParamsui("changeembed"),
        width = 250
    )
    body <- dashboardBody(
        use_waiter(),
        tabItems(
            tabItem(
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
            tabItem(
                tabName = "violinPlots",
                fluidRow(
                    plotHeatmapui("heatMap")
                ),
                fluidRow(
                    plotViolinui("violinPlot")
                )
            ),
            tabItem(
                tabName = "coveragePlots",
                fluidRow(
                    plotCoverage_UI("coverageplots")
                )
            ),
            tabItem(
                tabName = "reformatMetadata",
                fluidRow(
                    reformatMetadataDRui("reformatMetadataDR")
                )
            ),
            tabItem(
                tabName = "subsetObject",
                h2("Subset Object Input") %>%
                    default_helper(type = "markdown", content = "subsetObject"),
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
                    useShinyjs(),
                    # runcodeUI(code = "alert('Hello!')"),
                    textOutput("subsetMessages"),
                    width = 6
                ),
                chevreulBox(
                    title = "Selected Cells", tableSelectedui("subset"),
                    width = 6
                )
            ), tabItem(
                tabName = "findMarkers",
                h2("Find Markers"),
                chevreulMarkersui("findmarkers"),
                plotDimRedui("markerScatter")
            ), tabItem(
                tabName = "allTranscripts",
                h2("All Transcripts"),
                plotDimRedui("alltranscripts2"),
                allTranscriptsui("alltranscripts1")
            ),
            tabItem(
                tabName = "diffex",
                h2("Differential Expression") %>%
                    default_helper(type = "markdown", content = "diffex"),
                plotDimRedui("diffex"),
                diffexui("diffex")
            ),
            tabItem(
                tabName = "pathwayEnrichment",
                h2("Pathway Enrichment"),
                fluidRow(
                    pathwayEnrichmentui("pathwayEnrichment")
                )
            ),
            tabItem(
                tabName = "regressFeatures",
                fluidRow(
                    chevreulBox(
                        title = "Regress Features",
                        actionButton(
                            "regressAction",
                            "Regress Objects By Genes"
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
            ), tabItem(
                tabName = "techInfo",
                h2("Technical Information"),
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
        w <- Waiter$new()

        observe_helpers(help_dir = system.file("helpers", package = "chevreul", lib.loc = "/dataVolume/storage/rpkgs/devel_install/"))
        options(warn = -1)

        object <- reactiveVal(NULL)
        observe({
            object(single_cell_object)
        })

        organism_type <- reactive({
            "human"
        })

        loom_path <- reactive({
            loom_path
        })

        plot_types <- reactive({
            req(object())
            list_plot_types(object())
        })

        featureType <- reactive({
            "gene"
            # input$feature_type
        })


        observe({
            reformatted_object <- callModule(reformatMetadataDR, "reformatMetadataDR", object, featureType)
            object(reformatted_object())
        })

        reductions <- reactive({
            req(object())
          reducedDimNames(object())

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
        subset_selected_cells <- callModule(
            tableSelected, "subset",
            object
        )

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

        callModule(chevreulMarkers, "findmarkers", object, plot_types, featureType)

        callModule(pathwayEnrichment, "pathwayEnrichment", object)

        # # plot all transcripts
        # observe({
        #   # req(featureType())
        #   req(object())
        #   if (query_experiment(object(), "transcript")) {
        #     callModule(
        #       allTranscripts, "alltranscripts1", object, featureType,
        #       organism_type
        #     )
        #
        #     callModule(plotDimRed, "alltranscripts2", object, plot_types, featureType,
        #                organism_type = organism_type, reductions
        #     )
        #   }
        # })

        # plot all transcripts
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

        diffex_results <- callModule(
            diffex, "diffex", object, featureType,
            diffex_selected_cells
        )
        observeEvent(input$subsetAction, {
            req(input$subsetAction)
            req(subset_selected_cells())
            withCallingHandlers(
                {
                    html("subsetMessages", "")
                    message("Beginning")
                    subset_object <- object()[, colnames(object()) %in% subset_selected_cells()]
                    object(subset_object)
                    if (length(unique(object()$batch)) > 1) {
                        message(paste0("reintegrating gene expression"))
                        reintegrated_object <- reintegrate_object(object(),
                            resolution = seq(0.2, 2, by = 0.2),
                            organism = metadata(object())$experiment$organism
                        )
                        object(reintegrated_object)
                    } else {
                        processed_object <- object_pipeline(object(), resolution = seq(0.2, 2, by = 0.2))
                        object(processed_object)
                    }
                    message("Complete!")
                },
                message = function(m) {
                    html(id = "subsetMessages", html = paste0(
                        "Subsetting Object: ",
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
                    html("subsetMessages", "")
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
                            organism = metadata(object())$experiment$organism
                        )
                        object(reintegrated_object)
                    } else {
                        processed_object <- object_pipeline(object(), resolution = seq(0.2, 2, by = 0.2))
                        object(processed_object)
                    }
                    message("Complete!")
                },
                message = function(m) {
                    html(id = "subsetMessages", html = paste0(
                        "Subsetting Object: ",
                        m$message
                    ), add = FALSE)
                }
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
            req(object())
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
            regressed_object <- regress_by_features(object(),
                feature_set = list(input$geneSet), set_name = make_clean_names(input$geneSetName),
                regress = input$runRegression
            )
            object(regressed_object)
            removeModal()
        })

        callModule(techInfo, "techInfo", object)
    }
    shinyApp(ui, server, enableBookmarking = "server")
}
