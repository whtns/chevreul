#' Create a obect app on a shinyproxy docker instance
#'
#' @param object a obect object
#' @param appTitle a title for the app
#' @param organism_type human or mouse
#' @param futureMb the megabytes available for the future package
#' @param db_name  path to sqlite database listing bigwig files for cells 
#'
#' @return a dockerized shiny app
#' @export
dockerChevreulApp <- function(
        object = NULL,
        appTitle = NULL,
        organism_type = "human",
        futureMb = 13000,
        db_name = "single-cell-projects.db") {

  db_path <- file.path(rappdirs::user_cache_dir(appname="chevreul"), db_name)

    plan(strategy = "multicore", workers = 6)
    future_size <- futureMb * 1024^2
    options(future.globals.maxSize = future_size)
    options(shiny.maxRequestSize = 40 * 1024^2)
    options(DT.options = list(
        pageLength = 2000, paging = FALSE,
        info = TRUE, searching = TRUE, autoWidth = FALSE, ordering = TRUE,
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
                tabName = "comparePlots", icon = icon("chart-bar"), 
                selected = TRUE
            ), menuItem("Heatmap/Violin Plots",
                tabName = "violinPlots", icon = icon("sort")
            ), menuItem("Coverage Plots",
                tabName = "coveragePlots", icon = icon("mountain")
            ), menuItem("Differential Expression",
                tabName = "diffex", icon = icon("magnet")
                # ), menuItem("Pathway Enrichment Analysis",
                #   tabName = "pathwayEnrichment", icon = icon("sitemap")
            ), menuItem("Find Markers",
                tabName = "findMarkers", icon = icon("bullhorn")
            ), menuItem("Subset SingleCellExperiment Input",
                tabName = "subsetSingleCellExperiment", icon = icon("filter")
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
                    default_helper(type = "markdown", 
                                   content = "comparePlots"),
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
                tabName = "subsetSingleCellExperiment",
                h2("Subset SingleCellExperiment Input") %>%
                    default_helper(type = "markdown", 
                                   content = "subsetSingleCellExperiment"),
                plotDimRedui("subset"),
                chevreulBox(
                    title = "Subset Settings",
                    checkboxInput("legacySettingsSubset", 
                                  "Use Legacy Settings", value = FALSE),
                    actionButton("subsetAction", 
                                 "subset obect by selected cells"),
                    actionButton("subsetCsv", "subset obect by uploaded csv"),
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
                            "Regress SingleCellExperiment Objects By Genes"
                        ),
                        width = 12
                    ) %>%
                        default_helper(type = "markdown", 
                                       content = "regressFeatures")
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

        observe_helpers(
            help_dir = system.file(
                "helpers", package = "chevreul", 
                lib.loc = "/dataVolume/storage/rpkgs/devel_install/"))
        options(warn = -1)

        object <- reactiveVal(NULL)
        observe({
            object(object)
        })

        organism_type <- reactive({
            "human"
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
            reformatted_object <- callModule(reformatMetadataDR, 
                                             "reformatMetadataDR", 
                                             object, featureType)
            object(reformatted_object())
        })

        reductions <- reactive({
            req(object())
            reducedDimNames(object)
        })

        observe({
            req(object())

            callModule(plotDimRed, "plotdimred1", object, plot_types, 
                       featureType,
                organism_type = organism_type, reductions
            )
            callModule(plotDimRed, "plotdimred2", object, plot_types, 
                       featureType,
                organism_type = organism_type, reductions
            )
            callModule(plotDimRed, "diffex", object, plot_types, 
                       featureType,
                organism_type = organism_type, reductions
            )
            callModule(plotDimRed, "subset", object, plot_types, 
                       featureType,
                organism_type = organism_type, reductions
            )
            callModule(plotDimRed, "markerScatter", object, plot_types, 
                       featureType,
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
            plotCoverage, "coverageplots", object, plot_types, proj_dir, 
            organism_type, bigwig_db
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

        callModule(chevreulMarkers, "findmarkers", object, plot_types, 
                   featureType)

        callModule(pathwayEnrichment, "pathwayEnrichment", object)

        # plot all transcripts
        observe({
            req(featureType())
            req(object())
            if (query_experiment(object(), "transcript")) {
                callModule(
                    allTranscripts, "alltranscripts1", object, featureType,
                    organism_type
                )

                callModule(plotDimRed, "alltranscripts2", object, plot_types, 
                           featureType,
                    organism_type = organism_type, reductions
                )
            }
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
                    subset_object <- 
                        object()[, colnames(object()) %in% 
                                     subset_selected_cells()]
                    object(subset_object)
                    if (length(unique(object()$batch)) > 1) {
                        message("reintegrating gene expression")
                        reintegrated_object <- 
                            reintegrate_object(object(),
                            resolution = seq(0.2, 1, by = 0.2),
                            organism = metadata(object())[["experiment"]][["organism"]]
                        )
                        object(reintegrated_object)
                    } else {
                        processed_object <- 
                            obect_pipeline(object(), 
                                           resolution = seq(0.2, 1, by = 0.2))
                        object(processed_object)
                    }
                    message("Complete!")
                },
                message = function(m) {
                    html(id = "subsetMessages", html = paste0(
                        "Subsetting SingleCellExperiment Object: ",
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
                        message("reintegrating gene expression")
                        reintegrated_object <- reintegrate_object(object(),
                            resolution = seq(0.2, 1, by = 0.2),
                            organism = metadata(object())[["experiment"]][["organism"]]
                        )
                        object(reintegrated_object)
                    } else {
                        processed_object <- 
                            obect_pipeline(object(), 
                                           resolution = seq(0.2, 1, by = 0.2))
                        object(processed_object)
                    }
                    message("Complete!")
                },
                message = function(m) {
                    html(id = "subsetMessages", html = paste0(
                        "Subsetting SingleCellExperiment Object: ",
                        m$message
                    ), add = FALSE)
                }
            )
        })

        observeEvent(input$regressAction, {
            req(object())
            showModal(modalDialog(
                title = "Regressing out cell cycle effects",
                "This process may take a minute or two!"
            ))
            regressed_object <- regress_cell_cycle(object())
            object(regressed_object)
            removeModal()
        })

        callModule(techInfo, "techInfo", object)
    }
    shinyApp(ui, server, enableBookmarking = "server")
}
