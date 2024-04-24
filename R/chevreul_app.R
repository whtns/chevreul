#' Run Differential Expression
#'
#'
#'
#' @param object SingleCellExperiment object
#' @param cluster1 cluster 1
#' @param cluster2 cluster 2
#' @param resolution resolution
#' @param diffex_scheme scheme for differential expression
#' @param featureType gene or transcript
#' @param tests t, wilcox, or bimod
#'
#'
#' @return a dataframe with differential expression information
#' @export
#' @examples
#' chevreul_sce <- chevreuldata::human_gene_transcript_sce()
#' run_object_de(chevreul_sce,
#'     diffex_scheme = "louvain",
#'     cluster1 = 1, cluster2 = 2, tests = "t"
#' )
#'
#' cells1 <- colnames(chevreul_sce)[chevreul_sce$batch == "Zhong"]
#' cells2 <- colnames(chevreul_sce)[chevreul_sce$batch == "Kuwahara"]
#'
#' run_object_de(chevreul_sce,
#'     diffex_scheme = "custom",
#'     cluster1 = cells1, cluster2 = cells2, tests = "t"
#' )
#'
run_object_de <- function(object, cluster1, cluster2, resolution = 0.2,
                          diffex_scheme = "louvain", featureType = "gene",
                          tests = c("t", "wilcox", "bimod")) {

  data("grch38")
  data("grch38_tx2gene")

    match.arg(tests)

    if (featureType == "transcript") object <- altExp(object, "transcript")

    if (diffex_scheme == "louvain") {
        if (query_experiment(object, "integrated")) {
            active_experiment <- "integrated"
        } else {
            active_experiment <- "gene"
        }
        colLabels(object) <- colData(object)[[paste0(active_experiment,
                                                     "_snn_res.", resolution)]]
        object <- object[, colLabels(object) %in% c(cluster1, cluster2)]
        colLabels(object) <- factor(colLabels(object))
    } else if (diffex_scheme == "custom") {
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
        message(test)
        de <- findMarkers(object, test.type = test)
        if (featureType == "transcript") {
            de_cols <- c("enstxp", "ensgene", "symbol", "p_val" = "p.value",
                         "avg_log2FC", "pct.1", "pct.2", "p_val_adj" = "FDR")
            de <- de[[1]] %>%
                as.data.frame() %>%
                rownames_to_column("enstxp") %>%
                left_join(grch38_tx2gene, by = "enstxp") %>%
                left_join(grch38, by = "ensgene")
            if ("summary.logFC" %in% colnames(de)) {
                de <- mutate(de, avg_log2FC = log(exp(summary.logFC), 2))
            }
            de <- select(de, any_of(de_cols))
        } else if (featureType == "gene") {
            de_cols <- c("ensgene", "symbol", "p_val" = "p.value", "avg_log2FC",
                         "pct.1", "pct.2", "p_val_adj" = "FDR")
            de <- de[[1]] %>%
                as.data.frame() %>%
                rownames_to_column("symbol") %>%
                left_join(grch38, by = "symbol")
            if ("summary.logFC" %in% colnames(de)) {
                de <- mutate(de, avg_log2FC = log(exp(summary.logFC), 2))
            }
            de <- select(de, any_of(de_cols))
        }
        test_list[[match(test, tests)]] <- de
    }
    names(test_list) <- tests
    return(test_list)
}




#' Prep Slider Values
#'
#' @param default_val Provide a default value
#'
#' @noRd
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
#' @param organism_type human or mouse
#' @param futureMb amount of Mb allocated to future package
#' @param db_name sqlite database with list of saved
#' SingleCellExperiment objects
#' @return a shiny app
#' @export
#'
#' @examples
#' \donttest{
#' chevreulApp()
#' }
#'
chevreulApp <-
  function(preset_project,
           appTitle = "chevreul",
           organism_type = "human",
           futureMb = 13000,
           db_name = "single-cell-projects.db") {

    db_path = file.path(rappdirs::user_cache_dir(appname="chevreul"), db_name)

    message(packageVersion("chevreul"))
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
        uiOutput("projInput"),
        actionButton("loadProject", "Load Selected Project") %>%
            default_helper(type = "markdown", content = "overview"),
        textOutput("appTitle"),
        bookmarkButton(),
        prettyRadioButtons("organism_type",
            inline = TRUE,
            "Organism", choices = c("human", "mouse"), selected = organism_type
        ),
        shinyFilesButton("objectUpload", "Load a SingleCellExperiment Dataset",
            "Please select a .rds file",
            multiple = FALSE
        ),
        shinySaveButton("saveSCE", "Save Current Dataset",
            "Save file as...",
            filetype = list(rds = "rds")
        ),
        verbatimTextOutput("savefile"),
        sidebarMenu(
            menuItem("Integrate Projects",
                tabName = "integrateProjects", icon = icon("object-group")
            ), menuItem("Reformat Metadata",
                tabName = "reformatMetadataDR", icon = icon("columns")
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
            ), menuItem("RNA Velocity",
                tabName = "plotVelocity", icon = icon("tachometer-alt")
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
                tabName = "integrateProjects",
                fluidRow(
                    integrateProjui("integrateproj"),
                    width = 12
                )
            ),
            tabItem(
                tabName = "reformatMetadataDR",
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
                    checkboxInput("legacySettingsSubset", "Use Legacy Settings",
                                  value = FALSE),
                    actionButton("subsetAction",
                                 "subset object by selected cells"),
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
                tabName = "plotVelocity",
                h2("RNA Velocity"),
                fluidRow(
                    plotVelocityui("plotvelocity"),
                )
            ),
            tabItem(
                tabName = "diffex",
                h2("Differential Expression") %>%
                    default_helper(type = "markdown", content = "diffex"),
                plotDimRedui("diffex"),
                diffexui("diffex")
            ),
            # tabItem(
            #   tabName = "pathwayEnrichment",
            #   h2("Pathway Enrichment"),
            #   fluidRow(
            #     pathwayEnrichmentui("pathwayEnrichment")
            #   )
            # ),
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
        # runcodeServer()
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
            state$values$uploadSCEPath <- uploadSCEPath()
        })

        onRestore(function(state) {
            uploadSCEPath(state$values$uploadSCEPath[[1]])
        })

        w <- Waiter$new()

        observe_helpers(help_dir = system.file("helpers", package = "chevreul"))
        options(warn = -1)

        con <- reactive({
            dbConnect(
                SQLite(),
                db_path
            )
        })

        projList <- reactivePoll(4000, session, checkFunc = function() {
            if (file.exists(db_path)) {
                dbReadTable(con(), "projects_tbl") %>%
                    deframe()
            } else {
                ""
            }
        }, valueFunc = function() {
            dbReadTable(con(), "projects_tbl") %>%
                deframe()
        })

        output$projInput <- renderUI({
            selectizeInput("setProject", "Select Project to Load",
                choices = projList(), selected = preset_project,
                multiple = FALSE
            )
        })

        proj_matrices <- reactive({
            create_proj_matrix(projList())
        })

        object <- reactiveVal()
        proj_dir <- reactiveVal()
        uploadSCEPath <- reactiveVal()

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
            paste0("Loaded Project: ", path_file(proj_dir()))
        })
        dataset_volumes <- reactive({
            message(proj_dir())
            dataset_volumes <- c(
                Home = path(
                    proj_dir(),
                    "output", "singlecellexperiment"
                ), "R Installation" = R.home(),
                getVolumes()()
            )
        })
        observe({
            req(dataset_volumes())
            shinyFileChoose(input, "objectUpload",
                roots = dataset_volumes(), session = session
            )
        })
        observeEvent(input$objectUpload, {
            req(dataset_volumes())

            file <- parseFilePaths(
                dataset_volumes(),
                input$objectUpload
            )
            message(file)
            uploadSCEPath(file$datapath)
        })

        observe({
            req(uploadSCEPath())
            message("uploaded")
            withProgress(
                message = paste0("Uploading Data"),
                value = 0,
                {
                    incProgress(2 / 10)
                    message(uploadSCEPath())
                    updated_object <- readRDS(uploadSCEPath())
                    object(updated_object)
                    incProgress(6 / 10)

                    organism <- case_when(
                        str_detect(uploadSCEPath(), "Hs") ~ "human",
                        str_detect(uploadSCEPath(), "Mm") ~ "mouse"
                    )

                    message(uploadSCEPath())
                    message(names(object))
                    incProgress(8 / 10)
                }
            )
        })


        featureType <- reactive({
            "gene"
        })

        observe({
            req(dataset_volumes())
            shinyFileSave(input, "saveSCE",
                # roots = dataset_volumes(),
                roots = c(Home = path(
                    proj_dir(),
                    "output", "singlecellexperiment"
                )),
                session = session, restrictions = system.file(package = "base")
            )
        })
        subSingleCellExperimentPath <- eventReactive(input$saveSCE, {
            req(object())
            savefile <- parseSavePath(
                dataset_volumes(),
                input$saveSCE
            )
            return(savefile$datapath)
        })
        observeEvent(input$saveSCE, {
            req(object())
            req(subSingleCellExperimentPath())

            withProgress(
                message = paste0("Saving Data"),
                value = 0,
                {
                    Sys.sleep(6)
                    incProgress(2 / 10)
                    saveRDS(
                        object(),
                        subSingleCellExperimentPath()
                    )
                    incProgress(10 / 10)
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
            newintegrated_project <- set_names(
                integration_path,
                path_file(integration_path)
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
            reformatted_object <- callModule(reformatMetadataDR,
                                             "reformatMetadataDR", object,
                                             featureType)
            object(reformatted_object())
        })

        reductions <- reactive({
            req(object())
            reducedDimNames(object())
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
            callModule(plotDimRed, "diffex", object, plot_types, featureType,
                organism_type = organism_type, reductions
            )
            callModule(plotDimRed, "subset", object, plot_types, featureType,
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
            organism_type
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
                    html("subsetMessages", "")
                    message("Beginning")

                    subset_object <-
                      object()[, colnames(object()) %in% subset_selected_cells()]
                    object(subset_object)
                    if (length(unique(object()[["batch"]])) > 1) {
                        message("reintegrating gene expression")
                        reintegrated_object <- reintegrate_object(object(),
                            resolution = seq(0.2, 2, by = 0.2),
                            legacy_settings = input$legacySettingsSubset,
                            organism = metadata(object())$experiment$organism
                        )
                        object(reintegrated_object)
                    } else {
                        subset_object <- object_pipeline(
                          object(),
                          resolution = seq(0.2, 2, by = 0.2),
                          legacy_settings = input$legacySettingsSubset)
                        object(subset_object)
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
                            resolution = seq(0.2, 2, by = 0.2),
                            legacy_settings = input$legacySettingsSubset,
                            organism = metadata(object())$experiment$organism
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
                    html(id = "subsetMessages", html = paste0(
                        "Subsetting SingleCellExperiment Object: ",
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
        callModule(chevreulMarkers, "findmarkers", object, plot_types,
                   featureType)
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

            callModule(plotDimRed, "alltranscripts2", object, plot_types,
                       featureType,
                organism_type = organism_type, reductions
            )
        })

        observeEvent(input$regressAction, {
            req(object())
            showModal(modalDialog(
                title = "Regressing cycle effects",
                "This process may take a minute or two!"
            ))
            regressed_object <- regress_cell_cycle(object())
            object(regressed_object)
            removeModal()
        })

        observe({
            req(uploadSCEPath())
            req(object())

            proj_path <- str_replace(uploadSCEPath(), "output.*", "")

            proj_name <- path_file(proj_path)
            message(proj_name)

            loom_path <- path(proj_path, "output", "velocyto", paste0(proj_name,
                                                                      ".loom"))
            message(loom_path)
            # need to check if this file exists

            callModule(plotVelocity, "plotvelocity", object, loom_path)
        })

        callModule(techInfo, "techInfo", object)

        sessionId <- as.integer(runif(1, 1, 100000))
        output$sessionId <- renderText(paste0("Session id: ", sessionId))
        session$onSessionEnded(function() {
            cat(paste0("Ended: ", sessionId))
            observe(dbConnect(con()))
        })
    }

    # onStop(function() dbDisconnect(con))
    shinyApp(ui, server, enableBookmarking = "url")
}
