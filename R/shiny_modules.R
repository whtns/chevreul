plotClustree_UI <- function(id) {
    ns <- NS(id)
    tagList(
        chevreulBox(
            title = "Clustering Tree",
            plotOutput(ns("clustree"), height = "700px")
        )
    )
}

#' plot clustree server
#'
#' @noRd
plotClustree <- function(input, output, session, object) {
    # set appropriate experiment
    experiment <- reactive({
        ifelse(query_experiment(object(), "integrated"), "integrated", "gene")
    })

    output$clustree <- renderPlot({
        req(object())
        experiment <- ifelse("integrated" %in% get_feature_types(object()), 
                             "integrated", "gene")
        object <- set_feature_type(object(), experiment)
        clustree::clustree(object(), prefix = paste0(experiment, "_snn_res."))
    })
}


#' Title
#'
#' @param id
#'
#' @noRd
plotViolinui <- function(id) {
    ns <- NS(id)
    tagList(
        chevreulBox(
            title = "Violin Plots",
            uiOutput(ns("vln_group")),
            selectizeInput(ns("customFeature"),
                "Gene or transcript expression by which to color the plot 
                eg. 'NRL' or 'ENST00000488147'",
                choices = NULL, multiple = FALSE
            ),
            radioButtons(ns("slot"), "Data Type", choices = 
                             c("transformed" = "data", 
                               "raw counts" = "counts")),
            downloadButton(ns("downloadPlot")),
            plotlyOutput(ns("vplot"), height = 750),
            width = 11
        ) %>%
            default_helper(type = "markdown", content = "violinPlot", 
                           size = "l")
    )
}

#' Plot Violin Server
#'
#' Plots a Violin plot of a single data (gene expression, metrics, etc.) 
#' in the server SingleCellExperiment app.
#'
#' @param object SingleCellExperiment object
#' @param featureType Gene or Transcript
#' @param organism_type Organism
#'
#' @noRd
plotViolin <- function(input, output, session, object, featureType, 
                       organism_type) {
    ns <- session$ns
    prefill_feature <- reactive({
        req(featureType())
        if (featureType() == "transcript") {
            if (organism_type() == "human") {
                "ENST00000488147"
            } else if (organism_type() == "mouse") {
                "ENSG00000488147"
            }
        } else if (featureType() == "gene") {
            if (organism_type() == "human") {
                "NRL"
            } else if (organism_type() == "mouse") {
                "NRL"
            }
        }
    })
    observe({
        req(prefill_feature())
        req(object())
        updateSelectizeInput(session, "customFeature",
            choices = rownames(object()),
            selected = prefill_feature(), server = TRUE
        )
    })

    output$vln_group <- renderUI({
        req(object())
        selectizeInput(ns("vlnGroup"), "Grouping variable",
            choices = metadata_from_object(object()), selected = "batch"
        )
    })

    vln_plot <- reactive({
        req(input$customFeature)
        req(input$vlnGroup)

        vln_plot <-
            plot_violin(object(), plot_var = input$vlnGroup, 
                        features = input$customFeature, slot = input$slot)
    })

    output$downloadPlot <- downloadHandler(
        filename = function() {
            paste("violin", ".pdf", sep = "")
        },
        content = function(file) {
            ggsave(file, vln_plot() + theme_pubr(base_size = 20, 
                                                 x.text.angle = 45), 
                   width = 16, height = 12)
        }
    )

    output$vplot <- renderPlotly({
        req(object())
        req(input$vlnGroup)
        exclude_trace_number <- 
            length(unique(get_cell_metadata(object())[[input$vlnGroup]])) * 2

        vln_plot <- ggplotly(vln_plot(), height = 700) %>%
            style(opacity = 0.5) %>%
            style(hoverinfo = "skip", 
                  traces = c(seq_len(exclude_trace_number))) %>%
            plotly_settings(width = 1200) %>%
            toWebGL() %>%
            identity()
    })
}


#' Plot Heatmap ui
#'
#' @param id
#'
#' @noRd
plotHeatmapui <- function(id) {
    ns <- NS(id)
    tagList(
        chevreulBox(
            title = "Heatmap",
            uiOutput(ns("colAnnoVarui")),
            radioButtons(ns("assayName"), "Data Scaling", 
                         choices = c(scaled = "scaledata", 
                                     unscaled = "logcounts"), 
                         selected = "scaledata", inline = TRUE),
            selectizeInput(ns("dendroSelect"), 
                           "Clustering algorithm or metadata for column 
                           arrangement", choices = NULL, selected = NULL, 
                           multiple = TRUE),
            actionButton(ns("actionHeatmap"), "Plot Heatmap"),
            downloadButton(ns("downloadPlot"), "Download Heatmap"),
            selectizeInput(ns("customFeature"), "Gene or transcript 
                           expression by which to color the plot; eg. 
                           'NRL' or 'ENST00000488147'",
                choices = NULL, multiple = TRUE
            ),
            plotOutput(ns("heatmap"), height = 750),
            width = 12
        ) %>%
            default_helper(type = "markdown", content = "heatMap")
    )
}

#' Plot Heatmap
#'
#' @param object a SingleCellExperiment object
#' @param featureType gene or transcript
#' @param organism_type human or mouse
#'
#' @noRd
#'
plotHeatmap <- function(input, output, session, object, featureType, 
                        organism_type) {
    ns <- session$ns

    w <- Waiter$new(ns("heatmap"),
        html = spin_loaders(id = 1, color = "black", 
                            style = "position:relative;margin:auto;"),
        color = transparent(.5)
    )

    observe({
        req(object())
        if (query_experiment(object(), "integrated")) {
            experiment <- "integrated"
        } else {
            experiment <- "gene"
        }

        preset_features <- 
            get_variable_features(object(), experiment = experiment)[seq(50)]

        updateSelectizeInput(session, "customFeature",
            choices = get_features(object(), experiment = experiment),
            selected = preset_features, server = TRUE
        )
    })

    output$colAnnoVarui <- renderUI({
        req(object())

        formatted_col_names <- colnames(get_cell_metadata(object())) %>%
            make_chevreul_clean_names()

        selectizeInput(ns("colAnnoVar"), "Column Annotation(s)",
            choices = formatted_col_names, selected = "batch", multiple = TRUE
        )
    })

    observe({
        req(object())

        hclust_methods <- c("Ward" = "ward.D2", "single", 
                            "complete", "average")

        updateSelectizeInput(
            session, 
            "dendroSelect", 
            choices = c(hclust_methods, 
                        colnames(get_cell_metadata(object()))), 
            selected = "ward.D2")
    })

    heatmap_plot <- eventReactive(input$actionHeatmap, {
        req(input$customFeature)
        req(input$colAnnoVar)

        if (query_experiment(object(), "integrated")) {
            experiment <- "integrated"
        } else {
            experiment <- "gene"
        }

        hm <- make_complex_heatmap(object(), 
                                   features = input$customFeature, 
                                   experiment = experiment, 
                                   group.by = input$colAnnoVar, 
                                   assayName = input$assayName, 
                                   col_arrangement = input$dendroSelect)

        hm <- draw(hm)
        return(hm)
    })

    output$heatmap <- renderPlot({
        w$show()
        heatmap_plot()
    })

    output$downloadPlot <- downloadHandler(
        filename = function() {
            paste("heatmap", ".pdf", sep = "")
        },
        content = function(file) {
            ggsave(file, as.ggplot(heatmap_plot()) + 
                       theme_pubr(base_size = 20, x.text.angle = 45), 
                   width = 16, height = 12)
        }
    )
}

#' Integrate Project UI
#'
#'
#' @noRd
integrateProjui <- function(id) {
    ns <- NS(id)
    tagList(
        chevreulBox(
            title = "Integrate Projects",
            shinyFilesButton(ns('integratedProjDir'), 
                             label='Integrated Project Directory', 
                             title='Choose a directory to store integrated datasets', multiple=FALSE),
            actionButton(ns("integrateAction"), "Integrate Selected Projects"),
            # useShinyjs(),
            # runcodeUI(code = "alert('Hello!')", id = "subsetcode"),
            textOutput(ns("integrationMessages")),
            checkboxInput(ns("legacySettings"), "Use Legacy Settings", value = FALSE),
            textOutput(ns("integrationResult")),
            shinySaveButton(ns("saveIntegratedProject"), "Save Integrated Project", "Save project as..."),
            DTOutput(ns("myDatatable")),
            width = 12
        ) %>%
            default_helper(type = "markdown", content = "integrateProjects")
    )
}

#' Integrate Projects Server Function
#'
#' @param proj_matrices project matrices
#' @param object a SingleCellExperiment object
#' @param proj_dir project directory
#' @param con a connection
#'
#' @noRd
integrateProj <- function(input, output, session, proj_matrices, 
                          object, proj_dir, con) {
    ns <- session$ns

    dataset_volumes <- reactive({
      dataset_volumes <- c(
        Home =
          dirname(proj_dir()),
        "R Installation" = R.home(),
        getVolumes()()
      )
    })

    observe({
      req(dataset_volumes())
      shinyFileChoose(input, "integratedProjDir",
                      roots = dataset_volumes(), session = session
      )
    })

    proj_matrix <- reactive({
        proj_matrices()$primary_projects
    })

    clean_proj_matrix <- reactive({
        clean_proj_matrix <- proj_matrix() %>%
            select(-project_path) %>%
            identity()
    })

    output$myDatatable <- renderDT(clean_proj_matrix(),
        server = FALSE,
        rownames = TRUE
    )

    selectedRows <- eventReactive(input$integrateAction, {
        ids <- input$myDatatable_rows_selected
    })

    selectedProjects <- reactive({
        selectedProjects <- slice(proj_matrix(), selectedRows()) %>%
            pull(project_path) %>%
            identity()
    })

    mergedObjects <- reactiveVal()

    observeEvent(input$integrateAction, {
        req(selectedProjects())
        withCallingHandlers(
            {
                html("integrationMessages", "")
                message("Beginning")
                message(selectedProjects())
                batches <- path(selectedProjects(), "output", 
                                "singlecellexperiment", 
                                "unfiltered_object.rds") %>%
                    map(readRDS)

                names(batches) <- names(selectedProjects())
                message(names(batches))
                mergedObjects(
                    integration_workflow(
                        batches, legacy_settings = input$legacySettings))
                # mergedObjects(batches[[1]])

                message("Integration Complete!")
            },
            message = function(m) {
                html(id = "integrationMessages", 
                     html = paste0("Running Integration: ", m$message), 
                     add = FALSE)
            }
        )
    })

    newProjDir <- reactive({
        req(mergedObjects())

        newProjName <- paste0(map(path_file(selectedProjects()), ~ 
                                      gsub("_proj", "", .x)), collapse = "_")
        newProjDir <- path(input$integratedProjDir, newProjName)

        proj_dir(newProjDir)

        newProjDir
    })


    volumes <- reactive({
        volumes <- c(
            Home = input$integratedProjDir,
            "R Installation" = R.home(),
            getVolumes()()
        )
        # print(volumes)
        volumes
    })

    observe({
        shinyFileSave(input,
            "saveIntegratedProject",
            roots = volumes(),
            session = session,
            restrictions = system.file(package = "base")
        )
    })


    integratedProjectSavePath <- eventReactive(input$saveIntegratedProject, {
        savefile <- parseSavePath(volumes(), input$saveIntegratedProject)

        savefile$datapath
    })

    output$integrationResult <- renderText({
        integratedProjectSavePath()
    })

    observeEvent(input$saveIntegratedProject, {
        req(mergedObjects())
        req(integratedProjectSavePath())

        if (!is.null(integratedProjectSavePath())) {
            withProgress(
                message = paste0("Saving Integrated Dataset to ", 
                                 integratedProjectSavePath()),
                value = 0,
                {
                    # Sys.sleep(6)
                    incProgress(2 / 10)
                    save_object(mergedObjects(), 
                                proj_dir = integratedProjectSavePath())
                    writeLines(character(), path(integratedProjectSavePath(), 
                                                 ".here"))
                    # create_proj_db()
                    dbAppendTable(con, "projects_tbl", data.frame(
                        project_name = path_file(integratedProjectSavePath()),
                        project_path = integratedProjectSavePath(),
                        project_slug = str_remove(
                            path_file(integratedProjectSavePath()), "_proj$"),
                        project_type = "integrated_projects"
                    ))
                    incProgress(8 / 10)
                }
            )
        }
    })


    return(integratedProjectSavePath)
}


#' Change Embedding Parameters UI
#'
#' @noRd
changeEmbedParamsui <- function(id) {
    ns <- NS(id)

    minDist_vals <- prep_slider_values(0.3)
    negsamprate_vals <- prep_slider_values(5)

    tagList(
        selectizeInput(ns("dims"), label = "Dimensions from PCA", 
                       choices = seq(1, 99), multiple = TRUE, 
                       selected = seq(30)),
        sliderInput(ns("minDist"), label = "Minimum Distance", 
                    min = minDist_vals$min, max = minDist_vals$max, 
                    value = minDist_vals$value, step = minDist_vals$step),
        sliderInput(ns("negativeSampleRate"), label = "Negative Sample Rate",
                    min = negsamprate_vals$min, max = negsamprate_vals$max, 
                    value = negsamprate_vals$value, 
                    step = negsamprate_vals$step)
    )
}

#' Change Embedding Parameters
#'
#' @param object a SingleCellExperiment object
#'
#' @noRd
changeEmbedParams <- function(input, output, session, object) {
    ns <- session$ns

    object <- RunUMAP(object(), dims = as.numeric(input$dims), 
                      reduction = "PCA", min.dist = input$minDist, 
                      negative.sample.rate = input$negativeSampleRate)

    return(object)
}

#' Plot Dimensional Reduduction UI
#'
#' @noRd
plotDimRedui <- function(id) {
    ns <- NS(id)
    chevreulBox(
        title = "Embedding",
        chevreulDropDownButton(
            ns("dimPlotSettings"),
            selectizeInput(ns("embedding"), "Embedding", choices = NULL, 
                           selected = NULL),
            sliderInput(ns("dotSize"), "Size of Points in UMAP", min = 0.5, 
                        max = 2, step = 0.1, value = 1),
            selectizeInput(ns("dim1"), "Dimension 1", choices = seq(1, 99), 
                           selected = 1),
            selectizeInput(ns("dim2"), "Dimension 2", choices = seq(1, 99), 
                           selected = 2)
        ),
        selectizeInput(ns("plottype"), "Variable to Plot", choices = NULL, 
                       multiple = TRUE),
        selectizeInput(ns("customFeature"), 
                       "Gene or transcript expression by which to color 
                       the plot; eg. 'NRL' or 'ENST00000488147'", 
                       choices = NULL, multiple = FALSE),
        plotlyOutput(ns("dplot"), height = 500),
        width = 6
    )
}

#' Plot Dimensional Reduduction
#'
#' @param object a SingleCellExperiment object
#' @param plot_types plot types
#' @param featureType gene or transcript
#' @param organism_type human or mouse
#' @param reductions embeddings
#'
#' @noRd
plotDimRed <- function(
        input, output, session, object, plot_types, featureType,
        organism_type, reductions) {
    ns <- session$ns

    output$myPanel <- renderUI({
        req()
        lev <- sort(unique(input$select)) 
        cols <- gg_fill_hue(length(lev))

        # New IDs "colX1" so that it partly coincide with input$select...
        lapply(seq_along(lev), function(i) {
            colourInput(
                inputId = paste0("col", lev[i]),
                label = paste0("Choose colour for ", lev[i]),
                value = cols[i]
            )
        })
    })

    observe({
        req(object())
        updateSelectizeInput(session, "embedding",
            choices = reductions(),
            selected = rev(reductions())[1], server = TRUE
        )
    })

    selected_plot <- reactiveVal()
    observe({
        req(object())
        # selected_plot <- ifelse(is.null(selected_plot()), "louvain",
        #                         selected_plot())
        updateSelectizeInput(session, "plottype",
            choices = flatten_chr(plot_types()),
            selected = "batch"
        )
    })
    prefill_feature <- reactive({
        req(featureType())
        if (featureType() == "transcript") {
            if (organism_type() == "human") {
                "ENST00000488147"
            } else if (organism_type() == "mouse") {
                "ENSG00000488147"
            }
        } else if (featureType() == "gene") {
            if (organism_type() == "human") {
                "NRL"
            } else if (organism_type() == "mouse") {
                "NRL"
            }
        }
    })
    observe({
        req(prefill_feature())
        req(featureType())
        req(object())
        updateSelectizeInput(session, "customFeature",
            choices = get_features(object(), featureType()),
            selected = prefill_feature(), server = TRUE
        )
    })

    output$dplot <- renderPlotly({
        req(input$plottype)
        req(object())
        req(input$embedding)
        if (length(input$plottype) > 1) {
            cross_plot_object <- unite_metadata(object(), input$plottype)

            newcolname <- paste(input$plottype, collapse = "_")
            cross_plot_object[[newcolname]] <- Idents(cross_plot_object)

            selected_plot(newcolname)

            plot_var(cross_plot_object,
                dims = c(input$dim1, input$dim2),
                embedding = input$embedding, group = NULL, 
                pt.size = input$dotSize,
                return_plotly = TRUE
            )
        } else {
            if (input$plottype == "feature") {
                plot_feature(object(),
                    dims = c(
                        input$dim1,
                        input$dim2
                    ), embedding = input$embedding,
                    features = input$customFeature, pt.size = input$dotSize,
                    return_plotly = TRUE
                )
            } else if (input$plottype %in% plot_types()$continuous_vars) {
                plot_feature(object(),
                    dims = c(
                        input$dim1,
                        input$dim2
                    ), embedding = input$embedding,
                    features = input$plottype, pt.size = input$dotSize,
                    return_plotly = TRUE
                )
            } else if (input$plottype %in% plot_types()$category_vars) {
                plot_var(object(),
                    dims = c(input$dim1, input$dim2),
                    embedding = input$embedding, group = input$plottype, 
                    pt.size = input$dotSize,
                    return_plotly = TRUE
                )
            }
        }
    })
}


#' Create Table of Selected Cells UI
#'
#' @param id
#'
#' @noRd
tableSelectedui <- function(id) {
    ns <- NS(id)
    tagList(DTOutput(ns("brushtable")))
}

#' Create Table of Selected Cells
#'
#' @param object a SingleCellExperiment object
#'
#' @noRd
tableSelected <- function(input, output, session, object) {
    ns <- session$ns
    brush <- reactive({
        req(object())
        d <- event_data("plotly_selected")
        if (is.null(d)) {
            msg <- 
                "Click and drag events (i.e. select/lasso) appear here 
            (double-click to clear)"
            return(d)
        } else {
            # selected_cells <- colnames(object())[as.numeric(d$key)]
            d$key
        }
    })

    output$brushtable <- renderDT({
        req(object())
        req(brush())
        selected_meta <- data.frame(get_cell_metadata(object())[brush(), ])
        datatable(selected_meta,
            extensions = "Buttons",
            selection = list(mode = "multiple", 
                             selected = seq_len(nrow(selected_meta)), 
                             target = "row"),
            options = list(dom = "Bft", buttons = c("copy", "csv"), 
                           scrollX = "100px", scrollY = "800px")
        )
    })

    selected_cells <- reactive({
        selected_rows <- input$brushtable_rows_selected
        rownames(get_cell_metadata(object())[brush(), ])[selected_rows]
    })

    return(selected_cells)
}


#' Differential Expression UI
#'
#' @noRd
diffexui <- function(id) {
    ns <- NS(id)
    tagList(
        chevreulBox(
            title = "Differential Expression Settings",
            radioButtons(ns("diffex_scheme"),
                "Cells to Compare",
                choiceNames = 
                    c("SingleCellExperiment Cluster", "Custom Selection"), 
                choiceValues = c("louvain", "custom"),
                selected = "louvain",
                inline = TRUE
            ),
            conditionalPanel(
                ns = ns,
                condition = "input.diffex_scheme == 'louvain'",
                sliderInput(ns("objectResolution"), 
                            "Resolution of clustering algorithm (affects number of clusters)",
                    min = 0.2, max = 2, step = 0.2, value = 0.6
                ),
                numericInput(ns("cluster1"),
                    "first cluster to compare",
                    value = 0
                ), numericInput(ns("cluster2"),
                    "second cluster to compare",
                    value = 1
                )
            ),
            conditionalPanel(
                ns = ns,
                condition = "input.diffex_scheme == 'custom'",
                sliderInput(ns("customResolution"), "Resolution of clustering 
                            algorithm (affects number of clusters)",
                    min = 0.2, max = 2, step = 0.2, value = 0.6
                ),
                actionButton(
                    ns("saveClust1"),
                    "Save to Custom Cluster 1"
                ), actionButton(
                    ns("saveClust2"),
                    "Save to Custom Cluster 2"
                )
            ),
            uiOutput(ns("testChoices")),
            radioButtons(ns("featureType"), "Features to Compare", 
                         choices = c("gene", "transcript")),
            actionButton(
                ns("diffex"),
                "Run Differential Expression"
            ),
            downloadLink(ns("downloadData"), "Download Complete DE Results"),
            DTOutput(ns("DT1")),
            width = 6
        ),
        chevreulBox(
            title = "Volcano Plot",
            sliderInput(ns("FCcutoff"), "FC cutoff value (log2 fold change)",
                min = 0, max = 10, step = 0.5, value = 1
            ),
            sliderInput(ns("pCutoff"), "-log10 p adj value",
                min = 0, max = 5, step = 0.5, value = 1.5
            ),
            plotOutput(ns("volcano")),
            downloadButton(ns("downloadVolcanoPlot"), "Download Volcano Plot"),
            width = 6
        ),
        # chevreulBox(
        #   title = "Cells",
        #   tabsetPanel(type = "tabs",
        #               tabPanel("Selected Cells", tableSelectedui("diffex")),
        #               tabPanel("Custom Cluster 1", DTOutput(ns("cc1"))),
        #               tabPanel("Custom Cluster 2", DTOutput(ns("cc2")))
        #   ),
        #   width = 6
        # ),

        chevreulBox(
            title = "Selected Cells",
            tableSelectedui("diffex"),
            width = 12
        ),
        chevreulBox(
            title = "Custom Cluster 1", DTOutput(ns("cc1")),
            width = 6
        ), chevreulBox(
            title = "Custom Cluster 2", DTOutput(ns("cc2")),
            width = 6
        ),
    )
}


#' Title
#'
#' @noRd
cells_selected <- function(input) {
    if (identical(input, character(0))) {
        "Please selected desired cells by clicking on the table"
    } else {
        NULL
    }
}

#' Differential Expression
#'
#' @param object a SingleCellExperiment object
#' @param featureType gene or transcript
#' @param selected_cells selected cells
#' @param tests tests to use
#'
#' @noRd
diffex <- 
    function(input, output, session, object, featureType, selected_cells, 
             tests = c("t-test" = "t", "wilcoxon rank-sum test" = "wilcox", 
                       "Likelihood-ratio test (bimodal)" = "bimod", 
                       "MAST" = "MAST")) {
    ns <- session$ns

    experiment <- reactive({
        req(object())
        if (query_experiment(object(), "integrated")) {
            experiment <- "integrated"
        } else {
            experiment <- "gene"
        }
    })

    output$testChoices <- renderUI(
        selectizeInput(ns("diffex_method"),
            "Method of Differential Expression",
            choices = tests,
            selected = "t"
        )
    )

    brush <- reactive({
        req(object())
        d <- event_data("plotly_selected")
        if (is.null(d)) {
            msg <- "Click and drag events (i.e. select/lasso) appear here 
            (double-click to clear)"
            return(d)
        } else {
            # selected_cells <- colnames(object())[as.numeric(d$key)]
            d$key
        }
    })
    custom_cluster1 <- eventReactive(input$saveClust1, {
        validate(
            cells_selected(selected_cells())
        )
        isolate(selected_cells())
    })
    custom_cluster2 <- eventReactive(input$saveClust2, {
        validate(
            cells_selected(selected_cells())
        )
        isolate(selected_cells())
    })

    output$cc1 <- renderDT({
        req(custom_cluster1())
        selected_meta <- 
            data.frame(get_cell_metadata(object())[custom_cluster1(), ])
        datatable(selected_meta,
            extensions = "Buttons",
            options = list(dom = "Bft", buttons = c(
                "copy",
                "csv"
            ), scrollX = "100px", scrollY = "400px")
        )
    })
    output$cc2 <- renderDT({
        req(custom_cluster2())
        selected_meta <- 
            data.frame(get_cell_metadata(object())[custom_cluster2(), ])
        datatable(selected_meta,
            extensions = "Buttons",
            options = list(dom = "Bft", buttons = c(
                "copy",
                "csv"
            ), scrollX = "100px", scrollY = "400px")
        )
    })

    de_results <- eventReactive(input$diffex, {
        if (input$diffex_scheme == "louvain") {
            run_object_de(object(), input$cluster1, input$cluster2,
                resolution = input$objectResolution, 
                diffex_scheme = input$diffex_scheme, 
                input$featureType, tests = input$diffex_method
            )
        } else if (input$diffex_scheme == "custom") {
            # req(custom_cluster1())
            # req(custom_cluster2())
            cluster1 <- unlist(strsplit(
                custom_cluster1(),
                " "
            ))
            cluster2 <- unlist(strsplit(
                custom_cluster2(),
                " "
            ))
            run_object_de(object(), cluster1, cluster2,
                input$customResolution,
                diffex_scheme = input$diffex_scheme, input$featureType, 
                tests = input$diffex_method
            )
        }
    })

    output$DT1 <- renderDT(de_results()[[input$diffex_method]],
        extensions = "Buttons", options = list(
            dom = "Bfptr",
            buttons = c("copy", "csv"), scrollX = "100px", scrollY = "600px"
        ), class = "display"
    )


    Volcano <- reactive({
        de_results()[[input$diffex_method]] %>%
            distinct(symbol, .keep_all = TRUE) %>%
            column_to_rownames("symbol") %>%
            EnhancedVolcano(
                lab = rownames(.),
                x = "avg_log2FC",
                y = "p_val_adj",
                pCutoff = 1 / (10^as.numeric(input$pValCutoff)),
                FCcutoff = as.numeric(input$FCcutoff)
            )
    })

    output$volcano <- renderPlot({
        print(Volcano())
    })

    output$downloadVolcanoPlot <- downloadHandler(
        filename = function() {
            paste("DE_Volcano_plot", ".pdf", sep = "")
        },
        content = function(file) {
            ggsave(file, Volcano() + 
                       theme_pubr(base_size = 20, x.text.angle = 45), 
                   width = 16, height = 12)
        }
    )

    cluster_list <- reactive({
        if (input$diffex_scheme == "louvain") {
            object_meta <- 
                get_cell_metadata(object())[[paste0(DefaultAssay(object()), 
                                                    "_snn_res.", 
                                                    input$objectResolution)]]
            cluster1_cells <- rownames(
                object_meta[object_meta == input$cluster1, , drop = FALSE])
            cluster2_cells <- rownames(
                object_meta[object_meta == input$cluster2, , drop = FALSE])
            list(cluster1 = cluster1_cells, cluster2 = cluster2_cells)
        } else if (input$diffex_scheme == "feature") {
            list(cluster1 = custom_cluster1(), cluster2 = custom_cluster2())
        }
    })

    return(list(cluster_list = cluster_list, de_results = de_results))
}


#' Find Markers UI
#'
#' @noRd
chevreulMarkersui <- function(id) {
    ns <- NS(id)
    tagList(
        chevreulBox(
            title = "Find Markers",
            uiOutput(ns("dplottype")),
            sliderInput(ns("resolution2"), 
                        label = "Resolution of clustering algorithm (affects number of clusters)", 
                        min = 0.2, max = 2, step = 0.2, value = 0.6),
            numericInput(ns("num_markers"), 
                         "Select Number of Markers to Plot for Each Value", 
                         value = 5, min = 2, max = 20),
            uiOutput(ns("valueSelect")),
            radioButtons(ns("markerMethod"), "Method of Marker Selection", 
                         choices = c("wilcox"), selected = "wilcox", 
                         inline = TRUE),
            sliderInput(ns("pValCutoff"), "P Value cutoff", min = 0.01, 
                        max = 1, value = 1),
            selectizeInput(ns("dotFeature"), "Feature for Marker Plot", 
                           choices = NULL),
            actionButton(ns("plotDots"), "Plot Markers!"),
            downloadButton(ns("downloadMarkerTable"), "Download Markers!"),
            checkboxInput(ns("uniqueMarkers"), "Make Markers Unique", 
                          value = FALSE),
            checkboxInput(ns("hidePseudo"), "Hide Pseudogenes", value = TRUE),
            plotlyOutput(ns("markerplot"), height = 800),
            width = 6
        )
    ) %>%
        default_helper(type = "markdown", content = "findMarkers")
}

#' Find Markers
#'
#' @param object a SingleCellExperiment object
#'
#' @noRd
chevreulMarkers <- function(input, output, session, object, plot_types, 
                            featureType) {
    ns <- session$ns

    observe({
        req(object())
        updateSelectizeInput(session, "dotFeature", 
                             choices = c(mainExpName(object()), 
                                         altExpNames(object())), 
                             selected = "gene", server = TRUE)
    })

    output$dplottype <- renderUI({
        req(object())
        # selected_plot <- ifelse(is.null(selected_plot()), "louvain",
        #                         selected_plot())
        selectizeInput(ns("plottype"), "Variable to Plot",
            choices = flatten_chr(plot_types()),
            selected = "louvain", multiple = TRUE
        )
    })

    experiment <- reactive({
        req(object())
        if (query_experiment(object(), "integrated")) {
            experiment <- "integrated"
        } else {
            experiment <- "gene"
        }
    })

    group_by <- reactive({
        req(input$plottype)

        if (input$plottype == "louvain") {
            group_by <- paste0(experiment(), "_snn_res.", input$resolution2)
        } else {
            group_by <- input$plottype
        }
    })

    output$valueSelect <- renderUI({
        req(object())
        req(group_by())

        choices <- levels(get_cell_metadata(object())[[group_by()]])

        selectizeInput(ns("displayValues"), "Values to display", 
                       multiple = TRUE, choices = choices)
    })

    marker_plot_return <- eventReactive(input$plotDots, {
        plot_markers(object(), group_by = group_by(), 
                     num_markers = input$num_markers, 
                     selected_values = input$displayValues, 
                     marker_method = input$markerMethod, 
                     experiment = input$dotFeature, 
                     featureType = featureType(), 
                     hide_pseudo = input$hidePseudo, 
                     unique_markers = input$uniqueMarkers, 
                     p_val_cutoff = input$pValCutoff, return_plotly = TRUE)
    })

    output$markerplot <- renderPlotly({
        # req(input$displayClusters)
        marker_plot_return()$plot
    })

    output$downloadMarkerTable <- downloadHandler(
        filename = function() {
            paste(group_by(), "_markers.csv", sep = "")
        },
        content = function(file) {
            write_csv(marker_plot_return()$markers, file)
        }
    )
}

#' Plot Read Count UI
#'
#' @noRd
plotReadCountui <- function(id) {
    ns <- NS(id)
    chevreulBox(
        title = "Histogram (Read Counts, etc.)",
        uiOutput(ns("group_byui")),
        uiOutput(ns("colorbyui")),
        sliderInput(ns("resolution"), 
                    "Resolution of clustering algorithm 
                    (affects number of clusters)",
            min = 0.2, max = 2, step = 0.2, value = 0.6
        ),
        radioButtons(ns("yScale"), "Y axis transforamtion", 
                     choices = c("log", "linear"), selected = "linear", 
                     inline = TRUE),
        plotlyOutput(ns("rcplot"), height = 500),
        collapsed = TRUE
    )
}

#' Plot Read Count
#'
#' @param object a SingleCellExperiment object
#' @param plot_types plot types
#'
#' @noRd
plotReadCount <- function(input, output, session, object, plot_types) {
    ns <- session$ns

    output$colorbyui <- renderUI({
        req(object())
        selectInput(ns("colorby"), "Variable to Color the Plot by",
            choices = flatten_chr(plot_types()), selected = c("louvain"), 
            multiple = FALSE
        )
    })

    output$group_byui <- renderUI({
        req(object())
        selectInput(ns("group_by"), "Variable for x-axis",
            choices = flatten_chr(plot_types()), selected = c("nCount_RNA"), 
            multiple = FALSE
        )
    })

    output$rcplot <- renderPlotly({
        req(object())
        req(input$colorby)

        if (input$colorby == "louvain") {
            if (query_experiment(object(), "integrated")) {
                experiment <- "integrated"
            } else {
                experiment <- "gene"
            }

            louvain_resolution <- paste0(experiment, "_snn_res.", 
                                         input$resolution)
            plot_readcount(object(), group_by = input$group_by, 
                           fill_by = louvain_resolution, 
                           yscale = input$yScale, return_plotly = TRUE)
        } else if (input$colorby %in% flatten_chr(plot_types())) {
            plot_readcount(object(), group_by = input$group_by, 
                           fill_by = input$colorby, yscale = input$yScale, 
                           return_plotly = TRUE)
        }
    })
}

#' Cell Cycle Score UI
#'
#' @noRd
ccScoreui <- function(id) {
    ns <- NS(id)
    tagList()
}

#' Cell Cycle Score
#'
#'
#' @noRd
ccScore <- function(input, output, session) {
    ns <- session$ns
    output$rplot1 <- renderPlot({
        req(object())
        plot_ridge(object(), features = input$feature)
    })
    plotOutput("rplot1", height = 750)
}

#' Plot All Transcripts UI Module
#'
#' @noRd
allTranscriptsui <- function(id) {
    ns <- NS(id)
    tagList(
        default_helper(
            chevreulBox(
                title = "Transcript Expression per Gene",
                selectizeInput(ns("embeddingGene"), "Gene or transcript expression by which to color the plot; eg. 'NRL'", choices = NULL, selected = NULL),
                selectizeInput(ns("transcriptSelect"), "Transcript to Plot", choices = NULL),
                downloadButton(ns("downloadPlot"), "Download Transcript Plots"),
                selectizeInput(ns("embedding"), "Embedding", choices = NULL, selected = NULL),
                plotlyOutput(ns("transcriptPlot")),
                # uiOutput(ns("plotlys")),
                width = 6
            ),
            type = "markdown", content = "allTranscripts"
        ),
        default_helper(
            chevreulBox(
                title = "Transcript Expression per Gene",
                selectizeInput(ns("compositionGene"), "Gene or transcript expression by which to color the plot; eg. 'NRL'", choices = NULL, selected = NULL),
                selectizeInput(ns("groupby"), "Group by:", choices = NULL, selected = NULL),
                actionButton(ns("plotComposition"), "Plot transcript composition"),
                checkboxInput(ns("standardizeExpression"), "Standardize Expression", value = FALSE),
                checkboxInput(ns("dropZero"), "Drop Zero Values", value = FALSE),
                plotlyOutput(ns("compositionPlot")),
                DTOutput(ns("compositionDT")),
                width = 6
            ),
            type = "markdown", content = "allTranscripts"
        )
    )
}

#' Plot All Transcripts Server
#'
#' @param object a SingleCellExperiment object
#' @param featureType gene or transcript
#'
#' @noRd
allTranscripts <- function(
        input, output, session, object,
        featureType, organism_type) {
    ns <- session$ns

    observe({
        req(object())
        updateSelectizeInput(session, "compositionGene", 
                             choices = get_features(object(), 
                                                    experiment = "gene"), 
                             selected = "NRL", server = TRUE)
        updateSelectizeInput(session, "embeddingGene", 
                             choices = get_features(object(), 
                                                    experiment = "gene"), 
                             selected = "NRL", server = TRUE)

        formatted_col_names <- colnames(get_cell_metadata(object())) %>%
            make_chevreul_clean_names()

        updateSelectizeInput(session, "groupby", choices = formatted_col_names, 
                             selected = "batch", server = TRUE)
    })

    transcripts <- reactive({
        req(object())
        req(input$embeddingGene)
        if (query_experiment(object(), "transcript")) {
            get_transcripts_from_object(object(), input$embeddingGene, 
                                        organism = organism_type())
        }
    })

    observe({
        req(object())
        req(transcripts())
        updateSelectizeInput(session, "embedding", choices = c("PCA", "TSNE", "UMAP"), selected = "UMAP", server = TRUE)
        updateSelectizeInput(session, "transcriptSelect", choices = transcripts(), server = TRUE)
    })

    composition_plot <- eventReactive(input$plotComposition, {
        plot_transcript_composition(object(), 
                                    gene_symbol = input$compositionGene, 
                                    group.by = input$groupby, 
                                    standardize = input$standardizeExpression, 
                                    drop_zero = input$dropZero)
    })

    output$compositionPlot <- renderPlotly({
        composition_plot()$plot %>%
            ggplotly(height = 400) %>%
            plotly_settings() %>%
            toWebGL() %>%
            # partial_bundle() %>%
            identity()
    })

    output$compositionDT <- renderDT({
        datatable(composition_plot()$data,
            extensions = "Buttons",
            options = list(dom = "Bft", buttons = c(
                "copy",
                "csv"
            ), scrollX = "100px", scrollY = "400px")
        )
    })

    pList <- reactive({
        req(transcripts())
        req(input$embedding)
        pList <- plot_all_transcripts(object(), transcripts(), input$embedding, 
                                      from_gene = FALSE, combine = FALSE)
    })

    output$transcriptPlot <- renderPlotly({
        req(pList())
        pList()[[input$transcriptSelect]] %>%
            ggplotly(height = 400) %>%
            plotly_settings() %>%
            toWebGL()
    })

    output$downloadPlot <- downloadHandler(
        filename = function() {
            paste(input$embeddingGene, "_transcripts.pdf", sep = "")
        },
        content = function(file) {
            pdf(file)
            map(pList(), print)
            dev.off()
        }
    )
}

#' Title
#'
#' @noRd
pathwayEnrichmentui <- function(id) {
    ns <- NS(id)
    chevreulBox(
        title = "Enriched pathways by cluster",
        tagList(
            actionButton(ns("calcPathwayEnrichment"), 
                         "Calculate Pathway Enrichment"),
            uiOutput(ns("enriched_pathways_by_cluster_select_source_UI")),
            uiOutput(ns("enriched_pathways_by_cluster_UI"))
        ),
        width = 12
    )
}

#' pathway enrichment
#'
#' @noRd
pathwayEnrichment <- function(input, output, session, object, featureType) {
    ns <- session$ns

    enriched_pathways <- eventReactive(input$calcPathwayEnrichment, {
        req(object())
        if (featureType() == "gene") {
            enriched_object <- tryCatch(getEnrichedPathways(object()), 
                                        error = function(e) e)
            enrichr_available <- !any(is(enriched_object, "error"))
            if (enrichr_available) {
                object <- enriched_object
            }
        }

        metadata(object())$enriched_pathways
    })

    # UI element: choose source for pathway enrichement results (currently 
    # Enrichr or GSVA)
    output$enriched_pathways_by_cluster_select_source_UI <- renderUI({
        req(object())
        if (is.null(enriched_pathways())) {
            textOutput(ns("enriched_pathways_by_cluster_table_missing"))
        } else {
            selectInput(
                ns("enriched_pathways_by_cluster_select_source"),
                label = NULL,
                choices = names(enriched_pathways())
            )
        }
    })

    # UI element: display results or alternative text
    output$enriched_pathways_by_cluster_UI <- renderUI({
        req(object())
        req(input$enriched_pathways_by_cluster_select_source)
        if (input$enriched_pathways_by_cluster_select_source == "enrichr") {
            if (!is.null(enriched_pathways()$enrichr$by_cluster)) {
                if (is.list(enriched_pathways()$enrichr$by_cluster)) {
                    tagList(
                        fluidRow(
                            column(
                                4,
                                uiOutput(ns("enriched_pathways_by_cluster_select_cluster_UI"))
                            ),
                            column(
                                8,
                                uiOutput(ns("enriched_pathways_by_cluster_select_db_UI"))
                            )
                        ),
                        DTOutput(ns("enriched_pathways_by_cluster_table_present"))
                    )
                } else if (enriched_pathways()$enrichr$by_cluster == "no_markers_found") {
                    textOutput(ns("enriched_pathways_by_cluster_table_no_markers_found"))
                }
            } else {
                textOutput(ns("enriched_pathways_by_cluster_table_missing_enrichr"))
            }
        }
    })


    # UI element: choose cluster
    output$enriched_pathways_by_cluster_select_cluster_UI <- renderUI({
        req(object())
        req(input$enriched_pathways_by_cluster_select_source)
        if (input$enriched_pathways_by_cluster_select_source == "enrichr") {
            choices <- levels(enriched_pathways()$enrichr$by_cluster$cluster) %>%
                intersect(., unique(enriched_pathways()$enrichr$by_cluster$cluster))
        }
        selectInput(
            ns("enriched_pathways_by_cluster_select_cluster"),
            label = NULL,
            choices = choices
        )
    })

    # UI element: choose database
    output$enriched_pathways_by_cluster_select_db_UI <- renderUI({
        req(
            input$enriched_pathways_by_cluster_select_source,
            input$enriched_pathways_by_cluster_select_cluster
        )
        choices <- enriched_pathways()$enrichr$by_cluster %>%
            filter(cluster == input$enriched_pathways_by_cluster_select_cluster) %>%
            pull(db) %>%
            intersect(., levels(.))
        selectInput(
            ns("enriched_pathways_by_cluster_select_db"),
            label = NULL,
            choices = choices
        )
    })

    # table
    output$enriched_pathways_by_cluster_table_present <- renderDT(server = FALSE, {
        req(
            input$enriched_pathways_by_cluster_select_source,
            input$enriched_pathways_by_cluster_select_cluster,
            input$enriched_pathways_by_cluster_select_db
        )
        if (input$enriched_pathways_by_cluster_select_source == "enrichr" & is.data.frame(enriched_pathways()$enrichr$by_cluster)) {
            format_pathway_table(
                enriched_pathways()$enrichr$by_cluster,
                input$enriched_pathways_by_cluster_select_cluster,
                input$enriched_pathways_by_cluster_select_db
            )
        }
    })

    # # alternative text messages
    output$enriched_pathways_by_cluster_table_missing <- renderText({
        "Data not available. Possible reason: Data not generated."
    })

    output$enriched_pathways_by_cluster_table_no_markers_found <- renderText({
        "No marker genes identified to perform pathway enrichment analysis with."
    })

    output$enriched_pathways_by_cluster_table_missing_enrichr <- renderText({
        "Data not available. Possible reasons: Only 1 cluster in this data set, no marker genes found or data not generated."
    })

    output$enriched_pathways_by_cluster_table_no_gene_sets_enriched <- renderText({
        "Either the loaded data set consists of a single cluster (in which case GSVA cannot be applied) or no gene sets were found to be enriched (with the selected statistical thresholds) in any cluster."
    })

    output$enriched_pathways_by_cluster_table_only_one_cluster_in_data_set <- renderText({
        "The loaded data set consists of a single cluster which means GSVA cannot be applied."
    })

    output$enriched_pathways_by_cluster_table_missing_gsva <- renderText({
        "Data not available. Possible reason: Data not generated."
    })
    # info box
    observeEvent(input$enriched_pathways_by_cluster_info, {
        showModal(
            modalDialog(
                enriched_pathways_by_cluster_info$text,
                title = enriched_pathways_by_cluster_info$title,
                easyClose = TRUE,
                footer = NULL
            )
        )
    })
}

#' Title
#'
#' @noRd
techInfoui <- function(id) {
    ns <- NS(id)
    fluidRow(
        chevreulBox(
            title = "Information about samples and analysis",
            htmlOutput(ns("sample_info_general")),
            width = 12
        )
    )
}

#' Title
#'
#' @param object a SingleCellExperiment object
#'
#' @noRd
techInfo <- function(input, output, session, object) {
    ns <- session$ns

    object_metadata <- reactive({
        req(object())
        metadata(object())$experiment
    })

    observe({
        # general info
        output$sample_info_general <- renderText({
            info <- paste0(
                "<strong><u>General</u></strong>",
                "<ul>",
                "<li><b>Date of analysis:</b> ",
                object_metadata()$date_of_analysis,
                "<li><b>Date of export:</b> ",
                object_metadata()$date_of_export,
                "<li><b>Experiment name:</b> ",
                object_metadata()$experiment_name,
                "<li><b>Organism:</b> ",
                object_metadata()$organism,
                "</ul>",
                "<strong><u>Parameters</u></strong>",
                "<ul>",
                "<li><b>Discard genes in fewer than X cells:</b> ",
                object_metadata()$parameters$discard_genes_expressed_in_fewer_cells_than,
                "<li><b>Keep mitochondrial genes:</b> ",
                object_metadata()$parameters$keep_mitochondrial_genes,
                "<li><b>Min/max # of UMI:</b> ",
                paste0(
                    object_metadata()$filtering$UMI_min, " / ",
                    object_metadata()$filtering$UMI_max
                ),
                "<li><b>Min/max # of expressed genes:</b> ",
                paste0(
                    object_metadata()$filtering$genes_min, " / ",
                    object_metadata()$filtering$genes_max
                ),
                "<li><b>Cluster resolution: </b>",
                paste(object_metadata()$parameters$cluster_resolution, collapse = ","),
                "<li><b>Number of principal components: </b>",
                object_metadata()$parameters$number_PCs,
                "<li><b>Variables to regress: </b>",
                object_metadata()$parameters$variables_to_regress_out,
                "<li><b>tSNE perplexity: </b>",
                object_metadata()$parameters$tSNE_perplexity,
                "</ul>",
                "<strong><u>Marker genes</u></strong>",
                "<ul>",
                # "<li><b>Only positive:</b> ",
                # object_metadata()$marker_genes$parameters$only_positive,
                # "<li><b>Fraction of cells in group of interest that must express marker gene:</b> ",
                # object_metadata()$marker_genes$parameters$minimum_percentage,
                # "<li><b>LogFC threshold:</b> ",
                # object_metadata()$marker_genes$parameters$logFC_threshold,
                "<li><b>p-value threshold:</b> ",
                "0.05",
                # object_metadata()$marker_genes$parameters$p_value_threshold,
                "</ul>",
                "<strong><u>Pathway enrichment</u></strong>",
                "<ul>",
                # "<li><b>Enrichr:</b>",
                # "<ul>",
                # "<li><b>Databases:</b> ",
                # paste0(object_metadata()$enriched_pathways$enrichr$parameters$databases, collapse = ", "),
                # "<li><b>Adj. p-value cut-off:</b> ",
                # object_metadata()$enriched_pathways$enrichr$parameters$adj_p_cutoff,
                # "<li><b>Max. terms:</b> ",
                # object_metadata()$enriched_pathways$enrichr$parameters$max_terms,
                # "</ul>",
                "</ul>"
            )
            paste0(
              info,
                "<strong><u>Technical info (package versions)</u></strong>",
                "<ul>",
                "<li><strong>chevreul version:</strong> ",
                object_metadata()$chevreul_version,
                "<li><strong>SingleCellExperiment version:</strong> ",
                object_metadata()$SingleCellExperiment_version,
                "<li><strong>Session info:</strong> ",
                "</ul>",
              "<pre>",
              object_metadata()$sessionInfo,
              "</pre>"
            )
        })
    })
}


#' Title
#'
#' @noRd
plotCoverage_UI <- function(id) {
    ns <- NS(id)
    tagList(
        chevreulBox(
            title = "Plot Coverage",
            selectizeInput(ns("geneSelect"), "Select a Gene", choices = NULL, 
                           selected = "NRL", multiple = FALSE),
            selectizeInput(ns("varSelect"), "Color by Variable", 
                           choices = NULL, multiple = FALSE),
            actionButton(ns("plotCoverage"), "Plot Coverage"),
            downloadButton(ns("downloadPlot"), "Download Coverage Plot"),
            uiOutput(ns("displayvaluesui")),
            br(),
            chevreulDropDownButton(
                ns("coveragePlotSettings"),
                checkboxInput(ns("collapseIntrons"), "Collapse Introns", 
                              value = TRUE),
                checkboxInput(ns("meanCoverage"), 
                              "Summarize Coverage to Mean", value = TRUE),
                checkboxInput(ns("summarizeTranscripts"), 
                              "Summarize transcript models to gene", 
                              value = FALSE),
                radioButtons(ns("yScale"), "Scale Y Axis", 
                             choices = c("absolute", "log10"), 
                             selected = "log10"),
                numericInput(ns("start"), "start coordinate", value = NULL),
                numericInput(ns("end"), "end coordinate", value = NULL)
            ),
            DTOutput(ns("coverageTable")),
            plotOutput(ns("coveragePlot"), height = "1500px"),
            width = 12
        )
    )
}


#' Plot Coverage Module
#'
#' @param object a SingleCellExperiment object
#' @param plot_types plot types
#' @param bigwig_dir bigwig directory
#' @param organism_type human or mouse
#'
#' @noRd
plotCoverage <- function(input, output, session, object, plot_types, proj_dir,
                         organism_type = "human", 
                         bigwig_db = "~/.cache/chevreul/bw-files.db") {
    ns <- session$ns

    w <- Waiter$new(ns("coveragePlot"),
        html = spin_loaders(id = 1, color = "black", 
                            style = "position:relative;margin:auto;"),
        color = transparent(.5)
    )

    observe({
        req(object())
        updateSelectizeInput(session, "geneSelect", 
                             choices = get_features(object(), 
                                                    experiment = "gene"), 
                             server = TRUE)

        formatted_col_names <- colnames(get_cell_metadata(object())) %>%
            make_chevreul_clean_names()

        updateSelectizeInput(session, "varSelect", 
                             choices = formatted_col_names, selected = "batch")
    })

    displayvalues <- reactive({
        req(input$varSelect)
        req(object())
        unique(object()[][[input$varSelect]])
    })

    output$displayvaluesui <- renderUI({
        req(input$varSelect)
        selectizeInput(ns("displayvalues"), "groups to display", 
                       choices = displayvalues(), multiple = TRUE)
    })

    bigwig_tbl <- reactive({
        load_bigwigs(object(), bigwig_db)
    })

    coverage_return <- eventReactive(input$plotCoverage, {
        req(object())
        req(bigwig_tbl())

        plot_gene_coverage_by_var(
            genes_of_interest = input$geneSelect,
            cell_metadata = get_cell_metadata(object()),
            bigwig_tbl = bigwig_tbl(),
            group_by = input$varSelect,
            values_of_interest = input$displayvalues,
            organism = metadata(object())$experiment$organism,
            mean_only = input$meanCoverage,
            rescale_introns = input$collapseIntrons,
            scale_y = input$yScale,
            start = input$start,
            end = input$end,
            summarize_transcripts = input$summarizeTranscripts
        )
    })

    output$coveragePlot <- renderPlot({
        w$show()

        coverage_return()$plot
    })

    output$coverageTable <- renderDT({
        datatable(coverage_return()$table,
            extensions = "Buttons",
            options = list(dom = "Bft", buttons = c("copy", "csv"), 
                           scrollY = "400px")
        )
    })

    output$downloadPlot <- downloadHandler(
        filename = function() {
            paste("coverage", ".pdf", sep = "")
        },
        content = function(file) {
            ggsave(file, coverage_return()$plot, width = 16, height = 12)
        }
    )
}

#' Reformat SingleCellExperiment Object Metadata UI
#'
#' @noRd
reformatMetadataDRui <- function(id) {
    ns <- NS(id)
    tagList(
        chevreulBox(
            title = "Reformat Metadata",
            checkboxInput(ns("header"), "Header", TRUE),
            fileInput(
                ns("metaFile"), 
                "Choose CSV File of metadata with cell names in first column",
                accept = c(
                    "text/csv",
                    "text/comma-separated-values,text/plain",
                    ".csv"
                )
            ),
            actionButton(ns("updateMetadata"), "Update Metadata"),
            radioButtons(
                ns("updateMethod"), "Update By:", 
                         choices = c("table (below)" = "spreadsheet", "uploaded file" = "file"), inline = TRUE),
            width = 12,
            dataSelectUI(ns("select1")),
            dataFilterUI(ns("filter1")),
            hidden(actionButton(ns("sync"), label = NULL, icon = icon("sync"))),
            dataOutputUI(ns("output-active")),
            dataOutputUI(ns("output-update"), icon = "file-download"),
            hidden(actionButton(ns("cut"), label = NULL, icon = icon("cut"))),
            dataEditUI(ns("edit1"))
        ) %>%
            default_helper(type = "markdown", content = "reformatMetadata")
    )
}

#' Reformat SingleCellExperiment Object Metadata Server
#'
#' @param object a SingleCellExperiment object
#'
#' @noRd
reformatMetadataDR <- function(input, output, session, object, 
                               featureType = "gene",
    col_bind = NULL,
    col_edit = TRUE,
    col_options = NULL,
    col_stretch = FALSE,
    col_names = TRUE,
    col_readonly = NULL,
    col_factor = FALSE,
    row_bind = NULL,
    row_edit = TRUE,
    save_as = NULL,
    title = NULL,
    logo = NULL,
    logo_size = 30,
    logo_side = "left",
    viewer = "dialog",
    viewer_height = 800,
    viewer_width = 2000,
    theme = "yeti",
    read_fun = "read.csv",
    read_args = NULL,
    write_fun = "write.csv",
    write_args = NULL,
    quiet = FALSE,
    code = FALSE,
    hide = FALSE) {
    ns <- session$ns

    table_out <- reactive({
        req(object())
        get_cell_metadata(object())
    })

    values <- reactiveValues(
        data = NULL, data_active = NULL,
        rows = NULL, columns = NULL, cut = FALSE
    )

    observeEvent(table_out(), {
        values$rows <- NULL
        values$columns <- NULL

        values$data <- table_out() %>%
            data_bind_rows(row_bind = row_bind) %>%
            data_bind_cols(col_bind = col_bind) %>%
            identity()
    })

    data_select <- dataSelectServer("select1",
        data = reactive(values$data),
        hide = hide
    )
    data_filter <- dataFilterServer("filter1",
        data = reactive(values$data),
        hide = hide
    )
    observe({
        values$rows <- data_filter$rows()
        values$columns <- data_select$columns()
    })

    observe({
        if (length(values$rows) == 0 & length(values$columns) == 0) {
            values$data_active <- values$data
        } else {
            if (length(values$rows) != 0 & length(values$columns) == 0) {
                values$data_active <- 
                    values$data[values$rows, , drop = FALSE]
            } else if (length(values$rows) == 0 & length(values$columns) != 0) {
                values$data_active <- 
                    values$data[, values$columns, drop = FALSE]
            } else if (length(values$rows) != 0 & length(values$columns) != 0) {
                values$data_active <- 
                    values$data[values$rows, values$columns, drop = FALSE]
            }
        }
    })

    data_update <- dataEditServer("edit1",
        data = reactive({
            values$data_active
        }),
        col_bind = NULL, col_edit = col_edit, col_options = col_options,
        col_stretch = col_stretch, col_names = col_names,
        col_readonly = col_readonly, col_factor = col_factor,
        row_bind = NULL, row_edit = row_edit, quiet = quiet, 
        height = viewer_height, width = viewer_width
    )
    observe({
        values$data_active <- data_update()
    })

    observeEvent(input$updateMetadata, {
        if (input$updateMethod == "file") {
            inFile <- input$metaFile

            if (is.null(inFile)) {
                return(NULL)
            }

            object(set_cell_metadata(object(), read_csv(inFile$datapath)))
        } else if (input$updateMethod == "spreadsheet") {
            reformatted_object <- 
                propagate_spreadsheet_changes(values$data_active, object())
            object(reformatted_object)
        }
    })


    observeEvent(input$sync, {
        if (length(values$rows) == 0 & length(values$columns) == 0) {
            values$data <- values$data_active
        } else {
            if (length(values$rows) != 0 & length(values$columns) == 0) {
                values$data[values$rows, ] <- 
                    values$data_active
            } else if (length(values$rows) == 0 & length(values$columns) != 0) {
                values$data[, values$columns] <- 
                    values$data_active
            } else if (length(values$rows) != 0 & length(values$columns) != 0) {
                values$data[values$rows, values$columns] <- 
                    values$data_active
            }
            if (!is.null(values$data_active)) {
                if (!all(rownames(values$data_active) == rownames(values$data)[values$rows])) {
                    rownames(values$data)[values$rows] <- 
                        rownames(values$data_active)
                }
                if (!all(colnames(values$data_active) == colnames(values$data)[values$columns])) {
                    colnames(values$data)[values$columns] <- 
                        colnames(values$data_active)
                }
            }
        }
    })

    dataOutputServer("output-active",
        data = reactive({
            values$data_active
        }), save_as = "metadata.csv", write_fun = write_fun, 
        write_args = write_args,
        hide = hide
    )
    dataOutputServer("output-update",
        data = reactive({
            values$data
        }), save_as = "metadata.csv", write_fun = write_fun, 
        write_args = write_args,
        hide = hide
    )

    # SAVE AS
    if (!hide & !is.null(save_as)) {
        do.call(
            write_fun,
            c(list(x_edit, save_as), write_args)
        )
    }

    observeEvent(input$cut, {
        if (values$cut) {
            values$cut <- FALSE
            updateButton(session, "cut", NULL,
                block = FALSE,
                style = "danger"
            )
        } else {
            values$cut <- TRUE
            updateButton(session, "cut", NULL,
                block = FALSE,
                style = "success"
            )
        }
    })

    return(object)
}
