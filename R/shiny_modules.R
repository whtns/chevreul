#' plot clustree ui
#'
#' @param id
#'
#' @return
#' @export
#'
#' @examples
plotClustree_UI <- function(id) {
    ns <- NS(id)
    tagList(
        chevreulBox(
            title = "Clustering Tree",
            # textOutput(ns("checkSeu")),
            plotOutput(ns("clustree"), height = "700px")
        )
    )
}

#' plot clustree server
#'
#' @param input
#' @param output
#' @param session
#' @param seu
#'
#' @return
#' @export
#'
#' @examples
plotClustree <- function(input, output, session, seu) {
    # set appropriate assay
    # assay = reactive({
    #   ifelse("integrated" %in% names(seu()@assays), "integrated", "gene")
    # })

    output$checkSeu <- renderText({
        req(seu())
        "test"
    })

    output$clustree <- renderPlot({
        req(seu())

        assay <- ifelse("integrated" %in% names(seu()@assays), "integrated", "gene")
        # DefaultAssay(seu()) <- assay
        clustree::clustree(seu(), assay = assay)
    })
}


#' Title
#'
#' @param id
#'
#' @return
#' @export
#'
#' @examples
plotViolinui <- function(id) {
    ns <- NS(id)
    tagList(
        chevreulBox(
            title = "Violin Plots",
            uiOutput(ns("vln_group")),
            selectizeInput(ns("customFeature"),
                "Gene or transcript expression by which to color the plot eg. 'RXRG' or 'ENST00000488147'",
                choices = NULL, multiple = TRUE
            ),
            radioButtons(ns("slot"), "Data Type", choices = c("transformed" = "data", "raw counts" = "counts")),
            downloadButton(ns("downloadPlot")),
            plotly::plotlyOutput(ns("vplot"), height = 750),
            width = 11
        ) %>%
            default_helper(type = "markdown", content = "violinPlot", size = "l")
    )
}

#' Plot Violin Server
#'
#' Plots a Violin plot of a single data (gene expression, metrics, etc.) in the server Seurat app.
#'
#' @param input
#' @param output
#' @param session
#' @param seu Seurat object
#' @param featureType Gene or Transcript
#' @param organism_type Organism
#'
#' @return
#' @export
#'
#' @examples
plotViolin <- function(input, output, session, seu, featureType, organism_type) {
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
                "RXRG"
            } else if (organism_type() == "mouse") {
                "Rxrg"
            }
        }
    })
    observe({
        req(prefill_feature())
        req(seu())
        updateSelectizeInput(session, "customFeature",
            choices = rownames(seu()@assays[["gene"]]),
            selected = prefill_feature(), server = TRUE
        )
    })

    output$vln_group <- renderUI({
        req(seu())
        selectizeInput(ns("vlnGroup"), "Grouping variable",
            choices = colnames(seu()[[]]), selected = "batch"
        )
    })

    vln_plot <- reactive({
        req(input$customFeature)
        req(input$vlnGroup)

        vln_plot <-
            plot_violin(seu(), plot_var = input$vlnGroup, features = input$customFeature, slot = input$slot)
    })

    output$downloadPlot <- downloadHandler(
        filename = function() {
            paste("violin", ".pdf", sep = "")
        },
        content = function(file) {
            ggsave(file, vln_plot() + ggpubr::theme_pubr(base_size = 20, x.text.angle = 45), width = 16, height = 12)
        }
    )

    output$vplot <- plotly::renderPlotly({
        req(seu())
        req(input$vlnGroup)
        exclude_trace_number <- length(unique(seu()[[]][[input$vlnGroup]])) * 2

        vln_plot <- plotly::ggplotly(vln_plot(), height = 700) %>%
            plotly::style(opacity = 0.5) %>%
            plotly::style(hoverinfo = "skip", traces = c(1:exclude_trace_number)) %>%
            plotly_settings(width = 1200) %>%
            plotly::toWebGL() %>%
            identity()
    })
}


#' Plot Heatmap ui
#'
#' @param id
#'
#' @return
#' @export
#'
#' @examples
plotHeatmapui <- function(id) {
    ns <- NS(id)
    tagList(
        chevreulBox(
            title = "Heatmap",
            uiOutput(ns("colAnnoVarui")),
            radioButtons(ns("slot"), "Data Scaling", choices = c(scaled = "scale.data", unscaled = "data"), selected = "scale.data", inline = TRUE),
            selectizeInput(ns("dendroSelect"), "Clustering algorithm or metadata for column arrangement", choices = NULL, selected = NULL, multiple = TRUE),
            actionButton(ns("actionHeatmap"), "Plot Heatmap"),
            downloadButton(ns("downloadPlot"), "Download Heatmap"),
            selectizeInput(ns("customFeature"), "Gene or transcript expression by which to color the plot; eg. 'RXRG' or 'ENST00000488147'",
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
#' @param input
#' @param output
#' @param session
#' @param seu
#' @param featureType
#' @param organism_type
#'
#' @return
#' @export
#'
#' @examples
plotHeatmap <- function(input, output, session, seu, featureType, organism_type) {
    ns <- session$ns

    w <- waiter::Waiter$new(ns("heatmap"),
        html = waiter::spin_loaders(id = 1, color = "black", style = "position:relative;margin:auto;"),
        color = waiter::transparent(.5)
    )

    observe({
        req(seu())
        if ("integrated" %in% names(seu()@assays)) {
            assay <- "integrated"
        } else {
            assay <- "gene"
        }

        preset_features <- VariableFeatures(seu(), assay = assay)[1:50]

        updateSelectizeInput(session, "customFeature",
            choices = rownames(seu()@assays[["gene"]]),
            selected = preset_features, server = TRUE
        )
    })

    output$colAnnoVarui <- renderUI({
        req(seu())

        formatted_col_names <- colnames(seu()@meta.data) %>%
            make_chevreul_clean_names()

        selectizeInput(ns("colAnnoVar"), "Column Annotation(s)",
            choices = formatted_col_names, selected = "batch", multiple = TRUE
        )
    })

    observe({
        req(seu())

        hclust_methods <- c("Ward" = "ward.D2", "single", "complete", "average")

        updateSelectizeInput(session, "dendroSelect", choices = c(hclust_methods, colnames(seu()[[]])), selected = "ward.D2")
    })

    heatmap_plot <- eventReactive(input$actionHeatmap, {
        req(input$customFeature)
        req(input$colAnnoVar)

        if ("integrated" %in% names(seu()@assays)) {
            assay <- "integrated"
        } else {
            assay <- "gene"
        }

        hm <- seu_complex_heatmap(seu(), features = input$customFeature, assay = assay, group.by = input$colAnnoVar, slot = input$slot, col_arrangement = input$dendroSelect)

        hm <- ComplexHeatmap::draw(hm)
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
            ggsave(file, ggplotify::as.ggplot(heatmap_plot()) + ggpubr::theme_pubr(base_size = 20, x.text.angle = 45), width = 16, height = 12)
        }
    )
}

#' Integrate Project UI
#'
#' @param id
#'
#' @return
#' @export
#'
#' @examples
integrateProjui <- function(id) {
    ns <- NS(id)
    tagList(
        chevreulBox(
            title = "Integrate Projects",
            actionButton(ns("integrateAction"), "Integrate Selected Projects"),
            # shinyjs::useShinyjs(),
            # shinyjs::runcodeUI(code = "shinyjs::alert('Hello!')", id = "subsetcode"),
            textOutput(ns("integrationMessages")),
            checkboxInput(ns("legacySettings"), "Use Legacy Settings", value = FALSE),
            textOutput(ns("integrationResult")),
            shinyFiles::shinySaveButton(ns("saveIntegratedProject"), "Save Integrated Project", "Save project as..."),
            DT::dataTableOutput(ns("myDatatable")),
            width = 12
        ) %>%
            default_helper(type = "markdown", content = "integrateProjects")
    )
}

#' Integrate Projects Server Function
#'
#' @param input
#' @param output
#' @param proj_matrices
#' @param session
#'
#' @return
#' @export
#'
#' @examples
integrateProj <- function(input, output, session, proj_matrices, seu, proj_dir, con) {
    ns <- session$ns

    proj_matrix <- reactive({
        proj_matrices()$primary_projects
    })

    clean_proj_matrix <- reactive({
        clean_proj_matrix <- proj_matrix() %>%
            dplyr::select(-project_path) %>%
            identity()
    })

    output$myDatatable <- DT::renderDT(clean_proj_matrix(),
        server = FALSE,
        rownames = TRUE
    )

    selectedRows <- eventReactive(input$integrateAction, {
        ids <- input$myDatatable_rows_selected
    })

    selectedProjects <- reactive({
        selectedProjects <- dplyr::slice(proj_matrix(), selectedRows()) %>%
            dplyr::pull(project_path) %>%
            identity()
    })

    mergedSeus <- reactiveVal()

    observeEvent(input$integrateAction, {
        req(selectedProjects())
        withCallingHandlers(
            {
                shinyjs::html("integrationMessages", "")
                message("Beginning")
                message(selectedProjects())
                batches <- fs::path(selectedProjects(), "output", "seurat", "unfiltered_seu.rds") %>%
                    purrr::map(readRDS)

                names(batches) <- names(selectedProjects())
                print(names(batches))
                mergedSeus(integration_workflow(batches, legacy_settings = input$legacySettings))
                # mergedSeus(batches[[1]])

                message("Integration Complete!")
            },
            message = function(m) {
                shinyjs::html(id = "integrationMessages", html = paste0("Running Integration: ", m$message), add = FALSE)
            }
        )
    })

    newProjDir <- reactive({
        req(mergedSeus())
        print("foo created successfully")
        # print(names(mergedSeus()))
        #
        # for (i in names(mergedSeus())) {
        #   seu[[i]] <- mergedSeus()[[i]]
        # }
        #

        newProjName <- paste0(purrr::map(fs::path_file(selectedProjects()), ~ gsub("_proj", "", .x)), collapse = "_")
        integrated_proj_dir <- "/dataVolume/storage/single_cell_projects/integrated_projects/"
        newProjDir <- fs::path(integrated_proj_dir, newProjName)

        proj_dir(newProjDir)

        newProjDir
    })


    volumes <- reactive({
        volumes <- c(
            Home = "/dataVolume/storage/single_cell_projects/integrated_projects/",
            "R Installation" = R.home(),
            shinyFiles::getVolumes()()
        )
        # print(volumes)
        volumes
    })

    observe({
        shinyFiles::shinyFileSave(input,
            "saveIntegratedProject",
            roots = volumes(),
            session = session,
            restrictions = system.file(package = "base")
        )
    })


    integratedProjectSavePath <- eventReactive(input$saveIntegratedProject, {
        savefile <- shinyFiles::parseSavePath(volumes(), input$saveIntegratedProject)

        savefile$datapath
    })

    output$integrationResult <- renderText({
        integratedProjectSavePath()
    })

    observeEvent(input$saveIntegratedProject, {
        req(mergedSeus())
        req(integratedProjectSavePath())

        if (!is.null(integratedProjectSavePath())) {
            shiny::withProgress(
                message = paste0("Saving Integrated Dataset to ", integratedProjectSavePath()),
                value = 0,
                {
                    # Sys.sleep(6)
                    shiny::incProgress(2 / 10)
                    save_seurat(mergedSeus(), proj_dir = integratedProjectSavePath())
                    set_permissions_call <- paste0("chmod -R 775 ", integratedProjectSavePath())
                    system(set_permissions_call)
                    writeLines(character(), fs::path(integratedProjectSavePath(), ".here"))
                    # create_proj_db()
                    DBI::dbAppendTable(con, "projects_tbl", data.frame(
                        project_name = fs::path_file(integratedProjectSavePath()),
                        project_path = integratedProjectSavePath(),
                        project_slug = stringr::str_remove(fs::path_file(integratedProjectSavePath()), "_proj$"),
                        project_type = "integrated_projects"
                    ))
                    shiny::incProgress(8 / 10)

                    velocyto_dir <- fs::path(integratedProjectSavePath(), "output", "velocyto")
                    fs::dir_create(velocyto_dir)
                    new_loom_path <- fs::path(velocyto_dir, fs::path_file(integratedProjectSavePath()))
                    # need to configure conda for line below
                    combine_looms(selectedProjects(), new_loom_path)
                }
            )
        }
    })


    return(integratedProjectSavePath)
}


#' Change Embedding Parameters UI
#'
#' @param id
#'
#' @return
#' @export
#'
#' @examples
changeEmbedParamsui <- function(id) {
    ns <- NS(id)

    minDist_vals <- prep_slider_values(0.3)
    negsamprate_vals <- prep_slider_values(5)

    tagList(
        selectizeInput(ns("dims"), label = "Dimensions from PCA", choices = seq(1, 99), multiple = TRUE, selected = 1:30),
        sliderInput(ns("minDist"), label = "Minimum Distance", min = minDist_vals$min, max = minDist_vals$max, value = minDist_vals$value, step = minDist_vals$step),
        sliderInput(ns("negativeSampleRate"), label = "Negative Sample Rate", min = negsamprate_vals$min, max = negsamprate_vals$max, value = negsamprate_vals$value, step = negsamprate_vals$step)
    )
}

#' Change Embedding Parameters
#'
#' @param input
#' @param output
#' @param session
#' @param seu
#'
#' @return
#' @export
#'
#' @examples
changeEmbedParams <- function(input, output, session, seu) {
    ns <- session$ns

    seu <- RunUMAP(seu(), dims = as.numeric(input$dims), reduction = "pca", min.dist = input$minDist, negative.sample.rate = input$negativeSampleRate)

    return(seu)
}

#' Plot Dimensional Reduduction UI
#'
#' @param id
#'
#' @return
#' @export
#'
#' @examples
plotDimRedui <- function(id) {
    ns <- NS(id)
    chevreulBox(
        title = "Embedding",
        dropdownButton(
            ns("dimPlotSettings"),
            selectizeInput(ns("embedding"), "Embedding", choices = NULL, selected = NULL),
            sliderInput(ns("dotSize"), "Size of Points in UMAP", min = 0.5, max = 2, step = 0.1, value = 1),
            selectizeInput(ns("dim1"), "Dimension 1", choices = seq(1, 99), selected = 1),
            selectizeInput(ns("dim2"), "Dimension 2", choices = seq(1, 99), selected = 2)
        ),
        selectizeInput(ns("plottype"), "Variable to Plot", choices = NULL, multiple = TRUE),
        selectizeInput(ns("customFeature"), "Gene or transcript expression by which to color the plot; eg. 'RXRG' or 'ENST00000488147'", choices = NULL, multiple = TRUE),
        plotly::plotlyOutput(ns("dplot"), height = 500),
        width = 6
    )
}

#' Plot Dimensional Reduduction
#'
#' @param input
#' @param output
#' @param session
#' @param seu
#' @param plot_types
#' @param featureType
#' @param organism_type
#' @param reductions
#'
#' @return
#' @export
#'
#' @examples
plotDimRed <- function(input, output, session, seu, plot_types, featureType,
    organism_type, reductions) {
    ns <- session$ns

    output$myPanel <- renderUI({
        req()
        lev <- sort(unique(input$select)) # sorting so that "things" are unambigious
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
        req(seu())
        updateSelectizeInput(session, "embedding",
            choices = reductions(),
            selected = rev(reductions())[1], server = TRUE
        )
    })
    selected_plot <- reactiveVal()
    observe({
        req(seu())
        # selected_plot <- ifelse(is.null(selected_plot()), "louvain",
        #                         selected_plot())
        updateSelectizeInput(session, "plottype",
            choices = purrr::flatten_chr(plot_types()),
            selected = purrr::flatten_chr(plot_types())[[1]]
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
                "RXRG"
            } else if (organism_type() == "mouse") {
                "Rxrg"
            }
        }
    })
    observe({
        req(prefill_feature())
        req(seu())
        updateSelectizeInput(session, "customFeature",
            choices = rownames(seu()@assays[["gene"]]),
            selected = prefill_feature(), server = TRUE
        )
    })

    output$dplot <- plotly::renderPlotly({
        req(input$plottype)
        req(seu())
        req(input$embedding)
        if (length(input$plottype) > 1) {
            cross_plot_seu <- unite_metadata(seu(), input$plottype)

            newcolname <- paste(input$plottype, collapse = "_")
            cross_plot_seu[[newcolname]] <- Idents(cross_plot_seu)

            selected_plot(newcolname)

            plot_var(cross_plot_seu,
                dims = c(input$dim1, input$dim2),
                embedding = input$embedding, group = NULL, pt.size = input$dotSize,
                return_plotly = TRUE
            )
        } else {
            if (input$plottype == "feature") {
                plot_feature(seu(),
                    dims = c(
                        input$dim1,
                        input$dim2
                    ), embedding = input$embedding,
                    features = input$customFeature, pt.size = input$dotSize,
                    return_plotly = TRUE
                )
            } else if (input$plottype %in% plot_types()$continuous_vars) {
                plot_feature(seu(),
                    dims = c(
                        input$dim1,
                        input$dim2
                    ), embedding = input$embedding,
                    features = input$plottype, pt.size = input$dotSize,
                    return_plotly = TRUE
                )
            } else if (input$plottype %in% plot_types()$category_vars) {
                plot_var(seu(),
                    dims = c(input$dim1, input$dim2),
                    embedding = input$embedding, group = input$plottype, pt.size = input$dotSize,
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
#' @return
#' @export
#'
#' @examples
tableSelectedui <- function(id) {
    ns <- NS(id)
    tagList(DT::DTOutput(ns("brushtable")))
}

#' Create Table of Selected Cells
#'
#' @param input
#' @param output
#' @param session
#' @param seu
#'
#' @return
#' @export
#'
#' @examples
tableSelected <- function(input, output, session, seu) {
    ns <- session$ns
    brush <- reactive({
        req(seu())
        d <- plotly::event_data("plotly_selected")
        if (is.null(d)) {
            msg <- "Click and drag events (i.e. select/lasso) appear here (double-click to clear)"
            return(d)
        } else {
            # selected_cells <- colnames(seu())[as.numeric(d$key)]
            d$key
        }
    })

    output$brushtable <- DT::renderDT({
        req(seu())
        req(brush())
        selected_meta <- data.frame(seu()[[]][brush(), ])

        # selection = list(mode = 'multiple', selected = c(1, 3, 8), target = 'row'),
        DT::datatable(selected_meta,
            extensions = "Buttons",
            selection = list(mode = "multiple", selected = 1:nrow(selected_meta), target = "row"),
            options = list(dom = "Bft", buttons = c("copy", "csv"), scrollX = "100px", scrollY = "800px")
        )
    })

    selected_cells <- reactive({
        selected_rows <- input$brushtable_rows_selected
        rownames(seu()[[]][brush(), ])[selected_rows]
    })

    return(selected_cells)
}


#' Differential Expression UI
#'
#' @param id
#'
#' @return
#' @export
#'
#' @examples
diffexui <- function(id) {
    ns <- NS(id)
    tagList(
        chevreulBox(
            title = "Differential Expression Settings",
            radioButtons(ns("diffex_scheme"),
                "Cells to Compare",
                choiceNames = c("Seurat Cluster", "Custom Selection"), choiceValues = c("louvain", "custom"),
                selected = "louvain",
                inline = TRUE
            ),
            conditionalPanel(
                ns = ns,
                condition = "input.diffex_scheme == 'louvain'",
                sliderInput(ns("seuratResolution"), "Resolution of clustering algorithm (affects number of clusters)",
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
                sliderInput(ns("customResolution"), "Resolution of clustering algorithm (affects number of clusters)",
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
            radioButtons(ns("featureType"), "Features to Compare", choices = c("gene", "transcript")),
            actionButton(
                ns("diffex"),
                "Run Differential Expression"
            ),
            downloadLink(ns("downloadData"), "Download Complete DE Results"),
            DT::dataTableOutput(ns("DT1")),
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
        #               tabPanel("Custom Cluster 1", DT::DTOutput(ns("cc1"))),
        #               tabPanel("Custom Cluster 2", DT::DTOutput(ns("cc2")))
        #   ),
        #   width = 6
        # ),

        chevreulBox(
            title = "Selected Cells",
            tableSelectedui("diffex"),
            width = 12
        ),
        chevreulBox(
            title = "Custom Cluster 1", DT::DTOutput(ns("cc1")),
            width = 6
        ), chevreulBox(
            title = "Custom Cluster 2", DT::DTOutput(ns("cc2")),
            width = 6
        ),
    )
}


#' Title
#'
#' @param input
#'
#' @return
#' @export
#'
#' @examples
cells_selected <- function(input) {
    if (identical(input, character(0))) {
        "Please selected desired cells by clicking on the table"
    } else {
        NULL
    }
}

#' Differential Expression
#'
#' @param input
#' @param output
#' @param session
#' @param seu
#' @param featureType
#' @param selected_cells
#' @param tests
#'
#' @return
#' @export
#'
#' @examples
diffex <- function(input, output, session, seu, featureType, selected_cells, tests = c("t-test" = "t", "wilcoxon rank-sum test" = "wilcox", "Likelihood-ratio test (bimodal)" = "bimod", "MAST" = "MAST")) {
    ns <- session$ns

    assay <- reactive({
        req(seu())
        if ("integrated" %in% names(seu()@assays)) {
            assay <- "integrated"
        } else {
            assay <- "gene"
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
        req(seu())
        d <- plotly::event_data("plotly_selected")
        if (is.null(d)) {
            msg <- "Click and drag events (i.e. select/lasso) appear here (double-click to clear)"
            return(d)
        } else {
            # selected_cells <- colnames(seu())[as.numeric(d$key)]
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

    output$cc1 <- DT::renderDT({
        req(custom_cluster1())
        selected_meta <- data.frame(seu()[[]][custom_cluster1(), ])
        DT::datatable(selected_meta,
            extensions = "Buttons",
            options = list(dom = "Bft", buttons = c(
                "copy",
                "csv"
            ), scrollX = "100px", scrollY = "400px")
        )
    })
    output$cc2 <- DT::renderDT({
        req(custom_cluster2())
        selected_meta <- data.frame(seu()[[]][custom_cluster2(), ])
        DT::datatable(selected_meta,
            extensions = "Buttons",
            options = list(dom = "Bft", buttons = c(
                "copy",
                "csv"
            ), scrollX = "100px", scrollY = "400px")
        )
    })

    de_results <- eventReactive(input$diffex, {
        if (input$diffex_scheme == "louvain") {
            run_seurat_de(seu(), input$cluster1, input$cluster2,
                resolution = input$seuratResolution, diffex_scheme = "louvain", input$featureType, tests = input$diffex_method
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
            run_seurat_de(seu(), cluster1, cluster2,
                input$customResolution,
                diffex_scheme = "feature", input$featureType, tests = input$diffex_method
            )
        }
    })

    output$DT1 <- DT::renderDT(de_results()[[input$diffex_method]],
        extensions = "Buttons", options = list(
            dom = "Bfptr",
            buttons = c("copy", "csv"), scrollX = "100px", scrollY = "600px"
        ), class = "display"
    )


    Volcano <- reactive({
        de_results()[[input$diffex_method]] %>%
            dplyr::distinct(symbol, .keep_all = TRUE) %>%
            tibble::column_to_rownames("symbol") %>%
            EnhancedVolcano::EnhancedVolcano(
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
            ggplot2::ggsave(file, Volcano() + ggpubr::theme_pubr(base_size = 20, x.text.angle = 45), width = 16, height = 12)
        }
    )

    cluster_list <- reactive({
        if (input$diffex_scheme == "louvain") {
            seu_meta <- seu()[[paste0(DefaultAssay(seu()), "_snn_res.", input$seuratResolution)]]
            cluster1_cells <- rownames(seu_meta[seu_meta == input$cluster1, , drop = FALSE])
            cluster2_cells <- rownames(seu_meta[seu_meta == input$cluster2, , drop = FALSE])
            list(cluster1 = cluster1_cells, cluster2 = cluster2_cells)
        } else if (input$diffex_scheme == "feature") {
            list(cluster1 = custom_cluster1(), cluster2 = custom_cluster2())
        }
    })

    return(list(cluster_list = cluster_list, de_results = de_results))
}


#' Find Markers UI
#'
#' @param id
#'
#' @return
#' @export
#'
#' @examples
findMarkersui <- function(id) {
    ns <- NS(id)
    tagList(
        chevreulBox(
            title = "Find Markers",
            uiOutput(ns("dplottype")),
            sliderInput(ns("resolution2"), label = "Resolution of clustering algorithm (affects number of clusters)", min = 0.2, max = 2, step = 0.2, value = 0.6),
            numericInput(ns("num_markers"), "Select Number of Markers to Plot for Each Value", value = 5, min = 2, max = 20),
            uiOutput(ns("valueSelect")),
            radioButtons(ns("markerMethod"), "Method of Marker Selection", choices = c("presto", "genesorteR"), selected = "presto", inline = TRUE),
            sliderInput(ns("pValCutoff"), "P Value cutoff", min = 0.01, max = 1, value = 1),
            selectizeInput(ns("dotFeature"), "Feature for Marker Plot", choices = NULL),
            actionButton(ns("plotDots"), "Plot Markers!"),
            downloadButton(ns("downloadMarkerTable"), "Download Markers!"),
            checkboxInput(ns("uniqueMarkers"), "Make Markers Unique", value = FALSE),
            checkboxInput(ns("hidePseudo"), "Hide Pseudogenes", value = TRUE),
            plotly::plotlyOutput(ns("markerplot"), height = 800),
            width = 6
        )
    ) %>%
        default_helper(type = "markdown", content = "findMarkers")
}

#' Find Markers
#'
#' @param input
#' @param output
#' @param session
#' @param seu
#'
#' @return
#' @export
#'
#' @examples
findMarkers <- function(input, output, session, seu, plot_types, featureType) {
    ns <- session$ns

    observe({
        req(seu())
        updateSelectizeInput(session, "dotFeature", choices = names(seu()@assays), selected = "gene", server = TRUE)
    })

    output$dplottype <- renderUI({
        req(seu())
        # selected_plot <- ifelse(is.null(selected_plot()), "louvain",
        #                         selected_plot())
        selectizeInput(ns("plottype"), "Variable to Plot",
            choices = purrr::flatten_chr(plot_types()),
            selected = "louvain", multiple = TRUE
        )
    })

    assay <- reactive({
        req(seu())
        if ("integrated" %in% names(seu()@assays)) {
            assay <- "integrated"
        } else {
            assay <- "gene"
        }
    })

    metavar <- reactive({
        req(input$plottype)

        if (input$plottype == "louvain") {
            metavar <- paste0(assay(), "_snn_res.", input$resolution2)
        } else {
            metavar <- input$plottype
        }
    })

    output$valueSelect <- renderUI({
        req(seu())
        req(metavar())

        choices <- levels(seu()[[]][[metavar()]])

        selectizeInput(ns("displayValues"), "Values to display", multiple = TRUE, choices = choices)
    })

    # observe({
    #   req(assay())
    #   Seurat::DefaultAssay(seu()) <- "gene"
    # })

    marker_plot_return <- eventReactive(input$plotDots, {
        plot_markers(seu(), metavar = metavar(), num_markers = input$num_markers, selected_values = input$displayValues, marker_method = input$markerMethod, seurat_assay = input$dotFeature, featureType = featureType(), hide_pseudo = input$hidePseudo, unique_markers = input$uniqueMarkers, p_val_cutoff = input$pValCutoff, return_plotly = TRUE)
    })

    output$markerplot <- plotly::renderPlotly({
        # req(input$displayClusters)
        marker_plot_return()$plot
    })

    output$downloadMarkerTable <- downloadHandler(
        filename = function() {
            paste(metavar(), "_markers.csv", sep = "")
        },
        content = function(file) {
            write_csv(marker_plot_return()$markers, file)
        }
    )
}

#' Plot Read Count UI
#'
#' @param id
#' @param plot_types
#'
#' @return
#' @export
#'
#' @examples
plotReadCountui <- function(id) {
    ns <- NS(id)
    chevreulBox(
        title = "Histogram (Read Counts, etc.)",
        uiOutput(ns("metavarui")),
        uiOutput(ns("colorbyui")),
        sliderInput(ns("resolution"), "Resolution of clustering algorithm (affects number of clusters)",
            min = 0.2, max = 2, step = 0.2, value = 0.6
        ),
        radioButtons(ns("yScale"), "Y axis transforamtion", choices = c("log", "linear"), selected = "linear", inline = TRUE),
        plotly::plotlyOutput(ns("rcplot"), height = 500),
        collapsed = TRUE
    )
}

#' Plot Read Count
#'
#' @param input
#' @param output
#' @param session
#' @param seu
#' @param plot_types
#'
#' @return
#' @export
#'
#' @examples
plotReadCount <- function(input, output, session, seu, plot_types) {
    ns <- session$ns

    output$colorbyui <- renderUI({
        req(seu())
        shiny::selectInput(ns("colorby"), "Variable to Color the Plot by",
            choices = purrr::flatten_chr(plot_types()), selected = c("louvain"), multiple = FALSE
        )
    })

    output$metavarui <- renderUI({
        req(seu())
        shiny::selectInput(ns("metavar"), "Variable for x-axis",
            choices = purrr::flatten_chr(plot_types()), selected = c("nCount_RNA"), multiple = FALSE
        )
    })

    output$rcplot <- plotly::renderPlotly({
        req(seu())
        req(input$colorby)

        if (input$colorby == "louvain") {
            if ("integrated" %in% names(seu()@assays)) {
                assay <- "integrated"
            } else {
                assay <- "gene"
            }

            louvain_resolution <- paste0(assay, "_snn_res.", input$resolution)
            plot_readcount(seu(), metavar = input$metavar, color.by = louvain_resolution, yscale = input$yScale, return_plotly = TRUE)
        } else if (input$colorby %in% purrr::flatten_chr(plot_types())) {
            plot_readcount(seu(), metavar = input$metavar, color.by = input$colorby, yscale = input$yScale, return_plotly = TRUE)
        }
    })
}

#' Cell Cycle Score UI
#'
#' @param id
#'
#' @return
#' @export
#'
#' @examples
ccScoreui <- function(id) {
    ns <- NS(id)
    tagList()
}

#' Cell Cycle Score
#'
#' @param input
#' @param output
#' @param session
#'
#' @return
#' @export
#'
#' @examples
ccScore <- function(input, output, session) {
    ns <- session$ns
    output$rplot1 <- renderPlot({
        req(seu())
        plot_ridge(seu(), features = input$feature)
    })
    plotOutput("rplot1", height = 750)
}

#' Plot All Transcripts UI Module
#'
#' @param id
#'
#' @return
#' @export
#'
#' @examples
allTranscriptsui <- function(id) {
    ns <- NS(id)
    tagList(
        default_helper(
            chevreulBox(
                title = "Transcript Expression per Gene",
                selectizeInput(ns("embeddingGene"), "Gene or transcript expression by which to color the plot; eg. 'RXRG'", choices = NULL, selected = NULL),
                selectizeInput(ns("transcriptSelect"), "Transcript to Plot", choices = NULL),
                downloadButton(ns("downloadPlot"), "Download Transcript Plots"),
                selectizeInput(ns("embedding"), "Embedding", choices = NULL, selected = NULL),
                plotly::plotlyOutput(ns("transcriptPlot")),
                # uiOutput(ns("plotlys")),
                width = 6
            ),
            type = "markdown", content = "allTranscripts"
        ),
        default_helper(
            chevreulBox(
                title = "Transcript Expression per Gene",
                selectizeInput(ns("compositionGene"), "Gene or transcript expression by which to color the plot; eg. 'RXRG'", choices = NULL, selected = NULL),
                selectizeInput(ns("groupby"), "Group by:", choices = NULL, selected = NULL),
                actionButton(ns("plotComposition"), "Plot transcript composition"),
                checkboxInput(ns("standardizeExpression"), "Standardize Expression", value = FALSE),
                checkboxInput(ns("dropZero"), "Drop Zero Values", value = FALSE),
                plotly::plotlyOutput(ns("compositionPlot")),
                DT::DTOutput(ns("compositionDT")),
                width = 6
            ),
            type = "markdown", content = "allTranscripts"
        )
    )
}

#' Plot All Transcripts Server
#'
#' @param input
#' @param output
#' @param session
#' @param seu
#' @param featureType
#'
#' @return
#' @export
#'
#' @examples
allTranscripts <- function(input, output, session, seu,
    featureType, organism_type) {
    ns <- session$ns

    observe({
        req(seu())
        updateSelectizeInput(session, "compositionGene", choices = rownames(seu()[["gene"]]), selected = "RXRG", server = TRUE)
        updateSelectizeInput(session, "embeddingGene", choices = rownames(seu()[["gene"]]), selected = "RXRG", server = TRUE)

        formatted_col_names <- colnames(seu()@meta.data) %>%
            make_chevreul_clean_names()

        updateSelectizeInput(session, "groupby", choices = formatted_col_names, selected = "batch", server = TRUE)
    })

    transcripts <- reactive({
        req(seu())
        if ("transcript" %in% names(seu()@assays)) {
            get_transcripts_from_seu(seu(), input$embeddingGene, organism = organism_type())
        }
    })

    observe({
        req(seu())
        req(transcripts())
        updateSelectizeInput(session, "embedding", choices = c("pca", "tsne", "umap"), selected = "umap", server = TRUE)
        updateSelectizeInput(session, "transcriptSelect", choices = transcripts(), server = TRUE)
    })

    composition_plot <- eventReactive(input$plotComposition, {
        plot_transcript_composition(seu(), gene_symbol = input$compositionGene, group.by = input$groupby, standardize = input$standardizeExpression, drop_zero = input$dropZero)
    })

    output$compositionPlot <- plotly::renderPlotly({
        composition_plot()$plot %>%
            plotly::ggplotly(height = 400) %>%
            plotly_settings() %>%
            plotly::toWebGL() %>%
            # plotly::partial_bundle() %>%
            identity()
    })

    output$compositionDT <- DT::renderDT({
        DT::datatable(composition_plot()$data,
            extensions = "Buttons",
            options = list(dom = "Bft", buttons = c(
                "copy",
                "csv"
            ), scrollX = "100px", scrollY = "400px")
        )
    })

    pList <- reactive({
        req(transcripts())
        pList <- plot_all_transcripts(seu(), transcripts(), input$embedding, from_gene = FALSE, combine = FALSE)
    })

    output$transcriptPlot <- plotly::renderPlotly({
        pList()[[input$transcriptSelect]] %>%
            plotly::ggplotly(height = 400) %>%
            plotly_settings() %>%
            plotly::toWebGL()
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

#' RNA Velocity UI Module
#'
#' @param id
#'
#' @return
#' @export
#'
#' @examples
plotVelocityui <- function(id) {
    ns <- NS(id)
    tagList(
        chevreulBox(
            title = "Calculate Velocity",
            width = 4,
            textOutput(ns("velocityFlag")),
            radioButtons(ns("velocityMode"), "Velocity Mode", choices = c("deterministic (velocyto)" = "deterministic", "stochastic" = "stochastic", "dynamical" = "dynamical")),
            actionButton(ns("calc_velocity"), "calculate velocity"),
            textOutput(ns("scveloMessages")),
        ),
        chevreulBox(
            title = "Plot Veloctiy on Embedding",
            width = 12,
            selectizeInput(ns("embedding"), "dimensional reduction method",
                choices = c("pca", "tsne", "umap"),
                selected = "umap"
            ),
            selectizeInput(ns("varSelect"), "Color by Variable", choices = NULL, multiple = FALSE),
            sliderInput(ns("resolution"), "Resolution of clustering algorithm (affects number of clusters)", min = 0.2, max = 2, step = 0.2, value = 0.6),
            radioButtons(ns("plotFormat"), "velocity format", choices = c("arrow", "stream"), selected = "arrow", inline = TRUE),
            actionButton(ns("plot_velocity_embedding"), "plot velocity on embedding"),
            downloadButton(ns("downloadEmbeddingPlot"), label = "Download Plot"),
            imageOutput(ns("velocityEmbeddingPlot"), height = "800px")
        ),
        chevreulBox(
            title = "Plot Velocity and Expression",
            width = 12,
            selectizeInput(ns("geneSelect"), "Select a Gene", choices = NULL, selected = NULL, multiple = TRUE),
            actionButton(ns("plot_velocity_expression"), "plot velocity and expression"),
            downloadButton(ns("downloadExpressionPlot"), label = "Download Plot"),
            imageOutput(ns("velocityExpressionPlot"), height = "500px")
        ) %>%
            default_helper(type = "markdown", content = "plotVelocity")
    )
}

#' RNA Velocity Server Module
#'
#' @param input
#' @param output
#' @param session
#' @param seu
#' @param loom_path
#'
#' @return
#' @export
#'
#' @examples
plotVelocity <- function(input, output, session, seu, loom_path) {
    ns <- session$ns

    print("running scvelo")

    observe({
        req(seu())
        updateSelectizeInput(session, "varSelect", choices = colnames(seu()[[]]), selected = "batch", server = TRUE)
        updateSelectizeInput(session, "geneSelect", choices = rownames(seu()[["gene"]]), selected = "RXRG", server = TRUE)
    })

    # reactive val adata ------------------------------

    adata <- reactiveVal()

    observeEvent(input$calc_velocity, {
        req(seu())
        withCallingHandlers(
            {
                shinyjs::html("scveloMessages", "")
                message("Beginning")

                if ("integrated" %in% names(seu()@assays)) {
                    assay <- "integrated"
                } else {
                    assay <- "gene"
                }

                adata <- prep_scvelo(seu(), loom_path, velocity_mode = input$velocityMode)

                adata(adata)
                message("scvelo Complete!")
            },
            message = function(m) {
                shinyjs::html(id = "scveloMessages", html = paste0("Running scvelo: ", m$message), add = FALSE)
            }
        )
    })


    velocity_flag <- eventReactive(input$calc_velocity, {
        req(adata())
        "Velocity Calculated for this dataset"
    })

    output$velocityFlag <- renderText({
        req(adata())
        velocity_flag()
    })

    observe({
        req(adata())

        if ("integrated" %in% names(seu()@assays)) {
            assay <- "integrated"
        } else {
            assay <- "gene"
        }

        cluster_resolution <- paste0(assay, "_snn_res.", input$resolution)

        plot_scvelo(adata(), group.by = input$varSelect, plot_method = input$plotFormat)
        fig <- pyplot$gcf()
        # fig$savefig("velocity_embedding.pdf")
        fig$savefig("velocity_embedding.svg")
    })

    observe({
        req(adata())
        req(input$geneSelect)

        if ("integrated" %in% names(seu()@assays)) {
            assay <- "integrated"
        } else {
            assay <- "gene"
        }
        scvelo_expression(adata(), features = input$geneSelect)

        fig <- pyplot$gcf()
        # fig$savefig("velocity_expression.pdf")
        fig$savefig("velocity_expression.svg")
    })

    output$downloadEmbeddingPlot <- downloadHandler(
        filename = function() {
            paste("velocity_embedding", ".svg", sep = "")
        },
        content = function(file) {
            file.copy("velocity_embedding.svg", file, overwrite = TRUE)
        }
    )

    output$downloadExpressionPlot <- downloadHandler(
        filename = function() {
            paste("velocity_expression", ".svg", sep = "")
        },
        content = function(file) {
            file.copy("velocity_expression.svg", file, overwrite = TRUE)
        }
    )

    expression_path <- eventReactive(input$plot_velocity_expression, {
        "velocity_expression.svg"
    })

    embedding_path <- eventReactive(input$plot_velocity_embedding, {
        "velocity_embedding.svg"
    })

    output$velocityEmbeddingPlot <- renderImage(
        {
            # req(velocityExpressionPlot())
            # Get width and height of image output
            width <- session$clientData$output_image_width
            height <- session$clientData$output_image_height

            # Return a list containing information about the image
            list(
                src = embedding_path(),
                contentType = "image/svg+xml",
                width = 1200,
                height = 800,
                alt = "This is alternate text"
            )
        },
        deleteFile = FALSE
    )

    output$velocityExpressionPlot <- renderImage(
        {
            # req(velocityExpressionPlot())
            # Get width and height of image output
            width <- session$clientData$output_image_width
            height <- session$clientData$output_image_height

            # Return a list containing information about the image
            list(
                src = expression_path(),
                contentType = "image/svg+xml",
                width = 1500,
                height = 500,
                alt = "This is alternate text"
            )
        },
        deleteFile = FALSE
    )
}


#' Monocle UI Module
#'
#' @param id
#'
#' @return
#' @export
#'
#' @examples
monocleui <- function(id) {
    ns <- NS(id)
    tagList(
        fluidRow(
            chevreulBox(
                title = "Seurat Data",
                plotly::plotlyOutput(ns("seudimplot"), height = 500),
                width = 6
                # plotDimRedui(ns("plotdimred")
            ),
            chevreulBox(
                title = "Pseudotime Settings",
                actionButton(ns("subsetSeurat"), "Subset Seurat before Pseudotime Calculation"),
                actionButton(ns("calcCDS"), "Calculate Pseudotime"),
                sliderInput(ns("cdsResolution"), "Resolution of clustering algorithm (affects number of clusters)",
                    min = 0.2, max = 2, step = 0.2, value = 0.6
                ),
                actionButton(ns("subsetCells"), "Subset Monocle Object After Pseudotime Calculation"),
                uiOutput(ns("rootCellsui")),
                actionButton(ns("plotPseudotime"), "Calculate Pseudotime With Root Cells"),
                downloadButton(ns("downloadPT"), "Export Pseudotime"),
                checkboxInput(ns("flipPtime"), "Invert Pseudotime", value = TRUE),
                width = 6
            )
        ),
        chevreulBox(
            title = "Embedding Plot",
            selectizeInput(ns("plottype1"), "Variable to Plot", choices = c(Louvain = "louvain"), selected = "Louvain", multiple = TRUE),
            selectizeInput(ns("customFeature1"), "Gene or transcript expression by which to color the plot",
                choices = NULL, multiple = FALSE
            ),
            uiOutput(ns("moduleSelect1")),
            plotly::plotlyOutput(ns("monoclePlot1")),
            width = 6
        ),
        chevreulBox(
            title = "Embedding Plot",
            selectizeInput(ns("plottype2"), "Variable to Plot", choices = c(Louvain = "louvain"), selected = "Louvain", multiple = TRUE),
            selectizeInput(ns("customFeature2"), "gene or transcript on which to color the plot",
                choices = NULL, multiple = FALSE
            ),
            uiOutput(ns("moduleSelect2")),
            plotly::plotlyOutput(ns("monoclePlot2")),
            width = 6
        ),
        fluidRow(
            chevreulBox(
                title = "calculate pseudotime",
                radioButtons(ns("diffexFeature"), "Feature for differential expression", choices = c("gene", "transcript")),
                actionButton(ns("calcPtimeGenes"), "Find Pseudotime Correlated Genes"),
                sliderInput(ns("qvalThreshold"), "Set q value threshold for module calculation", min = 0.01, 0.1, value = 0.05, step = 0.01),
                textOutput("pseudotimeMessages"),
                uiOutput(ns("partitionSelect")),
                uiOutput(ns("genePlotQuery2")),
                DT::DTOutput(ns("ptimeGenesDT")),
                downloadButton(ns("downloadGenesDT"), "Download data as csv"),
                # uiOutput(ns("ptimeGenes")),
                width = 6
            ),
            chevreulBox(
                title = "Plot Feature Expression over Pseudotime",
                plotly::plotlyOutput(ns("ptimeGenesLinePlot")),
                width = 6,
                height = 650
            )
        ),
        chevreulBox(
            title = "Heatmap",
            uiOutput(ns("colAnnoVarui")),
            radioButtons(ns("heatmapRows"), "annotate heatmap rows by genes or modules?", choices = c("modules", "genes")),
            downloadButton(ns("downloadPlot"), "Download Heatmap"),
            downloadButton(ns("downloadCds"), "Download celldataset"),
            plotOutput(ns("monocleHeatmap"), width = "800px", height = "1200px")
        ),
        chevreulBox(
            title = "Modules",
            plotOutput(ns("modulePlot")),
            div(DT::dataTableOutput(ns("moduleTable")), style = "font-size: 75%")
        )
    )
}

#' Monocle Server Module
#'
#' @param input
#' @param output
#' @param session
#' @param cds
#' @param seu
#' @param plot_types
#' @param resolution
#'
#' @return
#' @export
#'
#' @examples
monocle <- function(input, output, session, seu, plot_types, featureType,
    organism_type, reductions) {
    ns <- session$ns

    # markermarker
    w <- waiter::Waiter$new(ns("monocleHeatmap"),
        html = waiter::spin_loaders(id = 1, color = "black", style = "position:relative;margin:auto;"),
        color = waiter::transparent(.5)
    )

    output$colAnnoVarui <- renderUI({
        req(seu())

        selectizeInput(ns("colAnnoVar"), "Column Annotation(s)",
            choices = colnames(seu()[[]]), selected = "batch", multiple = TRUE
        )
    })

    cds_rvs <- reactiveValues(selected = c(traj = TRUE, ptime = FALSE, diff_features = FALSE))
    cds_plot_types <- reactiveVal(c(Pseudotime = "pseudotime", Module = "module"))
    myplot_types <- reactive({
        c(purrr::flatten_chr(plot_types()), cds_plot_types())
    })

    # to be able to subset, create a new copy of the seurat object

    seu_monocle <- reactiveVal()

    observe({
        req(seu())
        seu_monocle(seu())
    })

    louvain_resolution <- reactive({
        if ("integrated" %in% names(seu()@assays)) {
            assay <- "integrated"
        } else {
            assay <- "gene"
        }

        paste0(assay, "_snn_res.", input$cdsResolution)
    })

    seudimplot <- reactive({
        req(seu_monocle())

        plot_var(seu_monocle(), embedding = "umap", group = louvain_resolution(), return_plotly = TRUE)
    })

    output$seudimplot <- plotly::renderPlotly({
        seudimplot()
    })

    # callModule(plotDimRed, "plotdimred", seu, plot_types, featureType,
    #            organism_type, reductions)


    observeEvent(input$subsetSeurat, {
        req(seu_monocle())

        d <- plotly::event_data("plotly_selected", priority = "event")
        if (is.null(d)) {
            msg <- "Click and drag events (i.e. select/lasso) appear here (double-click to clear)"
            print(d)
        } else {
            print(d$key)
            print(d)
            subset_monocle <- seu_monocle()[, d$key]
            seu_monocle(subset_monocle)
        }
    })

    observeEvent(input$calcCDS, {
        req(seu_monocle())
        cds_rvs$selected <- c(traj = TRUE, ptime = FALSE, diff_features = FALSE)
        cds <- convert_seu_to_cds(seu(), resolution = input$cdsResolution)
        # cds <- convert_seu_to_cds(seu_monocle(), resolution = input$cdsResolution)
        cds <- cds[, colnames(cds) %in% colnames(seu_monocle())]

        cds <- threshold_monocle_genes(seu_monocle(), cds)

        cds <- learn_graph_by_resolution(cds, seu_monocle(),
            resolution = input$cdsResolution
        )
        updateSelectizeInput(session, "plottype1", selected = "louvain", choices = myplot_types())
        updateSelectizeInput(session, "customFeature1", choices = rownames(cds), server = TRUE)
        updateSelectizeInput(session, "plottype2", selected = "louvain", choices = myplot_types())
        updateSelectizeInput(session, "customFeature2", choices = rownames(cds), server = TRUE)
        cds_rvs$traj <- cds
    })

    selected_plot <- reactiveVal()

    output$monoclePlot1 <- plotly::renderPlotly({
        req(input$plottype1)
        req(cds_rvs$traj)
        w$show()
        print(cds_rvs$selected)
        if (input$plottype1 == "louvain") {
            cluster_resolution <- reactive({
                if (any(stringr::str_detect(colnames(colData(cds_rvs$traj)), "integrated"))) {
                    paste0("integrated", "_snn_res.", input$cdsResolution)
                } else {
                    paste0("gene", "_snn_res.", input$cdsResolution)
                }
            })
            plot_cds(cds_rvs$traj, color_cells_by = cluster_resolution())
        } else if (input$plottype1 == "pseudotime") {
            plot_pseudotime(cds_rvs$traj, color_cells_by = "pseudotime", resolution = input$cdsResolution)
        } else if (input$plottype1 == "feature") {
            plot_monocle_features(cds_rvs$traj, genes = input$customFeature1, monocle_heatmap()$agg_mat)
        } else if (input$plottype1 == "module") {
            print(monocle_heatmap()$module_table)
            print(input$plotModule1)
            genes <- monocle_heatmap()$module_table %>%
                filter(module %in% input$plotModule1) %>%
                dplyr::mutate(module = factor(module))
            plot_monocle_features(cds_rvs$traj, genes = genes, monocle_heatmap()$agg_mat)
        } else {
            plot_cds(cds_rvs$traj, color_cells_by = input$plottype1)
        }
    })

    output$monoclePlot2 <- plotly::renderPlotly({
        req(input$plottype2)
        req(cds_rvs$traj)
        w$show()
        print(cds_rvs$selected)
        if (input$plottype2 == "louvain") {
            cluster_resolution <- reactive({
                if (any(stringr::str_detect(colnames(colData(cds_rvs$traj)), "integrated"))) {
                    paste0("integrated", "_snn_res.", input$cdsResolution)
                } else {
                    paste0("gene", "_snn_res.", input$cdsResolution)
                }
            })
            plot_cds(cds_rvs$traj, color_cells_by = cluster_resolution())
        } else if (input$plottype2 == "pseudotime") {
            plot_pseudotime(cds_rvs$traj, color_cells_by = "pseudotime", resolution = input$cdsResolution)
        } else if (input$plottype2 == "feature") {
            plot_monocle_features(cds_rvs$traj, genes = input$customFeature2, monocle_heatmap()$agg_mat)
        } else if (input$plottype2 == "module") {
            print(monocle_heatmap()$module_table)
            print(input$plotModule2)

            genes <- monocle_heatmap()$module_table %>%
                filter(module %in% input$plotModule2) %>%
                dplyr::mutate(module = factor(module))
            plot_monocle_features(cds_rvs$traj, genes = genes, monocle_heatmap()$agg_mat)
        } else {
            plot_cds(cds_rvs$traj, color_cells_by = input$plottype2)
        }
    })

    cdsbrush <- reactive({
        req(cds_rvs$traj)
        d <- plotly::event_data("plotly_selected")
        if (is.null(d)) {
            msg <- "Click and drag events (i.e. select/lasso) appear here (double-click to clear)"
            return(d)
        } else {
            # selected_cells <- colnames(cds_rvs$traj)[as.numeric(d$key)]
            d$key
        }
    })

    observeEvent(input$subsetCells, {
        req(cds_rvs$traj)
        print(cdsbrush())
        cds_rvs$traj <- cds_rvs$traj[, cdsbrush()]
    })

    output$rootCellsui <- renderUI({
        selectizeInput(ns("rootCells"), "Choose Root Cells", choices = c("Choose Root Cells" = "", colnames(cds_rvs$traj)), multiple = TRUE)
    })

    exported_pseudotime <- reactiveVal()

    observeEvent(input$plotPseudotime, {
        req(cds_rvs$traj)
        req(input$rootCells)
        cds_rvs$traj <- monocle3::order_cells(cds_rvs$traj, root_cells = input$rootCells)

        # # select only first partition
        # cds_rvs$traj <- cds_rvs$traj[, monocle3::partitions(cds_rvs$traj) == 1]


        if (input$flipPtime) {
            cds_rvs$traj <- flip_pseudotime(cds_rvs$traj)
        }
        updateSelectizeInput(session, "plottype1", selected = "pseudotime", choices = myplot_types())
        updateSelectizeInput(session, "plottype2", selected = "pseudotime", choices = myplot_types())
        cds_rvs$selected <- c(traj = FALSE, ptime = TRUE, diff_features = FALSE)
        # markermarker

        exported_pseudotime(export_pseudotime(cds_rvs$traj, input$rootCells))
    })


    output$downloadPT <- downloadHandler(
        filename = function() {
            paste("pseudotime-", Sys.Date(), ".csv", sep = "")
        },
        content = function(file) {
            readr::write_csv(exported_pseudotime(), file)
        }
    )

    # markermarker
    observeEvent(input$calcPtimeGenes, {
        req(input$diffexFeature)
        if (req(cds_rvs$selected["ptime"])) {
            # markermarker
            # cds_rvs$traj  <- swap_counts_from_feature(cds_rvs$traj, input$diffexFeature)

            showModal(modalDialog(
                title = "Calculating Pseudotime Correlated Features",
                "This may take a few minutes!"
            ))
            cds_rvs$traj@metadata[["diff_features"]] <- monocle3::graph_test(cds_rvs$traj, neighbor_graph = "principal_graph", cores = 4, expression_family = "negbinom")

            cds_rvs$selected <- c(traj = FALSE, ptime = FALSE, diff_features = TRUE)

            removeModal()
        }
    })

    cds_pr_test_res <- reactive({
        if (req(cds_rvs$selected["diff_features"])) {
            cds_rvs$traj@metadata$diff_features %>%
                # subset(q_value < 0.05) %>%
                dplyr::arrange(q_value) %>%
                dplyr::select(-status) %>%
                # dplyr::filter %>%
                identity()
        }
    })

    observe({
        req(cds_pr_test_res())
        if (req(cds_rvs$selected["diff_features"])) {
            output$genePlotQuery2 <- renderUI({
                selectizeInput(ns("genePlotQuery1"), "Pick Gene to Plot on Pseudotime", choices = rownames(cds_pr_test_res()), multiple = TRUE, selected = rownames(cds_pr_test_res())[1])
            })

            output$partitionSelect <- renderUI({
                selectizeInput(ns("partitions"), "Select a Partition to Plot", choices = levels(monocle3::partitions(cds_rvs$traj)), multiple = FALSE)
            })
        }
    })

    observe({
        req(cds_pr_test_res())
        req(input$genePlotQuery1)
        if (req(cds_rvs$selected["diff_features"])) {
            output$ptimeGenesLinePlot <- plotly::renderPlotly({
                genes_in_pseudotime <- prep_plot_genes_in_pseudotime(cds_rvs$traj, input$genePlotQuery1, input$cdsResolution)
                genes_in_pseudotime <-
                    genes_in_pseudotime %>%
                    plotly::ggplotly(height = 600) %>%
                    plotly_settings() %>%
                    plotly::toWebGL() %>%
                    # plotly::partial_bundle() %>%
                    identity()
            })

            output$ptimeGenesDT <- DT::renderDT({
                DT::datatable(cds_pr_test_res(),
                    extensions = "Buttons",
                    options = list(dom = "Bftp", buttons = c("copy", "csv"), scrollX = "100px", scrollY = "400px", pageLength = 200, paging = TRUE)
                )
            })

            output$downloadGenesDT <- downloadHandler(
                filename = function() {
                    paste("diffex_ptime-", Sys.Date(), ".csv", sep = "")
                },
                content = function(file) {
                    write.csv(cds_pr_test_res(), file)
                }
            )
        }
    })

    monocle_heatmap <- reactive({
        req(cds_rvs$traj)
        req(input$colAnnoVar)

        heatmap_genes <- cds_pr_test_res() %>%
            dplyr::filter(q_value < input$qvalThreshold)

        monocle_module_heatmap(cds_rvs$traj, rownames(heatmap_genes), input$cdsResolution, collapse_rows = input$heatmapRows, group.by = input$colAnnoVar)
    }) %>%
        bindCache(cds_rvs$traj, input$cdsResolution, input$heatmapRows)

    module_choices <- reactive({
        module_choices <- as.character(unique(monocle_heatmap()$module_table$module))
        # names(module_choices) <- paste("Module", module_choices)
    })

    output$moduleSelect1 <- renderUI({
        selectizeInput(ns("plotModule1"), "gene module to plot (if computed)", choices = module_choices(), multiple = TRUE)
    })
    output$moduleSelect2 <- renderUI({
        selectizeInput(ns("plotModule2"), "gene module to plot (if computed)", choices = module_choices(), multiple = TRUE)
    })

    observe({
        output$monocleHeatmap <- renderPlot({
            monocle_heatmap()$module_heatmap
        })

        output$moduleTable <- DT::renderDT({
            DT::datatable(monocle_heatmap()$module_table,
                extensions = "Buttons",
                options = list(dom = "Bft", buttons = c(
                    "copy",
                    "csv"
                ), scrollX = "100px", scrollY = "400px")
            )
        })

        output$modulePlot <- renderPlot({
            ggplot(monocle_heatmap()$module_table, aes(dim_1, dim_2, color = module)) +
                geom_point()
        })
    })

    output$downloadPlot <- downloadHandler(
        filename = function() {
            paste("heatmap", ".pdf", sep = "")
        },
        content = function(file) {
            ggsave(file, ggplotify::as.ggplot(monocle_heatmap()$module_heatmap), width = 16, height = 12)
        }
    )

    output$downloadCds <- downloadHandler(
        filename = function() {
            paste("cds", ".rds", sep = "")
        },
        content = function(file) {
            saveRDS(cds_rvs$traj, file)
        }
    )
}


#' Title
#'
#' @param id
#'
#' @return
#' @export
#'
#' @examples
pathwayEnrichmentui <- function(id) {
    ns <- NS(id)
    chevreulBox(
        title = "Enriched pathways by cluster",
        tagList(
            actionButton(ns("calcPathwayEnrichment"), "Calculate Pathway Enrichment"),
            uiOutput(ns("enriched_pathways_by_cluster_select_source_UI")),
            uiOutput(ns("enriched_pathways_by_cluster_UI"))
        ),
        width = 12
    )
}

#' pathway enrichment
#'
#' @param input
#' @param output
#' @param session
#'
#' @return
#' @export
#'
#' @examples
pathwayEnrichment <- function(input, output, session, seu, featureType) {
    ns <- session$ns

    ## ----------------------------------------------------------------------------##
    ## Tab: Enriched pathways
    ## ----------------------------------------------------------------------------##

    ## ----------------------------------------------------------------------------##
    ## Clusters.
    ## ----------------------------------------------------------------------------##

    enriched_pathways <- eventReactive(input$calcPathwayEnrichment, {
        req(seu())
        if (featureType() == "gene") {
            enriched_seu <- tryCatch(getEnrichedPathways(seu()), error = function(e) e)
            enrichr_available <- !any(class(enriched_seu) == "error")
            if (enrichr_available) {
                seu <- enriched_seu
            }
        }

        seu()@misc$enriched_pathways
    })

    # UI element: choose source for pathway enrichement results (currently Enrichr or GSVA)
    output$enriched_pathways_by_cluster_select_source_UI <- renderUI({
        req(seu())
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
        req(seu())
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
                        DT::dataTableOutput(ns("enriched_pathways_by_cluster_table_present"))
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
        req(seu())
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
            dplyr::filter(cluster == input$enriched_pathways_by_cluster_select_cluster) %>%
            dplyr::pull(db) %>%
            intersect(., levels(.))
        selectInput(
            ns("enriched_pathways_by_cluster_select_db"),
            label = NULL,
            choices = choices
        )
    })

    # table
    output$enriched_pathways_by_cluster_table_present <- DT::renderDataTable(server = FALSE, {
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
#' @param id
#'
#' @return
#' @export
#'
#' @examples
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
#' @param input
#' @param output
#' @param session
#' @param seu
#'
#' @return
#' @export
#'
#' @examples
techInfo <- function(input, output, session, seu) {
    ns <- session$ns
    ## ----------------------------------------------------------------------------##
    ## Tab: Analysis info.
    ## ----------------------------------------------------------------------------##

    misc <- reactive({
        req(seu())
        Seurat::Misc(seu())
    })

    observe({
        # general info
        output$sample_info_general <- renderText({
            info <- paste0(
                "<strong><u>General</u></strong>",
                "<ul>",
                "<li><b>Date of analysis:</b> ",
                misc()$experiment$date_of_analysis,
                "<li><b>Date of export:</b> ",
                misc()$experiment$date_of_export,
                "<li><b>Experiment name:</b> ",
                misc()$experiment$experiment_name,
                "<li><b>Organism:</b> ",
                misc()$experiment$organism,
                "</ul>",
                "<strong><u>Parameters</u></strong>",
                "<ul>",
                "<li><b>Discard genes in fewer than X cells:</b> ",
                misc()$experiment$parameters$discard_genes_expressed_in_fewer_cells_than,
                "<li><b>Keep mitochondrial genes:</b> ",
                misc()$experiment$parameters$keep_mitochondrial_genes,
                "<li><b>Min/max # of UMI:</b> ",
                paste0(
                    misc()$experiment$filtering$UMI_min, " / ",
                    misc()$experiment$filtering$UMI_max
                ),
                "<li><b>Min/max # of expressed genes:</b> ",
                paste0(
                    misc()$experiment$filtering$genes_min, " / ",
                    misc()$experiment$filtering$genes_max
                ),
                "<li><b>Cluster resolution: </b>",
                paste(misc()$experiment$parameters$cluster_resolution, collapse = ","),
                "<li><b>Number of principal components: </b>",
                misc()$experiment$parameters$number_PCs,
                "<li><b>Variables to regress: </b>",
                misc()$experiment$parameters$variables_to_regress_out,
                "<li><b>tSNE perplexity: </b>",
                misc()$experiment$parameters$tSNE_perplexity,
                "</ul>",
                "<strong><u>Gene lists</u></strong>",
                "<ul>",
                # "<li><b>Mitochondrial genes:</b> ",
                # paste0(mito_features[[misc()$experiment$organism]][["gene"]], collapse = ", "),
                # "<li><b>Ribosomal genes:</b> ",
                # paste0(ribo_features[[misc()$experiment$organism]][["gene"]], collapse = ", "),
                "<li><b>S phase genes:</b> ",
                paste0(cc.genes$s.genes, collapse = ", "),
                "<li><b>G2M phase genes:</b> ",
                paste0(cc.genes$g2m.genes, collapse = ", "),
                "</ul>",
                "<strong><u>Marker genes</u></strong>",
                "<ul>",
                # "<li><b>Only positive:</b> ",
                # misc()$marker_genes$parameters$only_positive,
                # "<li><b>Fraction of cells in group of interest that must express marker gene:</b> ",
                # misc()$marker_genes$parameters$minimum_percentage,
                # "<li><b>LogFC threshold:</b> ",
                # misc()$marker_genes$parameters$logFC_threshold,
                "<li><b>p-value threshold:</b> ",
                "0.05",
                # misc()$marker_genes$parameters$p_value_threshold,
                "</ul>",
                "<strong><u>Pathway enrichment</u></strong>",
                "<ul>",
                "<li><b>Enrichr:</b>",
                "<ul>",
                "<li><b>Databases:</b> ",
                paste0(misc()$enriched_pathways$enrichr$parameters$databases, collapse = ", "),
                "<li><b>Adj. p-value cut-off:</b> ",
                misc()$enriched_pathways$enrichr$parameters$adj_p_cutoff,
                "<li><b>Max. terms:</b> ",
                misc()$enriched_pathways$enrichr$parameters$max_terms,
                "</ul>",
                "</ul>"
            )
            info_R_raw <- misc()$experiment$technical_info$R
            info_R <- c()
            for (i in 1:length(info_R_raw)) {
                info_R <- paste(info_R, "<br>", info_R_raw[i])
            }
            paste0(
                info,
                "<strong><u>Technical info (package versions)</u></strong>",
                "<ul>",
                "<li><strong>chevreul version:</strong> ",
                misc()$experiment$technical_info$chevreul_version,
                "<li><strong>Seurat version:</strong> ",
                misc()$technical_info$seurat_version,
                "<li><strong>Session info:</strong> ",
                "</ul>",
                "<pre>",
                info_R,
                "</pre>"
            )
        })

        # R info
        output$sample_info_R <- renderPrint({
            if (!is.null(misc()$technical_info$R)) {
                capture.output(misc()$technical_info$R)
            } else {
                print("Not available")
            }
        })
    })
}

#' Title
#'
#' @param id
#'
#' @return
#' @export
#'
#' @examples
plotCoverage_UI <- function(id) {
    ns <- NS(id)
    tagList(
        chevreulBox(
            title = "Plot Coverage",
            selectizeInput(ns("geneSelect"), "Select a Gene", choices = NULL, selected = "RXRG", multiple = FALSE),
            selectizeInput(ns("varSelect"), "Color by Variable", choices = NULL, multiple = FALSE),
            actionButton(ns("plotCoverage"), "Plot Coverage"),
            downloadButton(ns("downloadPlot"), "Download Coverage Plot"),
            uiOutput(ns("displayvaluesui")),
            br(),
            dropdownButton(
                ns("coveragePlotSettings"),
                checkboxInput(ns("collapseIntrons"), "Collapse Introns", value = TRUE),
                checkboxInput(ns("meanCoverage"), "Summarize Coverage to Mean", value = TRUE),
                checkboxInput(ns("summarizeTranscripts"), "Summarize transcript models to gene", value = FALSE),
                radioButtons(ns("yScale"), "Scale Y Axis", choices = c("absolute", "log10"), selected = "log10"),
                numericInput(ns("start"), "start coordinate", value = NULL),
                numericInput(ns("end"), "end coordinate", value = NULL)
            ),
            DT::DTOutput(ns("coverageTable")),
            plotOutput(ns("coveragePlot"), height = "1500px"),
            width = 12
        )
    )
}


#' Plot Coverage Module
#'
#' @param input
#' @param output
#' @param session
#' @param seu
#' @param plot_types
#' @param bigwig_dir
#' @param organism_type
#'
#' @return
#' @export
#'
#' @examples
plotCoverage <- function(input, output, session, seu, plot_types, proj_dir, organism_type = "human", bigwig_db = "~/.cache/chevreul/bw-files.db") {
    ns <- session$ns

    w <- waiter::Waiter$new(ns("coveragePlot"),
        html = waiter::spin_loaders(id = 1, color = "black", style = "position:relative;margin:auto;"),
        color = waiter::transparent(.5)
    )

    observe({
        req(seu())
        updateSelectizeInput(session, "geneSelect", choices = rownames(seu()[["gene"]]), server = TRUE)

        formatted_col_names <- colnames(seu()@meta.data) %>%
            make_chevreul_clean_names()

        updateSelectizeInput(session, "varSelect", choices = formatted_col_names, selected = "batch")
    })

    displayvalues <- reactive({
        req(input$varSelect)
        req(seu())
        unique(seu()[][[input$varSelect]])
    })

    output$displayvaluesui <- renderUI({
        req(input$varSelect)
        selectizeInput(ns("displayvalues"), "groups to display", choices = displayvalues(), multiple = TRUE)
    })

    bigwig_tbl <- reactive({
        load_bigwigs(seu(), bigwig_db)
    })

    coverage_return <- eventReactive(input$plotCoverage, {
        req(seu())
        req(bigwig_tbl())

        plot_gene_coverage_by_var(
            genes_of_interest = input$geneSelect,
            cell_metadata = seu()@meta.data,
            bigwig_tbl = bigwig_tbl(),
            var_of_interest = input$varSelect,
            values_of_interest = input$displayvalues,
            organism = seu()@misc$experiment$organism,
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

    output$coverageTable <- DT::renderDT({
        DT::datatable(coverage_return()$table,
            extensions = "Buttons",
            options = list(dom = "Bft", buttons = c("copy", "csv"), scrollY = "400px")
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
