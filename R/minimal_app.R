#' Create a minimal seurat app
#'
#' @param object
#' @param loom_path
#' @param appTitle
#' @param feature_types
#' @param organism_type
#' @param futureMb
#'
#' @return
#' @export
#'
#' @examples
minimalSeuratApp <- function(object = seurat_pancreas_reduced, loom_path = NULL, appTitle = NULL, feature_types = "gene",
                             organism_type = "human", futureMb = 13000) {
  print(feature_types)
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
    textOutput("appTitle"),
    uiOutput("featureType"),
    shinyWidgets::prettyRadioButtons("organism_type",
      inline = TRUE,
      "Organism", choices = c("human", "mouse"), selected = organism_type
    ),
    shinydashboard::sidebarMenu(
      shinydashboard::menuItem("Integrate Projects",
        tabName = "integrateProjects", icon = icon("object-group")
      ), shinydashboard::menuItem("Reformat Metadata",
        tabName = "reformatMetadata", icon = icon("columns")
      ), shinydashboard::menuItem("Plot Data",
        tabName = "comparePlots", icon = icon("chart-bar"), selected = TRUE
      ), shinydashboard::menuItem("Heatmap/Violin Plots",
        tabName = "violinPlots", icon = icon("sort")
      ), shinydashboard::menuItem("Coverage Plots",
        tabName = "coveragePlots", icon = icon("mountain")
      ), shinydashboard::menuItem("Differential Expression",
        tabName = "diffex", icon = icon("magnet")
      ), shinydashboard::menuItem("Pathway Enrichment Analysis",
        tabName = "pathwayEnrichment", icon = icon("sitemap")
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
        seuratToolsBox(
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
        tabName = "reformatMetadata",
        fluidRow(
          reformatMetadataui("reformatmetadata")
        )
      ),
      shinydashboard::tabItem(
        tabName = "subsetSeurat",
        h2("Subset Seurat Input") %>%
          default_helper(type = "markdown", content = "subsetSeurat"),
        plotDimRedui("subset"),
        seuratToolsBox(
          title = "Subset Settings",
          checkboxInput("legacySettingsSubset", "Use Legacy Settings", value = FALSE),
          actionButton("subsetAction", "subset seurat by selected cells"),
          actionButton("subsetCsv", "subset seurat by uploaded csv"),
          fileInput("uploadCsv",
            "Upload .csv file with cells to include",
            accept = c(".csv")
          ),
          shinyjs::useShinyjs(),
          # shinyjs::runcodeUI(code = "shinyjs::alert('Hello!')"),
          textOutput("subsetMessages"),
          width = 6
        ),
        seuratToolsBox(
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
      shinydashboard::tabItem(
        tabName = "pathwayEnrichment",
        h2("Pathway Enrichment"),
        fluidRow(
          pathwayEnrichmentui("pathwayEnrichment")
        )
      ),
      shinydashboard::tabItem(
        tabName = "regressFeatures",
        fluidRow(
          seuratToolsBox(
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
    w <- waiter::Waiter$new()

    shinyhelper::observe_helpers(help_dir = system.file("helpers", package = "seuratTools", lib.loc = "/dataVolume/storage/rpkgs/devel_install/"))
    options(warn = -1)
    # shinylogs::track_usage(storage_mode = shinylogs::store_json(path = "logs/"))
    # projects_db <- "/dataVolume/storage/single_cell_projects/single_cell_projects.db"

    seu <- reactiveVal(NULL)
    observe({
      seu(object)
    })

    organism_type <- reactive({
      input$organism_type
    })

    loom_path <- reactive({
      loom_path
    })

    plot_types <- reactive({
      req(seu())
      list_plot_types(seu())
    })

    output$featureType <- renderUI({
      req(seu())
      shinyWidgets::prettyRadioButtons("feature_type",
        "Feature for Display",
        choices = c("gene", "transcript"),
        selected = "gene", inline = TRUE
      )
    })

    observeEvent(input$feature_type, {
      req(seu())
      # DefaultAssay(seu()) <- paste0(input$feature_type, ".", "integrated")
      # DefaultAssay(seu()) <- "gene"
    })

    featureType <- reactive({
      "gene"
      # input$feature_type
    })

    integrationResults <- callModule(
      integrateProj, "integrateproj",
      proj_matrices, seu, proj_dir, con
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
    #
    # # reformatted_seu <- reactive({
    # #   req(seu())
    # #   callModule(reformatMetadata, "reformatmetadata", seu)
    # # })
    # #
    # # observe({
    # #   req(reformatted_seu())
    # #   seu(reformatted_seu())
    # # })

    reductions <- reactive({
      req(seu())
      names(seu()@reductions)
      # c("pca", "tsne", "umap")
    })

    observe({
      req(seu())

      callModule(plotDimRed, "plotdimred1", seu, plot_types, featureType,
        organism_type = organism_type, reductions
      )
      callModule(plotDimRed, "plotdimred2", seu, plot_types, featureType,
        organism_type = organism_type, reductions
      )
      callModule(plotDimRed, "diffex", seu, plot_types, featureType,
        organism_type = organism_type, reductions
      )
      callModule(plotDimRed, "subset", seu, plot_types, featureType,
        organism_type = organism_type, reductions
      )
      callModule(plotDimRed, "markerScatter", seu, plot_types, featureType,
        organism_type = organism_type, reductions
      )
    })

    callModule(plotReadCount, "plotreadcount1", seu, plot_types)
    callModule(plotReadCount, "plotreadcount2", seu, plot_types)

    callModule(
      plotViolin, "violinPlot", seu, featureType,
      organism_type
    )
    callModule(
      plotHeatmap, "heatMap", seu, featureType,
      organism_type
    )
    callModule(plotClustree, "clustreePlot", seu)
    callModule(tableSelected, "tableselected", seu)
    diffex_selected_cells <- callModule(
      tableSelected, "diffex",
      seu
    )
    subset_selected_cells <- callModule(
      tableSelected, "subset",
      seu
    )

    observeEvent(input$changeEmbedAction, {
      showModal(modalDialog(
        title = "Recalculating Embedding",
        "This process may take a minute or two!"
      ))
      seu <- callModule(
        changeEmbedParams, "changeembed",
        seu
      )
      removeModal()
    })

    callModule(findMarkers, "findmarkers", seu, plot_types, featureType)

    callModule(pathwayEnrichment, "pathwayEnrichment", seu)

    # plot all transcripts
    observe({
      req(featureType())
      req(seu())
      if ("transcript" %in% names(seu()@assays)) {
        callModule(
          allTranscripts, "alltranscripts1", seu, featureType,
          organism_type
        )

        callModule(plotDimRed, "alltranscripts2", seu, plot_types, featureType,
          organism_type = organism_type, reductions
        )
      }
    })

    diffex_results <- callModule(
      diffex, "diffex", seu, featureType,
      diffex_selected_cells
    )
    observeEvent(input$subsetAction, {
      req(input$subsetAction)
      req(subset_selected_cells())
      withCallingHandlers(
        {
          shinyjs::html("subsetMessages", "")
          message("Beginning")
          subset_seu <- seu()[, colnames(seu()) %in% subset_selected_cells()]
          seu(subset_seu)
          if (length(unique(seu()$batch)) > 1) {
            message(paste0("reintegrating ", i, " expression"))
            reintegrated_seu <- reintegrate_seu(seu(),
              feature = i,
              resolution = seq(0.2, 2, by = 0.2),
              organism = seu()@misc$experiment$organism
            )
            seu(reintegrated_seu)
          }
          else {
            processed_seu <- seurat_pipeline(seu(), resolution = seq(0.2, 2, by = 0.2))
            seu(processed_seu)
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
          subset_seu <- subset_by_meta(
            input$uploadCsv$datapath,
            seu()
          )
          seu(subset_seu)

          if (length(unique(seu()[[]]$batch)) > 1) {
            message(paste0("reintegrating ", i, " expression"))
            reintegrated_seu <- reintegrate_seu(seu(),
              feature = i,
              resolution = seq(0.2, 2, by = 0.2),
              organism = seu()@misc$experiment$organism
            )
            seu(reintegrated_seu)
          }
          else {
            processed_seu <- seurat_pipeline(seu(), resolution = seq(0.2, 2, by = 0.2))
            seu(processed_seu)
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

    prior_gene_set <- reactive({
      # req(input$priorGeneSet)
      req(seu())

      if (is.null(input$priorGeneSet)) {
        ""
      } else if (input$priorGeneSet == "Apoptosis") {
        marker_genes <- c(
          "CASP3", "CASP7", "BAX", "BAK1", "BID", "BBC3",
          "BCL2", "MCL1"
        )

        marker_genes[marker_genes %in% rownames(seu())]
      }
      else if (input$priorGeneSet == "Cell Cycle") {
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

        marker_genes[marker_genes %in% rownames(seu())]
      } else if (input$priorGeneSet == "Mitochondrial") {
        marker_genes <- mito_features[[organism_type()]][["gene"]]

        marker_genes[marker_genes %in% rownames(seu())]
      } else if (input$priorGeneSet == "Ribosomal") {
        marker_genes <- ribo_features[[organism_type()]][["gene"]]

        marker_genes[marker_genes %in% rownames(seu())]
      }
    })

    observe({
      req(seu())
      updateSelectizeInput(session, "geneSet",
        choices = rownames(seu()),
        selected = prior_gene_set(),
        server = TRUE
      )
    })
    observeEvent(input$regressAction, {
      req(seu())
      showModal(modalDialog(
        title = "Regressing out provided list of features",
        "This process may take a minute or two!"
      ))
      regressed_seu <- seuratTools::regress_by_features(seu(),
        feature_set = list(input$geneSet), set_name = janitor::make_clean_names(input$geneSetName),
        regress = input$runRegression
      )
      seu(regressed_seu)
      removeModal()
    })

    observe({
      req(reductions())
      callModule(
        monocle, "monocle", seu, plot_types, featureType,
        organism_type, reductions
      )
    })

    observe({
      req(seu())
      req(loom_path())

      callModule(plotVelocity, "plotvelocity", seu, loom_path())
    })
  }
  shinyApp(ui, server, enableBookmarking = "server")
}
