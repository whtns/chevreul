#' Create Seurat App
#'
#' @param seurat_object a seurat object
#' @param feature_type "gene" or "transcript"
#' @param organism_type "human" or "mouse"
#' @param futureMb amount of Mb allocated to future package
#'
#' @return
#' @export
#'
#' @examples
viewSeurat <- function(seurat_object, feature_type = "gene",
                       organism_type = "human", futureMb = 13000) {
  print(feature_type)
  future::plan(strategy = "multicore", workers = 6)
  future_size <- futureMb * 1024^2
  options(future.globals.maxSize = future_size)
  options(shiny.maxRequestSize = 40 * 1024^2)
  options(DT.options = list(
    pageLength = 2000, paging = FALSE,
    info = TRUE, searching = TRUE, autoWidth = F, ordering = TRUE,
    scrollX = TRUE, language = list(search = "Filter:")
  ))
  header <- shinydashboard::dashboardHeader(title = "Seurat Tool")
  sidebar <- shinydashboard::dashboardSidebar(textOutput("appTitle"),
                                              shinyWidgets::prettyRadioButtons("organism_type",
                                                                               "Organism",
                                                                               choices = c("human", "mouse"), selected = organism_type),
  actionButton("stopApp", "Return to Rstudio"),
  verbatimTextOutput("savefile"), actionButton("changeEmbedAction",
    label = "Change Embedding Parameters"
  ), changeEmbedParamsui("changeembed"),
  shinydashboard::sidebarMenu(
    shinydashboard::menuItem("Reformat Metadata",
    tabName = "reformatMetadata"
  ), shinydashboard::menuItem("Compare Scatter Plots",
    tabName = "comparePlots"
  ), shinydashboard::menuItem("Compare Read Counts",
    tabName = "compareReadCount"
  ), shinydashboard::menuItem("Violin/Heatmap Plots",
    tabName = "violinPlots"
  ), shinydashboard::menuItem("Differential Expression",
    tabName = "diffex"
  ), shinydashboard::menuItem("Find Markers",
    tabName = "findMarkers"
  ), shinydashboard::menuItem("Subset Seurat Input",
    tabName = "subsetSeurat"
  ), shinydashboard::menuItem("Monocle",
    tabName = "monocle"
  ), shinydashboard::menuItem("Regress Features",
    tabName = "regressFeatures"
  )),
  width = 450
  )
  body <- shinydashboard::dashboardBody(shinydashboard::tabItems(
    shinydashboard::tabItem(
      tabName = "violinPlots",
      h2("Violin Plots"), fluidRow(
        box(plotViolinui("violinPlot")),
        box(plotHeatmapui("heatMap"))
      )
    ), shinydashboard::tabItem(
      tabName = "comparePlots",
      h2("Compare Plots"), fluidRow(
        box(plotDimRedui("hello")),
        box(plotDimRedui("howdy"))
      ), fluidRow(box(
        title = "Selected Cells",
        tableSelectedui("hello"), width = 6
      ), box(plotClustree_UI("clustreePlot")))
    ),
    shinydashboard::tabItem(
      tabName = "reformatMetadata",
      h2("Reformat Metadata"), fluidRow((reformatMetadataui("hello")))
    ),
    shinydashboard::tabItem(
      tabName = "compareReadCount",
      h2("Compare Read Counts"), fluidRow(
        box(plotReadCountui("hello2")),
        box(plotReadCountui("howdy2"))
      )
    ), shinydashboard::tabItem(
      tabName = "subsetSeurat",
      h2("Subset Seurat Input"), column(box(plotDimRedui("subset"),
        width = 12
      ), width = 6), column(box(shinyWidgets::actionBttn(
        "subsetAction",
        "subset seurat by selected cells"
      ), shinyWidgets::actionBttn(
        "subsetCsv",
        "subset seurat by uploaded csv"
      ), fileInput("uploadCsv",
        "Upload .csv file with cells to include",
        accept = c(".csv")
      ),
      shinyjs::useShinyjs(), textOutput("subsetMessages"),
      width = 12
      ), box(
        title = "Selected Cells", tableSelectedui("subset"),
        width = 12
      ), width = 6)
    ), shinydashboard::tabItem(
      tabName = "findMarkers",
      h2("Find Markers"), fluidRow(
        box(findMarkersui("hello")),
        box(plotDimRedui("markerScatter"))
      )
    ), shinydashboard::tabItem(
      tabName = "diffex",
      h2("Differential Expression"), column(box(plotDimRedui("diffex"),
        width = 12
      ), box(tableSelectedui("diffex"),
        width = 12
      ), width = 6), column(box(diffexui("hello"),
        width = 12
      ), width = 6)
    ), shinydashboard::tabItem(
      tabName = "regressFeatures",
      h2("Regress Features"), fluidRow(actionButton(
        "regressAction",
        "Regress Seurat Objects By Genes"
      ), box(checkboxInput("runRegression",
        "Run Regression?",
        value = FALSE
      ), checkboxGroupInput("priorGeneSet",
        "Choose a marker gene set:",
        choices = c(
          "Apoptosis",
          "Cell Cycle"
        )
      ), selectizeInput("geneSet",
        "List of genes",
        choices = NULL, multiple = TRUE
      ),
      textInput("geneSetName", "Name for Gene Set"),
      width = 12
      ))
    ), shinydashboard::tabItem(
      tabName = "monocle",
      h2("Monocle"), fluidRow(
        box(actionButton(
          "calcCDS",
          "Calculate Pseudotime"
        ), shinyFiles::shinySaveButton("saveCDS",
          "Save Existing Pseudotime to File", "Save file as...",
          filetype = list(rds = "rds")
        ), shinyFiles::shinyFilesButton("loadCDS",
          "Load Pseudotime from File", "Load Pseudotime File",
          multiple = FALSE
        ), sliderInput("cdsResolution",
          "Resolution of clustering algorithm (affects number of clusters)",
          min = 0.2, max = 2, step = 0.2, value = 0.6
        )),
        fluidRow(box(monocleui("arrow"), width = 12))
      )
    )
  ))
  ui <- function(request) {
    ui <- dashboardPage(
      header = header, sidebar = sidebar,
      body = body
    )
  }
  myserver <- function(input, output, session) {
    # options(warn = -1)

    observe({
      if(input$stopApp > 0){
        stopApp(seu[[feature_type]])
      }
    })

    seu <- reactiveValues()
    observe({
      seu[[feature_type]] <- seurat_object
      seu$active <- seurat_object
      print(names(seu))
    })
    featureType <- reactive({
      feature_type
    })
    organismType <- reactive({
      organism_type
    })
    plot_types <- reactive({
      list_plot_types(seu$active)
    })
    # observe({
    #   shinyFiles::shinyFileSave(input, "saveSeurat",
    #     roots = dataset_volumes(),
    #     session = session, restrictions = system.file(package = "base")
    #   )
    # })
    # subSeuratPath <- eventReactive(input$saveSeurat, {
    #   req(seu$active)
    #   savefile <- shinyFiles::parseSavePath(
    #     dataset_volumes(),
    #     input$saveSeurat
    #   )
    #   return(savefile$datapath)
    # })
    # observeEvent(input$saveSeurat, {
    #   req(seu$active)
    #   req(subSeuratPath())
    #   shiny::withProgress(
    #     message = paste0("Saving Data"),
    #     value = 0,
    #     {
    #       Sys.sleep(6)
    #       shiny::incProgress(2 / 10)
    #       saveRDS(
    #         shiny::reactiveValuesToList(seu),
    #         subSeuratPath()
    #       )
    #       shiny::incProgress(10 / 10)
    #     }
    #   )
    # })


    seu <- callModule(reformatMetadata, "hello", seu)
    callModule(plotDimRed, "hello", seu, plot_types, featureType,
      organism_type = organismType
    )
    callModule(plotDimRed, "howdy", seu, plot_types, featureType,
      organism_type = organismType
    )
    callModule(plotDimRed, "diffex", seu, plot_types, featureType,
      organism_type = organismType
    )
    callModule(plotDimRed, "subset", seu, plot_types, featureType,
      organism_type = organismType
    )
    callModule(plotDimRed, "markerScatter", seu, plot_types,
      featureType,
      organism_type = organismType
    )
    callModule(plotReadCount, "hello2", seu, plot_types)
    callModule(plotReadCount, "howdy2", seu, plot_types)
    callModule(
      plotViolin, "violinPlot", seu, featureType,
      organismType
    )
    callModule(
      plotHeatmap, "heatMap", seu, featureType,
      organismType
    )
    callModule(plotClustree, "clustreePlot", seu)
    callModule(tableSelected, "hello", seu)
    diffex_selected_cells <- callModule(
      tableSelected, "diffex",
      seu
    )
    subset_selected_cells <- callModule(
      tableSelected, "subset",
      seu
    )
    observeEvent(input$subsetAction, {
      req(subset_selected_cells())
      withCallingHandlers(
        {
          shinyjs::html("subsetMessages", "")
          message("Beginning")
            seu[[feature_type]] <- seu[[feature_type]][, subset_selected_cells()]
          if (length(unique(seu$gene[[]]$batch)) > 1) {
            print(names(seu))
              message(paste0("reintegrating ", feature_type, " expression"))
              seu[[feature_type]] <- reintegrate_seu(seu[[feature_type]],
                feature = feature_type,
                resolution = seq(0.2, 2, by = 0.2)
              )
          }
          else {
              seu[[feature_tyupe]] <- seurat_pipeline(seu[[feature_type]], resolution = seq(0.2,
                2,
                by = 0.2
              ))
          }
          seu$active <- seu[[feature_type]]
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
            seu[[feature_type]] <- subset_by_meta(
              input$uploadCsv$datapath,
              seu[[feature_type]]
            )
          if (length(unique(seu$gene[[]]$batch)) > 1) {
              message(paste0("reintegrating ", feature_type, " expression"))
              seu[[feature_type]] <- reintegrate_seu(seu[[feature_type]],
                feature = feature_type,
                resolution = seq(0.2, 2, by = 0.2)
              )
          }
          else {
              seu[[feature_type]] <- seurat_pipeline(seu[[feature_type]], resolution = seq(0.2,
                2,
                by = 0.2
              ))
          }
          seu$active <- seu[[feature_type]]
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
      seu <- callModule(
        changeEmbedParams, "changeembed",
        seu
      )
      removeModal()
    })
    callModule(findMarkers, "hello", seu)
    diffex_results <- callModule(
      diffex, "hello", seu, featureType,
      diffex_selected_cells
    )

    prior_gene_set <- reactive({
      req(input$priorGeneSet)
      if (input$priorGeneSet == "Apoptosis") {
        if (organism_type == "human"){
          c(
            "CASP3", "CASP7", "BAX", "BAK1", "BID", "BBC3",
            "BCL2", "MCL1"
          )
        } else if (organism_type == "mouse"){
          c("Casp3", "Casp7", "Bax", "Bak1", "Bid", "Bbc3", "Bcl2",
            "Mcl1")

        }
      }
      else if (input$priorGeneSet == "Cell Cycle") {
        if (organism_type == "human"){
          c(
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
        } else if (organism_type == "mouse") {
          c("Mcm5", "Pcna", "Tyms", "Fen1", "Mcm2", "Mcm4", "Rrm1",
            "Ung", "Gins2", "Mcm6", "Cdca7", "Dtl", "Prim1", "Uhrf1",
            "Mlf1ip", "Hells", "Rfc2", "Rpa2", "Nasp", "Rad51ap1",
            "Gmnn", "Wdr76", "Slbp", "Ccne2", "Ubr7", "Pold3", "Msh2",
            "Atad2", "Rad51", "Rrm2", "Cdc45", "Cdc6", "Exo1", "Tipin",
            "Dscc1", "Blm", "Casp8ap2", "Usp1", "Clspn", "Pola1", "Chaf1b",
            "Brip1", "E2f8")

        }

      }
      else if (is.null(input$priorGeneSet)) {
        c("")
      }
    })
    observe({
      updateSelectizeInput(session, "geneSet",
        choices = rownames(seu$active),
        selected = prior_gene_set(), server = TRUE
      )
    })
    observeEvent(input$regressAction, {
      req(seu$active)
      showModal(modalDialog(
        title = "Regressing out provided list of features",
        "This process may take a minute or two!"
      ))
      seu[[feature_type]] <- seuratTools::regress_by_features(seu[[feature_type]],
        feature_set = list(input$geneSet), set_name = janitor::make_clean_names(input$geneSetName),
        regress = input$runRegression
      )
      seu$active <- seu[[feature_type]]
      removeModal()
    })
    cds <- reactiveValues()
    observeEvent(input$calcCDS, {
      req(seu$active)
      cds$traj <- convert_seu_to_cds(seu$active, resolution = input$cdsResolution)
    })
    observeEvent(input$calcCDS, {
      req(cds$traj)
      cds$traj <- learn_graph_by_resolution(cds$traj,
        seu$active,
        resolution = input$cdsResolution
      )
    })
    observe({
      shinyFiles::shinyFileChoose(input, "loadCDS",
        roots = dataset_volumes(),
        session = session
      )
    })
    cdsLoadPath <- eventReactive(input$loadCDS, {
      file <- shinyFiles::parseFilePaths(
        dataset_volumes(),
        input$loadCDS
      )
      file$datapath
    })
    observeEvent(input$loadCDS, {
      req(cdsLoadPath())
      shiny::withProgress(
        message = paste0("Uploading Data"),
        value = 0,
        {
          Sys.sleep(6)
          shiny::incProgress(2 / 10)
          dataset <- readRDS(cdsLoadPath())
          shiny::incProgress(10 / 10)
          for (i in names(dataset)) {
            cds[[i]] <- dataset[[i]]
          }
        }
      )
    })

    callModule(monocle, "arrow", cds, seu, featureType,
      resolution = reactive(input$cdsResolution)
    )

  }
  runApp(shinyApp(ui, myserver))
}
