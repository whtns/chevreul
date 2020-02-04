#' Run Seurat Differential Expression
#'
#' @param seu
#' @param cluster1
#' @param cluster2
#' @param resolution
#' @param diffex_scheme
#'
#' @return
#' @export
#'
#' @examples
run_seurat_de <- function(seu, cluster1, cluster2, resolution, diffex_scheme = "seurat", featureType, tests = c("t", "wilcox", "bimod")) {

  if (diffex_scheme == "seurat"){

    if ("integrated" %in% names(seu@assays)){
      active_assay <- "integrated"
    } else {
      active_assay <- "RNA"
    }

    Idents(seu) <- paste0(active_assay, "_snn_res.", resolution)
    seu <- subset(seu, idents = c(cluster1, cluster2))
  } else if (diffex_scheme == "custom"){
    # subset by supplied cell ids
    #
    seu <- seu[,c(cluster1, cluster2)]

    keep_cells <- c(cluster1, cluster2)
    new_idents <- c(rep(1, length(cluster1)), rep(2, length(cluster2)))
    names(new_idents) <- keep_cells
    new_idents <- new_idents[colnames(seu)]
    Idents(seu) <- new_idents
    cluster1 = 1
    cluster2 = 2

  }

  test_list <- vector("list", length(tests))

  for (test in tests){
    print(test)
    de <- FindMarkers(seu,
                      ident.1 = cluster1,
                      ident.2 = cluster2,
                      test.use = test)


    if (featureType() == "transcript"){
      de_cols <- c("enstxp", "ensgene", "symbol", "p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj")

      de <-
        de %>%
        tibble::rownames_to_column("enstxp") %>%
        dplyr::left_join(annotables::grch38_tx2gene, by = "enstxp") %>%
        dplyr::left_join(annotables::grch38, by = "ensgene") %>%
        dplyr::select(one_of(de_cols))

    } else if (featureType() == "gene"){
      de_cols <- c("ensgene", "symbol", "p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj")

      de <-
        de %>%
        tibble::rownames_to_column("symbol") %>%
        dplyr::left_join(annotables::grch38, by = "symbol") %>%
        dplyr::select(one_of(de_cols))
    }

    test_list[[match(test, tests)]] = de

  }
  names(test_list) <- tests
  return(test_list)
}


#' Run Enrichment Browser on Differentially Expressed Genes
#'
#' @param seu
#' @param cluster1_cells
#' @param cluster2_cells
#'
#' @return
#' @export
#'
#' @examples
run_enrichmentbrowser <- function(seu, cluster1_cells, cluster2_cells, ...){

    # subset by supplied cell ids
    #
  seu <- seu[,c(cluster1_cells, cluster2_cells)]

  keep_cells <- c(cluster1_cells, cluster2_cells)
  new_idents <- c(rep(0, length(cluster1_cells)), rep(1, length(cluster2_cells)))
  names(new_idents) <- keep_cells
  new_idents <- new_idents[colnames(seu)]
  Idents(seu) <- new_idents

  counts <- GetAssayData(seu, slot = "counts")
  counts <- as.matrix(counts)
  mode(counts) <- "integer"

  rowData <- data.frame(rownames(seu), row.names = rownames(seu))

  colData <- as.data.frame(seu[[]])

  se <- SummarizedExperiment::SummarizedExperiment(assays=list(counts=counts),
                                                   rowData=rowData, colData=colData)

  se$GROUP <- forcats::fct_inseq(Idents(seu))

  se <- EnrichmentBrowser::deAna(se, grp = se$GROUP, de.method = "edgeR")
  se <- EnrichmentBrowser::idMap(se, org = "hsa", from = "SYMBOL", to = "ENTREZID")

  outdir <- fs::path("www", "enrichmentbrowser")
  report.name = "mainpage.html"

  hsa.grn <- EnrichmentBrowser::compileGRN(org="hsa", db="kegg")

  # EnrichmentBrowser::ebrowser( meth=c("ora", "ggea"), perm=0, comb=TRUE,
  #           exprs=se, gs=go.gs, grn=hsa.grn, org="hsa", nr.show=3,
  #           out.dir=outdir, report.name=report.name, browse = FALSE)

  sbea.res <- EnrichmentBrowser::sbea(method = "ora", se = se, gs = go.gs, perm = 0,
                                      alpha = 0.1)

  nbea.res <- EnrichmentBrowser::nbea(method="ggea", se=se, gs=go.gs, grn=hsa.grn)

  res <- EnrichmentBrowser::combResults(list(sbea.res, nbea.res))

  EnrichmentBrowser::eaBrowse(res, html.only = TRUE, out.dir = outdir, graph.view=hsa.grn,
                              report.name = report.name)

  return(fs::path("enrichmentbrowser", "mainpage.html"))
}

#' Prep Slider Values
#'
#' @param default_val
#'
#' @return
#' @export
#'
#' @examples
prep_slider_values <- function(default_val){
  min <- round(default_val*0.25, digits = 1)
  max <- round(default_val*2.0, digits = 1)
  step = 10^((ceiling(log10(default_val)))-1)
  value = default_val
  return(list(min = min, max = max, value = value, step = step))
}


#' Create Seurat App
#'
#' @param preset_project A preloaded project to start the app with
#' @param filterTypes A named vector of file suffixes corresponding to subsets of the data
#' @param appTitle A title of the App
#' @param futureMb amount of Mb allocated to future package
#' @param featureTypes
#' @param preset_project
#'
#' @return
#' @export
#'
#' @examples
seuratApp <- function(preset_project, filterTypes, appTitle, feature_types = "gene", futureMb = 849){

  print(feature_types)

  # projList <- system("locate -d /dataVolume/storage/single_cell_projects/single_cell_projects.db '*_proj'", intern = TRUE) %>%
  #   purrr::set_names(fs::path_file(.)) %>%
  #   identity()
  #
  # pc_projList <- fs::path(fs::path_file(fs::path_dir(projList)), fs::path_file(projList))
  # pc_preset <- fs::path(fs::path_file(fs::path_dir(preset_project)), fs::path_file(preset_project))
  #
  # preset_project <- projList[match(pc_preset, pc_projList)]
  #
  # print(projList)
  # print(preset_project)

  future::plan(strategy = "multicore", workers = 6)
  future_size = futureMb*1024^2
  options(future.globals.maxSize= future_size)
  options(DT.options = list(pageLength = 2000, paging = FALSE,
                            info = TRUE, searching = TRUE, autoWidth = F, ordering = TRUE, scrollX = TRUE,
                            language = list(search = "Filter:")))

  # header ------------------------
  header <- shinydashboard::dashboardHeader()

  # sidebar ------------------------
 sidebar <- shinydashboard::dashboardSidebar(
  # selectizeInput("setProject", "Select Project to Load", choices = projList, selected = preset_project, multiple = F),
  uiOutput("projInput"),
  actionButton("loadProject", "Load Selected Project"),
  textOutput("appTitle"),
  uiOutput("featureType"),
  # shinyWidgets::prettyRadioButtons("featureType", "Feature for Display", choices = featureTypes, selected = "gene"),
  shinyWidgets::prettyRadioButtons("organism_type", "Organism", choices = c("human", "mouse"), selected = "human"),
  shinyFiles::shinyFilesButton("seuratUpload", "Load a Seurat Dataset", "Please select a .rds file", multiple = FALSE),
  shinyFiles::shinySaveButton("saveSeurat", "Save current Dataset", "Save file as...", filetype = list(rds = "rds")),
  verbatimTextOutput("savefile"),
  actionButton("changeEmbedAction", label = "Change Embedding Parameters"),
  changeEmbedParamsui("changeembed"),
  # bookmarkButton(),
  shinydashboard::sidebarMenu(
    shinydashboard::menuItem("Integrate Projects", tabName = "integrateProjects"),
    shinydashboard::menuItem("Reformat Metadata", tabName = "reformatMetadata"),
    shinydashboard::menuItem("Compare Scatter Plots", tabName = "comparePlots"),
    shinydashboard::menuItem("Compare Read Counts", tabName = "compareReadCount"),
    shinydashboard::menuItem("Violin Plots", tabName = "violinPlots"),
    shinydashboard::menuItem("Differential Expression", tabName = "diffex"),
    shinydashboard::menuItem("Gene Enrichment Analysis", tabName = "geneEnrichment"),
    shinydashboard::menuItem("Find Markers", tabName = "findMarkers"),
    shinydashboard::menuItem("Subset Seurat Input", tabName = "subsetSeurat"),
    shinydashboard::menuItem("All Transcripts", tabName = "allTranscripts"),
    # shinydashboard::menuItem("RNA Velocity", tabName = "rnaVelocity"),
    shinydashboard::menuItem("Monocle", tabName = "monocle"),
    # shinydashboard::menuItem("cellAlign", tabName = "cellAlign"),
    shinydashboard::menuItem("Regress Features", tabName = "regressFeatures")
  ),
  width = 450
)

# body ------------------------
  body <- shinydashboard::dashboardBody(shinydashboard::tabItems(
    shinydashboard::tabItem(
      tabName = "violinPlots",
      h2("Violin Plots"), fluidRow(
        plotViolinui("violinPlot")
      )
    ),
  shinydashboard::tabItem(
    tabName = "comparePlots",
    h2("Compare Plots"), fluidRow(
      box(plotDimRedui(
        "hello"
      )), box(plotDimRedui(
        "howdy"
      ))
    ), fluidRow(box(
      title = "Selected Cells",
      tableSelectedui("hello"), width = 12
    ))
  ),
  shinydashboard::tabItem(
    tabName = "integrateProjects",
    h2("Integrate Projects"), fluidRow(
      (integrateProjui(
        "hello"
      )))
  ), shinydashboard::tabItem(
    tabName = "reformatMetadata",
    h2("Reformat Metadata"), fluidRow(
      (reformatMetadataui(
        "hello"
      )))
  ), shinydashboard::tabItem(
    tabName = "compareReadCount",
    h2("Compare Read Counts"), fluidRow(box(plotReadCountui(
      "hello2"
    )), box(plotReadCountui(
      "howdy2"
    )))
  ), shinydashboard::tabItem(
    tabName = "subsetSeurat",
    h2("Subset Seurat Input"), column(box(plotDimRedui(
      "subset"
    ), width = 12), width = 6),
    column(
      box(shinyWidgets::actionBttn(
        "subsetAction",
        "subset seurat by selected cells"
      ),
      shinyWidgets::actionBttn(
        "subsetCsv",
        "subset seurat by uploaded csv"
      ),
      fileInput("uploadCsv", "Upload .csv file with cells to include", accept = c(".csv")),
      shinyjs::useShinyjs(),
      textOutput("subsetMessages"),
      width = 12
      ),
      box(
        title = "Selected Cells",
        tableSelectedui("subset"), width = 12
      ),
      width = 6
    )
  ),
  shinydashboard::tabItem(
    tabName = "findMarkers", h2("Find Markers"),
    fluidRow(
      box(findMarkersui("hello")),
      box(plotDimRedui("markerScatter"))
      )
  ), shinydashboard::tabItem(
    tabName = "allTranscripts",
    h2("All Transcripts"), fluidRow(shinyWidgets::actionBttn(
      "plotTrx",
      "Plot all transcripts"
    )), fluidRow(column(allTranscriptsui("hello"),
      width = 6
    ), column(allTranscriptsui("howdy"),
      width = 6
    ))
  ), shinydashboard::tabItem(
    tabName = "diffex",
    h2("Differential Expression"), column(box(plotDimRedui(
      "diffex"
    ), width = 12),
      box(
        tableSelectedui("diffex"), width = 12
      ), width = 6),
    column(
      box(
        diffexui("hello"),
        width = 12
      ),
      width = 6
    )
  ),
  shinydashboard::tabItem(
    tabName = "geneEnrichment",
    h2("Gene Enrichment"),
    geneEnrichmentui("hello")
  ),
  shinydashboard::tabItem(
    tabName = "regressFeatures",
    h2("Regress Features"),
    fluidRow(
      actionButton("regressAction", "Regress Seurat Objects By Genes"),
      box(
        selectizeInput("geneSet", "List of genes", choices = NULL, multiple = TRUE),
        textInput("geneSetName", "Name for Gene Set"),
        width = 12
      )
    )
  ),
  shinydashboard::tabItem(
    tabName = "monocle",
    h2("Monocle"),
    fluidRow(
      box(
      actionButton("calcCDS", "Calculate Pseudotime"),
      shinyFiles::shinySaveButton("saveCDS", "Save Existing Pseudotime to File", "Save file as...", filetype = list(rds = "rds")),
      shinyFiles::shinyFilesButton("loadCDS", "Load Pseudotime from File", "Load Pseudotime File", multiple = FALSE),
      sliderInput("cdsResolution", "Resolution of clustering algorithm (affects number of clusters)", min = 0.2, max = 2, step = 0.2, value = 0.6)
      ),
      fluidRow(
        box(monocleui("arrow"),
            width = 12
        )
      )
    )
  )
))

  # shinydashboard ui ------------------------
  ui <- function(request){
    ui <- dashboardPage(
      header = header,
      sidebar = sidebar,
      body = body
    )
  }

# shinydashboard server ------------------------
  server <- function(input, output, session) {
    options(warn = -1)

    projList <- reactive({
      system("locate -d /dataVolume/storage/single_cell_projects/single_cell_projects.db '*_proj'", intern = TRUE) %>%
      purrr::set_names(fs::path_file(.))
    })

    # presetProject <- reactive({
    #   pc_projList <- fs::path(fs::path_file(fs::path_dir(projList())), fs::path_file(projList()))
    #   pc_preset <- fs::path(fs::path_file(fs::path_dir(preset_project)), fs::path_file(preset_project))
    #
    #   projList()[match(pc_preset, pc_projList())]
    # })

    output$projInput = renderUI({
      selectizeInput("setProject", "Select Project to Load", choices = projList(), selected = preset_project, multiple = F)
    })

    proj_matrices <- reactive({
      create_proj_matrix(projList())
    })

# sidebar ------------------------
    seu <- reactiveValues()
    proj_dir <- reactiveVal()

    if(!is.null(preset_project)){
      proj_dir(preset_project)
    }

    organism_type <- reactive({
      input$organism_type
    })

    plot_types <- reactive({
      list_plot_types(seu$active)
    })

    # upload seurat object

    observeEvent(input$loadProject, {
      proj_dir(input$setProject)
    })

    output$appTitle <- renderText({
      req(proj_dir())
      paste0("Loaded Project: ", fs::path_file(proj_dir()))
    })

    # list volumes
    volumes <- reactive({
      print(proj_dir())
      volumes <- c(Home = fs::path(proj_dir(), "output", "seurat"), "R Installation" = R.home(), shinyFiles::getVolumes())
      # print(volumes)
      })


    observe({
      req(volumes())
      shinyFiles::shinyFileChoose(input, "seuratUpload", roots = volumes(), session = session)
    })

    uploadSeuratPath <- eventReactive(input$seuratUpload, {
      req(volumes())
      file <- shinyFiles::parseFilePaths(volumes(), input$seuratUpload)
      file$datapath
    })

    observeEvent(input$seuratUpload, {
      req(uploadSeuratPath())

      shiny::withProgress(
        message = paste0("Uploading Data"),
        value = 0,
        {
          # Sys.sleep(6)
          shiny::incProgress(2/10)
          # Sys.sleep(12)
          # shiny::incProgress(4/10)
          # Sys.sleep(18)
          # shiny::incProgress(6/10)
          print(uploadSeuratPath())
          dataset <- readRDS(uploadSeuratPath())

          if(!typeof(dataset[[1]]@misc$markers[[1]]) == "list"){
            dataset <- purrr::map(dataset, find_all_markers)
            saveRDS(dataset, uploadSeuratPath())
          }

          shiny::incProgress(6/10)
          # dataset$gene <- find_all_markers(dataset$gene)
          # dataset$transcript <- find_all_markers(dataset$transcript)

          seu_names <- names(dataset)[!names(dataset) == "active"]

          for (i in seu_names){
            seu[[i]] <- dataset[[i]]
          }

          print(names(seu))

          shiny::incProgress(8/10)
        }
      )

    })

    output$featureType <- renderUI({
      req(seu)

      seu_names <- names(seu)[!(names(seu) == "active")]

      shinyWidgets::prettyRadioButtons("feature_type", "Feature for Display", choices = seu_names, selected = "gene", inline = T)
    })

    observeEvent(input$feature_type, {
      seu$active <- seu[[input$feature_type]]
    })

    featureType <- reactive({
      featureType <- input$feature_type
      # "gene"
    })

    # save seurat object
    observe({
      shinyFiles::shinyFileSave(input, "saveSeurat", roots = volumes(), session = session, restrictions = system.file(package = "base"))
    })


    subSeuratPath <- eventReactive(input$saveSeurat, {

      req(seu$active)
      savefile <- shinyFiles::parseSavePath(volumes(), input$saveSeurat)

      return(savefile$datapath)


    })

    observe({
      req(seu$active)
      req(subSeuratPath())
      shiny::withProgress(
        message = paste0("Saving Data"),
        value = 0,
        {
          Sys.sleep(6)
          shiny::incProgress(2/10)
          saveRDS(shiny::reactiveValuesToList(seu), subSeuratPath())
          shiny::incProgress(10/10)

        }
      )
    })

# body ------------------------

    integrationResults <- callModule(integrateProj, "hello", proj_matrices, seu, proj_dir)

    observe({
      req(integrationResults())
      create_proj_db()
      integration_path <- paste0(integrationResults(), "_proj")
      proj_dir(integration_path)

      newintegrated_project <- purrr::set_names(integration_path, fs::path_file(integration_path))

      projList <- c(projList, newintegrated_project)

      updateSelectizeInput(session, "setProject",
                        label = "Select input label",
                        choices = projList,
      )
    })

    seu <- callModule(reformatMetadata, "hello", seu)
    callModule(plotDimRed, "hello", seu, plot_types, featureType, organism_type = organism_type)
    callModule(plotDimRed, "howdy", seu, plot_types, featureType, organism_type = organism_type)
    callModule(plotDimRed, "diffex", seu, plot_types, featureType, organism_type = organism_type)
    callModule(plotDimRed, "subset", seu, plot_types, featureType, organism_type = organism_type)
    callModule(plotDimRed, "markerScatter", seu, plot_types, featureType, organism_type = organism_type)
    callModule(plotReadCount, "hello2", seu, plot_types)
    callModule(plotReadCount, "howdy2", seu, plot_types)
    callModule(plotViolin, "violinPlot", seu, featureType)
    callModule(tableSelected, "hello", seu)
    diffex_selected_cells <- callModule(tableSelected, "diffex", seu)

    subset_selected_cells <- callModule(tableSelected, "subset", seu)

    upload_cells <- reactive({
      req(input$uploadCsv)
      upload_cells <- readr::read_csv(input$uploadCsv$datapath) %>%
        dplyr::pull(X1)
    })

    observeEvent(input$subsetAction, {

      withCallingHandlers({
        shinyjs::html("subsetMessages", "")
        message("Beginning")

        for (i in names(seu)[!(names(seu) == "active")]){
          seu[[i]] <- seu[[i]][, subset_selected_cells()]
        }

        if(length(unique(seu$gene[[]]$batch)) > 1){
          print(names(seu)[!(names(seu) == "active")])
          for (i in names(seu)[!(names(seu) == "active")]){
            message(paste0("reintegrating ", i, " expression"))
            # harmony
            # seu[[i]] <- seurat_pipeline(seu[[i]], reduction = "harmony", resolution = seq(0.2, 2.0, by = 0.2))
            # seurat cca
            seu[[i]] <- reintegrate_seu(seu[[i]], feature = i, resolution = seq(0.2, 2.0, by = 0.2))
          }

        } else {

          for (i in names(seu)[!(names(seu) == "active")]){
            seu[[i]] <- seurat_pipeline(seu[[i]], resolution = seq(0.2, 2.0, by = 0.2))
          }

        }
        seu$active <- seu[[input$feature_type]]

        message("Complete!")

      },
      message = function(m) {
        shinyjs::html(id = "subsetMessages", html = paste0("Subsetting Seurat Object: ", m$message), add = FALSE)
      })
    })

    observeEvent(input$subsetCsv, {

      withCallingHandlers({
        shinyjs::html("subsetMessages", "")
        message("Beginning")

        for (i in names(seu)[!(names(seu) == "active")]){
          seu[[i]] <- seu[[i]][, upload_cells()]
        }

        if(length(unique(seu$gene[[]]$batch)) > 1){

          for (i in names(seu)[!(names(seu) == "active")]){
            message(paste0("reintegrating ", i, " expression"))
            # harmony
            # seu[[i]] <- seurat_pipeline(seu[[i]], reduction = "harmony", resolution = seq(0.2, 2.0, by = 0.2))
            # seurat cca
            seu[[i]] <- reintegrate_seu(seu[[i]], feature = i, resolution = seq(0.2, 2.0, by = 0.2))
          }

        } else {

          for (i in names(seu)[!(names(seu) == "active")]){
            seu[[i]] <- seurat_pipeline(seu[[i]], resolution = seq(0.2, 2.0, by = 0.2))
          }

        }
        seu$active <- seu[[input$feature_type]]

        message("Complete!")

      },
      message = function(m) {
        shinyjs::html(id = "subsetMessages", html = paste0("Subsetting Seurat Object: ", m$message), add = FALSE)
      })
    })

    observeEvent(input$changeEmbedAction,{

      showModal(modalDialog(
        title = "Recalculating Embedding",
        "This process may take a minute or two!"
      ))

      seu <- callModule(changeEmbedParams, "changeembed", seu)

      # seu$active <- callModule(embedParam, "minDist", seu$active)
      # seu$active <- callModule(embedParam, "negativeSampleRate", seu$active)
      removeModal()
    })


    callModule(findMarkers, "hello", seu)

    diffex_results <- callModule(diffex, "hello", seu, featureType, diffex_selected_cells)

    callModule(geneEnrichment, "hello", seu, diffex_results)
    observeEvent(input$plotTrx, {
      showModal(modalDialog(
        title = "Plotting Transcripts",
        "This process may take a minute or two!"
      ))
      callModule(allTranscripts, "hello", seu, featureType, organism_type)
      callModule(allTranscripts, "howdy", seu, featureType, organism_type)
      removeModal()
    })

    # callModule(rnaVelocity, "arrow", seu, featureType, "arrow")
    # callModule(rnaVelocity, "grid", seu, featureType, "grid")

    observe({
      updateSelectizeInput(session,
                           'geneSet',
                           choices = annotables::grch38$symbol,
                           server = TRUE)
    })

    observeEvent(input$regressAction,{
      req(seu$active)
      showModal(modalDialog(
        title = "Regressing out provided list of features",
        "This process may take a minute or two!"
      ))
      seu$gene <- seuratTools::regress_by_features(seu$gene, feature_set = list(input$geneSet), set_name = janitor::make_clean_names(input$geneSetName))
      # seu$gene <- regressed_seu$gene
      # seu$transcript <- regressed_seu$transcript
      seu$active <- seu[[input$featureType]]
      removeModal()
    })

    cds <- reactiveValues()

    observeEvent(input$calcCDS, {
      req(seu$active)
      cds$traj <- convert_seu_to_cds(seu$active, resolution = input$cdsResolution)
    })

    observeEvent(input$calCDS, {
      req(cds$traj)
      cds$traj <- learn_graph_by_resolution(cds$traj, seu$active, resolution = input$cdsResolution)
    })

    observe({
      shinyFiles::shinyFileChoose(input, "loadCDS", roots = volumes(), session = session)
    })


    cdsLoadPath <- eventReactive(input$loadCDS, {
      file <- shinyFiles::parseFilePaths(volumes(), input$loadCDS)
      file$datapath
    })

      observeEvent(input$loadCDS, {
        req(cdsLoadPath())
        shiny::withProgress(
          message = paste0("Uploading Data"),
          value = 0,
          {
            Sys.sleep(6)
            shiny::incProgress(2/10)
            # Sys.sleep(12)
            # shiny::incProgress(4/10)
            # Sys.sleep(18)
            # shiny::incProgress(6/10)
            # Sys.sleep(24)
            # shiny::incProgress(8/10)
            # browser()
            dataset <- readRDS(cdsLoadPath())
            shiny::incProgress(10/10)

            for (i in names(dataset)){
              cds[[i]] <- dataset[[i]]
            }

            # cds$traj <- dataset$traj
          }
        )
    })

    callModule(monocle, "arrow", cds, seu, featureType, resolution = reactive(input$cdsResolution))

    observe({
      shinyFiles::shinyFileSave(input, "saveCDS", roots = volumes(), session = session, restrictions = system.file(package = "base"))
    })


    cdsSavePath <- eventReactive(input$saveCDS, {
      savefile <- shinyFiles::parseSavePath(volumes(), input$saveCDS)

      savefile$datapath
    })

    observeEvent(input$saveCDS, {

      req(cds$traj)
      req(cdsSavePath())

      if (!is.null(cdsSavePath())){
        shiny::withProgress(
          message = paste0("Saving Data"),
          value = 0,
          {
            Sys.sleep(6)
            shiny::incProgress(2/10)
            # Sys.sleep(12)
            # shiny::incProgress(4/10)
            # Sys.sleep(18)
            # shiny::incProgress(6/10)
            # Sys.sleep(24)
            # shiny::incProgress(8/10)

            cdsList <- reactiveValuesToList(cds)
            names(cdsList) <- names(cds)


            saveRDS(cdsList, cdsSavePath())
            shiny::incProgress(10/10)

          }
        )
      }

    })


  }
  shinyApp(ui, server, enableBookmarking = "server")

}
