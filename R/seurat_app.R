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
  if (diffex_scheme == "seurat") {
    if ("integrated" %in% names(seu@assays)) {
      active_assay <- "integrated"
    } else {
      active_assay <- "RNA"
    }


    Idents(seu) <- paste0(active_assay, "_snn_res.", resolution)
    seu <- subset(seu, idents = c(cluster1, cluster2))
  } else if (diffex_scheme == "custom") {
    # subset by supplied cell ids
    #
    seu <- seu[, c(cluster1, cluster2)]

    keep_cells <- c(cluster1, cluster2)
    new_idents <- c(rep(1, length(cluster1)), rep(2, length(cluster2)))
    names(new_idents) <- keep_cells
    new_idents <- new_idents[colnames(seu)]
    Idents(seu) <- new_idents
    cluster1 <- 1
    cluster2 <- 2
  }

  test_list <- vector("list", length(tests))

  for (test in tests) {
    print(test)
    de <- FindMarkers(seu,
      ident.1 = cluster1,
      ident.2 = cluster2,
      test.use = test
    )


    if (featureType() == "transcript") {
      de_cols <- c("enstxp", "ensgene", "symbol", "p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj")

      de <-
        de %>%
        tibble::rownames_to_column("enstxp") %>%
        dplyr::left_join(annotables::grch38_tx2gene, by = "enstxp") %>%
        dplyr::left_join(annotables::grch38, by = "ensgene") %>%
        dplyr::select(one_of(de_cols))
    } else if (featureType() == "gene") {
      de_cols <- c("ensgene", "symbol", "p_val", "avg_logFC", "pct.1", "pct.2", "p_val_adj")

      de <-
        de %>%
        tibble::rownames_to_column("symbol") %>%
        dplyr::left_join(annotables::grch38, by = "symbol") %>%
        dplyr::select(one_of(de_cols))
    }

    test_list[[match(test, tests)]] <- de
  }
  names(test_list) <- tests
  return(test_list)
}


#' Run Enrichment Browser on Differentially Expressed Genes
#'
#' @param seu
#' @param enrichment_method
#' @param ...
#' @param cluster_list
#' @param de_results
#'
#' @return
#' @export
#'
#' @examples
run_enrichmentbrowser <- function(seu, cluster_list, de_results, enrichment_method = c("ora"), ...) {
  cluster1_cells <- cluster_list$cluster1
  cluster2_cells <- cluster_list$cluster2

  test_diffex_results <- de_results$t %>%
    dplyr::mutate(FC = log2(exp(avg_logFC))) %>%
    dplyr::mutate(ADJ.PVAL = p_val_adj) %>%
    dplyr::distinct(symbol, .keep_all = TRUE) %>%
    tibble::column_to_rownames("symbol") %>%
    identity()

  # subset by supplied cell ids
  #
  seu <- seu[, c(cluster1_cells, cluster2_cells)]
  seu <- seu[rownames(seu) %in% de_results$t$symbol, ]

  seu[["RNA"]]@meta.features <- test_diffex_results

  keep_cells <- c(cluster1_cells, cluster2_cells)
  new_idents <- c(rep(0, length(cluster1_cells)), rep(1, length(cluster2_cells)))
  names(new_idents) <- keep_cells
  new_idents <- new_idents[colnames(seu)]
  Idents(seu) <- new_idents

  counts <- GetAssayData(seu, assay = "RNA", slot = "counts")
  counts <- as.matrix(counts)
  mode(counts) <- "integer"

  rowData <- data.frame(
    FC = seu[["RNA"]][[]]$FC,
    ADJ.PVAL = seu[["RNA"]][[]]$ADJ.PVAL,
    row.names = rownames(seu@assays$RNA)
  )

  colData <- as.data.frame(seu[[]])

  se <- SummarizedExperiment::SummarizedExperiment(
    assays = list(counts = counts),
    rowData = rowData, colData = colData
  )

  se$GROUP <- forcats::fct_inseq(Idents(seu))

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


#' Create Seurat App
#'
#' @param preset_project A preloaded project to start the app with
#' @param filterTypes A named vector of file suffixes corresponding to subsets of the data
#' @param appTitle A title of the App
#' @param futureMb amount of Mb allocated to future package
#' @param preset_project
#' @param feature_types
#' @param organism_type
#'
#' @return
#' @export
#'
#' @examples
seuratApp <- function(preset_project, filterTypes, appTitle = NULL, feature_types = "gene",
                      organism_type = "human", db_path = "/dataVolume/storage/single_cell_projects/rsqlite/single-cell-projects.db", futureMb = 13000) {
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
  sidebar <- shinydashboard::dashboardSidebar(uiOutput("projInput"),
    actionButton("loadProject", "Load Selected Project"),
    shinyFiles::shinyDirButton(
      "deleteProject", "Delete an Integrated Project or Dataset",
      "Please select a file or directory to delete"
    ),
    textOutput("appTitle"), uiOutput("featureType"), shinyWidgets::prettyRadioButtons("organism_type",
      "Organism",
      choices = c("human", "mouse"), selected = organism_type
    ),
    shinyFiles::shinyFilesButton("seuratUpload", "Load a Seurat Dataset",
      "Please select a .rds file",
      multiple = FALSE
    ),
    shinyFiles::shinySaveButton("saveSeurat", "Save Current Dataset",
      "Save file as...",
      filetype = list(rds = "rds")
    ),
    verbatimTextOutput("savefile"), actionButton("changeEmbedAction",
      label = "Change Embedding Parameters"
    ), changeEmbedParamsui("changeembed"),
    shinydashboard::sidebarMenu(shinydashboard::menuItem("Integrate Projects",
      tabName = "integrateProjects"
    ), shinydashboard::menuItem("Reformat Metadata",
      tabName = "reformatMetadata"
    ), shinydashboard::menuItem("Compare Scatter Plots",
      tabName = "comparePlots"
    ), shinydashboard::menuItem("Compare Read Counts",
      tabName = "compareReadCount"
    ), shinydashboard::menuItem("Violin/Heatmap Plots",
      tabName = "violinPlots"
    ), shinydashboard::menuItem("Differential Expression",
      tabName = "diffex"
    ), shinydashboard::menuItem("Gene Enrichment Analysis",
      tabName = "geneEnrichment"
    ), shinydashboard::menuItem("Find Markers",
      tabName = "findMarkers"
    ), shinydashboard::menuItem("Subset Seurat Input",
      tabName = "subsetSeurat"
    ), shinydashboard::menuItem("All Transcripts",
      tabName = "allTranscripts"
    ), shinydashboard::menuItem("RNA Velocity",
      tabName = "plotVelocity"
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
        box(plotDimRedui("plotdimred1")),
        box(plotDimRedui("plotdimred2"))
      ), fluidRow(box(
        title = "Selected Cells",
        tableSelectedui("tableselected"), width = 6
      ), box(plotClustree_UI("clustreePlot")))
    ),
    shinydashboard::tabItem(
      tabName = "integrateProjects",
      h2("Integrate Projects"), fluidRow((integrateProjui("integrateproj")))
    ),
    shinydashboard::tabItem(
      tabName = "reformatMetadata",
      h2("Reformat Metadata"), fluidRow((reformatMetadataui("reformatmetadata")))
    ),
    shinydashboard::tabItem(
      tabName = "compareReadCount",
      h2("Compare Read Counts"), fluidRow(
        box(plotReadCountui("plotreadcount1")),
        box(plotReadCountui("plotreadcount2"))
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
        box(findMarkersui("findmarkers")),
        box(plotDimRedui("markerScatter"))
      )
    ), shinydashboard::tabItem(
      tabName = "allTranscripts",
      h2("All Transcripts"), fluidRow(shinyWidgets::actionBttn(
        "plotTrx",
        "Plot all transcripts"
      )), fluidRow(column(allTranscriptsui("alltranscripts1"),
        width = 6
      ), column(allTranscriptsui("alltranscripts2"),
        width = 6
      ))
    ),
    shinydashboard::tabItem(
      tabName = "plotVelocity",
      h2("RNA Velocity"),
      fluidRow(
        box(
          plotVelocityui("plotvelocity")
        )
      )
    ),
    shinydashboard::tabItem(
      tabName = "diffex",
      h2("Differential Expression"), column(box(plotDimRedui("diffex"),
        width = 12
      ), box(tableSelectedui("diffex"),
        width = 12
      ), width = 6), column(box(diffexui("diffex"),
        width = 12
      ), width = 6)
    ), shinydashboard::tabItem(
      tabName = "geneEnrichment",
      h2("Gene Enrichment"), geneEnrichmentui("geneenrichment"),
      downloadTable_UI("downloadtable")
    ), shinydashboard::tabItem(
      tabName = "regressFeatures",
      h2("Regress Features"), fluidRow(actionButton(
        "regressAction",
        "Regress Seurat Objects By Genes"
      ), box(checkboxInput("runRegression",
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
  server <- function(input, output, session) {
    options(warn = -1)
    shinylogs::track_usage(storage_mode = shinylogs::store_json(path = "logs/"))
    # projects_db <- "/dataVolume/storage/single_cell_projects/single_cell_projects.db"
    rsqlite_db <- db_path

    con <- DBI::dbConnect(
      RSQLite::SQLite(),
      rsqlite_db
    )

    projList <- reactivePoll(4000, session, checkFunc = function() {
      if (file.exists(rsqlite_db)) {
        # system("locate -d /dataVolume/storage/single_cell_projects/single_cell_projects.db '*.here'",
        #   intern = TRUE
        # ) %>%
        #   fs::path_dir() %>%
        #   purrr::set_names(fs::path_file(.))
        DBI::dbReadTable(con, "projects") %>%
          tibble::deframe()
      } else {
        ""
      }
    }, valueFunc = function() {
      # system("locate -d /dataVolume/storage/single_cell_projects/single_cell_projects.db '*.here'",
      #   intern = TRUE
      # ) %>%
      #   fs::path_dir() %>%
      #   purrr::set_names(fs::path_file(.))
      DBI::dbReadTable(con, "projects") %>%
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
    seu <- reactiveValues()
    proj_dir <- reactiveVal()
    if (!is.null(preset_project)) {
      proj_dir(preset_project)
    }
    organism_type <- reactive({
      input$organism_type
    })
    plot_types <- reactive({
      list_plot_types(seu$active)
    })
    observeEvent(input$loadProject, {
      proj_dir(input$setProject)
    })
    output$appTitle <- renderText({
      req(proj_dir())
      paste0("Loaded Project: ", fs::path_file(proj_dir()))
    })
    project_volumes <- reactive({
      print(proj_dir())
      project_volumes <- c(
        Home = fs::path("/dataVolume/storage/single_cell_projects/integrated_projects"),
        `R Installation` = R.home(), shinyFiles::getVolumes()
      )
    })
    dataset_volumes <- reactive({
      print(proj_dir())
      dataset_volumes <- c(
        Home = fs::path(
          proj_dir(),
          "output", "seurat"
        ), `R Installation` = R.home(),
        shinyFiles::getVolumes()
      )
    })
    observe({
      req(dataset_volumes())
      shinyFiles::shinyFileChoose(input, "seuratUpload",
        roots = dataset_volumes(), session = session
      )
    })
    observe({
      req(project_volumes())
      shinyFiles::shinyDirChoose(input, "deleteProject",
        roots = project_volumes(), session = session,
        restrictions = system.file(package = "base")
      )
    })
    uploadSeuratPath <- eventReactive(input$seuratUpload, {
      req(dataset_volumes())
      file <- shinyFiles::parseFilePaths(
        dataset_volumes(),
        input$seuratUpload
      )
      file$datapath
    })
    deleteSeuratPath <- eventReactive(input$deleteProject, {
      req(project_volumes())
      file <- shinyFiles::parseDirPath(
        project_volumes(),
        input$deleteProject
      )
    })
    observeEvent(input$seuratUpload, {
      req(uploadSeuratPath())
      shiny::withProgress(
        message = paste0("Uploading Data"),
        value = 0,
        {
          shiny::incProgress(2 / 10)
          print(uploadSeuratPath())
          dataset <- readRDS(uploadSeuratPath())
          shiny::incProgress(6 / 10)
          seu_names <- names(dataset)[!names(dataset) ==
            "active"]
          for (i in seu_names) {
            seu[[i]] <- dataset[[i]]
          }
          seu$active <- seu[["gene"]]
          print(uploadSeuratPath())
          print(names(seu))
          shiny::incProgress(8 / 10)
        }
      )
    })
    observeEvent(input$deleteProject, {
      req(deleteSeuratPath())
      message <- paste0("Deleting Project")
      print(deleteSeuratPath())
      fs::file_delete(deleteSeuratPath())
      showModal(modalDialog(
        title = "Project Deleted",
        paste0("You successfully deleted: ", deleteSeuratPath()),
        easyClose = TRUE, footer = NULL
      ))
    })
    output$featureType <- renderUI({
      req(seu)
      seu_names <- names(seu)[!(names(seu) == "active")]
      shinyWidgets::prettyRadioButtons("feature_type",
        "Feature for Display",
        choices = seu_names,
        selected = "gene", inline = T
      )
    })
    observeEvent(input$feature_type, {
      seu$active <- seu[[input$feature_type]]
    })
    featureType <- reactive({
      featureType <- input$feature_type
    })
    observe({
      shinyFiles::shinyFileSave(input, "saveSeurat",
        roots = dataset_volumes(),
        session = session, restrictions = system.file(package = "base")
      )
    })
    subSeuratPath <- eventReactive(input$saveSeurat, {
      req(seu$active)
      savefile <- shinyFiles::parseSavePath(
        dataset_volumes(),
        input$saveSeurat
      )
      return(savefile$datapath)
    })
    observeEvent(input$saveSeurat, {
      req(seu$active)
      req(subSeuratPath())
      shiny::withProgress(
        message = paste0("Saving Data"),
        value = 0,
        {
          Sys.sleep(6)
          shiny::incProgress(2 / 10)
          saveRDS(
            shiny::reactiveValuesToList(seu),
            subSeuratPath()
          )
          shiny::incProgress(10 / 10)
        }
      )
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
    observe({
      print(newprojList())
      updateSelectizeInput(session, "setProject",
        label = "Select input label",
        choices = newprojList(),
        server = TRUE
      )
    })
    seu <- callModule(reformatMetadata, "reformatmetadata", seu)

    observe({
      req(seu$active)

      reductions <- reactive({
        # names(seu$active@reductions)
        c("pca", "tsne", "umap")
      })

      callModule(plotDimRed, "plotdimred1", seu, plot_types, featureType,
                 organism_type = organism_type, reductions)
      callModule(plotDimRed, "plotdimred2", seu, plot_types, featureType,
                 organism_type = organism_type, reductions)
      # callModule(plotDimRed, "plotdimred", seu, plot_types, featureType,
      #            organism_type = organism_type, reductions)
      callModule(plotDimRed, "diffex", seu, plot_types, featureType,
                 organism_type = organism_type, reductions)
      callModule(plotDimRed, "subset", seu, plot_types, featureType,
                 organism_type = organism_type, reductions)
      callModule(plotDimRed, "markerScatter", seu, plot_types, featureType,
                 organism_type = organism_type, reductions)
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
    observeEvent(input$subsetAction, {
      req(subset_selected_cells())
      withCallingHandlers(
        {
          shinyjs::html("subsetMessages", "")
          message("Beginning")
          for (i in names(seu)[!(names(seu) == "active")]) {
            seu[[i]] <- seu[[i]][, colnames(seu[[i]]) %in% subset_selected_cells()]
          }
          if (length(unique(seu$gene[[]]$batch)) > 1) {
            print(names(seu)[!(names(seu) == "active")])
            for (i in names(seu)[!(names(seu) == "active")]) {
              message(paste0("reintegrating ", i, " expression"))
              seu[[i]] <- reintegrate_seu(seu[[i]],
                feature = i,
                resolution = seq(0.2, 2, by = 0.2)
              )
            }
          }
          else {
            for (i in names(seu)[!(names(seu) == "active")]) {
              seu[[i]] <- seurat_pipeline(seu[[i]], resolution = seq(0.2,
                2,
                by = 0.2
              ))
            }
          }
          seu$active <- seu[[input$feature_type]]
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
          for (i in names(seu)[!(names(seu) == "active")]) {
            seu[[i]] <- subset_by_meta(
              input$uploadCsv$datapath,
              seu[[i]]
            )
          }
          if (length(unique(seu$gene[[]]$batch)) > 1) {
            for (i in names(seu)[!(names(seu) == "active")]) {
              message(paste0("reintegrating ", i, " expression"))
              seu[[i]] <- reintegrate_seu(seu[[i]],
                feature = i,
                resolution = seq(0.2, 2, by = 0.2)
              )
            }
          }
          else {
            for (i in names(seu)[!(names(seu) == "active")]) {
              seu[[i]] <- seurat_pipeline(seu[[i]], resolution = seq(0.2,
                2,
                by = 0.2
              ))
            }
          }
          seu$active <- seu[[input$feature_type]]
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
    callModule(findMarkers, "findmarkers", seu)
    diffex_results <- callModule(
      diffex, "diffex", seu, featureType,
      diffex_selected_cells
    )
    enrichment_report <- callModule(
      geneEnrichment, "geneenrichment",
      seu, diffex_results
    )
    observe({
      req(enrichment_report())
      callModule(downloadTable, "downloadtable", enrichment_report)
    })

    observeEvent(input$plotTrx, {
      showModal(modalDialog(
        title = "Plotting Transcripts",
        "This process may take a minute or two!"
      ))
      callModule(
        allTranscripts, "alltranscripts1", seu, featureType,
        organism_type
      )
      callModule(
        allTranscripts, "alltranscripts2", seu, featureType,
        organism_type
      )
      removeModal()
    })

    prior_gene_set <- reactive({
      # req(input$priorGeneSet)
      req(seu$active)

      if (is.null(input$priorGeneSet)){
        ""
      } else if (input$priorGeneSet == "Apoptosis") {
        marker_genes <- c(
          "CASP3", "CASP7", "BAX", "BAK1", "BID", "BBC3",
          "BCL2", "MCL1"
        )

        marker_genes[marker_genes %in% rownames(seu$active)]
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
          "E2F8")

          marker_genes[marker_genes %in% rownames(seu$active)]
      }  else if (input$priorGeneSet == "Mitochondrial") {
        marker_genes <- mito_features[[organism_type()]][["gene"]]

        marker_genes[marker_genes %in% rownames(seu$active)]
      } else if (input$priorGeneSet == "Ribosomal") {
        marker_genes <- ribo_features[[organism_type()]][["gene"]]

        marker_genes[marker_genes %in% rownames(seu$active)]
      }
    })

    observe({

      updateSelectizeInput(session, "geneSet",
        choices = rownames(seu$active),
        selected = prior_gene_set(),
        server = TRUE
      )
    })
    observeEvent(input$regressAction, {
      req(seu$active)
      showModal(modalDialog(
        title = "Regressing out provided list of features",
        "This process may take a minute or two!"
      ))
      seu$gene <- seuratTools::regress_by_features(seu$gene,
        feature_set = list(input$geneSet), set_name = janitor::make_clean_names(input$geneSetName),
        regress = input$runRegression
      )
      seu$active <- seu[[input$feature_type]]
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

    observe({
      req(uploadSeuratPath())
      req(seu)

      proj_path <- stringr::str_replace(uploadSeuratPath(), "output.*", "")

      proj_name <- fs::path_file(proj_path)
      print(proj_name)

      loom_path <- fs::path(proj_path, "output", "velocyto", paste0(proj_name, ".loom"))

      print(loom_path)

      callModule(plotVelocity, "plotvelocity", seu, loom_path, featureType)
    })


    observe({
      shinyFiles::shinyFileSave(input, "saveCDS",
        roots = dataset_volumes(),
        session = session, restrictions = system.file(package = "base")
      )
    })
    cdsSavePath <- eventReactive(input$saveCDS, {
      savefile <- shinyFiles::parseSavePath(
        volumes(),
        input$saveCDS
      )
      savefile$datapath
    })
    observeEvent(input$saveCDS, {
      req(cds$traj)
      req(cdsSavePath())
      if (!is.null(cdsSavePath())) {
        shiny::withProgress(
          message = paste0("Saving Data"),
          value = 0,
          {
            Sys.sleep(6)
            shiny::incProgress(2 / 10)
            cdsList <- reactiveValuesToList(cds)
            names(cdsList) <- names(cds)
            saveRDS(cdsList, cdsSavePath())
            shiny::incProgress(10 / 10)
          }
        )
      }
    })
  }
  shinyApp(ui, server, enableBookmarking = "server")
}
