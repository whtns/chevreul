
#' Plot RNA velocity Computed by Velocyto.R
#'
#' @param vel velocity from velocyto.R
#' @param emb embedding -- either tnse or umap
#' @param cell.colors A vector of colors to use for the plot
#' @param format either an arrow or grid format
#'
#' @return
#' @export
#'
#' @examples
plot_velocity <- function(vel, emb, cell.colors, format = "arrow") {

  arrow.scale=3; cell.alpha=1.0; cell.cex=1; fig.height=4; fig.width=4.5;

  if (format == "arrow") {
    velocyto.R::show.velocity.on.embedding.cor(emb, vel, n=100, scale='sqrt',
                                   cell.colors=velocyto.R::ac(cell.colors, alpha=cell.alpha),
                                   cex=cell.cex, arrow.scale=arrow.scale, arrow.lwd=1)
  } else if (format == "grid"){
    #Alternatively, the same function can be used to calculate a velocity vector field:
    velocyto.R::show.velocity.on.embedding.cor(emb, vel, n=100, scale='sqrt',
                                   cell.colors=velocyto.R::ac(cell.colors, alpha=cell.alpha),
                                   cex=cell.cex, arrow.scale=arrow.scale,
                                   show.grid.flow=TRUE, min.grid.cell.mass=0.5,
                                   grid.n=20, arrow.lwd=2)
  }


}

#' Plot Metadata Variables
#'
#' @param seu
#' @param embedding
#' @param group
#'
#' @return
#' @export
#'
#' @examples
plot_var <- function(seu, embedding = "umap", group = "batch", dims = c(1,2)){
  #
  metadata <- as_tibble(seu[[]][Seurat::Cells(seu),], rownames = "sID")
  cellid <- metadata[["sID"]]
  key <- rownames(metadata)

  if (embedding == "umap"){
    dims = c(1,2)
  } else if (embedding == "tsne"){
    dims = c(1,2)
  }

  dims <- as.numeric(dims)

  d <- Seurat::DimPlot(object = seu, dims = dims, reduction = embedding, group.by = group) +
    aes(key = key, cellid = cellid)

  plotly::ggplotly(d, tooltip = "cellid", height  = 750) %>%
    plotly::layout(dragmode = "lasso") %>%
    identity()

}

#' Plot Features
#'
#' @param seu
#' @param embedding
#' @param features
#' @param dims
#'
#' @return
#' @export
#'
#' @examples
plot_feature <- function(seu, embedding, features, dims = c(1,2)){

  metadata <- as_tibble(seu[[]][Seurat::Cells(seu),], rownames = "sID")

  cellid <- metadata[["sID"]]
  key <- rownames(metadata)

  if (embedding == "umap"){
    dims = c(1,2)
  } else if (embedding == "tsne"){
    dims = c(1,2)
  }

  dims <- as.numeric(dims)

  fp <- Seurat::FeaturePlot(object = seu, features = features, dims = dims, reduction = embedding)	+
    aes(key = key, cellid = cellid)

  plotly::ggplotly(fp, tooltip = "cellid", height = 750) %>%
    plotly::layout(dragmode = "lasso") %>%
    identity()

}

#' Plot Rides
#'
#' @param seu
#' @param features
#'
#' @return
#' @export
#'
#' @examples
plot_ridge <- function(seu, features){

  cc_genes_path <- "~/single_cell_projects/resources/regev_lab_cell_cycle_genes.txt"
  cc.genes <- readLines(con = cc_genes_path)
  s.genes <- cc.genes[1:43]
  g2m.genes <- cc.genes[44:97]

  seu <- CellCycleScoring(object = seu, s.genes, g2m.genes,
                          set.ident = TRUE)

  RidgePlot(object = seu, features = features)

  # plotly::ggplotly(r, height = 750)
  #
}


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
run_seurat_de <- function(seu, cluster1, cluster2, resolution, diffex_scheme = "seurat") {

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

  tests <- c("t", "wilcox", "bimod")
  test_list <- vector("list", length(tests))

  for (test in tests){
    print(test)
    de <- FindMarkers(seu,
                      ident.1 = cluster1,
                      ident.2 = cluster2,
                      test.use = test)
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

#' Plot Cluster Marker Genes
#'
#' @param seu
#' @param resolution
#'
#' @return
#' @export
#'
#' @examples
plot_markers <- function(seu, resolution){

  markerplot <- DotPlot(seu, features = unique(seu@misc$markers[[resolution]]), group.by = resolution) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))

  plotly::ggplotly(markerplot, height = 800) %>%
    plotly::layout(dragmode = "lasso")

}

#' Plot Read Count
#'
#' @param seu
#' @param plot_type
#'
#' @return
#' @export
#'
#' @examples
plot_readcount <- function(seu, plot_type){

  seu_tbl <- tibble::rownames_to_column(seu[[]], "SID")

  rc_plot <- ggplot(seu_tbl, aes(x=reorder(SID, -nCount_RNA), y = nCount_RNA, fill = !!as.symbol(plot_type))) +
    # scale_y_continuous(breaks = seq(0, 8e7, by = 5e5)) +
    scale_y_log10() +
    geom_bar(position = "identity", stat = "identity") +
    # geom_text(data=subset(agg_qc_wo_na, Sample %in% thresholded_cells & align_type == "paired_total"),
    #   aes(Sample, count, label=Sample)) +
    # geom_text(data=subset(agg_qc_wo_na, Sample %in% low_read_count_cells & align_type == "paired_aligned_one"),
    #   aes(Sample, count, label=Sample)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    # scale_fill_manual(values = c( "low_read_count"="tomato", "keep"="gray" ), guide = FALSE ) +
    labs(title = "Paired Unique Reads", x = "Sample") +
    NULL

  rc_plot <- plotly::ggplotly(rc_plot)
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
#' @param proj_dir The project directory of the base dataset
#' @param filterTypes A named vector of file suffixes corresponding to subsets of the data
#' @param appTitle A title of the App
#' @param futureMb amount of Mb allocated to future package
#'
#' @return
#' @export
#'
#' @examples
seuratApp <- function(proj_dir, filterTypes, appTitle, feature_types = "gene", futureMb = 849){

  futureMb = 850
  future::plan(strategy = "multicore", workers = 6)
  future_size = futureMb*1024^2
  options(future.globals.maxSize= future_size)
  options(DT.options = list(pageLength = 2000, paging = FALSE,
                            info = TRUE, searching = TRUE, autoWidth = F, ordering = TRUE,
                            language = list(search = "Filter:")))

  header <- shinydashboard::dashboardHeader(title = appTitle)

  sidebar <- shinydashboard::dashboardSidebar(shinyWidgets::prettyRadioButtons("feature_type", "Feature for Display", choices = feature_types, selected = "gene"),
                                              shinyWidgets::prettyRadioButtons("organism_type", "Organism", choices = c("human", "mouse"), selected = "human"),
                                              shinyFiles::shinyFilesButton("seuratUpload", "Load a Seurat Dataset", "Please select a .rds file", multiple = FALSE),
                                              shinyFiles::shinySaveButton("saveSeurat", "Save current Dataset", "Save file as...", filetype = list(rds = "rds")),
                                              verbatimTextOutput("savefile"),
                                              actionButton("changeEmbedAction", label = "Change Embedding Parameters"),
                                              changeEmbedParamsui("changeembed"),
                                              shinydashboard::sidebarMenu(shinydashboard::menuItem("Compare Plots", tabName = "comparePlots"),
                                                                          shinydashboard::menuItem("Compare Read Counts", tabName = "compareReadCount"),
                                                                          shinydashboard::menuItem("Differential Expression", tabName = "diffex"),
                                                                          shinydashboard::menuItem("Gene Enrichment Analysis", tabName = "geneEnrichment"),
                                                                          shinydashboard::menuItem("Find Markers", tabName = "findMarkers"),
                                                                          shinydashboard::menuItem("Subset Seurat Input", tabName = "subsetSeurat"),
                                                                          shinydashboard::menuItem("All Transcripts", tabName = "allTranscripts"),
                                                                          # shinydashboard::menuItem("RNA Velocity", tabName = "rnaVelocity"),
                                                                          shinydashboard::menuItem("Monocle", tabName = "monocle"),
                                                                          shinydashboard::menuItem("cellAlign", tabName = "cellAlign"),
                                                                          shinydashboard::menuItem("Regress Features", tabName = "regressFeatures")),

                                              width = 450)



  body <- shinydashboard::dashboardBody(shinydashboard::tabItems(
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
        "subset seurat"
      ),
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
    fluidRow(box(findMarkersui("hello")))
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
    h2("Differential Expression"), column(diffexui("hello"),
      width = 6
    ), column(box(plotDimRedui(
      "diffex"
    ), width = 12),
    box(
      title = "Selected Cells", tableSelectedui("diffex"),
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
  # shinydashboard::tabItem(
  #   tabName = "rnaVelocity",
  #   h2("RNA Velocity"),
  #   fluidRow(
  #     box(rnaVelocityui("arrow"),
  #         width = 12)),
  #   fluidRow(
  #     box(rnaVelocityui("grid"),
  #         width = 12))
  # ),
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
  ),
  shinydashboard::tabItem(
    tabName = "cellAlign",
    h2("cellalign")
  )
))

  ui <- dashboardPage(
    header = header,
    sidebar = sidebar,
    body = body
  )

  server <- function(input, output, session) {
    options(warn = -1)

    # sidebar
    # seu <- callModule(loadData, "loadDataui", proj_dir, feature_types = feature_types, selected_feature = input$feature_type)

    seu <- reactiveValues()

    observeEvent(input$feature_type, {
      seu$active <- seu[[input$feature_type]]
    })

    feature_type <- reactive({
      input$feature_type
    })

    plot_types <- reactive({
      list_plot_types(seu$active)
    })

    # upload seurat object

    # list volumes
    volumes <- c(Home = fs::path(proj_dir, "output", "seurat"), "R Installation" = R.home(), shinyFiles::getVolumes())

    shinyFiles::shinyFileChoose(input, "seuratUpload", roots = volumes, session = session)
    uploadSeuratPath <- eventReactive(input$seuratUpload, {
      file <- shinyFiles::parseFilePaths(volumes, input$seuratUpload)
      file$datapath
    })

    observe({
      req(uploadSeuratPath())

      shiny::withProgress(
        message = paste0("Uploading Data"),
        value = 0,
        {
          # Sys.sleep(6)
          # shiny::incProgress(2/10)
          # Sys.sleep(12)
          # shiny::incProgress(4/10)
          # Sys.sleep(18)
          # shiny::incProgress(6/10)
          # Sys.sleep(24)
          # shiny::incProgress(8/10)
          # browser()
          dataset <- readRDS(uploadSeuratPath())
          seu$gene <- dataset$gene
          seu$transcript <- dataset$transcript
          seu$active <- dataset$gene
        }
      )

    })

    # save seurat object
    shinyFiles::shinyFileSave(input, "saveSeurat", roots = volumes, session = session, restrictions = system.file(package = "base"))

    subSeuratPath <- eventReactive(input$saveSeurat, {

      req(seu$active)
      savefile <- shinyFiles::parseSavePath(volumes, input$saveSeurat)

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
          saveRDS(list(gene = seu$gene, transcript = seu$transcript), subSeuratPath())
          shiny::incProgress(10/10)

        }
      )
    })


    # observeEvent(input$subSeuratPath, { print(parseSavePath(volumes, input$subSeuratPath)$datapath) })
    #
    # shiny::withProgress(
    #   message = paste0("Downloading Data"),
    #   value = 0,
    #   {
    #     shiny::incProgress(1/10)
    #     Sys.sleep(1)
    #     shiny::incProgress(5/10)
    #     saveRDS(list(gene = seu$gene, transcript = seu$transcript), fileinfo$datapath)
    #
    #   }
    # )


    # body
    callModule(plotDimRed, "hello", seu, plot_types, feature_type, organism_type = input$organism_type)
    callModule(plotDimRed, "howdy", seu, plot_types, feature_type, organism_type = input$organism_type)
    callModule(plotDimRed, "diffex", seu, plot_types, feature_type, organism_type = input$organism_type)
    callModule(plotDimRed, "subset", seu, plot_types, feature_type, organism_type = input$organism_type)
    callModule(plotReadCount, "hello2", seu, plot_types)
    callModule(plotReadCount, "howdy2", seu, plot_types)
    callModule(tableSelected, "hello", seu)
    callModule(tableSelected, "diffex", seu)
    selected_cells <- callModule(tableSelected, "subset", seu)
    observeEvent(input$subsetAction, {

      withCallingHandlers({
        shinyjs::html("subsetMessages", "")
        message("Beginning")

        seu$gene <- seu$gene[, selected_cells()]
        seu$transcript <- seu$transcript[, selected_cells()]

        if(length(unique(seu$gene[[]]$batch)) > 1){

          for (i in names(seu)[!names(seu) == "active"]){
            message(paste0("reintegrating ", i, " expression"))
            # harmony
            # seu[[i]] <- seurat_pipeline(seu[[i]], reduction = "harmony", resolution = seq(0.2, 2.0, by = 0.2))
            # seurat cca
            seu[[i]] <- reintegrate_seu(seu[[i]], feature = i, resolution = seq(0.2, 2.0, by = 0.2))
          }

        } else {

          for (i in names(seu)[!names(seu) == "active"]){
            seu[[i]] <- seuratTools::seurat_pipeline(seu[[i]], resolution = seq(0.2, 2.0, by = 0.2))
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
    diffex_results <- callModule(diffex, "hello", seu)
    callModule(geneEnrichment, "hello", seu, diffex_results)
    observeEvent(input$plotTrx, {
      showModal(modalDialog(
        title = "Plotting Transcripts",
        "This process may take a minute or two!"
      ))
      callModule(allTranscripts, "hello", seu, feature_type)
      callModule(allTranscripts, "howdy", seu, feature_type)
      removeModal()
    })

    # callModule(rnaVelocity, "arrow", seu, feature_type, "arrow")
    # callModule(rnaVelocity, "grid", seu, feature_type, "grid")

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
      regressed_seu <- seuratTools::regress_by_gene_set(list(gene = seu$gene, transcript = seu$transcript), gene_set = list(input$geneSet), set_name = janitor::make_clean_names(input$geneSetName))
      seu$gene <- regressed_seu$gene
      seu$transcript <- regressed_seu$transcript
      seu$active <- seu$gene
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

    shinyFiles::shinyFileChoose(input, "loadCDS", roots = volumes, session = session)

    cdsLoadPath <- eventReactive(input$loadCDS, {
      file <- shinyFiles::parseFilePaths(volumes, input$loadCDS)
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

    callModule(monocle, "arrow", cds, seu, feature_type, resolution = reactive(input$cdsResolution))

    shinyFiles::shinyFileSave(input, "saveCDS", roots = volumes, session = session, restrictions = system.file(package = "base"))

    cdsSavePath <- eventReactive(input$saveCDS, {
      savefile <- shinyFiles::parseSavePath(volumes, input$saveCDS)

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

    callModule(cellAlign, "hello")


  }
  shinyApp(ui, server)
}
