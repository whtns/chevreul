
#' Plot RNA velocity Computed by Velocyto.R
#'
#' @param vel
#' @param emb
#' @param cell.colors
#' @param format
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
plot_var <- function(seu, embedding = "umap", group = "batch"){
  #
  metadata <- as_tibble(seu[[]][Seurat::Cells(seu),], rownames = "sID")
  cellid <- metadata[["sID"]]
  key <- rownames(metadata)

  d <- Seurat::DimPlot(object = seu, reduction = embedding, group.by = group) +
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
#'
#' @return
#' @export
#'
#' @examples
plot_feature <- function(seu, embedding, features){

  metadata <- as_tibble(seu[[]][Seurat::Cells(seu),], rownames = "sID")

  cellid <- metadata[["sID"]]
  key <- rownames(metadata)

  fp <- Seurat::FeaturePlot(object = seu, reduction = embedding, features = features)	+
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
    Idents(seu) <- paste0("clusters_", resolution)
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

#' TPlot Cluster Marker Genes
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
  #
  rc_plot <- ggplot(data.frame(seu[[]]), aes(x=reorder(Sample_ID, -nCount_RNA), y = nCount_RNA, fill = !!as.symbol(plot_type))) +
    # scale_y_continuous(breaks = seq(0, 8e7, by = 5e5)) +
    scale_y_log10() +
    geom_bar(position = "identity", stat = "identity") +
    # geom_text(data=subset(agg_qc_wo_na, Sample %in% thresholded_cells & align_type == "paired_total"),
    #   aes(Sample, count, label=Sample)) +
    # geom_text(data=subset(agg_qc_wo_na, Sample %in% low_read_count_cells & align_type == "paired_aligned_one"),
    #   aes(Sample, count, label=Sample)) +
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1)) +
    # scale_fill_manual(values = c( "low_read_count"="tomato", "keep"="gray" ), guide = FALSE ) +
    labs(title = "Paired Aligned One Reads", x = "Sample") +
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
#' @param proj_dir The project directory of the base dataset ex. "~/single_cell_projects/sc_cone_devel/proj"
#' @param plot_types The types of plots to be shown as a named list containing two vectors: 1) category_vars and 2) continuous_vars
#' @param filterTypes A named vector of file suffixes corresponding to subsets of the data, ex. filterTypes <- c("", "remove_lowrc") %>% set_names(c("Unfiltered", "low read count cells"))
#' @param appTitle A title of the app
#' @param futureMb amount of Mb allocated to future package
#'
#' @return
#' @export
#'
#' @examples
seuratApp <- function(proj_dir, plot_types, filterTypes, appTitle, futureMb = 850, ...){

  futureMb = 850
  future::plan(strategy = "multicore", workers = 6)
  future_size = futureMb*1024^2
  options(future.globals.maxSize= future_size)
  options(DT.options = list(pageLength = 2000, paging = FALSE,
                            info = TRUE, searching = TRUE, autoWidth = F, ordering = TRUE,
                            language = list(search = "Filter:")))

  header <- shinydashboard::dashboardHeader(title = appTitle)


  loadDataui <- function(id, label = "Load Data", filterTypes) {
    ns <- NS(id)
    tagList(
      shinyWidgets::prettyRadioButtons(ns("filterType"), "dataset to include", choices = filterTypes, selected = ""),
      shinyWidgets::actionBttn(ns("loadButton"), "Load Default Dataset")
      # fileInput(ns("seuratUpload"), "Upload .rds file")
    )
  }

  loadData <- function(input, output, session, proj_dir, feature_type) {
    ns <- session$ns
    seu <- reactiveValues()

    observeEvent(input$loadButton, {
      showModal(modalDialog("Loading Data", footer = NULL))
      if (!input$filterType == "") {
        filterType = paste0("_", input$filterType)
      }
      else {
        filterType = input$filterType
      }
      seu_paths <- rprojroot::find_root_file("output/sce", criterion = rprojroot::has_file_pattern("*.Rproj"),
                                             path = proj_dir) %>% fs::dir_ls() %>% fs::path_filter(paste0("*_seu",
                                                                                                          filterType, ".rds")) %>% identity()
      feature_seus <- purrr::map(seu_paths, readRDS) %>%
        purrr::set_names(c("gene", "transcript"))
      Seurat::DefaultAssay(feature_seus$gene) <- "RNA"
      Seurat::DefaultAssay(feature_seus$transcript) <- "RNA"
      removeModal()
      seu$gene <- feature_seus$gene
      seu$transcript <- feature_seus$transcript
      seu$active <- feature_seus[[feature_type]]
    })

    return(seu)

  }

  customFeatureui <- function(id) {
    ns <- NS(id)
    tagList(choice <- reactive({
      if (input$feature_type == "transcript") {
        def_text <- "ENST00000488147"
      }
      else if (input$feature_type == "gene") {
        def_text <- "RXRG"
      }
    }))
  }
  customFeature <- function(input, output, session) {
    ns <- session$ns
    output$featuretext <- renderUI({
      textInput("feature", "gene or transcript on which to color the plot; eg. 'RXRG' or 'ENST00000488147'",
                value = choice())
    })
    uiOutput("featuretext")
  }

  changeEmbedParamsui <- function(id){
    ns <- NS(id)

    minDist_vals <- prep_slider_values(0.3)
    negsamprate_vals <- prep_slider_values(5)

    tagList(
      sliderInput(ns("minDist"), label = "Minimum Distance", min = minDist_vals$min, max = minDist_vals$max, value = minDist_vals$value, step = minDist_vals$step),
      sliderInput(ns("negativeSampleRate"), label = "NegativeSampleRate", min = negsamprate_vals$min, max = negsamprate_vals$max, value = negsamprate_vals$value, step = negsamprate_vals$step)
    )
  }

  changeEmbedParams <- function(input, output, session, seu){
    ns <- session$ns
    #
    # output$embedControls <- renderUI({
    #   tagList(
    #     sliderInput(ns("minDist"), label = "Minimum Distance", min = minDist_vals$min, max = minDist_vals$max, value = minDist_vals$value, step = minDist_vals$step),
    #     sliderInput(ns("negativeSampleRate"), label = "NegativeSampleRate", min = minDist_vals$min, max = minDist_vals$max, value = minDist_vals$value, step = minDist_vals$step)
    #   )
    # })

    seu$gene <- RunUMAP(seu$gene, dims = 1:30, reduction = "pca", min.dist = input$minDist, negative.sample.rate = input$negativeSampleRate)
    seu$transcript <- RunUMAP(seu$transcript, dims = 1:30, reduction = "pca", min.dist = input$minDist, negative.sample.rate = input$negativeSampleRate)
    seu$active <- seu$gene


    return(seu)

  }

  sidebar <- shinydashboard::dashboardSidebar(loadDataui("loadDataui", filterTypes = filterTypes),
                                              shinyWidgets::prettyRadioButtons("feature_type", "cluster on genes or transcripts?", choices = c("gene", "transcript"), selected = "gene"),
                                              shinyWidgets::prettyRadioButtons("organism_type", "What Organism was Sequenced?", choices = c("human", "mouse"), selected = "human"),
                                              shinyFiles::shinyFilesButton("seuratUpload", "Load a Custom Dataset", "Please select a .rds file", multiple = FALSE),
                                              shinyFiles::shinySaveButton("saveSeurat", "Save current Dataset", "Save file as...", filetype = list(rds = "rds")),
                                              verbatimTextOutput("savefile"),
                                              actionButton("changeEmbedAction", label = "Change Embedding Parameters"),
                                              changeEmbedParamsui("changeembed"),
                                              shinydashboard::sidebarMenu(shinydashboard::menuItem("Compare Plots", tabName = "comparePlots"),
                                                                          shinydashboard::menuItem("Compare Read Counts", tabName = "compareReadCount"),
                                                                          shinydashboard::menuItem("Differential Expression", tabName = "diffex"),
                                                                          shinydashboard::menuItem("Find Markers", tabName = "findMarkers"),
                                                                          shinydashboard::menuItem("Subset Seurat Input", tabName = "subsetSeurat"),
                                                                          shinydashboard::menuItem("All Transcripts", tabName = "allTranscripts"),
                                                                          shinydashboard::menuItem("RNA Velocity", tabName = "rnaVelocity")),

                                              width = 450)

  plotDimRedui <- function(id, plot_types) {
    ns <- NS(id)
    tagList(selectizeInput(ns("dplottype"), "Variable to Plot",
                           choices = plot_types, selected = c("custom"), multiple = TRUE),
            shinyWidgets::prettyRadioButtons(ns("embedding"),
                                             "dimensional reduction method", choices = c("pca", "harmony", "tsne", "umap"), selected = "umap", inline = TRUE),
            uiOutput(ns("featuretext")),
            sliderInput(ns("resolution"), "Resolution of clustering algorithm (affects number of clusters)", min = 0.2, max = 2, step = 0.2, value = 0.6),
            plotly::plotlyOutput(ns("dplot"), height = 750) %>%
              shinycssloaders::withSpinner())
  }
  plotDimRed <- function(input, output, session, seu, plot_types, feature_type, organism_type) {
    ns <- session$ns
    continuous_vars <- plot_types$continuous_vars
    category_vars <- plot_types$category_vars
    plot_types <- purrr::flatten_chr(plot_types)
    prefill_feature <- reactive({
      if (feature_type() == "transcript") {
        "ENST00000488147"
      }
      else if (feature_type() == "gene") {
        "RXRG"
      }
    })
    output$featuretext <- renderUI({
      textInput(ns("customFeature"), "gene or transcript on which to color the plot; eg. 'RXRG' or 'ENST00000488147'",
                value = prefill_feature())
    })
    output$dplot <- plotly::renderPlotly({
      req(seu$active)
      req(input$customFeature)
      if (length(input$dplottype) > 1) {
        mycols = input$dplottype
        louvain_resolution = paste0("clusters_", input$resolution)
        leiden_resolution = paste0("leiden_clusters_", input$resolution)
        mycols <- gsub("^seurat$", louvain_resolution,
                       mycols)
        newcolname = paste(mycols, collapse = "_")
        newdata = as_tibble(seu$active[[mycols]], rownames = "Sample_ID") %>%
          tidyr::unite(!!newcolname, mycols) %>% deframe() %>%
          identity()
        seu$active <- AddMetaData(seu$active, metadata = newdata,
                                  col.name = newcolname)
        plot_var(seu$active, embedding = input$embedding,
                 group = newcolname)
      }
      else {
        if (input$dplottype == "custom") {
          plot_feature(seu$active, embedding = input$embedding,
                       features = input$customFeature)
        }
        else if (input$dplottype %in% continuous_vars) {
          plot_feature(seu$active, embedding = input$embedding,
                       features = input$dplottype)
        }
        else if (input$dplottype == "seurat") {
          louvain_resolution = paste0("clusters_", input$resolution)
          plot_var(seu$active, embedding = input$embedding,
                   group = louvain_resolution)
        }
        else if (input$dplottype == "leiden") {
          leiden_resolution = paste0("leiden_clusters_", input$resolution)
          plot_var(seu$active, embedding = input$embedding,
                   group = leiden_resolution)
        }
        else if (input$dplottype %in% plot_types) {
          plot_var(seu$active, embedding = input$embedding,
                   group = input$dplottype)
        }
      }
    })
  }
  tableSelectedui <- function(id) {
    ns <- NS(id)
    tagList(DT::DTOutput(ns("brushtable")))
  }
  tableSelected <- function(input, output, session, seu) {
    ns <- session$ns
    brush <- reactive({
      req(seu$active)
      d <- plotly::event_data("plotly_selected")
      if (is.null(d)) {
        msg <- "Click and drag events (i.e. select/lasso) appear here (double-click to clear)"
        return(d)
      }
      else {
        selected_cells <- colnames(seu$active)[as.numeric(d$key)]
      }
    })
    output$brushtable <- DT::renderDT({
      req(seu$active)
      req(brush())
      selected_meta <- data.frame(seu$active[[]][brush(),
                                                 ])
      DT::datatable(selected_meta, extensions = "Buttons",
                    options = list(dom = "Bft", buttons = c("copy",
                                                            "csv"), scrollX = "100px", scrollY = "400px"))
    })
    selected_cells <- brush
    return(selected_cells)
  }
  displaySubsetui <- function(id) {
    ns <- NS(id)
    tagList(fluidRow(box(DT::DTOutput(ns("seu_meta")), width = 6),
                     box(DT::DTOutput(ns("sub_seu_meta")), width = 6)))
  }
  displaySubset <- function(input, output, session, seu, plot_types) {
    ns <- session$ns
    plot_types <- plot_types[-which(plot_types %in% c("seurat",
                                                      "custom"))]
    output$seu_meta <- DT::renderDT({
      req(seu$active)
      seu_meta <- data.frame(seu$active[[]]) %>% dplyr::select(Sample_ID,
                                                               batch, one_of(plot_types), everything())
      DT::datatable(seu_meta, extensions = "Buttons",
                    options = list(dom = "Bft", buttons = c("copy",
                                                            "csv"), scrollX = "100px", scrollY = "1000px"))
    }, server = TRUE)
    output$sub_seu_meta <- DT::renderDT({
      req(seu$active)
      req(input$seu_meta_rows_selected)
      sub_seu_meta <- data.frame(seu$active[[]][input$seu_meta_rows_selected,
                                                ]) %>% dplyr::select(Sample_ID, batch, one_of(plot_types),
                                                                     everything())
      brushtable <- DT::datatable(sub_seu_meta, extensions = "Buttons",
                                  options = list(dom = "Bft", buttons = c("copy",
                                                                          "csv"), scrollX = "100px", scrollY = "1000px"))
    }, server = TRUE)
    selected_rows <- reactive({
      input$seu_meta_rows_selected
    })
    return(selected_rows)
  }
  subsetSeuratui <- function(id) {
    ns <- NS(id)
    tagList()
  }

  subsetSeurat <- function(input, output, session, seu, selected_rows) {
    ns <- session$ns
    sub_seu <- reactive({
      showModal(modalDialog(title = "Subsetting and Recalculating Embeddings",
                            "This process may take a minute or two!"))
      seu$gene <- seu$gene[, selected_rows()]
      seu$gene <- seuratTools::seurat_pipeline(seu$gene,
                                               resolution = seq(0.6, 2, by = 0.2))
      seu$transcript <- seu$transcript[, selected_rows()]

      seu$transcript <- seuratTools::seurat_pipeline(seu$transcript,
                                                     resolution = seq(0.6, 2, by = 0.2))
      seu$active <- seu$gene
      removeModal()
    })
    return(sub_seu)
  }


  diffexui <- function(id) {
    ns <- NS(id)
    tagList(box(sliderInput(ns("resolution"), "Resolution of clustering algorithm (affects number of clusters)",
                            min = 0.2, max = 2, step = 0.2, value = 0.6
    ), shinyWidgets::prettyRadioButtons(ns("diffex_scheme"),
                                        "Cells to Compare",
                                        choiceNames = c(
                                          "Seurat Cluster",
                                          "Custom"
                                        ), choiceValues = c("seurat", "custom"),
                                        selected = "seurat"
    ), conditionalPanel(
      ns = ns,
      condition = "input.diffex_scheme == 'seurat'", numericInput(ns("cluster1"),
                                                                  "first cluster to compare",
                                                                  value = 0
      ), numericInput(ns("cluster2"),
                      "second cluster to compare",
                      value = 1
      )
    ), conditionalPanel(
      ns = ns,
      condition = "input.diffex_scheme == 'custom'", shinyWidgets::actionBttn(
        ns("saveClust1"),
        "Save to Custom Cluster 1"
      ), shinyWidgets::actionBttn(
        ns("saveClust2"),
        "Save to Custom Cluster 2"
      )
    ), shinyWidgets::prettyRadioButtons(ns("diffex_method"),
                                        "Method of Differential Expression",
                                        choiceNames = c(
                                          "t-test",
                                          "wilcoxon rank-sum test", "Likelihood-ratio test (bimodal)"
                                        ),
                                        choiceValues = c("t", "wilcox", "bimod")
    ), shinyWidgets::actionBttn(
      ns("diffex"),
      "Run Differential Expression"
    ), DT::dataTableOutput(ns("DT1")),
    downloadLink(ns("downloadData"), "Download Complete DE Results"),
    width = 12
    ), box(
      title = "Custom Cluster 1", DT::DTOutput(ns("cc1")),
      width = 12
    ), box(
      title = "Custom Cluster 2", DT::DTOutput(ns("cc2")),
      width = 12
    ))
  }
  diffex <- function(input, output, session, seu) {
    ns <- session$ns
    brush <- reactive({
      req(seu$active)
      d <- plotly::event_data("plotly_selected")
      if (is.null(d)) {
        msg <- "Click and drag events (i.e. select/lasso) appear here (double-click to clear)"
        return(d)
      }
      else {
        selected_cells <- colnames(seu$active)[as.numeric(d$key)]
      }
    })
    custom_cluster1 <- eventReactive(input$saveClust1,
                                     {
                                       isolate(brush())
                                     })
    custom_cluster2 <- eventReactive(input$saveClust2,
                                     {
                                       isolate(brush())
                                     })
    output$cc1 <- DT::renderDT({
      req(custom_cluster1())
      selected_meta <- data.frame(seu$active[[]][custom_cluster1(),
                                                 ])
      DT::datatable(selected_meta, extensions = "Buttons",
                    options = list(dom = "Bft", buttons = c("copy",
                                                            "csv"), scrollX = "100px", scrollY = "400px"))
    })
    output$cc2 <- DT::renderDT({
      req(custom_cluster2())
      selected_meta <- data.frame(seu$active[[]][custom_cluster2(),
                                                 ])
      DT::datatable(selected_meta, extensions = "Buttons",
                    options = list(dom = "Bft", buttons = c("copy",
                                                            "csv"), scrollX = "100px", scrollY = "400px"))
    })
    de_results <- eventReactive(input$diffex, {
      if (input$diffex_scheme == "seurat") {
        run_seurat_de(seu$active, input$cluster1, input$cluster2,
                      input$resolution, diffex_scheme = "seurat")
      }
      else if (input$diffex_scheme == "custom") {
        cluster1 <- unlist(strsplit(custom_cluster1(),
                                    " "))
        cluster2 <- unlist(strsplit(custom_cluster2(),
                                    " "))
        run_seurat_de(seu$active, cluster1, cluster2,
                      input$resolution, diffex_scheme = "custom")
      }
    })
    output$DT1 <- DT::renderDT(de_results()[[input$diffex_method]],
                               extensions = "Buttons", options = list(dom = "Bftpr",
                                                                      buttons = c("copy", "csv")), class = "display")
  }

  findMarkersui <- function(id) {
    ns <- NS(id)
    tagList(
      sliderInput(ns("resolution2"), label = "Resolution of clustering algorithm (affects number of clusters)", min = 0.2, max = 2, step = 0.2, value = 0.6),
      plotly::plotlyOutput(ns("markerplot"), height = 800)
    )
  }

  findMarkers <- function(input, output, session, seu) {
    ns <- session$ns

    resolution <- reactive({paste0("clusters_", input$resolution2)})

    output$markerplot <- plotly::renderPlotly({
      req(seu)
      #
      plot_markers(seu$active, resolution())
    })
  }

  plotReadCountui <- function(id, plot_types) {
    ns <- NS(id)
    tagList(selectInput(ns("rcplottype"), "Variable to Plot",
                        choices = plot_types, selected = c("custom"), multiple = TRUE),
            sliderInput(ns("resolution"), "Resolution of clustering algorithm (affects number of clusters)",
                        min = 0.2, max = 2, step = 0.2, value = 0.6),
            plotly::plotlyOutput(ns("rcplot"), height = 750) %>%
              shinycssloaders::withSpinner())
  }

  plotReadCount <- function(input, output, session, seu, plot_types) {
    ns <- session$ns
    output$rcplot <- plotly::renderPlotly({
      req(seu$active)
      if (input$rcplottype == "custom") {
        plot_readcount(seu$active, input$rcplottype)
      }
      else if (input$rcplottype == "seurat") {
        louvain_resolution = paste0("clusters_", input$resolution)
        plot_readcount(seu$active, louvain_resolution)
      }
      else if (input$rcplottype %in% plot_types) {
        plot_readcount(seu$active, input$rcplottype)
      }
    })
  }
  ccScoreui <- function(id) {
    ns <- NS(id)
    tagList()
  }
  ccScore <- function(input, output, session) {
    ns <- session$ns
    output$rplot1 <- renderPlot({
      req(seu$active)
      plot_ridge(seu$active, features = input$feature)
    })
    plotOutput("rplot1", height = 750)
  }

  allTranscriptsui <- function(id) {
    ns <- NS(id)
    tagList(fluidRow(box(shinyWidgets::prettyRadioButtons(ns("color_feature"), "cluster on genes or transcripts?", choices = c("gene", "transcript"), selected = "gene"),
                         textInput(ns("feature"), "gene or transcript on which to color the plot; eg. 'RXRG' or 'ENST00000488147'"),
                         uiOutput(ns("outfile")), uiOutput(ns("downloadPlot")),
                         width = 12)), fluidRow(uiOutput(ns("plotlys"))))
  }

  allTranscripts <- function(input, output, session, seu,
                             feature_type) {
    ns <- session$ns
    transcripts <- reactiveValues()
    transcripts <- reactive({
      req(feature_type())
      req(input$feature)
      req(seu)
      if (feature_type() == "gene") {
        transcripts <- dplyr::filter(annotables::grch38,
                                     symbol == input$feature) %>% dplyr::inner_join(annotables::grch38_tx2gene,
                                                                                    by = "ensgene") %>% dplyr::pull(enstxp)
        transcripts <- transcripts[transcripts %in%
                                     rownames(seu$transcript)]
      }
      else if (feature_type() == "transcript") {
        transcripts <- input$feature
        transcripts <- transcripts[transcripts %in%
                                     rownames(seu$transcript)]
      }
    })

    pList <- reactive({
      req(seu$active)

      if(input$color_feature == "gene"){
        # browser()
        transcript_cols <- as.data.frame(t(as.matrix(seu$transcript[["RNA"]][transcripts(),])))

        cells <- rownames(transcript_cols)
        transcript_cols <- as.list(transcript_cols) %>%
          purrr::map(~purrr::set_names(.x, cells))

        seu$gene[[transcripts()]] <- transcript_cols

        pList <- purrr::map(transcripts(), ~plot_feature(seu$gene,
                                                         embedding = input$embedding, features = .x))
        names(pList) <- transcripts()

      } else if (input$color_feature == "transcript"){
        pList <- purrr::map(transcripts(), ~plot_feature(seu$transcript,
                                                         embedding = input$embedding, features = .x))
        names(pList) <- transcripts()
      }
      return(pList)
    })

    output$plotlys <- renderUI({

      # result <- vector("list", n)

      plot_output_list <- purrr::map(names(pList()), ~plotly::plotlyOutput(ns(.x), height = 750))


      #   plot_output_list <- lapply(1:length(pList()), function(i) {
      #     plotname <- transcripts()[[i]]
      #     plotly::plotlyOutput(ns(plotname), height = 750)
      #   })
      do.call(tagList, plot_output_list)
    })

    observe({
      # browser()
      for (i in 1:length(pList())) {
        local({
          my_i <- i
          plotname <- transcripts()[[my_i]]
          output[[plotname]] <- plotly::renderPlotly({
            pList()[[my_i]]
          })
        })
      }

    })

    output$outfile <- renderUI({
      req(pList())
      textInput(ns("outfile"), "a descriptive name for the output file",
                value = paste0(input$feature, "_transcripts", "_clustered_by_", input$color_feature, ".pdf"))
    })

    output$plots <- downloadHandler(filename = function() {
      paste(input$outfile, ".pdf", sep = "")
    }, content = function(file) {
      pdf(file)
      lapply(pList(), print)
      dev.off()
    })

    output$downloadPlot <- renderUI({
      req(pList())
      downloadButton(ns("plots"), label = "Download plots")
    })
  }

  rnaVelocityui <- function(id){
    ns <- NS(id)
    tagList(
      shinyWidgets::prettyRadioButtons(ns("embedding"), "dimensional reduction method", choices = c("pca", "tsne", "umap"), selected = "umap", inline = TRUE),
      sliderInput(ns("resolution"), "Resolution of clustering algorithm (affects number of clusters)", min = 0.2, max = 2, step = 0.2, value = 0.6),
      plotOutput(ns("vel_plot"), height = "800px")
    )
  }

  rnaVelocity <- function(input, output, session, seu, feature_type, format = "grid"){
    ns <- session$ns
    output$vel_plot <- renderPlot({
      req(seu$active)

      showModal(modalDialog("Loading Plots", footer=NULL))
      vel <- seu$active@misc$vel
      emb <- Embeddings(seu$active, reduction = input$embedding)

      louvain_resolution = paste0("clusters_", input$resolution)

      cell.colors <- as_tibble(seu$active[[louvain_resolution]], rownames = "cellid") %>%
        tibble::deframe() %>%
        as.factor()

      levels(cell.colors) <- scales::hue_pal()(length(levels(cell.colors)))

      plot_velocity(vel, emb, cell.colors, format = format)
      removeModal()
    })
  }

  body <- shinydashboard::dashboardBody(shinydashboard::tabItems(
    shinydashboard::tabItem(
      tabName = "comparePlots",
      h2("Compare Plots"), fluidRow(box(plotDimRedui(
        "hello",
        purrr::flatten_chr(plot_types)
      )), box(plotDimRedui(
        "howdy",
        purrr::flatten_chr(plot_types)
      ))), fluidRow(box(
        title = "Selected Cells",
        tableSelectedui("hello"), width = 12
      ))
    ), shinydashboard::tabItem(
      tabName = "compareReadCount",
      h2("Compare Read Counts"), fluidRow(box(plotReadCountui(
        "hello",
        purrr::flatten_chr(plot_types)
      )), box(plotReadCountui(
        "howdy",
        purrr::flatten_chr(plot_types)
      )))
    ), shinydashboard::tabItem(
      tabName = "subsetSeurat",
      h2("Subset Seurat Input"), column(box(plotDimRedui(
        "subset",
        purrr::flatten_chr(plot_types)
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
        ), width = 6)
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
        "diffex",
        purrr::flatten_chr(plot_types)
      ), width = 12),
      box(
        title = "Selected Cells", tableSelectedui("diffex"),
        width = 12
      ),
      width = 6
      )
    ),
    shinydashboard::tabItem(
      tabName = "rnaVelocity",
      h2("RNA Velocity"),
      fluidRow(
        box(rnaVelocityui("arrow"),
            width = 12)),
      fluidRow(
        box(rnaVelocityui("grid"),
            width = 12))
    )
  ))

  ui <- dashboardPage(
    header = header, sidebar = sidebar,
    body = body
  )

  server <- function(input, output, session) {
    options(warn = -1)
    # sidebar
    seu <- callModule(loadData, "loadDataui", proj_dir,
                      input$feature_type)
    observeEvent(input$feature_type, {
      seu$active <- seu[[input$feature_type]]
    })
    feature_type <- reactive({
      input$feature_type
    })

    # list volumes
    volumes <- c(Home = fs::path_home(), "R Installation" = R.home(), shinyFiles::getVolumes()())

    # upload seurat object

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
          Sys.sleep(12)
          shiny::incProgress(4/10)
          Sys.sleep(18)
          shiny::incProgress(6/10)
          Sys.sleep(24)
          shiny::incProgress(8/10)
          saveRDS(list(gene = seu$gene, transcript = seu$transcript), subSeuratPath())

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
    callModule(plotReadCount, "hello", seu, purrr::flatten_chr(plot_types))
    callModule(plotReadCount, "howdy", seu, purrr::flatten_chr(plot_types))
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

          message("reintegrating transcript expression")
          seu$transcript <- seurat_pipeline(seu$transcript, reduction = "harmony", resolution = seq(0.2, 2.0, by = 0.2))

          message("reintegrating gene expression")
          seu$gene <- seurat_pipeline(seu$gene, reduction = "harmony", resolution = seq(0.2, 2.0, by = 0.2))

          # subset_seus <- seuratTools::reintegrate_seus(list(gene = seu$gene,
          #                                                   transcript = seu$transcript), reduction = "harmony", temp = TRUE)

        } else {

          seu$transcript <- seuratTools::seurat_pipeline(seu$transcript, resolution = seq(0.2, 2.0, by = 0.2))
          seu$gene <- seuratTools::seurat_pipeline(seu$gene, resolution = seq(0.2, 2.0, by = 0.2))

        }
        seu$active <- seu$gene

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
    callModule(diffex, "hello", seu)
    observeEvent(input$plotTrx, {
      showModal(modalDialog(
        title = "Plotting Transcripts",
        "This process may take a minute or two!"
      ))
      callModule(allTranscripts, "hello", seu, feature_type)
      callModule(allTranscripts, "howdy", seu, feature_type)
      removeModal()
    })

    callModule(rnaVelocity, "arrow", seu, feature_type, "arrow")
    callModule(rnaVelocity, "grid", seu, feature_type, "grid")

  }
  shinyApp(ui, server)
}
