#' Load Data UI
#'
#' @param id
#' @param label
#' @param filterTypes
#'
#' @return
#' @export
#'
#' @examples
loadDataui <- function(id, label = "Load Data", filterTypes) {
  ns <- NS(id)
  tagList(
    shinyWidgets::prettyRadioButtons(ns("filterType"), "dataset to include", choices = filterTypes, selected = ""),
    shinyWidgets::actionBttn(ns("loadButton"), "Load Default Dataset")
    # fileInput(ns("seuratUpload"), "Upload .rds file")
  )
}

#' Load Data
#'
#' @param input
#' @param output
#' @param session
#' @param feature_types
#' @param selected_feature
#' @param proj_dir
#'
#' @return
#' @export
#'
#' @examples
loadData <- function(input, output, session, proj_dir, feature_types = c("gene"), selected_feature) {
  ns <- session$ns
  seu <- reactiveValues()

  observeEvent(input$loadButton, {
    showModal(modalDialog("Loading Data", footer = NULL))
    # if (!input$filterType == "") {
    #   filterType = paste0("_", input$filterType)
    # }
    # else {
    #   filterType = input$filterType
    # }

    # seu_paths <- paste0("*", feature_types, "_seu", filterType, ".rds")
    #
    # seu_paths <- fs::path(proj_dir, "output", "seurat") %>%
    #   fs::dir_ls() %>%
    #   fs::path_filter(seu_paths) %>%
    #   purrr::set_names(feature_types) %>%
    #   identity()

    seu_paths <- load_seurat_path(proj_dir, features = feature_types, suffix = input$filterType)

    feature_seus <- purrr::map(seu_paths, readRDS) %>%
      purrr::set_names(feature_types)

    feature_seus <- purrr::map(feature_seus, SetDefaultAssay, "RNA")
    removeModal()

    for (i in names(feature_seus)){
      seu[[i]] <- feature_seus[[i]]
    }
    seu$active <- feature_seus[[selected_feature]]

    # seu$transcript <- feature_seus[["transcript"]]
    # seu$gene <- feature_seus[["gene"]]
    # seu$active <- feature_seus[[selected_feature]]


  })

  return(seu)

}

#' Plot Custom Feature UI
#'
#' @param id
#'
#' @return
#' @export
#'
#' @examples
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

#' Plot Custom Feature
#'
#' @param input
#' @param output
#' @param session
#'
#' @return
#' @export
#'
#' @examples
customFeature <- function(input, output, session) {
  ns <- session$ns
  output$featuretext <- renderUI({
    textInput("feature", "gene or transcript on which to color the plot; eg. 'RXRG' or 'ENST00000488147'",
              value = choice())
  })
  uiOutput("featuretext")
}

#' Change Embedding Parameters UI
#'
#' @param id
#'
#' @return
#' @export
#'
#' @examples
changeEmbedParamsui <- function(id){
  ns <- NS(id)

  minDist_vals <- prep_slider_values(0.3)
  negsamprate_vals <- prep_slider_values(5)

  tagList(
    selectizeInput(ns("dims"), label = "Dimensions from PCA", choices = seq(1,99), multiple = TRUE, selected = 1:30),
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
changeEmbedParams <- function(input, output, session, seu){
  ns <- session$ns
  #
  # output$embedControls <- renderUI({
  #   tagList(
  #     sliderInput(ns("minDist"), label = "Minimum Distance", min = minDist_vals$min, max = minDist_vals$max, value = minDist_vals$value, step = minDist_vals$step),
  #     sliderInput(ns("negativeSampleRate"), label = "NegativeSampleRate", min = minDist_vals$min, max = minDist_vals$max, value = minDist_vals$value, step = minDist_vals$step)
  #   )
  # })

  seu$gene <- RunUMAP(seu$gene, dims = as.numeric(input$dims), reduction = "pca", min.dist = input$minDist, negative.sample.rate = input$negativeSampleRate)
  seu$transcript <- RunUMAP(seu$transcript, dims = as.numeric(input$dims), reduction = "pca", min.dist = input$minDist, negative.sample.rate = input$negativeSampleRate)
  seu$active <- seu$gene


  return(seu)

}

#' Plot Dimensionally Reduced Data UI
#'
#' @param id
#'
#' @return
#' @export
#'
#' @examples
plotDimRedui <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("dplottype")),
    shinyWidgets::prettyRadioButtons(ns("embedding"),
      "dimensional reduction method",
      choices = c("pca", "harmony", "tsne", "umap"), selected = "umap", inline = TRUE
    ),
    fluidRow(
      column(2, selectizeInput(ns("dim1"), "Dimension 1", choices = seq(1,99), selected = 1)),
      column(2, selectizeInput(ns("dim2"), "Dimension 2", choices = seq(1,99), selected = 2))
    ),
    uiOutput(ns("featuretext")),
    sliderInput(ns("resolution"), "Resolution of clustering algorithm (affects number of clusters)", min = 0.2, max = 2, step = 0.2, value = 0.6),
    plotly::plotlyOutput(ns("dplot"), height = 750)
  )
}

#' Plot Dimensionally Reduced Data
#'
#' @param input
#' @param output
#' @param session
#' @param seu
#' @param plot_types
#' @param feature_type
#' @param organism_type
#'
#' @return
#' @export
#'
#' @examples
plotDimRed <- function(input, output, session, seu, plot_types, feature_type, organism_type) {
  ns <- session$ns

  output$dplottype <- renderUI({
    req(seu$active)
    selectizeInput(ns("plottype"), "Variable to Plot",
                   choices = purrr::flatten_chr(plot_types()), selected = c("custom"), multiple = TRUE)
  })

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
    req(input$plottype)
    req(seu$active)

    req(input$customFeature)
    if (length(input$dplottype) > 1) {
      mycols = input$dplottype

      louvain_resolution = paste0(DefaultAssay(seu$active), "_snn_res.", input$resolution)
      leiden_resolution = paste0("leiden_clusters_", input$resolution)
      mycols <- gsub("^seurat$", louvain_resolution,
                     mycols)
      newcolname = paste(mycols, collapse = "_")
      newdata = as_tibble(seu$active[[mycols]], rownames = "Sample_ID") %>%
        tidyr::unite(!!newcolname, mycols) %>% deframe() %>%
        identity()
      seu$active <- AddMetaData(seu$active, metadata = newdata,
                                col.name = newcolname)
      plot_var(seu$active, dims = c(input$dim1, input$dim2), embedding = input$embedding,
               group = newcolname)
    }
    else {
      if (input$plottype == "custom") {
        plot_feature(seu$active, dims = c(input$dim1, input$dim2), embedding = input$embedding,
                     features = input$customFeature)
      }
      else if (input$plottype %in% plot_types()$continuous_vars) {
        plot_feature(seu$active, dims = c(input$dim1, input$dim2), embedding = input$embedding,
                     features = input$plottype)
      }
      else if (input$plottype == "seurat") {

        if ("integrated" %in% names(seu$active@assays)){
          active_assay <- "integrated"
        } else {
          active_assay <- "RNA"
        }

        louvain_resolution = paste0(active_assay, "_snn_res.", input$resolution)
        plot_var(seu$active, dims = c(input$dim1, input$dim2), embedding = input$embedding,
                 group = louvain_resolution)
      }
      else if (input$plottype == "leiden") {
        leiden_resolution = paste0("leiden_clusters_", input$resolution)
        plot_var(seu$active, dims = c(input$dim1, input$dim2), embedding = input$embedding,
                 group = leiden_resolution)
      }
      else if (input$plottype %in% plot_types()$category_vars) {
        plot_var(seu$active, dims = c(input$dim1, input$dim2), embedding = input$embedding,
                 group = input$plottype)
      }
    }

    # if ("integrated" %in% names(seu$active@assays)){
    #   DefaultAssay(seu$active) <- "RNA"
    # }

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

#' Subset Seruat UI
#'
#' @param id
#'
#'
#' @return
#' @export
#'
#' @examples
subsetSeuratui <- function(id) {
  ns <- NS(id)
  tagList()
}

#' Subset Seruat
#'
#' @param input
#' @param output
#' @param session
#' @param seu
#' @param selected_rows
#'
#' @return
#' @export
#'
#' @examples
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
  ),
  shinyWidgets::actionBttn(ns("diffex"),
                              "Run Differential Expression"),
  downloadLink(ns("downloadData"), "Download Complete DE Results"),
  DT::dataTableOutput(ns("DT1")),
  width = 12
  ), box(
    title = "Custom Cluster 1", DT::DTOutput(ns("cc1")),
    width = 12
  ), box(
    title = "Custom Cluster 2", DT::DTOutput(ns("cc2")),
    width = 12
  ))
}



#' Differential Expression
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

  cluster_list <- reactive({
    if (input$diffex_scheme == "seurat"){
      seu_meta <- seu$active[[paste0(DefaultAssay(seu$active), "_snn_res.", input$resolution)]]
      cluster1_cells <- rownames(seu_meta[seu_meta == input$cluster1, , drop = FALSE])
      cluster2_cells <- rownames(seu_meta[seu_meta == input$cluster2, , drop = FALSE])
      list(cluster1 = cluster1_cells, cluster2 = cluster2_cells)
    } else if (input$diffex_scheme == "custom"){
      list(cluster1 = custom_cluster1(), cluster2 = custom_cluster2())
    }

  })

  return(cluster_list)

}

#' Gene Enrichment UI
#'
#' @param id
#'
#' @return
#' @export
#'
#' @examples
geneEnrichmentui <- function(id){
  ns <- NS(id)
  tagList(
    actionButton(ns("enrichmentAction"), "Run Enrichment Analysis"),
    textOutput(ns("enrichmentMessages")),
    uiOutput(ns("reportLink"))
    # tags$a("Results of Functional Enrichment Analysis", target = "_blank", href = "enrichmentbrowser/mainpage.html")
    # tags$iframe(style = "height:1400px; width:100%", src = "enrichmentbrowser/mainpage.html")
  )

}


#' Gene Enrichment
#'
#' @param input
#' @param output
#' @param session
#' @param seu
#' @param diffex_clusters
#'
#' @return
#' @export
#'
#' @examples
geneEnrichment <- function(input, output, session, seu, diffex_results){
    ns <- session$ns
    enrichmentReport <- eventReactive(input$enrichmentAction, {
      withCallingHandlers({
        shinyjs::html("enrichmentMessages", "")
        message("Beginning")

        # showModal(modalDialog("Calculating Functional Enrichment", footer=NULL))
        enrichmentReport <- run_enrichmentbrowser(seu = seu$active,
                              cluster1_cells = diffex_results()$cluster1,
                              cluster2_cells = diffex_results()$cluster2)

        # enrichmentReport <- "enrichmentbrowser2/mainpage.html"
        # removeModal()

        return(enrichmentReport)

    },
      message = function(m) {
        shinyjs::html(id = "enrichmentMessages", html = paste0("Running Functional Enrichment Analysis: ", m$message), add = FALSE)
      })
    })

    output$reportLink <- renderUI({
      tags$a("Results of Functional Enrichment Analysis", target = "_blank", href = enrichmentReport())
    })

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
    sliderInput(ns("resolution2"), label = "Resolution of clustering algorithm (affects number of clusters)", min = 0.2, max = 2, step = 0.2, value = 0.6),
    plotly::plotlyOutput(ns("markerplot"), height = 800)
  )
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
findMarkers <- function(input, output, session, seu) {
  ns <- session$ns

  resolution <- reactive({
    if ("integrated" %in% names(seu$active@assays)){
      active_assay <- "integrated"
    } else {
      active_assay <- "RNA"
    }
    paste0(active_assay, "_snn_res.", input$resolution2)

    })

  output$markerplot <- plotly::renderPlotly({
    req(seu)
    #
    plot_markers(seu$active, resolution())
  })
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
  tagList(
    uiOutput(ns("rcplottype")),
    sliderInput(ns("resolution"), "Resolution of clustering algorithm (affects number of clusters)",
      min = 0.2, max = 2, step = 0.2, value = 0.6
    ),
    plotly::plotlyOutput(ns("rcplot"), height = 750)
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

  output$rcplottype <- renderUI({
    req(seu$active)
    shiny::selectInput(ns("plottype"), "Variable to Plot",
                       choices = purrr::flatten_chr(plot_types()), selected = c("seurat"), multiple = TRUE
    )
  })

  output$rcplot <- plotly::renderPlotly({
    req(seu$active)
    req(input$plottype)

    if (input$plottype == "seurat") {

      if ("integrated" %in% names(seu$active@assays)){
        active_assay <- "integrated"
      } else {
        active_assay <- "RNA"
      }

      louvain_resolution = paste0(active_assay, "_snn_res.", input$resolution)
      plot_readcount(seu$active, louvain_resolution)
    }
    else if (input$plottype %in% purrr::flatten_chr(plot_types())) {
      plot_readcount(seu$active, input$plottype)
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
    req(seu$active)
    plot_ridge(seu$active, features = input$feature)
  })
  plotOutput("rplot1", height = 750)
}

#' Plot All Transcripts UI
#'
#' @param id
#'
#' @return
#' @export
#'
#' @examples
allTranscriptsui <- function(id) {
  ns <- NS(id)
  tagList(fluidRow(box(shinyWidgets::prettyRadioButtons(ns("color_feature"), "cluster on genes or transcripts?", choices = c("gene", "transcript"), selected = "gene"),
                       textInput(ns("feature"), "gene or transcript on which to color the plot; eg. 'RXRG' or 'ENST00000488147'"),
                       uiOutput(ns("outfile")), uiOutput(ns("downloadPlot")),
                       width = 12)), fluidRow(uiOutput(ns("plotlys"))))
}

#' Plot All Transcripts
#'
#' @param input
#' @param output
#' @param session
#' @param seu
#' @param feature_type
#'
#' @return
#' @export
#'
#' @examples
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

#' RNA Velocity UI
#'
#' @param id
#'
#' @return
#' @export
#'
#' @examples
rnaVelocityui <- function(id){
  ns <- NS(id)
  tagList(
    shinyWidgets::prettyRadioButtons(ns("embedding"), "dimensional reduction method", choices = c("pca", "tsne", "umap"), selected = "umap", inline = TRUE),
    sliderInput(ns("resolution"), "Resolution of clustering algorithm (affects number of clusters)", min = 0.2, max = 2, step = 0.2, value = 0.6),
    plotOutput(ns("vel_plot"), height = "800px")
  )
}

#' RNA Velocity
#'
#' @param input
#' @param output
#' @param session
#' @param seu
#' @param feature_type
#' @param format
#'
#' @return
#' @export
#'
#' @examples
rnaVelocity <- function(input, output, session, seu, feature_type, format = "grid"){
  ns <- session$ns
  output$vel_plot <- renderPlot({
    req(seu$active)

    showModal(modalDialog("Loading Plots", footer=NULL))
    vel <- seu$active@misc$vel
    emb <- Embeddings(seu$active, reduction = input$embedding)

    louvain_resolution = paste0(DefaultAssay(seu$active), "_snn_res.", input$resolution)

    cell.colors <- as_tibble(seu$active[[louvain_resolution]], rownames = "cellid") %>%
      tibble::deframe() %>%
      as.factor()

    levels(cell.colors) <- scales::hue_pal()(length(levels(cell.colors)))

    plot_velocity(vel, emb, cell.colors, format = format)
    removeModal()
  })
}


monocleui <- function(id){
    ns <- NS(id)
    tagList(
      sliderInput(ns("resolution"), "Resolution of clustering algorithm (affects number of clusters)", min = 0.2, max = 2, step = 0.2, value = 0.6),
      # actionButton(ns("plotMonocle"), "plot trajectory"),
      plotlyOutput(ns("monoclePlot"))
      # box(textOutput(ns("monocleText")))

        )
}

monocle <- function(input, output, session, seu, input_type){
    ns <- session$ns

    cds <- reactive({
      req(seu$active)
      convert_seu_to_cds(seu$active)
    })

    output$monoclePlot <- renderPlotly({
      req(seu$active)

      plot_cds(cds(), input$resolution)
    })

    root_cells <- reactive({
      d <- plotly::event_data("plotly_selected")
      if (is.null(d)) {
        msg <- "Click and drag events (i.e. select/lasso) appear here (double-click to clear)"
        return(d)
      }
      else {
        selected_cells <- colnames(cds())[as.numeric(d$key)]
      }
    })

    output$monocleText <- renderText({
      req(root_cells())
      root_cells()
    })

}
