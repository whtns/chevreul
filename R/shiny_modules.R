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
    seuratToolsBox(
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

  output$checkSeu <- renderText({
    req(seu$active)
    "test"
  })

  output$clustree <- renderPlot({
    req(seu$active)
    clustree::clustree(seu$active)
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
plotViolinui <- function(id){
  ns <- NS(id)
  tagList(
    seuratToolsBox(
      title = "Violin Plots",
      uiOutput(ns("vln_group")),
      selectizeInput(ns("customFeature"),
                     "Gene or transcript expression by which to color the plot eg. 'RXRG' or 'ENST00000488147'",
                     choices = NULL, multiple = TRUE),
      plotly::plotlyOutput(ns("vplot"), height = 750),
      width = 12)
    )
}

#' Plot Violin Server
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
plotViolin <- function(input, output, session, seu, featureType, organism_type){
  ns <- session$ns
  prefill_feature <- reactive({
    req(featureType())
    if (featureType() == "transcript") {
      if (organism_type() == "human") {
        "ENST00000488147"
      }
      else if (organism_type() == "mouse") {
        "ENSG00000488147"
      }
    }
    else if (featureType() == "gene") {
      if (organism_type() == "human") {
        "RXRG"
      }
      else if (organism_type() == "mouse") {
        "Rxrg"
      }
    }
  })
  observe({
    req(prefill_feature())
    req(seu$active)
    updateSelectizeInput(session, "customFeature", choices = rownames(seu$active@assays$RNA),
                         selected = prefill_feature(), server = TRUE)
  })

  output$vln_group <- renderUI({
    req(seu$active)
    selectizeInput(ns("vlnGroup"), "Grouping variable",
                   choices = colnames(seu$active[[]]), selected = "batch")
  })
  output$vplot <- plotly::renderPlotly({
    req(input$customFeature)
    req(input$vlnGroup)

    plot_violin(seu$active, plot_var = input$vlnGroup, features = input$customFeature)
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
plotHeatmapui <- function(id){
  ns <- NS(id)
  tagList(
    seuratToolsBox(
      title = "Heatmap",
      uiOutput(ns("colAnnoVarui")),
      radioButtons(ns("slot"), "Data Scaling", choices = c(scaled = "scale.data", unscaled = "data"), selected = "scale.data", inline = TRUE),
      selectizeInput(ns("dendroSelect"), "Clustering algorigthm for column dendrogram or metadata for annotation", choices = NULL, selected = NULL),
      actionButton(ns("actionHeatmap"), "Plot Heatmap"),
      downloadButton(ns("downloadPlot"), "Download Heatmap"),
      selectizeInput(ns("customFeature"), "Gene or transcript expression by which to color the plot; eg. 'RXRG' or 'ENST00000488147'",
                     choices = NULL, multiple = TRUE),
      plotOutput(ns("heatmap"), height = 750),
      width = 12) %>%
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
plotHeatmap <- function(input, output, session, seu, featureType, organism_type){
  ns <- session$ns

  observe({
    req(seu$active)

    if ("integrated" %in% names(seu$active@assays)) {
      default_assay = "integrated"
    } else {
      default_assay = "RNA"
    }

    preset_features <- VariableFeatures(seu$active, assay = default_assay)[1:50]

    updateSelectizeInput(session, "customFeature", choices = rownames(seu$active@assays$RNA),
                         selected = preset_features, server = TRUE)
  })

  output$colAnnoVarui <- renderUI({
    req(seu$active)

    selectizeInput(ns("colAnnoVar"), "Column Annotation(s)",
                   choices = colnames(seu$active[[]]), selected = "batch", multiple = TRUE)
  })

  observe({
    req(seu$active)

    hclust_methods <- c("Ward" = "ward.D2", "single", "complete", "average")

    updateSelectizeInput(session, "dendroSelect", choices = c(hclust_methods, colnames(seu$active[[]])), selected = "ward.D2")
  })

  heatmap_plot <- eventReactive(input$actionHeatmap, {
    req(input$customFeature)
    req(input$colAnnoVar)

    if ("integrated" %in% names(seu$active@assays)) {
      default_assay = "integrated"
    } else {
      default_assay = "RNA"
    }

    hm <- seu_complex_heatmap(seu$active, features = input$customFeature, assay = default_assay, group.by = input$colAnnoVar, slot = input$slot, col_dendrogram = input$dendroSelect)

    return(hm)
  })

  output$heatmap <- renderPlot({
    heatmap_plot()
  })

  output$downloadPlot <- downloadHandler(
    filename = function() { paste("heatmap", '.pdf', sep='') },
    content = function(file) {
      ggsave(file, ggplotify::as.ggplot(heatmap_plot()), width = 16, height = 12)
    })



}





#' Reformat Seurat Object Metadata UI
#'
#' @param id
#'
#' @return
#' @export
#'
#' @examples
reformatMetadataui <- function(id) {
  ns <- NS(id)
  tagList(
    seuratToolsBox(
      title = "Reformat Metadata",
    checkboxInput(ns("header"), "Header", TRUE),
    fileInput(ns("metaFile"), "Choose CSV File of metadata with cell names in first column",
              accept = c(
                "text/csv",
                "text/comma-separated-values,text/plain",
                ".csv")
    ),
    actionButton(ns("updateMetadata"), "Update Metadata"),
    radioButtons(ns("updateMethod"), "Update By:", choices = c("table (below)" = "spreadsheet", "uploaded file" = "file"), inline = TRUE),
    # downloadLink(ns("downloadMetadata"), "Download Metadata"),
    rhandsontable::rHandsontableOutput(ns("seuTable")),
    width = 12
    ) %>%
      default_helper(type = "markdown", content = "reformatMetadata")

  )
}

#' Reformat Seurat Object Metadata Server
#'
#' @param input
#' @param output
#' @param sessionk
#' @param seu
#'
#' @return
#' @export
#'
#' @examples
reformatMetadata <- function(input, output, session, seu) {
  ns <- session$ns

  meta <- reactiveValues()

  observe({
    meta$old <- seu$active@meta.data
  })

  observeEvent(input$updateMetadata, {

    if (input$updateMethod == "file"){
      inFile <- input$metaFile

      if (is.null(inFile))
        return(NULL)

      meta$new <- format_new_metadata(inFile$datapath)

    } else if (input$updateMethod == "spreadsheet"){
      meta$new <- propagate_spreadsheet_changes(input$seuTable)
    }

    for (i in names(seu)){
      seu[[i]] <- Seurat::AddMetaData(seu[[i]], meta$new)
    }

  })


  table_out <- reactive({
    meta$new %||% meta$old
  })

  output$seuTable <- rhandsontable::renderRHandsontable({

    rhandsontable::rhandsontable(table_out(), rowHeaderWidth = 200, height = 700) %>%
      rhandsontable::hot_context_menu(
        customOpts = list(
          csv = list(name = "Download to CSV",
                     callback = htmlwidgets::JS(
                       "function (key, options) {
                         var csv = csvString(this, sep=',', dec='.');

                         var link = document.createElement('a');
                         link.setAttribute('href', 'data:text/plain;charset=utf-8,' +
                           encodeURIComponent(csv));
                         link.setAttribute('download', 'data.csv');

                         document.body.appendChild(link);
                         link.click();
                         document.body.removeChild(link);
                       }"))))

  })

  # output$downloadMetadata <- downloadHandler(
  #   filename = function() {
  #     paste("data-", Sys.Date(), ".csv", sep="")
  #   },
  #   content = function(file) {
  #     write.csv(table_out(), file)
  #   })

  return(seu)

}


#' Integrate Project UI
#'
#' @param id
#'
#' @return
#' @export
#'
#' @examples
integrateProjui <- function(id){
    ns <- NS(id)
    tagList(
      seuratToolsBox(
        title = "Integrate Projects",
        actionButton(ns("integrateAction"), "Integrate Selected Projects"),
        textOutput(ns("integrationComplete")),
        shinyjs::useShinyjs(),
        textOutput(ns("integrationMessages")),
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
integrateProj <- function(input, output, session, proj_matrices, seu, proj_dir, con){
    ns <- session$ns

    proj_matrix <- reactive({
      proj_matrices()$primary_projects
    })

    clean_proj_matrix <- reactive({


      clean_proj_matrix <- proj_matrix() %>%
        dplyr::select(-project_path) %>%
        identity()
    })



    output$myDatatable <- DT::renderDataTable(clean_proj_matrix(),
                                              server = FALSE,
                                              rownames=TRUE)

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
          withCallingHandlers({
            shinyjs::html("integrationMessages", "")
            message("Beginning")

            # check if seurat paths exist

            # validate(
            #   need(input$data != "", "Please select a data set")
            # )

            batches <- fs::path(selectedProjects(), "output", "seurat", "unfiltered_seu.rds") %>%
              purrr::map(readRDS)

            names(batches) <- names(selectedProjects())

            mergedSeus(integration_workflow(batches))

            message("Integration Complete!")

          },
          message = function(m) {
            shinyjs::html(id = "integrationMessages", html = paste0("Running Integration: ", m$message), add = FALSE)
          })
      })

    newProjDir <- reactive({
      req(mergedSeus())

      print(names(mergedSeus()))

      for (i in names(mergedSeus())){
        seu[[i]] <- mergedSeus()[[i]]
      }

      newProjName <- paste0(purrr::map(fs::path_file(selectedProjects()), ~gsub("_proj", "", .x)), collapse = "_")
      integrated_proj_dir <- "/dataVolume/storage/single_cell_projects/integrated_projects/"
      newProjDir <- fs::path(integrated_proj_dir, newProjName)

      proj_dir(newProjDir)

      newProjDir

    })

    output$integrationComplete <- renderText({
      req(mergedSeus())
      # print("integration complete!")
      print("")
    })


    volumes <- reactive({
      volumes <- c(Home = fs::path("/dataVolume/storage/single_cell_projects/integrated_projects"), "R Installation" = R.home(), shinyFiles::getVolumes())
      # print(volumes)
      volumes
    })

    observe({
      shinyFiles::shinyFileSave(input, "saveIntegratedProject", roots = volumes(), session = session, restrictions = system.file(package = "base"))
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

      if (!is.null(integratedProjectSavePath())){

        shiny::withProgress(
          message = paste0("Saving Integrated Dataset to ", integratedProjectSavePath()),
          value = 0,
          {
            # Sys.sleep(6)
            shiny::incProgress(2/10)
            # myseuratdir <- fs::path(paste0(integratedProjectSavePath(), "_proj"), "output", "seurat")
            # dir.create(myseuratdir)
            # myseuratpath <- fs::path(myseuratdir, "unfiltered_seu.rds")
            # saveRDS(mergedSeus(), myseuratpath)
            # Sys.chmod(myseuratpath)
            save_seurat(mergedSeus(), proj_dir = integratedProjectSavePath())
            set_permissions_call <- paste0("chmod -R 775 ", integratedProjectSavePath())
            system(set_permissions_call)
            writeLines(character(), fs::path(integratedProjectSavePath(), ".here"))
            # create_proj_db()
            DBI::dbAppendTable(con, "projects", data.frame(project_name = fs::path_file(integratedProjectSavePath()), project_path = integratedProjectSavePath()))
            # system("updatedb -l 0 -U /dataVolume/storage/single_cell_projects/ -o /dataVolume/storage/single_cell_projects/single_cell_projects.db", wait = TRUE)
            # # print(set_permissions_call)
            shiny::incProgress(8/10)

            velocyto_dir <- fs::path(integratedProjectSavePath(), "output", "velocyto")
            fs::dir_create(velocyto_dir)
            new_loom_path <- fs::path(velocyto_dir, fs::path_file(integratedProjectSavePath()))
            combine_looms(selectedProjects(), new_loom_path)

          })

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

#' Plot Dimensional Reduduction UI
#'
#' @param id
#'
#' @return
#' @export
#'
#' @examples
plotDimRedui <- function(id){
  ns <- NS(id)
  seuratToolsBox(title = "Embedding",
                 uiOutput(ns("dplottype")),
          uiOutput(ns("embeddings")),
          sliderInput(ns("dotSize"), "Size of Points in UMAP", min = 1, max = 3, step = 0.5, value = 1),
          fluidRow(
            column(2,
                   selectizeInput(ns("dim1"), "Dimension 1", choices = seq(1, 99), selected = 1)
                   ),
            column(2,
                   selectizeInput(ns("dim2"), "Dimension 2", choices = seq(1, 99), selected = 2)
                   )
            ),
          selectizeInput(ns("customFeature"), "Gene or transcript expression by which to color the plot; eg. 'RXRG' or 'ENST00000488147'", choices = NULL, multiple = TRUE),
          sliderInput(ns("resolution"), "Resolution of clustering algorithm (affects number of clusters)", min = 0.2, max = 2, step = 0.2, value = 0.6),
          plotly::plotlyOutput(ns("dplot"), height = 500),
          width = 12)
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
                       organism_type, reductions){
  ns <- session$ns

  output$embeddings <- renderUI({
    req(seu$active)
    selectizeInput(ns("embedding"), "dimensional reduction method", choices = reductions(), selected = rev(reductions())[1])
  })

  selected_plot <- reactiveVal()
  output$dplottype <- renderUI({
    req(seu$active)
    # selected_plot <- ifelse(is.null(selected_plot()), "seurat",
    #                         selected_plot())
    selectizeInput(ns("plottype"), "Variable to Plot", choices = purrr::flatten_chr(plot_types()),
                   selected = "seurat", multiple = TRUE)
  })
  prefill_feature <- reactive({
    req(featureType())
    if (featureType() == "transcript") {
      if (organism_type() == "human") {
        "ENST00000488147"
      }
      else if (organism_type() == "mouse") {

        "ENSG00000488147"
      }
    }
    else if (featureType() == "gene") {
      if (organism_type() == "human") {
        "RXRG"
      }
      else if (organism_type() == "mouse") {
        "Rxrg"
      }
    }
  })
  observe({
    req(prefill_feature())
    req(seu$active)
    updateSelectizeInput(session, "customFeature", choices = rownames(seu$active@assays$RNA),
                         selected = prefill_feature(), server = TRUE)
  })

  output$dplot <- plotly::renderPlotly({
    req(input$plottype)
    req(seu$active)
    req(input$embedding)
    if (length(input$plottype) > 1) {

      seu$active <- cross_plot_vars(seu$active, input$resolution, input$plottype)

      newcolname <- paste(input$plottype, collapse = "_")
      seu$active[[newcolname]] <- Idents(seu$active)

      selected_plot(newcolname)

      plot_var(seu$active, dims = c(input$dim1, input$dim2),
               embedding = input$embedding, group = NULL, pt.size = input$dotSize)
    }
    else {
      if (input$plottype == "custom") {
        plot_feature(seu$active, dims = c(input$dim1,
                                          input$dim2), embedding = input$embedding,
                     features = input$customFeature, pt.size = input$dotSize)
      }
      else if (input$plottype %in% plot_types()$continuous_vars) {
        plot_feature(seu$active, dims = c(input$dim1,
                                          input$dim2), embedding = input$embedding,
                     features = input$plottype, pt.size = input$dotSize)
      }
      else if (input$plottype == "seurat") {
        if ("integrated" %in% names(seu$active@assays)) {
          active_assay <- "integrated"
        }
        else {
          active_assay <- "RNA"
        }

        louvain_resolution = reactive({
          paste0(active_assay, "_snn_res.", input$resolution)
        })

        plot_var(seu$active, dims = c(input$dim1, input$dim2),
                 embedding = input$embedding, group = louvain_resolution(), pt.size = input$dotSize)
      }
      else if (input$plottype %in% plot_types()$category_vars) {
        plot_var(seu$active, dims = c(input$dim1, input$dim2),
                 embedding = input$embedding, group = input$plottype, pt.size = input$dotSize)
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
    selected_meta <- data.frame(seu$active[[]][brush(),])

    # selection = list(mode = 'multiple', selected = c(1, 3, 8), target = 'row'),
    DT::datatable(selected_meta, extensions = "Buttons",
                  selection = list(mode = 'multiple', selected = 1:nrow(selected_meta), target = 'row'),
                  options = list(dom = "Bft", buttons = c("copy", "csv"), scrollX = "100px", scrollY = "800px"))
  })

  selected_cells <- reactive({
    selected_rows <- input$brushtable_rows_selected
    rownames(seu$active[[]][brush(),])[selected_rows]
  })

  return(selected_cells)
}


# subsetSeuratui <- function(id) {
#   ns <- NS(id)
#   tagList()
# }
#
#
# subsetSeurat <- function(input, output, session, seu, selected_rows) {
#   ns <- session$ns
#   sub_seu <- reactive({
#     showModal(modalDialog(title = "Subsetting and Recalculating Embeddings",
#                           "This process may take a minute or two!"))
#     seu$gene <- seu$gene[, selected_rows()]
#     seu$gene <- seuratTools::seurat_pipeline(seu$gene, resolution = seq(0.6, 2, by = 0.2))
#     seu$transcript <- seu$transcript[, selected_rows()]
#
#     seu$transcript <- seuratTools::seurat_pipeline(seu$transcript, resolution = seq(0.6, 2, by = 0.2))
#     seu$active <- seu$gene
#     removeModal()
#   })
#   return(sub_seu)
# }


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
  tagList(seuratToolsBox(title = "Differential Expression Settings",
                         shinyWidgets::prettyRadioButtons(ns("diffex_scheme"),
    "Cells to Compare",
    choiceNames = c(
      "Seurat Cluster",
      "Custom"
    ), choiceValues = c("seurat", "custom"),
    selected = "seurat",
    inline = TRUE
  ), conditionalPanel(
    ns = ns,
    condition = "input.diffex_scheme == 'seurat'",
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
  ), conditionalPanel(
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
  ), uiOutput(ns("testChoices")),
  actionButton(
    ns("diffex"),
    "Run Differential Expression"
  ),
  downloadLink(ns("downloadData"), "Download Complete DE Results"),
  DT::dataTableOutput(ns("DT1")),
  width = 12
  ), seuratToolsBox(
    title = "Custom Cluster 1", DT::DTOutput(ns("cc1")),
    width = 12
  ), seuratToolsBox(
    title = "Custom Cluster 2", DT::DTOutput(ns("cc2")),
    width = 12
  ))

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

  default_assay <- reactive({
    req(seu$active)
    if ("integrated" %in% names(seu$active@assays)){
      active_assay <- "integrated"
    } else {
      active_assay <- "RNA"
    }
  })

  observe({
    req(seu$active)
    Seurat::DefaultAssay(seu$active) <- "RNA"
  })

  output$testChoices <- renderUI(
    selectizeInput(ns("diffex_method"),
                                     "Method of Differential Expression",
                                     choices = tests
    )

  )

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
                                     validate(
                                       cells_selected(selected_cells())
                                     )
                                     isolate(selected_cells())
                                   })
  custom_cluster2 <- eventReactive(input$saveClust2,
                                   {
                                     validate(
                                       cells_selected(selected_cells())
                                     )
                                     isolate(selected_cells())
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
                    resolution = input$seuratResolution, diffex_scheme = "seurat", featureType, tests = tests)
    }

    else if (input$diffex_scheme == "custom") {
      cluster1 <- unlist(strsplit(custom_cluster1(),
                                  " "))
      cluster2 <- unlist(strsplit(custom_cluster2(),
                                  " "))
      run_seurat_de(seu$active, cluster1, cluster2,
                    input$customResolution, diffex_scheme = "custom", featureType, tests = tests)
    }
  })

  output$DT1 <- DT::renderDT(de_results()[[input$diffex_method]],
                             extensions = "Buttons", options = list(dom = "Bfptr",
                                                                    buttons = c("copy", "csv"), scrollX = "100px", scrollY = "600px"), class = "display")

  cluster_list <- reactive({
    if (input$diffex_scheme == "seurat"){
      seu_meta <- seu$active[[paste0(DefaultAssay(seu$active), "_snn_res.", input$seuratResolution)]]
      cluster1_cells <- rownames(seu_meta[seu_meta == input$cluster1, , drop = FALSE])
      cluster2_cells <- rownames(seu_meta[seu_meta == input$cluster2, , drop = FALSE])
      list(cluster1 = cluster1_cells, cluster2 = cluster2_cells)
    } else if (input$diffex_scheme == "custom"){
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
    seuratToolsBox(
      title = "Find Markers",
      uiOutput(ns("dplottype")),
      sliderInput(ns("resolution2"), label = "Resolution of clustering algorithm (affects number of clusters)", min = 0.2, max = 2, step = 0.2, value = 0.6),
      numericInput(ns("num_markers"), "Select Number of Markers to Plot for Each Value", value = 5, min = 2, max = 20),
      uiOutput(ns("valueSelect")),
      radioButtons(ns("markerMethod"), "Method of Marker Selection", choices = c("presto", "genesorteR"), selected = "presto", inline = TRUE),
      actionButton(ns("plotDots"), "Plot Markers!"),
      checkboxInput(ns("hidePseudo"), "Hide Pseudogenes", value = TRUE),
      plotly::plotlyOutput(ns("markerplot"), height = 800),
      width = 12
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
#'
findMarkers <- function(input, output, session, seu, plot_types, featureType) {
  ns <- session$ns

  output$dplottype <- renderUI({
    req(seu$active)
    # selected_plot <- ifelse(is.null(selected_plot()), "seurat",
    #                         selected_plot())
    selectizeInput(ns("plottype"), "Variable to Plot", choices = purrr::flatten_chr(plot_types()),
                   selected = "seurat", multiple = TRUE)
  })

  default_assay <- reactive({
    req(seu$active)
    if ("integrated" %in% names(seu$active@assays)){
      active_assay <- "integrated"
    } else {
      active_assay <- "RNA"
    }
  })

  metavar <- reactive({
    req(input$plottype)

    if (input$plottype == "seurat"){
      metavar <- paste0(default_assay(), "_snn_res.", input$resolution2)
    } else {
      metavar <- input$plottype
    }

  })

  output$valueSelect <- renderUI({
    req(seu$active)
    req(metavar())

    choices = levels(seu$active[[]][[metavar()]])

    selectizeInput(ns("displayValues"), "Values to display", multiple = TRUE, choices = choices)
  })

  observe({
    req(default_assay())
    Seurat::DefaultAssay(seu$active) <- "RNA"
  })

  marker_plot <- eventReactive(input$plotDots, {
    plot_markers(seu = seu$active, metavar = metavar(), num_markers = input$num_markers, selected_values = input$displayValues, marker_method = input$markerMethod, featureType = featureType(), hide_pseudo = input$hidePseudo)
  })

  output$markerplot <- plotly::renderPlotly({
    # req(input$displayClusters)
    marker_plot()
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
  seuratToolsBox(
    title = "Histogram (Read Counts, etc.)",
    uiOutput(ns("metavarui")),
    uiOutput(ns("colorbyui")),
    sliderInput(ns("resolution"), "Resolution of clustering algorithm (affects number of clusters)",
      min = 0.2, max = 2, step = 0.2, value = 0.6
    ),
    plotly::plotlyOutput(ns("rcplot"), height = 500),
    collapsed = TRUE)
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
    req(seu$active)
    shiny::selectInput(ns("colorby"), "Variable to Color the Plot by",
                       choices = purrr::flatten_chr(plot_types()), selected = c("seurat"), multiple = FALSE
    )
  })

  output$metavarui <- renderUI({
    req(seu$active)
    shiny::selectInput(ns("metavar"), "Variable for x-axis",
                       choices = purrr::flatten_chr(plot_types()), selected = c("nCount_RNA"), multiple = FALSE
    )
  })

  output$rcplot <- plotly::renderPlotly({
    req(seu$active)
    req(input$colorby)

    if (input$colorby == "seurat") {

      if ("integrated" %in% names(seu$active@assays)){
        active_assay <- "integrated"
      } else {
        active_assay <- "RNA"
      }

      louvain_resolution = paste0(active_assay, "_snn_res.", input$resolution)
      plot_readcount(seu$active, metavar = input$metavar, color.by = louvain_resolution)
    }
    else if (input$colorby %in% purrr::flatten_chr(plot_types())) {
      plot_readcount(seu$active, metavar = input$metavar, color.by = input$colorby)
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
  tagList(fluidRow(seuratToolsBox(
    title = "Transcript Expression per Gene",
    actionButton(ns("plotTrx"), "Plot all transcripts"),
    uiOutput(ns("embeddings")),
   selectizeInput(ns("feature"), "Gene or transcript expression by which to color the plot; eg. 'RXRG'", choices = NULL, selected = NULL),
   uiOutput(ns("plotlys")),
   width = 12) %>%
     default_helper(type = "markdown", content = "allTranscripts"),
   ))
}

#' Plot All Transcripts
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
    req(seu$gene)
    updateSelectizeInput(session, "feature", choices = rownames(seu$gene@assays$RNA), selected = "RXRG", server = TRUE)
  })

  output$embeddings <- renderUI({
    req(seu$active)
    selectizeInput(ns("embedding"), "dimensional reduction method", choices = c("pca", "tsne", "umap"), selected = "umap")
  })

  transcripts <- reactive({

    get_transcripts_from_seu(seu$transcript, input$feature, organism = organism_type())

  })

  pList <- eventReactive(input$plotTrx, {


    pList <- plot_all_transcripts(seu$transcript, seu$gene, transcripts(), input$embedding)

  })

  observe({

    for (i in 1:length(pList())) {
      local({
        my_i <- i
        plotname <- transcripts()[[my_i]]
        output[[plotname]] <- renderPlot({
          pList()[[my_i]]
        })
      })
    }

  })

  output$plotlys <- renderUI({
    req(pList())

    plot_output_list <- purrr::map(names(pList()), ~plotOutput(ns(.x), height = 500))

    do.call(tagList, plot_output_list)
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
plotVelocityui <- function(id){
  ns <- NS(id)
  tagList(
    selectizeInput(ns("embedding"), "dimensional reduction method", choices = c("pca", "tsne", "umap"),
                                     selected = "umap"),
    sliderInput(ns("resolution"), "Resolution of clustering algorithm (affects number of clusters)", min = 0.2, max = 2, step = 0.2, value = 0.6),
    actionButton(ns("calc_velocity"), "calculate velocity"),
    textOutput(ns("velocityFlag")),
    actionButton(ns("plot_velocity"), "plot velocity"),
    radioButtons(ns("plotFormat"), "velocity format", choices = c("arrow", "grid"), selected = "grid"),
    downloadButton(ns("downloadPlot"), label = "Download plots"),
    plotOutput(ns("velocityPlot"), height = "800px")
  ) %>%
    default_helper(type = "markdown", content = "plotVelocity")
}

#' RNA Velocity
#'
#' @param input
#' @param output
#' @param session
#' @param seu
#' @param loom_path
#' @param featureType
#'
#' @return
#' @export
#'
#' @examples
plotVelocity <- function(input, output, session, seu, loom_path, featureType){
  ns <- session$ns

  observeEvent(input$calc_velocity, {
    req(seu$active)
    print(loom_path)

    if(is.null(seu$active@misc$vel)){
      if(is.null(seu[[featureType()]]@misc$vel)){
        showModal(modalDialog("calculating velocity", footer=NULL))
        seu[[featureType()]] <- velocyto_assay(seu[[featureType()]], loom_path)
        seu$active@misc$vel <- seu[[featureType()]]@misc$vel
        removeModal()
      } else if (!is.null(seu[[featureType()]]@misc$vel)) {
        seu$active@misc$vel <- seu[[featureType()]]@misc$vel
      }

    }
  })

  velocity_flag <- eventReactive(input$calc_velocity, {
    req(seu$active)
    if(!is.null(seu$active@misc$vel)){
      "Velocity Calculated for this dataset"
    }
  })

  output$velocityFlag <- renderText({
    velocity_flag()
  })

  velocity <- reactive({
    req(seu$active)

    if(!is.null(seu$active@misc$vel)){
      return(seu$active@misc$vel)
    } else {
      FALSE
    }
  })


  plotInput <- function(plot_format = "grid"){

    showModal(modalDialog("Loading Plots", footer=NULL))

    if ("integrated" %in% names(seu$active@assays)) {
      default_assay = "integrated"
    } else {
      default_assay = "RNA"
    }

    cluster_resolution = paste0(default_assay, "_snn_res.", input$resolution)

    cell.colors <- tibble::as_tibble(seu$active[[cluster_resolution]], rownames = "cellid") %>%
      tibble::deframe() %>%
      as.factor()

    levels(cell.colors) <- scales::hue_pal()(length(levels(cell.colors)))

    if (!is.null(seu$active@misc$cc)){
      plot_velocity_arrows(seu$active, velocity(), reduction = input$embedding, cell.colors, plot_format = plot_format, cc = seu$active@misc$cc)
    } else {
      plot_velocity_arrows(seu$active, velocity(), reduction = input$embedding, cell.colors, plot_format = plot_format)
    }

    removeModal()
  }

  velocityPlot <- eventReactive(input$plot_velocity, {
    req(seu$active)

    validate(
      need(velocity(), "Please calculate velocity")
    )
    plotInput(plot_format = input$plotFormat)
  })

  output$downloadPlot <- downloadHandler(
    filename = function() { paste("velocity", '.pdf', sep='') },
    content = function(file) {
      pdf(file)
      plotInput(plot_format = input$plotFormat)
      dev.off()
    })

  output$velocityPlot <- renderPlot({
    velocityPlot()
  })

  return(seu[[featureType()]])
}


#' Title
#'
#' @param id
#'
#' @return
#' @export
#'
#' @examples
monocleui <- function(id){
    ns <- NS(id)
    tagList(
      fluidRow(
        seuratToolsBox(
          plotly::plotlyOutput(ns("seudimplot"), height = 500),
          width = 6
          # plotDimRedui(ns("plotdimred"))
        ),
        seuratToolsBox(
          title = "Pseudotime Settings",
          actionButton(ns("subsetSeurat"), "Subset Seurat before Pseudotime Calculation"),
          actionButton(ns("calcCDS"), "Calculate Pseudotime"),
          sliderInput(ns("cdsResolution"), "Resolution of clustering algorithm (affects number of clusters)",
                      min = 0.2, max = 2, step = 0.2, value = 0.6),
          actionButton(ns("subsetCells"), "Subset Monocle Object After Pseudotime Calculation"),
          uiOutput(ns("rootCellsui")),
          actionButton(ns("plotPseudotime"), "Calculate Pseudotime With Root Cells"),
          checkboxInput(ns("flipPtime"), "Invert Pseudotime", value = TRUE),
          width = 6
        )
      ),
      fluidRow(
        seuratToolsBox(
          title = "Embedding Plot",
          selectizeInput(ns("plottype1"), "Variable to Plot", choices = c(Seurat = "seurat"), selected = "Seurat", multiple = TRUE),
          selectizeInput(ns("customFeature1"), "Gene or transcript expression by which to color the plot",
                         choices = NULL, multiple = FALSE),
          uiOutput(ns("moduleSelect1")),
          plotly::plotlyOutput(ns("monoclePlot1")),
          width = 6
        ),
        seuratToolsBox(
          title = "Embedding Plot",
          selectizeInput(ns("plottype2"), "Variable to Plot", choices = c(Seurat = "seurat"), selected = "Seurat", multiple = TRUE),
          selectizeInput(ns("customFeature2"), "gene or transcript on which to color the plot",
                         choices = NULL, multiple = FALSE),
          uiOutput(ns("moduleSelect2")),
          plotly::plotlyOutput(ns("monoclePlot2")),
          width = 6
        )
      ),
      fluidRow(
        seuratToolsBox(title = "calculate pseudotime",
                       actionButton(ns("calcPtimeGenes"), "Find Pseudotime Correlated Genes"),
            uiOutput(ns("partitionSelect")),
            uiOutput(ns("genePlotQuery2")),
            uiOutput(ns("ptimeGenes")),
            width = 12
            )
      ),
      fluidRow(
        seuratToolsBox(
          title = "Heatmap",
          radioButtons(ns("pickHeatmap"), "heatmap on clusters or cells?", choices = c(clusters = TRUE, cells = FALSE), selected = TRUE),
          iheatmapr::iheatmaprOutput(ns("monocleHeatmap"), width = "800px", height = "600px")
        ),
        seuratToolsBox(
          div(DT::dataTableOutput(ns("moduleTable")), style = "font-size: 75%")
        )
      )
      )
}

#' Title
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
#'
#'
#'
monocle <- function(input, output, session, seu, plot_types, featureType,
                    organism_type, reductions){
    ns <- session$ns

    cds <- reactiveValues(selected = "traj")
    cds_plot_types <- reactiveVal(c(Pseudotime = "pseudotime", Module = "module"))
    myplot_types <- reactive({
      c(purrr::flatten_chr(plot_types()), cds_plot_types())
    })

    observe({
      seu$monocle <- seu$active
    })

    seudimplot <- reactive({
      req(seu$monocle)
      if ("integrated" %in% names(seu$active@assays)) {
        active_assay <- "integrated"
      }
      else {
        active_assay <- "RNA"
      }

      louvain_resolution = reactive({
        paste0(active_assay, "_snn_res.", input$cdsResolution)
      })

      plot_var(seu$monocle, embedding = "umap", group = louvain_resolution())

    })

    output$seudimplot <- plotly::renderPlotly({
      seudimplot()
    })

    # callModule(plotDimRed, "plotdimred", seu, plot_types, featureType,
    #            organism_type, reductions)

    seubrush <- reactive({
      req(seu$monocle)
      d <- plotly::event_data("plotly_selected")
      if (is.null(d)) {
        msg <- "Click and drag events (i.e. select/lasso) appear here (double-click to clear)"
        return(d)
      }
      else {
        selected_cells <- colnames(seu$monocle)[as.numeric(d$key)]
      }
    })

    observeEvent(input$subsetSeurat, {
      req(seu$monocle)
      print(seubrush())
      seu$monocle <- seu$monocle[,seubrush()]
    })

    observeEvent(input$calcCDS, {
      req(seu$monocle)
        cds$traj <- convert_seu_to_cds(seu$monocle, resolution = input$cdsResolution)
        cds$traj <- learn_graph_by_resolution(cds$traj,
                                              seu$monocle,
                                              resolution = input$cdsResolution)
        updateSelectizeInput(session, "plottype1", selected = "seurat", choices = myplot_types())
        updateSelectizeInput(session, "customFeature1", choices = rownames(cds$traj), server = TRUE)
        updateSelectizeInput(session, "plottype2", selected = "seurat", choices = myplot_types())
        updateSelectizeInput(session, "customFeature2", choices = rownames(cds$traj), server = TRUE)
    })

    selected_plot <- reactiveVal()

    output$monoclePlot1 <- plotly::renderPlotly({
      req(input$plottype1)
      req(cds$traj)
      print(cds$selected)
      if (input$plottype1 == "seurat") {
        cluster_resolution = reactive({
          paste0("integrated", "_snn_res.", input$cdsResolution)
        })
        plot_cds(cds$traj, color_cells_by = cluster_resolution())
      } else if (input$plottype1 == "pseudotime"){
        plot_pseudotime(cds$traj, color_cells_by = "pseudotime", resolution = input$cdsResolution)
      } else if (input$plottype1 == "custom") {
        plot_monocle_features(cds$traj, genes = input$customFeature1, monocle_heatmap()$agg_mat)
      } else if (input$plottype1 == "module") {
        print(monocle_heatmap()$module_table)
        print(input$plotModule1)
        genes = monocle_heatmap()$module_table %>%
          filter(module %in% input$plotModule1) %>%
          dplyr::mutate(module = factor(module))
        plot_monocle_features(cds$traj, genes = genes, monocle_heatmap()$agg_mat)
      } else {
        plot_cds(cds$traj, color_cells_by = input$plottype1)
      }
    })

    output$monoclePlot2 <- plotly::renderPlotly({
      req(input$plottype2)
      req(cds$traj)
      print(cds$selected)
      if (input$plottype2 == "seurat") {
        cluster_resolution = reactive({
          paste0("integrated", "_snn_res.", input$cdsResolution)
        })
        plot_cds(cds$traj, color_cells_by = cluster_resolution())
      } else if (input$plottype2 == "pseudotime"){
        plot_pseudotime(cds$traj, color_cells_by = "pseudotime", resolution = input$cdsResolution)
      } else if (input$plottype2 == "custom") {
        plot_monocle_features(cds$traj, genes = input$customFeature2, monocle_heatmap()$agg_mat)
      } else if (input$plottype2 == "module") {
        print(monocle_heatmap()$module_table)
        print(input$plotModule2)

        genes = monocle_heatmap()$module_table %>%
          filter(module %in% input$plotModule2) %>%
          dplyr::mutate(module = factor(module))
        plot_monocle_features(cds$traj, genes = genes, monocle_heatmap()$agg_mat)
      } else {
        plot_cds(cds$traj, color_cells_by = input$plottype2)
      }
    })

    cdsbrush <- reactive({
      req(cds$traj)
      d <- plotly::event_data("plotly_selected")
      if (is.null(d)) {
        msg <- "Click and drag events (i.e. select/lasso) appear here (double-click to clear)"
        return(d)
      }
      else {
        selected_cells <- colnames(cds$traj)[as.numeric(d$key)]
      }
    })

    observeEvent(input$subsetCells, {
      req(cds$traj)
      print(cdsbrush())
      cds$traj <- cds$traj[,cdsbrush()]
    })

    output$rootCellsui <- renderUI({
      selectizeInput(ns("rootCells"), "Choose Root Cells", choices = c("Choose Root Cells" = "", colnames(cds$traj)), multiple = TRUE)
    })

    observeEvent(input$plotPseudotime, {
      req(cds$traj)
      req(input$rootCells)
      cds$traj <- monocle3::order_cells(cds$traj, root_cells = input$rootCells)
      if(input$flipPtime){
        cds$traj <- flip_pseudotime(cds$traj)
      }
      updateSelectizeInput(session, "plottype1", selected = "pseudotime", choices = myplot_types())
      updateSelectizeInput(session, "plottype2", selected = "pseudotime", choices = myplot_types())
      cds$selected = "ptime"
    })

    observeEvent(input$calcPtimeGenes, {
      if (req(cds$selected) == "ptime"){
        showModal(modalDialog(
          title = "Calculating features that vary over pseudotime",
          "This process may take a minute or two!"
        ))

        cds_pr_test_res = monocle3::graph_test(cds$traj, neighbor_graph="principal_graph", cores=4, expression_family = "negbinom")
        removeModal()

        cds$traj@metadata[["diff_genes"]] <- cds_pr_test_res
        cds$selected = "diff_genes"

      }
    })
    cds_pr_test_res <- reactive({
      if (req(cds$selected) == "diff_genes"){
        cds_pr_test_res <- cds$traj@metadata$diff_genes

        cds_pr_test_res <-
          cds_pr_test_res %>%
          subset(q_value < 0.05) %>%
          dplyr::arrange(q_value) %>%
          dplyr::select(-status)

      }
    })

    observe({
      if (req(cds$selected) == "diff_genes"){

      output$genePlotQuery2 <- renderUI({
        selectizeInput(ns("genePlotQuery1"), "Pick Gene to Plot on Pseudotime", choices = rownames(cds_pr_test_res()), multiple = TRUE, selected = rownames(cds_pr_test_res())[1])
      })

      output$partitionSelect <- renderUI({
        selectizeInput(ns("partitions"), "Select a Partition to Plot", choices = levels(monocle3::partitions(cds$traj)), multiple = FALSE)
      })

      output$ptimeGenesLinePlot <- plotly::renderPlotly({

        genes_in_pseudotime <- prep_plot_genes_in_pseudotime(cds$traj, input$genePlotQuery1, input$cdsResolution)

        genes_in_pseudotime <-
          genes_in_pseudotime %>%
          plotly::ggplotly(height = 400) %>%
          plotly_settings() %>%
          plotly::toWebGL() %>%
          # plotly::partial_bundle() %>%
          identity()

      })

      # mymarker

      output$ptimeGenesDT <- DT::renderDT({

        DT::datatable(cds_pr_test_res(), extensions = 'Buttons',
                      options = list(dom = "Bft", buttons = c("copy", "csv"), scrollX = "100px", scrollY = "400px"))

      })

      output$ptimeGenes <- renderUI({
        tagList(
          seuratToolsBox(
            div(DT::DTOutput(ns("ptimeGenesDT")), style = "font-size: 75%"),
              width = 4),
          seuratToolsBox(plotly::plotlyOutput(ns("ptimeGenesLinePlot")),
              width = 8)
                )
      })
      }
    })

      monocle_heatmap <- reactive({
        req(cds$traj)
        req(cds_pr_test_res())
        monocle_module_heatmap(cds$traj, rownames(cds_pr_test_res()), input$cdsResolution, collapse_cols = input$pickHeatmap)
      })

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
        output$monocleHeatmap <- iheatmapr::renderIheatmap({
          monocle_heatmap()$module_heatmap
        })

        output$moduleTable <- DT::renderDataTable({
          DT::datatable(monocle_heatmap()$module_table,
                        extensions = "Buttons",
                        options = list(dom = "Bft", buttons = c("copy",
                                                                "csv"), scrollX = "100px", scrollY = "400px"))
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
pathwayEnrichmentui <- function(id){
    ns <- NS(id)
        seuratToolsBox(
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
pathwayEnrichment <- function(input, output, session, seu, featureType){
    ns <- session$ns

    ##----------------------------------------------------------------------------##
    ## Tab: Enriched pathways
    ##----------------------------------------------------------------------------##

    ##----------------------------------------------------------------------------##
    ## Clusters.
    ##----------------------------------------------------------------------------##

    enriched_pathways <- eventReactive(input$calcPathwayEnrichment, {
      req(seu$active)
      if (featureType() == "gene"){
        enriched_seu <- tryCatch(getEnrichedPathways(seu$active), error = function(e) e)
        enrichr_available <- !any(class(enriched_seu) == "error")
        if(enrichr_available){
          seu$active <- enriched_seu
        }
      }

      seu$active@misc$enriched_pathways
    })

    # UI element: choose source for pathway enrichement results (currently Enrichr or GSVA)
    output$enriched_pathways_by_cluster_select_source_UI <- renderUI({
      req(seu$active)
      if (is.null(enriched_pathways()) ) {
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
      req(seu$active)
      req(input$enriched_pathways_by_cluster_select_source)
      if ( input$enriched_pathways_by_cluster_select_source == "enrichr" ) {
        if ( !is.null(enriched_pathways()$enrichr$by_cluster) ) {
          if ( is.list(enriched_pathways()$enrichr$by_cluster) ) {
            tagList(
              fluidRow(
                column(4,
                       uiOutput(ns("enriched_pathways_by_cluster_select_cluster_UI"))
                ),
                column(8,
                       uiOutput(ns("enriched_pathways_by_cluster_select_db_UI"))
                )
              ),
              DT::dataTableOutput(ns("enriched_pathways_by_cluster_table_present"))
            )
          } else if ( enriched_pathways()$enrichr$by_cluster == "no_markers_found" ) {
            textOutput(ns("enriched_pathways_by_cluster_table_no_markers_found"))
          }
        } else {
          textOutput(ns("enriched_pathways_by_cluster_table_missing_enrichr"))
        }
      }
    })

    # UI element: choose cluster
    output$enriched_pathways_by_cluster_select_cluster_UI <- renderUI({
      req(seu$active)
      req(input$enriched_pathways_by_cluster_select_source)
      if ( input$enriched_pathways_by_cluster_select_source == 'enrichr' ) {
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
      if ( input$enriched_pathways_by_cluster_select_source == "enrichr" & is.data.frame(enriched_pathways()$enrichr$by_cluster) ) {

        format_pathway_table(enriched_pathways()$enrichr$by_cluster,
                             input$enriched_pathways_by_cluster_select_cluster,
                             input$enriched_pathways_by_cluster_select_db)

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

techInfoui <- function(id){
    ns <- NS(id)
    fluidRow(
      seuratToolsBox(
        title = "Information about samples and analysis",
        htmlOutput(ns("sample_info_general")),
        width = 12
      )
    )
    }

techInfo <- function(input, output, session, seu){
    ns <- session$ns
    ##----------------------------------------------------------------------------##
    ## Tab: Analysis info.
    ##----------------------------------------------------------------------------##

    misc <- reactive({
      req(seu$active)
      seu$active@misc
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
        for ( i in 1:length(info_R_raw) ) {
          info_R <- paste(info_R, "<br>", info_R_raw[i])
        }
        paste0(
          info,
          "<strong><u>Technical info (package versions)</u></strong>",
          "<ul>",
          "<li><strong>seuratTools version:</strong> ",
          misc()$experiment$technical_info$seuratTools_version,
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
        if ( !is.null(misc()$technical_info$R) ) {
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
    seuratToolsBox(
      title = "Plot Coverage",
      selectizeInput(ns("geneSelect"), "Select a Gene", choices = NULL, selected = "RXRG", multiple = FALSE),
      selectizeInput(ns("varSelect"), "Color by Variable", choices = NULL, multiple = FALSE),
      actionButton(ns("plotCoverage"), "Plot Coverage"),
      downloadButton(ns("downloadPlot"), "Download Coverage Plot"),
      plotOutput(ns("coveragePlot"), height = "1000px"),
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
#' @param plot_types
#' @param bigwig_dir
#' @param organism_type
#'
#' @return
#' @export
#'
#' @examples
plotCoverage <- function(input, output, session, seu, plot_types, proj_dir, organism_type = "human") {

    observe({
      req(seu$active)
      updateSelectizeInput(session, "geneSelect", choices = rownames(seu$active), server = TRUE)
      updateSelectizeInput(session, "varSelect", choices = plot_types())
    })

    bigwig_tbl <- reactive({
      load_bigwigs(seu$active, proj_dir())
    })

    coverage_plot <- eventReactive(input$plotCoverage, {
      req(seu$active)
      req(bigwig_tbl())

      plot_gene_coverage_by_var(genes_of_interest = input$geneSelect,
                                cell_metadata = seu$active@meta.data,
                                bigwig_tbl = bigwig_tbl(),
                                var_of_interest = input$varSelect,
                                edb = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86,
                                mean_only = FALSE)
    })

    output$coveragePlot <- renderPlot({
      coverage_plot()
    })

    output$downloadPlot <- downloadHandler(
      filename = function() { paste("coverage", '.pdf', sep='') },
      content = function(file) {
        ggsave(file, coverage_plot(), width = 16, height = 12)
      })

}
