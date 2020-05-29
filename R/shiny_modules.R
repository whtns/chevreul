#' Title
#'
#' @param id
#'
#' @return
#' @export
#'
#' @examples
downloadTable_UI <- function(id) {
  ns <- NS(id)
  tagList(
    uiOutput(ns("geaDownloadButton"))
    # downloadButton(ns("mydata"), "my data")
  )
}

#' Title
#'
#' @param input
#' @param output
#' @param session
#' @param mytable
#'
#' @return
#' @export
#'
#' @examples
downloadTable <- function(input, output, session, mytable) {
  ns <- session$ns
  results <- reactive({
    req(mytable())
    mytable()$results
  })

  reportLink <- reactive({
    req(mytable())
    mytable()$report
  })

  output$geaDownloadButton <- renderUI({
    req(mytable())
    downloadButton(ns("mydata"), "Gene Enrichment Results")
  })

  # Downloadable csv of selected dataset ----
  output$mydata <- downloadHandler(
    filename = function(){
      paste0(fs::path_file(reportLink()), ".csv")
      },
    content = function(myfile) {
      write.csv(results(), myfile, row.names = FALSE)
    }
  )
}

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
    # textOutput(ns("checkSeu")),
    plotOutput(ns("clustree"))

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
  tagList(uiOutput(ns("vln_split")), uiOutput(ns("split_val")),
          uiOutput(ns("vln_group")), selectizeInput(ns("customFeature"),
                                                    "gene or transcript on which to color the plot; eg. 'RXRG' or 'ENST00000488147'",
                                                    choices = NULL, multiple = TRUE), plotOutput(ns("vplot"),
                                                                                                 height = 750))
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
    updateSelectizeInput(session, "customFeature", choices = rownames(seu$active),
                         selected = prefill_feature(), server = TRUE)
  })

  output$featuretext <- renderUI({
    textInput(ns("customFeature"), "gene or transcript on which to color the plot; eg. 'RXRG' or 'ENST00000488147'",
              value = prefill_feature())
  })

  output$vln_split <- renderUI({
    req(seu$active)
    selectizeInput(ns("vlnSplit"), "choose variable filter by",
                   choices = c("None", colnames(seu$active[[]])), selected = "batch")
  })
  output$split_val <- renderUI({
    req(seu$active)
    req(input$vlnSplit)
    if (input$vlnSplit == "None") {
      choices = "None"
    }
    else {
      choices = c("None", unique(seu$active[[input$vlnSplit]][,
                                                              1]))
    }
    selectizeInput(ns("splitVal"), "choose value to filter by",
                   choices = choices)
  })
  output$vln_group <- renderUI({
    req(seu$active)
    selectizeInput(ns("vlnGroup"), "choose variable to group by",
                   choices = colnames(seu$active[[]]), selected = "batch")
  })
  output$vplot <- renderPlot({
    req(input$customFeature)
    req(input$vlnGroup)
    req(input$vlnSplit)
    if (!("None" %in% c(input$vlnSplit, input$splitVal))) {
      selected_cells <- tibble::as_tibble(seu$active[[input$vlnSplit]],
                                          rownames = "sample_id") %>% dplyr::filter(!!sym(input$vlnSplit) ==
                                                                                      input$splitVal) %>% dplyr::pull(sample_id)
      sub_seu <- seu$active[, selected_cells]
    }
    else {
      sub_seu <- seu$active
    }
    plot_violin(sub_seu, plot_var = input$vlnGroup, features = input$customFeature)
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
  tagList(uiOutput(ns("hm_group")), selectizeInput(ns("customFeature"),
                                                   "gene or transcript on which to color the plot; eg. 'RXRG' or 'ENST00000488147'",
                                                   choices = NULL, multiple = TRUE), plotOutput(ns("heatmap"),
                                                                                                height = 750))
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

    updateSelectizeInput(session, "customFeature", choices = rownames(seu$active),
                         selected = preset_features, server = TRUE)
  })
  output$hm_group <- renderUI({
    req(seu$active)
    selectizeInput(ns("hmGroup"), "choose variable to group by",
                   choices = colnames(seu$active[[]]), selected = "batch")
  })
  output$heatmap <- renderPlot({
    req(input$customFeature)
    req(input$hmGroup)

    plot_heatmap(seu$active, group.by = input$hmGroup, features = input$customFeature)
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
    box(
    uiOutput(ns("colNames")),
    textInput(ns("newCol"), "provide a name for the new column"),
    actionButton(ns("mergeCol"), "Merge Selected Columns"),
    checkboxInput(ns("header"), "Header", TRUE),
    fileInput(ns("addCols"), "Choose CSV File of metadata with cell names in first column",
              accept = c(
                "text/csv",
                "text/comma-separated-values,text/plain",
                ".csv")
    ),
    downloadLink(ns("downloadMetadata"), "Download Metadata"),
    DT::DTOutput(ns("seuTable")),
    width = 12
    )

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
    # req(seu$active)
    meta$old <- data.frame(seu$active[[]]) %>%
      identity()

    # na_cols <- purrr::map_lgl(meta$old, ~all(is.na(.x)))
    # cluster_cols <- grepl("^cluster|snn_res", colnames(meta$old))
    #
    # keep_cols <- !(na_cols | cluster_cols)
    #
    # meta$old <- meta$old[,keep_cols]
  })

  seuColNames <- reactive({
    seuColNames <- colnames(seu$gene[[]]) %>%
      purrr::set_names(.)
  })

  output$colNames <- renderUI({
    selectizeInput(ns("col_names"), "Seurat Object Metadata Columns", choices = seuColNames(), multiple = TRUE)
  })

  observeEvent(input$mergeCol, {

    combined_cols <- combine_cols(seu, input$col_names, input$newCol)
    meta$new <- combined_cols

    for (i in names(seu)){
      seu[[i]] <- Seurat::AddMetaData(seu[[i]], meta$new)
    }

    meta$old <- cbind(meta$old, meta$new)

  })

  observeEvent(input$addCols, {

    inFile <- input$addCols

    if (is.null(inFile))
      return(NULL)

    meta$new <- format_new_metadata(inFile$datapath, header = input$header, row.names = 1)

    # meta$new <- read.csv(inFile$datapath, header = input$header, row.names = 1)

    for (i in names(seu)){
      seu[[i]] <- Seurat::AddMetaData(seu[[i]], meta$new)
      print(colnames(seu[[i]][[]]))
    }

    meta$old <- cbind(meta$old, meta$new)

  })

  output$seuTable <- DT::renderDT({
    # req(meta$old)

    DT::datatable(meta$old, extensions = 'Buttons', options = list(buttons = c("copy", "csv"),
                                                                   paging = TRUE, pageLength = 15, scrollX = "100px"))


  })

  output$downloadMetadata <- downloadHandler(
    filename = function() {
      paste("data-", Sys.Date(), ".csv", sep="")
    },
    content = function(file) {
      write.csv(meta$old, file)
    })

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
      box(
        actionButton(ns("integrateAction"), "Integrate Selected Projects"),
        DT::dataTableOutput(ns("myDatatable"))
        ),
      box(
        textOutput(ns("integrationComplete")),
        shinyjs::useShinyjs(),
        textOutput(ns("integrationMessages")),
        textOutput(ns("integrationResult")),
        shinyFiles::shinySaveButton(ns("saveIntegratedProject"), "Save Integrated Project", "Save project as...")
      )
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

            mergedSeus(integration_workflow(selectedProjects()))

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
  tagList(uiOutput(ns("dplottype")),
          uiOutput(ns("embeddings")),
          fluidRow(
            column(2,
                   selectizeInput(ns("dim1"), "Dimension 1", choices = seq(1, 99), selected = 1)
                   ),
            column(2,
                   selectizeInput(ns("dim2"), "Dimension 2", choices = seq(1, 99), selected = 2)
                   )
            ),
          selectizeInput(ns("customFeature"), "gene or transcript on which to color the plot; eg. 'RXRG' or 'ENST00000488147'",
                                                                                                                                                                                      choices = NULL, multiple = TRUE), sliderInput(ns("resolution"),
                                                                                                                                                                                                                                    "Resolution of clustering algorithm (affects number of clusters)", min = 0.2, max = 2, step = 0.2, value = 0.6),
          plotly::plotlyOutput(ns("dplot"), height = 750))
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
    radioButtons(ns("embedding"), "dimensional reduction method", choices = reductions(), inline = TRUE)
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
    updateSelectizeInput(session, "customFeature", choices = rownames(seu$active),
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
               embedding = input$embedding, group = NULL)
    }
    else {
      if (input$plottype == "custom") {
        plot_feature(seu$active, dims = c(input$dim1,
                                          input$dim2), embedding = input$embedding,
                     features = input$customFeature)
      }
      else if (input$plottype %in% plot_types()$continuous_vars) {
        plot_feature(seu$active, dims = c(input$dim1,
                                          input$dim2), embedding = input$embedding,
                     features = input$plottype)
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
                 embedding = input$embedding, group = louvain_resolution())
      }
      else if (input$plottype %in% plot_types()$category_vars) {
        plot_var(seu$active, dims = c(input$dim1, input$dim2),
                 embedding = input$embedding, group = input$plottype)
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
  tagList(box(shinyWidgets::prettyRadioButtons(ns("diffex_scheme"),
    "Cells to Compare",
    choiceNames = c(
      "Seurat Cluster",
      "Custom"
    ), choiceValues = c("seurat", "custom"),
    selected = "seurat"
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
    shinyWidgets::actionBttn(
      ns("saveClust1"),
      "Save to Custom Cluster 1"
    ), shinyWidgets::actionBttn(
      ns("saveClust2"),
      "Save to Custom Cluster 2"
    )
  ), uiOutput(ns("testChoices")),
  shinyWidgets::actionBttn(
    ns("diffex"),
    "Run Differential Expression"
  ),
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

  output$testChoices <- renderUI(
    shinyWidgets::prettyRadioButtons(ns("diffex_method"),
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
                                                                    buttons = c("copy", "csv"), scrollX = "100px", pageLength = 20, paging = FALSE), class = "display")

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
    fluidRow(
      box(
        actionButton(ns("enrichmentAction"), "Run Enrichment Analysis"),
        textOutput(ns("enrichmentMessages")),
        radioButtons(ns("enrichmentMethod"), "Enrichment Method to Use:",
                     c("Gene Set Enrichment Analysis" = "gsea",
                       "GO Over-representation Analysis" = "ora",
                       "GO Network Analysis" = "nbea"),
                     selected = c("ora")),
      ),
      box(
        fileInput(ns("uploadDiffex"), "Choose CSV File Differential Expression Results",
                  accept = c(
                    "text/csv",
                    "text/comma-separated-values,text/plain",
                    ".csv")
        )
      )
    ),
    htmlOutput(ns("map"))
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

    observeEvent(input$uploadDiffex, {

      inFile <- input$addCols

      if (is.null(inFile))
        return(NULL)

      # meta$new <- read.csv(inFile$datapath, header = input$header, row.names = 1)

    })

    enrichmentReport <- eventReactive(input$enrichmentAction, {
      withCallingHandlers({
        shinyjs::html("enrichmentMessages", "")
        message("Beginning")

        # showModal(modalDialog("Calculating Functional Enrichment", footer=NULL))
        enrichmentReport <- run_enrichmentbrowser(seu = seu$active,
                              cluster_list = diffex_results$cluster_list(),
                              de_results = diffex_results$de_results(),
                              enrichment_method = input$enrichmentMethod)

        # enrichmentReport <- "enrichmentbrowser2/mainpage.html"
        # removeModal()

        # zip::zipr("test.zip", "enrichmentreport")

        return(enrichmentReport)

    },
      message = function(m) {
        shinyjs::html(id = "enrichmentMessages", html = paste0("Running Functional Enrichment Analysis: ", m$message), add = FALSE)
      })
    })

    output$reportLink <- renderUI({
      tags$a("Results of Functional Enrichment Analysis", target = "_blank", href = enrichmentReport()$report)
    })

    output$map <- renderUI({
      tags$iframe(seamless="seamless", src= enrichmentReport()$report, width=1000, height=800)
    })

    return(enrichmentReport)


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
      numericInput(ns("num_markers"), "Select Number of Markers to Plot for Each Cluster", value = 5, min = 2, max = 20),
      uiOutput(ns("clusterSelect")),
      actionButton(ns("plotDots"), "Plot Markers!"),
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
#'
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

  output$clusterSelect <- renderUI({
    req(seu$active)
    req(resolution())

    choices = levels(seu$active[[]][[resolution()]])

    selectizeInput(ns("displayClusters"), "Clusters to display", multiple = TRUE, choices = choices)
  })

  marker_plot <- eventReactive(input$plotDots, {
    plot_markers(seu = seu$active, resolution = resolution(), num_markers = input$num_markers, selected_clusters = input$displayClusters)
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
  tagList(fluidRow(box(radioButtons(ns("embedding"), "dimensional reduction method", choices = c("pca", "tsne", "umap"), inline = TRUE),
                       textInput(ns("feature"), "gene on which to colour the plot; eg. 'RXRG'"),
                       # uiOutput(ns("outfile")),
                       # uiOutput(ns("downloadPlot")),
                       width = 12)), fluidRow(uiOutput(ns("plotlys"))))
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

  output$embeddings <- renderUI({
    req(seu$active)
    radioButtons(ns("embedding"), "dimensional reduction method", choices = c("pca", "tsne", "umap"), inline = TRUE)
  })

  transcripts <- reactiveValues()
  transcripts <- reactive({
    req(featureType())
    req(organism_type())
    req(input$feature)
    req(seu)
    transcripts <- get_transcripts_from_seu(seu, input$feature, organism = organism_type())

  })

  pList <- reactive({
    req(seu$active)

    if(featureType() == "gene"){
      # browser()
      transcript_cols <- as.data.frame(t(as.matrix(seu$transcript[["RNA"]][transcripts(),])))

      cells <- rownames(transcript_cols)
      transcript_cols <- as.list(transcript_cols) %>%
        purrr::map(~purrr::set_names(.x, cells))

      seu$gene[[transcripts()]] <- transcript_cols

      pList <- purrr::map(transcripts(), ~plot_feature(seu$gene,
                                                       embedding = input$embedding, features = .x))
      names(pList) <- transcripts()

    } else if (featureType() == "transcript"){
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
              value = paste0(input$feature, "_transcripts", "_clustered_by_", featureType(), ".pdf"))
  })

  output$plots <-
    downloadHandler(filename = function() {
      input$outfile
    },
    content = function(file) {
      pdf(file)
      lapply(pList(), print)
      dev.off()
  })

  output$downloadPlot <- downloadHandler(
    filename = function() { paste(input$dataset, '.png', sep='') },
    content = function(file) {
      png(file)
      print(plotInput())
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
plotVelocityui <- function(id){
  ns <- NS(id)
  tagList(
    shinyWidgets::prettyRadioButtons(ns("embedding"),
                                     "dimensional reduction method", choices = c("pca", "tsne", "umap"),
                                     selected = "umap", inline = TRUE),
    sliderInput(ns("resolution"), "Resolution of clustering algorithm (affects number of clusters)", min = 0.2, max = 2, step = 0.2, value = 0.6),
    actionButton(ns("calc_velocity"), "calculate velocity"),
    textOutput(ns("velocityFlag")),
    actionButton(ns("plot_velocity"), "plot velocity"),
    radioButtons(ns("plotFormat"), "velocity format", choices = c("arrow", "grid"), selected = "grid"),
    downloadButton(ns("downloadPlot"), label = "Download plots"),
    plotOutput(ns("velocityPlot"), height = "800px")
  )
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

    plot_velocity_arrows(seu$active, velocity(), reduction = input$embedding, cell.colors, plot_format = plot_format)
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
    filename = function() { paste("velocity", '.png', sep='') },
    content = function(file) {
      png(file)
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
        box(
          plotly::plotlyOutput(ns("seudimplot")),
          width = 6
          # plotDimRedui(ns("plotdimred"))
        ),
        box(
          actionButton(ns("subsetSeurat"), "Subset Seurat before Pseudotime Calculation"),
          actionButton(ns("calcCDS"), "Calculate Pseudotime"),
          sliderInput(ns("cdsResolution"), "Resolution of clustering algorithm (affects number of clusters)",
                      min = 0.2, max = 2, step = 0.2, value = 0.6),
          actionButton(ns("subsetCells"), "Subset Monocle Object After Pseudotime Calculation"),
          uiOutput(ns("rootCellsui")),
          actionButton(ns("plotPseudotime"), "Calculate Pseudotime With Root Cells"),
          width = 6
        )
      ),
      fluidRow(
        box(
          selectizeInput(ns("plottype1"), "Variable to Plot", choices = c(Seurat = "seurat"), selected = "Seurat", multiple = TRUE),
          selectizeInput(ns("customFeature1"), "gene or transcript on which to color the plot",
                         choices = NULL, multiple = FALSE),
          uiOutput(ns("moduleSelect1")),
          plotly::plotlyOutput(ns("monoclePlot1")),
          width = 6
        ),
        box(
          selectizeInput(ns("plottype2"), "Variable to Plot", choices = c(Seurat = "seurat"), selected = "Seurat", multiple = TRUE),
          selectizeInput(ns("customFeature2"), "gene or transcript on which to color the plot",
                         choices = NULL, multiple = FALSE),
          uiOutput(ns("moduleSelect2")),
          plotly::plotlyOutput(ns("monoclePlot2")),
          width = 6
        )
      ),
      fluidRow(
        box(actionButton(ns("calcPtimeGenes"), "Find Pseudotime Correlated Genes"),
            uiOutput(ns("partitionSelect")),
            uiOutput(ns("genePlotQuery2")),
            uiOutput(ns("ptimeGenes")),
            width = 12
            )
      ),
      fluidRow(
        box(
          radioButtons(ns("pickHeatmap"), "heatmap on clusters or cells?", choices = c(clusters = TRUE, cells = FALSE), selected = TRUE),
          iheatmapr::iheatmaprOutput(ns("monocleHeatmap"), width = "600px", height = "600px")
        ),
        box(
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
      cds$traj <- monocle3::order_cells(cds$traj, root_cells = input$rootCells)
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

        cds_pr_test_res = monocle3::graph_test(cds$traj, neighbor_graph="principal_graph", cores=4)
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
          box(
            div(DT::DTOutput(ns("ptimeGenesDT")), style = "font-size: 75%"),
              width = 4),
          box(plotly::plotlyOutput(ns("ptimeGenesLinePlot")),
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
