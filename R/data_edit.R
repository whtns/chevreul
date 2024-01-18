#' Reformat Seurat Object Metadata UI
#'
#' @param id
#'
#' @return
#' @export
#'
#' @examples
#' @importFrom DataEditR dataEditUI dataOutputUI dataSelectUI dataFilterUI
reformatMetadataDRui <- function(id) {
    ns <- NS(id)
    tagList(
        chevreulBox(
            title = "Reformat Metadata",
            checkboxInput(ns("header"), "Header", TRUE),
            fileInput(ns("metaFile"), "Choose CSV File of metadata with cell names in first column",
                accept = c(
                    "text/csv",
                    "text/comma-separated-values,text/plain",
                    ".csv"
                )
            ),
            actionButton(ns("updateMetadata"), "Update Metadata"),
            radioButtons(ns("updateMethod"), "Update By:", choices = c("table (below)" = "spreadsheet", "uploaded file" = "file"), inline = TRUE),
            # rhandsontable::rHandsontableOutput(ns("seuTable")),
            width = 12,
            dataSelectUI(ns("select1")),
            dataFilterUI(ns("filter1")),
            shinyjs::hidden(actionButton(ns("sync"), label = NULL, icon = icon("sync"))),
            dataOutputUI(ns("output-active")),
            dataOutputUI(ns("output-update"), icon = "file-download"),
            shinyjs::hidden(actionButton(ns("cut"), label = NULL, icon = icon("cut"))),
            dataEditUI(ns("edit1"))
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
#' @importFrom DataEditR dataSelectServer dataFilterServer dataOutputServer dataEditServer
reformatMetadataDR <- function(input, output, session, seu, featureType = "gene",
    col_bind = NULL,
    col_edit = TRUE,
    col_options = NULL,
    col_stretch = FALSE,
    col_names = TRUE,
    col_readonly = NULL,
    col_factor = FALSE,
    row_bind = NULL,
    row_edit = TRUE,
    save_as = NULL,
    title = NULL,
    logo = NULL,
    logo_size = 30,
    logo_side = "left",
    viewer = "dialog",
    viewer_height = 800,
    viewer_width = 2000,
    theme = "yeti",
    read_fun = "read.csv",
    read_args = NULL,
    write_fun = "write.csv",
    write_args = NULL,
    quiet = FALSE,
    code = FALSE,
    hide = FALSE) {
    ns <- session$ns

    # meta <- reactiveValues()
    #
    # observe({
    #   meta$old <- seu()@meta.data
    # })
    #
    # observeEvent(input$updateMetadata, {
    #   req(featureType())
    #
    #   if (input$updateMethod == "file"){
    #     inFile <- input$metaFile
    #
    #     if (is.null(inFile)){
    #       return(NULL)
    #     }
    #
    #     for (i in names(seu)){
    #       print(i)
    #       seu[[i]] <- format_new_metadata(seu[[i]], inFile$datapath)
    #     }
    #
    #   } else if (input$updateMethod == "spreadsheet"){
    #     meta$new <- propagate_spreadsheet_changes(input$seuTable)
    #   }
    #
    #   meta$new <- seu[[featureType()]]@meta.data
    #
    # })
    #
    #
    # table_out <- reactive({
    #   meta$new %||% meta$old
    # })

    # new section ------------------------------

    table_out <- reactive({
        req(seu())
        seu()@meta.data
        # seu$gene@meta.data
        # mtcars
    })

    values <- reactiveValues(
        data = NULL, data_active = NULL,
        rows = NULL, columns = NULL, cut = FALSE
    )

    observeEvent(table_out(), {
        values$rows <- NULL
        values$columns <- NULL

        values$data <- table_out() %>%
            DataEditR:::data_bind_rows(row_bind = row_bind) %>%
            DataEditR:::data_bind_cols(col_bind = col_bind) %>%
            identity()
    })

    data_select <- dataSelectServer("select1",
        data = reactive(values$data),
        hide = hide
    )
    data_filter <- dataFilterServer("filter1",
        data = reactive(values$data),
        hide = hide
    )
    observe({
        values$rows <- data_filter$rows()
        values$columns <- data_select$columns()
    })

    observe({
        if (length(values$rows) == 0 & length(values$columns) == 0) {
            values$data_active <- values$data
        } else {
            if (length(values$rows) != 0 & length(values$columns) == 0) {
                values$data_active <- values$data[values$rows, , drop = FALSE]
            } else if (length(values$rows) == 0 & length(values$columns) != 0) {
                values$data_active <- values$data[, values$columns, drop = FALSE]
            } else if (length(values$rows) != 0 & length(values$columns) != 0) {
                values$data_active <- values$data[values$rows, values$columns, drop = FALSE]
            }
        }
    })

    data_update <- dataEditServer("edit1",
        data = reactive({
            values$data_active
        }),
        col_bind = NULL, col_edit = col_edit, col_options = col_options,
        col_stretch = col_stretch, col_names = col_names,
        col_readonly = col_readonly, col_factor = col_factor,
        row_bind = NULL, row_edit = row_edit, quiet = quiet, height = viewer_height, width = viewer_width
    )
    observe({
        values$data_active <- data_update()
    })

    observeEvent(input$updateMetadata, {
        if (input$updateMethod == "file") {
            inFile <- input$metaFile

            if (is.null(inFile)) {
                return(NULL)
            }

            reformatted_seu <- format_new_metadata(seu(), inFile$datapath)
            seu(reformatted_seu)
        } else if (input$updateMethod == "spreadsheet") {
            reformatted_seu <- propagate_spreadsheet_changes(values$data_active, seu())
            seu(reformatted_seu)
        }
    })


    observeEvent(input$sync, {
        if (length(values$rows) == 0 & length(values$columns) == 0) {
            values$data <- values$data_active
        } else {
            if (length(values$rows) != 0 & length(values$columns) == 0) {
                values$data[values$rows, ] <- values$data_active
            } else if (length(values$rows) == 0 & length(values$columns) != 0) {
                values$data[, values$columns] <- values$data_active
            } else if (length(values$rows) != 0 & length(values$columns) != 0) {
                values$data[values$rows, values$columns] <- values$data_active
            }
            if (!is.null(values$data_active)) {
                if (!all(rownames(values$data_active) == rownames(values$data)[values$rows])) {
                    rownames(values$data)[values$rows] <- rownames(values$data_active)
                }
                if (!all(colnames(values$data_active) == colnames(values$data)[values$columns])) {
                    colnames(values$data)[values$columns] <- colnames(values$data_active)
                }
            }
        }
    })

    dataOutputServer("output-active",
        data = reactive({
            values$data_active
        }), save_as = "metadata.csv", write_fun = write_fun, write_args = write_args,
        hide = hide
    )
    dataOutputServer("output-update",
        data = reactive({
            values$data
        }), save_as = "metadata.csv", write_fun = write_fun, write_args = write_args,
        hide = hide
    )

    # SAVE AS
    if (!hide & !is.null(save_as)) {
        do.call(
            write_fun,
            c(list(x_edit, save_as), write_args)
        )
    }

    observeEvent(input$cut, {
        if (values$cut) {
            values$cut <- FALSE
            updateButton(session, "cut", NULL,
                block = FALSE,
                style = "danger"
            )
        } else {
            values$cut <- TRUE
            updateButton(session, "cut", NULL,
                block = FALSE,
                style = "success"
            )
        }
    })

    return(seu)
}
