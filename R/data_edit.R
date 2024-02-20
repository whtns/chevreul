#' Reformat Seurat Object Metadata UI
#'
#' @param id id
#'
#' @return UI to reformat the metadata of a single cell object
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
            width = 12,
            dataSelectUI(ns("select1")),
            dataFilterUI(ns("filter1")),
            hidden(actionButton(ns("sync"), label = NULL, icon = icon("sync"))),
            dataOutputUI(ns("output-active")),
            dataOutputUI(ns("output-update"), icon = "file-download"),
            hidden(actionButton(ns("cut"), label = NULL, icon = icon("cut"))),
            dataEditUI(ns("edit1"))
        ) %>%
            default_helper(type = "markdown", content = "reformatMetadata")
    )
}

#' Reformat Seurat Object Metadata Server
#'
#' @param object single cell object
#'
#'@noRd
reformatMetadataDR <- function(input, output, session, object, featureType = "gene",
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

    table_out <- reactive({
        req(object())
        get_cell_metadata(object())
    })

    values <- reactiveValues(
        data = NULL, data_active = NULL,
        rows = NULL, columns = NULL, cut = FALSE
    )

    observeEvent(table_out(), {
        values$rows <- NULL
        values$columns <- NULL

        values$data <- table_out() %>%
            data_bind_rows(row_bind = row_bind) %>%
            data_bind_cols(col_bind = col_bind) %>%
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

            reformatted_object <- format_new_metadata(object(), inFile$datapath)
            object(reformatted_object)
        } else if (input$updateMethod == "spreadsheet") {
            reformatted_object <- propagate_spreadsheet_changes(values$data_active, object())
            object(reformatted_object)
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

    return(object)
}

#' Add new columns to data
#'
#' from https://github.com/DillonHammill/DataEditR/blob/850feb2e29a2e94151c84300d1564ba19dcba415/R/helpers.R#L148
#'
#' @noRd
data_bind_cols <- function(data = NULL,
                           col_bind = NULL) {

  # BIND COLUMNS
  if(!is.null(data)) {
    # COLUMNS
    if (!is.null(col_bind)) {
      # NEW COLUMNS
      if (is.null(dim(col_bind))) {
        # COLUMNS AS LIST
        if (inherits(col_bind, "list")) {
          # NAMES
          if (is.null(names(col_bind))) {
            names(col_bind) <- paste0("V", length(col_bind))
          }
          # LENGTHS
          ind <- which(!unlist(lapply(col_bind, length)) == nrow(data))
          if (length(ind) > 0) {
            for (z in ind) {
              col_bind[[z]] <- rep(col_bind[[z]], nrow(data))
            }
          }
          # MATRIX
          col_bind <- do.call("cbind", col_bind)
          # COLUMN NAMES
        } else {
          col_bind <- matrix(
            rep("", nrow(data) * length(col_bind)),
            ncol = length(col_bind),
            dimnames = list(
              rownames(data),
              col_bind
            )
          )
        }
      }
      # BIND NEW COLUMNS
      data <- cbind(
        data,
        col_bind[1:nrow(data), , drop = FALSE]
      )
    }
  }

  return(data)

}

## BIND ROWS -------------------------------------------------------------------

#' Add new rows to data
#'
#' from https://github.com/DillonHammill/DataEditR/blob/850feb2e29a2e94151c84300d1564ba19dcba415/R/helpers.R#L148
#'
#' @noRd
data_bind_rows <- function(data = NULL,
                           row_bind = NULL) {

  # BIND ROWS
  if(!is.null(data)) {
    # ROWS
    if (!is.null(row_bind)) {
      # NEW ROWS
      if (is.null(dim(row_bind))) {
        # ROWS AS LIST
        if (inherits(row_bind, "list")) {
          # NAMES NOT NECESSARY
          # LENGTHS
          ind <- which(!unlist(lapply(row_bind, length)) == ncol(data))
          if (length(ind) > 0) {
            for (z in ind) {
              row_bind[[z]] <- rep(row_bind[[z]], ncol(data))
            }
          }
          # MATRIX
          row_bind <- do.call("rbind", row_bind)
          # ROW NAMES
        } else {
          row_bind <- matrix(
            rep("", ncol(data) * length(row_bind)),
            nrow = length(row_bind),
            dimnames = list(
              row_bind,
              colnames(data)
            )
          )
        }
      }
      # BIND NEW ROWS
      data <- rbind(data, row_bind[, 1:ncol(data)])
    }
  }

  return(data)

}
