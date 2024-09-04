#' Default Shiny Helper
#'
#' @param ui_element a shiny ui element
#' @param title a title for the helper
#' @param content content from an rmarkdown file
#' @param type style of the helper
#' @param ... extra params to helper
#'
#' @return a shinyhelper object
#'
#' @noRd
default_helper <- function(ui_element, title = "", 
                           content = "test", type = "inline", ...) {
    helper(ui_element,
        type = type,
        title = title,
        content = content,
        ...
    )
}

#' custom collapsible box
#'
#' @param content interactive content
#' @param title a title for the UI box
#'
#' @return a ui box with standard parameters
#'
#' @noRd
chevreulBox <- function(title = "", ...) {
    box(
        title = title,
        status = "primary",
        solidHeader = TRUE,
        collapsible = TRUE,
        ...
    )
}

#' Custom Dropdown Button
#'
#' @param content interactive content
#' @param title a title
#' @param ...
#'
#' @return a dropdown button
#'
#' @importFrom shinyWidgets dropdownButton
#'
#' @noRd
chevreulDropDownButton <- function(inputId, ...) {
    tagList(
        dropdownButton(
            inputId = inputId,
            ...,
            circle = TRUE, status = "info",
            icon = icon("bars"), width = "300px"
        ),
        br()
    )
}
