#' Default Shiny Helper
#'
#' @param ui_element
#' @param title
#' @param content
#' @param type
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
default_helper <- function(ui_element, title = "", content = "test", type = "inline", ...){
  shinyhelper::helper(ui_element,
                      type = type,
                      title = title,
                      content = content,
                      ...)
}

#' custom collapsible box
#'
#' @param content
#' @param title
#'
#' @return
#' @export
#'
#' @examples
seuratToolsBox <- function(content, title = "", ...) {
  shinydashboard::box(
    title = title,
    status = "primary",
    solidHeader = TRUE,
    collapsible = TRUE,
    content,
    ...
  )
}
