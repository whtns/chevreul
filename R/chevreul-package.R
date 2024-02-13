#' @keywords internal
"_PACKAGE"

## usethis namespace: start
#' @importFrom AnnotationDbi mapIds
#' @importFrom clustree clustree
#' @importFrom ComplexHeatmap draw
#' @importFrom DBI dbConnect
#' @importFrom DBI dbDisconnect
#' @importFrom DBI dbReadTable
#' @importFrom DBI dbWriteTable
#' @importFrom dplyr arrange
#' @importFrom dplyr bind_rows
#' @importFrom dplyr case_when
#' @importFrom dplyr distinct
#' @importFrom dplyr filter
#' @importFrom dplyr group_by
#' @importFrom dplyr left_join
#' @importFrom dplyr mutate
#' @importFrom dplyr pull
#' @importFrom dplyr rename
#' @importFrom dplyr select
#' @importFrom dplyr summarize
#' @importFrom EnhancedVolcano EnhancedVolcano
#' @importFrom EnrichmentBrowser eaBrowse
#' @importFrom EnrichmentBrowser idMap
#' @importFrom EnrichmentBrowser nbea
#' @importFrom EnrichmentBrowser sbea
#' @importFrom fs dir_ls
#' @importFrom fs path
#' @importFrom fs path_abs
#' @importFrom fs path_dir
#' @importFrom fs path_expand
#' @importFrom fs path_ext_remove
#' @importFrom fs path_file
#' @importFrom graphics points
#' @importFrom graphics text
#' @importFrom grDevices dev.off
#' @importFrom grDevices pdf
#' @importFrom harmony RunHarmony
#' @importFrom methods slot
#' @importFrom methods slot<-
#' @importFrom plotly config
#' @importFrom plotly event_data
#' @importFrom plotly ggplotly
#' @importFrom plotly highlight
#' @importFrom plotly layout
#' @importFrom plotly partial_bundle
#' @importFrom plotly plotlyOutput
#' @importFrom plotly renderPlotly
#' @importFrom plotly style
#' @importFrom plotly toWebGL
#' @importFrom purrr map
#' @importFrom purrr map_if
#' @importFrom purrr set_names
#' @importFrom scran cyclone
#' @importFrom scuttle addPerCellQCMetrics
#' @importFrom shinydashboard box
#' @importFrom shinydashboard dashboardBody
#' @importFrom shinydashboard dashboardHeader
#' @importFrom shinydashboard dashboardSidebar
#' @importFrom shinydashboard menuItem
#' @importFrom shinydashboard sidebarMenu
#' @importFrom shinydashboard tabItem
#' @importFrom shinydashboard tabItems
#' @importFrom shinyhelper observe_helpers
#' @importFrom shinyjs alert
#' @importFrom shinyjs hidden
#' @importFrom shinyjs html
#' @importFrom shinyjs runcodeServer
#' @importFrom shinyjs runcodeUI
#' @importFrom shinyjs useShinyjs
#' @importFrom stats as.dendrogram
#' @importFrom stats dist
#' @importFrom stats hclust
#' @importFrom stats reorder
#' @importFrom stats runif
#' @importFrom stats setNames
#' @importFrom stringr str_remove
#' @importFrom stringr str_replace
#' @importFrom stringr str_subset
#' @importFrom tibble column_to_rownames
#' @importFrom tibble deframe
#' @importFrom tibble enframe
#' @importFrom tibble rownames_to_column
#' @importFrom tidyr drop_na
#' @importFrom utils capture.output
#' @importFrom utils head
#' @importFrom utils packageVersion
#' @importFrom utils sessionInfo
## usethis namespace: end
#' @import ggplot2
#' @import Seurat
#' @import SingleCellExperiment
#' @import Seurat
#' @import shiny
NULL

# 1. change to ~/rpkgs
# 2. run /opt/R/devel/bin/R CMD build chevreul
# 3. run /opt/R/devel/bin/R CMD check chevreul_0.5.0.tar.gz
# 4. run /opt/R/devel/bin/R
# 5. run BiocCheck::BiocCheck()
