test_that("heatmap is ordered by pseudotime", {

  # load packages------------------------------
  library(tidyverse)
  library(shiny)
  library(shinydashboard)
  library(SingleCellExperiment)
  library(ggraph)
  library(formattable)
  library(clustree)
  # library(seuratTools)
  # library(seuratTools, lib.loc = "/dataVolume/storage/rpkgs/devel_install/")
  devtools::load_all()
  library(Seurat)

  seu <- readRDS("~/single_cell_projects/integrated_projects/7_seq_04142020/output/seurat/unfiltered_seu.rds")

  # debug(convert_seu_to_cds)
  # debug(monocle3::preprocess_cds)
  # debug(monocle_module_heatmap)

  resolution <- 1.6

  cds <- convert_seu_to_cds(seu$gene, resolution = resolution)

  # cds <- convert_seu_to_cds(seu$transcript, resolution = resolution)

  cds <- learn_graph_by_resolution(cds,
    seu$gene,
    resolution = resolution
  )

  cds <- monocle3::order_cells(cds, root_cells = c("ds20170407_S442"))

  cds_pr_test_res <- monocle3::graph_test(cds, neighbor_graph = "principal_graph", cores = 4, expression_family = "negbinom")

  cds_pr_test_res %>%
    subset(q_value < 0.05) %>%
    dplyr::arrange(q_value) %>%
    dplyr::select(-status)

  test0 <- monocle_module_heatmap(cds, rownames(cds_pr_test_res), 1.6, collapse_rows = TRUE, group.by = "batch")

  test1 <- monocle_module_heatmap(cds, rownames(cds_pr_test_res), 1.6, collapse_rows = TRUE, group.by = "batch")

  test0$module_heatmap
})
