# load packages------------------------------
library(tidyverse)
library(shiny)
# devmode()
library(shinydashboard)
library(SingleCellExperiment)
library(ggraph)
library(formattable)
library(clustree)
library(fs)
library(seuratTools)
# library(seuratTools, lib.loc = "/dataVolume/storage/rpkgs/devel_install/")
library(InteractiveComplexHeatmap)
#library(monocle3)
# devtools::load_all()
# library(Seurat)

`%notin%` <- function(x,y) !(x %in% y)
# project directory----------------is--------------
proj_dir <-  fs::path("~/single_cell_projects/sc_cone_devel/sc_cone_devel_organoid/")

# apptitle------------------------------


# call app------------------------------
loom_path <- fs::path(proj_dir, "output", "velocyto", "20181001-DS-organoid-Hs_proj2.loom")

# debug(parseFilePaths)
#debug(monocle3::graph_test)
# seuratApp(proj_dir)

retained_cells <- read_csv("~/single_cell_projects/sc_cone_devel/sc_cone_devel_organoid/20181001-DS-organoid-Hs_proj/results/101322_orgSeq_noSG_coneCellList.csv")
colnames(retained_cells)[[1]] <- "sample_id"

mgCells <- c("ds20181001-0018","ds20181001-0061","ds20181001-0411","ds20181001-0417","ds20181001-0797","ds20181001-0955")

retained_cells<- retained_cells[retained_cells$sample_id %notin% mgCells,]

seu <- readRDS("~/single_cell_projects/sc_cone_devel/sc_cone_devel_organoid/20181001-DS-organoid-Hs_proj/output/seurat/011221_noLR_noiPSC_nonPRgrp_noSideGrps_byUmapPosition_seu.rds")

subset_seu <- seu[,retained_cells$sample_id]

cds <- convert_seu_to_cds(seu, resolution = 2)
# cds <- convert_seu_to_cds(seu_monocle(), resolution = input$cdsResolution)

# cds <- cds[, colnames(subset_seu)] # this is wrong and results in a scrambled cds!
cds <- cds[, colnames(cds) %in% colnames(subset_seu)] # this is the correct way to subset

cds <- threshold_monocle_genes(subset_seu, cds)

cds2 <- learn_graph_by_resolution(cds, 
                                  subset_seu,
                                 resolution = "2"
                                 )

saveRDS(cds2, "cds_script.rds")

cds3 <- monocle3::order_cells(cds2, root_cells = "ds20181001-0562") #0897, 0936 to root in cone1. 562 is new choice

#cds3 <- flip_pseudotime(cds3)

test0 <- plot_cds(cds3, color_cells_by = "gene_snn_res.2")

test1 <- plot_pseudotime(cds3, color_cells_by = "pseudotime", resolution = "gene_snn_res.2")

cds3@metadata[["diff_features"]] <- monocle3::graph_test(cds3, neighbor_graph = "principal_graph", cores = 4, expression_family = "negbinom")

heatmap_genes <- dplyr::filter(cds3@metadata[["diff_features"]], q_value < 0.05)


#SAVE pt calc cds output for use
monocle3::save_monocle_objects(cds3,"src/pseudotime/Org_Res1-8_PT_root_0936_allGenes_112522.cds")

module_heatmap <- monocle_module_heatmap(cds3, 
                                rownames(heatmap_genes), 
                                "integrated_snn_res.1.6", 
                                collapse_rows = "genes", 
                                cluster_rows = FALSE,
                                resolution = 3e-3, # this the key parameter controlling the number of modules output
                                group.by = "batch")

gene_heatmap <- monocle_module_heatmap(cds3, 
                                rownames(heatmap_genes), 
                                "integrated_snn_res.1.6", 
                                collapse_rows = "genes", 
                                resolution = 3e-3, # this the key parameter controlling the number of modules output
                                group.by = "batch")

write_csv(test1$module_table, "results/monocle_modules.csv")

ggsave("results/module_plot.pdf", ggplotify::as.ggplot(test1$module_heatmap))

ggplot(test1$module_table, aes(dim_1, dim_2, color = module)) +
    geom_point()

test3 <- 
    t(module_heatmap$agg_mat) %>% 
    as.data.frame() %>% 
    # tibble::rownames_to_column("sample_id") %>% 
    identity()

colnames(test3) <- paste0("monocle_module_", colnames(test3))

test4 <- Seurat::AddMetaData(seu, test3)

FeaturePlot(test4, cols = c('grey','red'), features = colnames(test3))
