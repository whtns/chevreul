test_that("cluster cell metadata created", {
    resolution <- seq(0.2, 1, by = 0.2)
    clustered_object <- clustering_workflow(chevreul_sce, resolution = resolution)
    cluster_cols <- glue("{mainExpName(chevreul_sce)}_snn_res.{resolution}")
    expect_contains(colnames(colData(clustered_object)), cluster_cols)
})
