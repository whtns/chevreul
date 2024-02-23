test_that("cluster cell metadata created", {
  resolution = seq(0.2, 2, by = 0.2)
  clustered_object <- clustering_workflow(human_gene_transcript_sce, resolution = resolution)
  cluster_cols <- glue("{mainExpName(human_gene_transcript_sce)}_snn_res.{resolution}")
  expect_contains(colnames(colData(clustered_object)), cluster_cols)
})
