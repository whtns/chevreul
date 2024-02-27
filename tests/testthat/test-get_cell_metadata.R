test_that("Metadata pulled", {
  meta<-get_cell_metadata(chevreul_sce)
  expect_contains(colnames(meta),colnames(colData(chevreul_sce)))
})
