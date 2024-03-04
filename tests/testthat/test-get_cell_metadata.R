test_that("Metadata pulled", {
  chevreul_sce <- chevreuldata::human_gene_transcript_sce()
  meta<-get_cell_metadata(chevreul_sce)
  expect_contains(colnames(meta),colnames(colData(chevreul_sce)))
})
